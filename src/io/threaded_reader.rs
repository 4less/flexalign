//! A reader whose underlying source is decoded on its own thread.
//!
//! The FASTQ workers take batches from one shared reader under a lock, and with gzip input that
//! lock held a single-threaded inflate: on a 2 M-pair sample the workers' CPU time fell from 36 s
//! (plain FASTQ) to 22 s (gzip) for a LONGER wall time, because they were waiting on it. Here a
//! thread decodes ahead into a bounded queue of chunks, so each mate file inflates on a thread of
//! its own, in parallel with each other and with the workers. Chunk buffers are recycled through
//! a return queue: a fresh allocation per chunk cost page faults and zeroing that showed up as
//! +20% wall time. Plain files are read directly; this is for decoders.

use std::io::{self, Read};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread;

/// Bytes per chunk handed over the queue, and how many chunks may wait in it.
const CHUNK: usize = 4 << 20;
const QUEUE: usize = 4;

pub struct ThreadedReader {
    rx: Receiver<io::Result<Vec<u8>>>,
    back: SyncSender<Vec<u8>>,
    current: Vec<u8>,
    pos: usize,
    done: bool,
}

impl ThreadedReader {
    pub fn new<R: Read + Send + 'static>(mut inner: R) -> Self {
        let (tx, rx) = sync_channel::<io::Result<Vec<u8>>>(QUEUE);
        let (back, reuse) = sync_channel::<Vec<u8>>(QUEUE + 2);
        thread::spawn(move || {
            let mut spare = 0usize;
            loop {
                // A recycled buffer when one is back, a new one while the pool is still filling.
                let mut buf = match reuse.try_recv() {
                    Ok(b) => b,
                    Err(_) if spare < QUEUE + 2 => { spare += 1; Vec::with_capacity(CHUNK) }
                    Err(_) => match reuse.recv() { Ok(b) => b, Err(_) => return },
                };
                buf.resize(CHUNK, 0);
                let mut filled = 0;
                let err = loop {
                    match inner.read(&mut buf[filled..]) {
                        Ok(0) => break None,
                        Ok(n) => {
                            filled += n;
                            if filled == CHUNK {
                                break None;
                            }
                        }
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                        Err(e) => break Some(e),
                    }
                };
                let at_end = filled < CHUNK;
                if filled > 0 {
                    buf.truncate(filled);
                    if tx.send(Ok(buf)).is_err() {
                        return;
                    }
                }
                if let Some(e) = err {
                    let _ = tx.send(Err(e));
                    return;
                }
                if at_end {
                    return;
                }
            }
        });
        Self { rx, back, current: Vec::new(), pos: 0, done: false }
    }
}

impl Read for ThreadedReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.pos == self.current.len() {
            if self.done {
                return Ok(0);
            }
            // Hand the drained buffer back for reuse (the thread may already be gone).
            let drained = std::mem::take(&mut self.current);
            if drained.capacity() > 0 {
                let _ = self.back.try_send(drained);
            }
            self.pos = 0;
            match self.rx.recv() {
                Ok(Ok(chunk)) => self.current = chunk,
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => {
                    self.done = true;
                    return Ok(0);
                }
            }
        }
        let n = out.len().min(self.current.len() - self.pos);
        out[..n].copy_from_slice(&self.current[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}
