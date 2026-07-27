use std::{fs::File, io::{self, Write}, sync::{Arc, Mutex}};



#[derive(Debug)]
pub enum OutputTarget {
    Stdout(io::Stdout),
    File(File),
}


impl Write for OutputTarget {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutputTarget::Stdout(ref mut stdout) => stdout.write(buf),
            OutputTarget::File(ref mut file) => file.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputTarget::Stdout(ref mut stdout) => stdout.flush(),
            OutputTarget::File(ref mut file) => file.flush(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct OutputBuffer {
    writer: Arc<Mutex<OutputTarget>>,
    pub buffer: Vec<u8>,
    pub threshold: usize,
}

impl OutputBuffer {
    pub fn new(writer: Arc<Mutex<OutputTarget>>, threshold: usize) -> Self {
        OutputBuffer {
            writer,
            buffer: Vec::new(),
            threshold: threshold,
        }
    }

    /// Pushes whatever is buffered to the shared writer immediately. Needed when write ordering
    /// matters across buffers -- a SAM header must reach the target before any worker's records.
    pub fn flush(&mut self) {
        let mut wr = self.writer.lock().expect("Cannot lock writer");
        let _ = wr.write_all(&self.buffer);
        self.buffer.clear();
    }

    pub fn write(&mut self, str: String) {
        let _ = write!(self.buffer, "{}", str);

        if self.buffer.len() > self.threshold {
            let mut wr: std::sync::MutexGuard<'_, OutputTarget> = self.writer.lock().expect("Cannot lock writer");
            let _ = wr.write_all(&self.buffer);
            self.buffer.clear();
        }
    }
}

impl Drop for OutputBuffer {
    fn drop(&mut self) {
        let mut wr = self.writer.lock().expect("Cannot lock writer");
        let _ = wr.write_all(&self.buffer);
    }
}
