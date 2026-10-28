//! The `@SQ` block of the SAM header, kept next to the index as `<bundle>.samheader`.
//!
//! A marker-gene database has one `@SQ` line per reference -- 14.5 M lines, 350 MB on the full
//! protal database -- and formatting them from the names and lengths took ~1 s per input, on one
//! thread, before the first read was processed: more than a third of a 5-sample run. Copying the
//! finished block into the output takes ~0.06 s, so the block is formatted once per index and
//! copied into every SAM after that.
//!
//! The file is a cache, never a source: everything in it is derived from the bundle's names and
//! sequence lengths. Its first line is an `@CO` stamp naming the bundle's size and modification
//! time; when the stamp does not match the bundle (or the file is missing, or unreadable) the block
//! is formatted from the index again and the file rewritten. Deleting it is always safe.
//!
//! The stamp line itself is not copied, so the output is byte-identical to a header formatted
//! in place. A blob index has no bundle and a read-only index directory cannot take the file;
//! both keep the formatted block in memory and reuse it for the run's later inputs.

use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom, Write};
use std::sync::OnceLock;
use std::time::UNIX_EPOCH;

use crate::database::common::{DBPaths, FlexalignDatabase, IndexFormat};
use crate::io::output_buffer::OutputBuffer;
use crate::options::Options;

const STAMP_PREFIX: &str = "@CO\tflexalign samheader v1";

/// The block formatted in this run, when it could not be (or need not be) written to disk.
static IN_MEMORY: OnceLock<Vec<u8>> = OnceLock::new();

pub fn path_for(bundle: &str) -> String {
    format!("{bundle}.samheader")
}

fn stamp_of(bundle: &str) -> std::io::Result<String> {
    let meta = std::fs::metadata(bundle)?;
    let modified = meta.modified()?.duration_since(UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    Ok(format!("{STAMP_PREFIX}\tbundle_len={}\tbundle_modified_ns={modified}\n", meta.len()))
}

/// The `@SQ` lines of every reference with a length, in id order. Id 0 is flexmap's zero-length
/// sentinel and is skipped, as is any other empty entry (SAM requires LN >= 1).
pub fn format_sq_block<D: FlexalignDatabase>(db: &D) -> Vec<u8> {
    let mut out = Vec::with_capacity(64 << 20);
    let mut digits = itoa::Buffer::new();
    let mut id = 0usize;
    while let Some(name) = db.get_rname(id) {
        let length = db.get_reference_len(id).unwrap_or(0);
        if length > 0 {
            out.extend_from_slice(b"@SQ\tSN:");
            out.extend_from_slice(name.as_bytes());
            out.extend_from_slice(b"\tLN:");
            out.extend_from_slice(digits.format(length).as_bytes());
            out.push(b'\n');
        }
        id += 1;
    }
    out
}

/// The cached block for `bundle`, opened and positioned just past its stamp, if it matches the
/// bundle as it is now.
fn open_valid(bundle: &str) -> Option<File> {
    let stamp = stamp_of(bundle).ok()?;
    let mut file = File::open(path_for(bundle)).ok()?;
    let mut first = String::new();
    BufReader::new(&file).read_line(&mut first).ok()?;
    if first != stamp {
        return None;
    }
    file.seek(SeekFrom::Start(first.len() as u64)).ok()?;
    Some(file)
}

/// Writes the stamp and `block` to `<bundle>.samheader` through a temporary file and a rename,
/// so a concurrent run sees either the old file or the whole new one.
fn save(bundle: &str, block: &[u8]) -> std::io::Result<()> {
    let stamp = stamp_of(bundle)?;
    let path = path_for(bundle);
    let tmp = format!("{path}.tmp.{}", std::process::id());
    let result = (|| {
        let mut f = std::io::BufWriter::with_capacity(8 << 20, File::create(&tmp)?);
        f.write_all(stamp.as_bytes())?;
        f.write_all(block)?;
        f.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        std::fs::rename(&tmp, &path)
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    result
}

/// Writes the whole SAM header -- `@HD`, the `@SQ` block, `@PG` -- to `buffer`'s target and flushes
/// it, so it reaches the output before any worker's records.
pub fn write<D: FlexalignDatabase>(options: &Options, db: &D, buffer: &mut OutputBuffer) {
    // SO:unsorted -- workers write records as reads stream past, in no reference order.
    buffer.write_with(|b| b.extend_from_slice(b"@HD\tVN:1.6\tSO:unsorted\n"));
    buffer.flush();

    let format = crate::database::common::index_format();
    let bundle = (format != IndexFormat::Blob)
        .then(|| DBPaths::new(&options.reference).bundle_file(format));

    if let Some(block) = IN_MEMORY.get() {
        buffer.write_with(|b| b.extend_from_slice(block));
    } else if let Some(mut cached) = bundle.as_deref().and_then(open_valid) {
        if let Err(e) = buffer.copy_from(&mut cached) {
            panic!("could not copy the SAM header from {}: {e}", path_for(bundle.as_deref().unwrap()));
        }
    } else {
        let block = format_sq_block(db);
        let saved = match bundle.as_deref() {
            Some(bundle) => match save(bundle, &block) {
                Ok(()) => {
                    log::info!("SAM header formatted from the index and saved to {}", path_for(bundle));
                    true
                }
                Err(e) => {
                    log::warn!("could not save the SAM header to {}: {e}; keeping it in memory for this run", path_for(bundle));
                    false
                }
            },
            None => false,
        };
        buffer.write_with(|b| b.extend_from_slice(&block));
        if !saved {
            let _ = IN_MEMORY.set(block);
        }
    }

    buffer.write_with(|b| {
        b.extend_from_slice(b"@PG\tID:flexalign\tPN:flexalign\tVN:");
        b.extend_from_slice(env!("CARGO_PKG_VERSION").as_bytes());
        b.push(b'\n');
    });
    buffer.flush();
}
