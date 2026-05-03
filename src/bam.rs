//! BAM format output with BGZF compression.
//!
//! BAM is binary compressed SAM using BGZF (Blocked GNU Zip Format).
//! BGZF adds a virtual file pointer system on top of standard gzip,
//! allowing efficient random access to bgzip-compressed files.

use crate::error::BwaError;
use crate::reference::Reference;
use crate::sam::SAMRecord;
use flate2::write::DeflateEncoder;
use flate2::Compression;
use std::io::{BufWriter, Write};
use std::path::Path;

const BGZF_MAX_BLOCK_SIZE: usize = 65536;
const BGZF_FEXTRA_LEN: u16 = 6;
const BGZF_FEXTRA_ID: [u8; 2] = [b'B', b'C'];
const BAM_HEADER_MAGIC: [u8; 4] = [b'B', b'A', b'M', b'\x01'];
const EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
];

fn generate_header_text(reference: &Reference) -> String {
    let mut text = String::new();
    text.push_str("@HD\tVN:1.6\tSO:coordinate\n");
    for contig in &reference.contigs {
        text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", contig.name, contig.len()));
    }
    text.push_str("@PG\tID:bwa-rs\tPN:bwa-rs\tVN:0.1.0\tCL:bwa-rs\n");
    text
}

pub struct BAMWriter<W: Write> {
    writer: W,
    reference: Reference,
    header_text: Vec<u8>,
}

impl BAMWriter<BufWriter<std::fs::File>> {
    pub fn from_path(path: &Path, reference: Reference) -> Result<Self, BwaError> {
        let file = std::fs::File::create(path)?;
        let writer = BufWriter::new(file);
        let mut bam_writer = Self::new(writer, reference)?;
        bam_writer.write_header()?;
        Ok(bam_writer)
    }
}

impl<W: Write> BAMWriter<W> {
    pub fn new(writer: W, reference: Reference) -> Result<Self, BwaError> {
        let header_text = generate_header_text(&reference).into_bytes();
        Ok(Self {
            writer,
            reference,
            header_text,
        })
    }

    fn write_header(&mut self) -> Result<(), BwaError> {
        let header_len = self.header_text.len() as u32;
        let n_ref = self.reference.contigs.len() as u32;

        let mut header_block = Vec::new();
        header_block.extend_from_slice(&BAM_HEADER_MAGIC);
        header_block.extend_from_slice(&header_len.to_le_bytes());
        header_block.extend_from_slice(&self.header_text);
        header_block.extend_from_slice(&n_ref.to_le_bytes());

        for contig in &self.reference.contigs {
            let name_bytes = contig.name.as_bytes();
            let name_len = (name_bytes.len() + 1) as u32;
            let ref_len = contig.len() as u32;

            header_block.extend_from_slice(&name_len.to_le_bytes());
            header_block.extend_from_slice(name_bytes);
            header_block.push(0);
            header_block.extend_from_slice(&ref_len.to_le_bytes());
        }

        let compressed = self.bgzf_compress(&header_block)?;
        self.writer.write_all(&compressed).map_err(BwaError::Io)?;
        Ok(())
    }

    fn bgzf_compress(&self, data: &[u8]) -> Result<Vec<u8>, BwaError> {
        let mut result = Vec::new();
        let mut offset = 0;

        while offset < data.len() {
            let chunk_size = (BGZF_MAX_BLOCK_SIZE - 18).min(data.len() - offset);
            let block_data = &data[offset..offset + chunk_size];

            let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(block_data).map_err(BwaError::Io)?;
            let compressed_block = encoder.finish().map_err(BwaError::Io)?;

            let total_block_size = block_data.len() + 18;
            let bs = (total_block_size - 1) as u16;

            let mut block_header = Vec::with_capacity(18);
            block_header.push(0x1f);
            block_header.push(0x8b);
            block_header.push(0x08);
            block_header.push(0x04);
            block_header.push(0x00);
            block_header.push(0x00);
            block_header.push(0x00);
            block_header.push(0x00);
            block_header.push(0x00);
            block_header.push(0xff);
            block_header.extend_from_slice(&BGZF_FEXTRA_LEN.to_le_bytes());
            block_header.extend_from_slice(&BGZF_FEXTRA_ID);
            block_header.push((bs & 0xff) as u8);
            block_header.push(((bs >> 8) & 0xff) as u8);
            block_header.push(0x00);
            block_header.push(0x00);

            let crc = crc32_slice(block_data);
            let data_len = block_data.len() as u32;

            let mut trailer = Vec::new();
            trailer.extend_from_slice(&data_len.to_le_bytes());
            trailer.extend_from_slice(&crc.to_le_bytes());

            result.extend_from_slice(&block_header);
            result.extend_from_slice(&compressed_block);
            result.extend_from_slice(&trailer);

            offset += chunk_size;
        }

        Ok(result)
    }

    pub fn write_record(&mut self, record: &SAMRecord) -> Result<(), BwaError> {
        let encoded = self.encode_record(record)?;
        let compressed = self.bgzf_compress(&encoded)?;
        self.writer.write_all(&compressed).map_err(BwaError::Io)?;
        Ok(())
    }

    fn encode_record(&self, record: &SAMRecord) -> Result<Vec<u8>, BwaError> {
        let mut data = Vec::new();

        let ref_id = self
            .reference
            .contigs
            .iter()
            .position(|c| c.name == record.rname)
            .map(|i| i as i32)
            .unwrap_or(-1);

        let ref_len = Self::cigar_to_ref_len(&record.cigar);
        let bin = Self::calculate_bin(record.pos, record.pos.saturating_add(ref_len));

        let mut l_qname = record.qname.len() as u32 + 1;
        while !l_qname.is_multiple_of(4) {
            l_qname += 1;
        }

        let seq_bytes = Self::encode_seq(&record.seq);
        let l_seq = record.seq.len().div_ceil(2) as u32;

        let qual_bytes: Vec<u8> = if record.qual == "*" {
            vec![0xff; record.seq.len()]
        } else {
            record.qual.as_bytes().to_vec()
        };

        let n_cigar = Self::count_cigar_ops(&record.cigar) as u32;

        let unrolled_data_len = 4
            + 4
            + 2
            + 4
            + 4
            + 4
            + 4
            + 4
            + 4
            + l_qname as usize
            + (l_seq as usize * 2)
            + qual_bytes.len()
            + (n_cigar as usize * 4);

        data.extend_from_slice(&(unrolled_data_len as u32).to_le_bytes());
        data.extend_from_slice(&l_qname.to_le_bytes());
        data.extend_from_slice(&bin.to_le_bytes());
        data.extend_from_slice(&((n_cigar << 16) | (l_seq & 0xffff)).to_le_bytes());
        data.extend_from_slice(&record.flag.to_le_bytes());
        data.extend_from_slice(&ref_id.to_le_bytes());
        data.extend_from_slice(&record.pos.to_le_bytes());
        data.extend_from_slice(&l_seq.to_le_bytes());
        data.extend_from_slice(&(record.mapq as u32).to_le_bytes());

        data.extend_from_slice(record.qname.as_bytes());
        while data.len() % 4 != 0 {
            data.push(0);
        }

        data.extend_from_slice(&seq_bytes);
        data.extend_from_slice(&qual_bytes);

        let cigar_values = Self::parse_cigar(&record.cigar);
        for (op, len) in cigar_values {
            data.extend_from_slice(&len.to_le_bytes());
            data.push(op);
        }

        Ok(data)
    }

    fn encode_seq(seq: &str) -> Vec<u8> {
        let len = seq.len();
        let mut result = vec![0; len.div_ceil(2)];
        for (i, c) in seq.bytes().enumerate() {
            let nibble = match c {
                b'A' | b'a' => 1,
                b'C' | b'c' => 2,
                b'G' | b'g' => 4,
                b'T' | b't' => 8,
                _ => 15,
            };
            if i % 2 == 0 {
                result[i / 2] = nibble << 4;
            } else {
                result[i / 2] |= nibble;
            }
        }
        result
    }

    fn cigar_to_ref_len(cigar: &str) -> u32 {
        let mut len = 0u32;
        let mut current = 0u32;
        for c in cigar.bytes() {
            if c.is_ascii_digit() {
                current = current * 10 + (c - b'0') as u32;
            } else if c.is_ascii_alphabetic() {
                if matches!(c, b'M' | b'D' | b'N' | b'=' | b'X') {
                    len += current;
                }
                current = 0;
            }
        }
        len
    }

    fn count_cigar_ops(cigar: &str) -> usize {
        cigar.chars().filter(|c| c.is_ascii_alphabetic()).count()
    }

    fn parse_cigar(cigar: &str) -> Vec<(u8, u32)> {
        let mut result = Vec::new();
        let mut current = 0u32;
        for c in cigar.chars() {
            if c.is_ascii_digit() {
                current = current * 10 + c.to_digit(10).unwrap();
            } else if c.is_ascii_alphabetic() {
                let op = match c {
                    'M' => 0,
                    'I' => 1,
                    'D' => 2,
                    'N' => 3,
                    'S' => 4,
                    'H' => 5,
                    'P' => 6,
                    '=' => 7,
                    'X' => 8,
                    _ => continue,
                };
                result.push((op, current));
                current = 0;
            }
        }
        result
    }

    fn calculate_bin(start: u32, end: u32) -> u16 {
        let start = start.saturating_sub(1);
        let end = end.max(start + 1);
        let mut bin: u16 = 0;
        for l in (0..15).rev().step_by(3) {
            if start >> l == end >> l {
                bin = (bin << 3) | ((start >> l) & 0x7) as u16;
            }
        }
        bin
    }

    pub fn flush(&mut self) -> Result<(), BwaError> {
        self.writer.flush().map_err(BwaError::Io)
    }

    pub fn finish(&mut self) -> Result<(), BwaError> {
        self.writer.write_all(&EOF_MARKER).map_err(BwaError::Io)?;
        self.flush()
    }
}

fn crc32_slice(data: &[u8]) -> u32 {
    let mut crc = 0xffffffffu32;
    for &byte in data {
        let idx = ((crc ^ byte as u32) & 0xff) as usize;
        crc = CRC_TABLE[idx] ^ (crc >> 8);
    }
    crc ^ 0xffffffff
}

static CRC_TABLE: [u32; 256] = [
    0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3,
    0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988, 0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
    0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
    0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
    0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172, 0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
    0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940, 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
    0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f,
    0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924, 0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
    0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
    0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
    0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e, 0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457,
    0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
    0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb,
    0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0, 0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
    0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
    0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
    0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a, 0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683,
    0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
    0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7,
    0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc, 0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
    0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
    0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
    0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236, 0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f,
    0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
    0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713,
    0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38, 0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21,
    0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
    0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
    0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2, 0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db,
    0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
    0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf,
    0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94, 0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_seq() {
        assert_eq!(
            BAMWriter::<std::io::Sink>::encode_seq("ACGT"),
            vec![0x12, 0x48]
        );
        assert_eq!(BAMWriter::<std::io::Sink>::encode_seq("A"), vec![0x10]);
        assert_eq!(BAMWriter::<std::io::Sink>::encode_seq("AC"), vec![0x12]);
    }

    #[test]
    fn test_cigar_to_ref_len() {
        assert_eq!(BAMWriter::<std::io::Sink>::cigar_to_ref_len("10M"), 10);
        assert_eq!(BAMWriter::<std::io::Sink>::cigar_to_ref_len("5M2D3M"), 10);
        assert_eq!(BAMWriter::<std::io::Sink>::cigar_to_ref_len("5M2I3M"), 8);
    }

    #[test]
    fn test_count_cigar_ops() {
        assert_eq!(BAMWriter::<std::io::Sink>::count_cigar_ops("10M"), 1);
        assert_eq!(BAMWriter::<std::io::Sink>::count_cigar_ops("5M2D3M"), 3);
    }

    #[test]
    fn test_parse_cigar() {
        let ops = BAMWriter::<std::io::Sink>::parse_cigar("5M2D3M");
        assert_eq!(ops, vec![(0, 5), (2, 2), (0, 3)]);
    }

    #[test]
    fn test_calculate_bin() {
        assert_eq!(BAMWriter::<std::io::Sink>::calculate_bin(1, 10), 0);
        assert_eq!(BAMWriter::<std::io::Sink>::calculate_bin(1, 2), 0);
    }

    #[test]
    fn test_crc32() {
        let data = b"test";
        let result = crc32_slice(data);
        assert_ne!(result, 0);
        let data2 = b"";
        assert_eq!(crc32_slice(data2), 0);
    }
}
