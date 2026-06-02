use bwa_mem::paired::{is_proper_pair, mate_fields, InsertSizeDistribution};
use bwa_mem::sam::write_paired_record;
use bwa_mem::types::{AlignmentResult, Cigar, CigarOp};

fn mapped_read(position: usize, ref_len: u32, reverse: bool) -> AlignmentResult {
    let mut cigar = Cigar::new();
    cigar.push(CigarOp::M, ref_len);
    let mut r = AlignmentResult::new(position, cigar);
    r.reverse_strand = reverse;
    r
}

fn unmapped_read() -> AlignmentResult {
    let mut r = AlignmentResult::new(0, Cigar::new());
    r.flag = 0x4;
    r
}

fn parse_sam_fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

#[test]
fn proper_fr_pair_flags_and_rnext() {
    // r1 at 100 (fwd, 100M), r2 at 300 (rev, 100M) — classic FR proper pair
    let r1 = mapped_read(100, 100, false);
    let r2 = mapped_read(300, 100, true);

    let mut dist = InsertSizeDistribution::new();
    dist.add(300); // one observation so bounds are computed

    let proper = is_proper_pair(&r1, &r2, &dist);
    assert!(proper, "FR pair within distribution should be proper");

    let mf1 = mate_fields(&r1, &r2, true, proper);
    let mf2 = mate_fields(&r2, &r1, false, proper);

    // --- flag checks for read 1 ---
    assert_ne!(mf1.flag & 0x1, 0, "r1: paired");
    assert_ne!(mf1.flag & 0x2, 0, "r1: proper pair");
    assert_eq!(mf1.flag & 0x4, 0, "r1: mapped");
    assert_eq!(mf1.flag & 0x8, 0, "r1: mate mapped");
    assert_eq!(mf1.flag & 0x10, 0, "r1: forward strand");
    assert_ne!(mf1.flag & 0x20, 0, "r1: mate on reverse strand");
    assert_ne!(mf1.flag & 0x40, 0, "r1: first in pair");
    assert_eq!(mf1.flag & 0x80, 0, "r1: not second in pair");

    // --- flag checks for read 2 ---
    assert_ne!(mf2.flag & 0x1, 0, "r2: paired");
    assert_ne!(mf2.flag & 0x2, 0, "r2: proper pair");
    assert_eq!(mf2.flag & 0x4, 0, "r2: mapped");
    assert_eq!(mf2.flag & 0x8, 0, "r2: mate mapped");
    assert_ne!(mf2.flag & 0x10, 0, "r2: reverse strand");
    assert_eq!(mf2.flag & 0x20, 0, "r2: mate on forward strand");
    assert_eq!(mf2.flag & 0x40, 0, "r2: not first in pair");
    assert_ne!(mf2.flag & 0x80, 0, "r2: second in pair");

    // RNEXT = "=" for both (mates are mapped)
    assert_eq!(mf1.rnext, "=");
    assert_eq!(mf2.rnext, "=");

    // PNEXT: r1's pnext = r2's 1-based pos, and vice versa
    assert_eq!(mf1.pnext, 301, "r1 pnext = r2 position + 1");
    assert_eq!(mf2.pnext, 101, "r2 pnext = r1 position + 1");

    // TLEN: non-zero, and r1/r2 are negatives of each other
    assert_ne!(mf1.tlen, 0, "tlen nonzero for mapped pair");
    assert_eq!(mf1.tlen, -mf2.tlen, "tlen antisymmetric");
    assert!(mf1.tlen > 0, "left-anchor read tlen is positive");
    assert!(mf2.tlen < 0, "right-anchor read tlen is negative");

    // placed_pos: None for both (both are mapped)
    assert_eq!(mf1.placed_pos, None);
    assert_eq!(mf2.placed_pos, None);
}

#[test]
fn improper_pair_flags() {
    // Same-strand pair — cannot be proper
    let r1 = mapped_read(100, 100, false);
    let r2 = mapped_read(300, 100, false);
    let dist = InsertSizeDistribution::with_params(300.0, 50.0);

    let proper = is_proper_pair(&r1, &r2, &dist);
    assert!(!proper);

    let mf1 = mate_fields(&r1, &r2, true, proper);
    assert_ne!(mf1.flag & 0x1, 0, "still paired");
    assert_eq!(mf1.flag & 0x2, 0, "not proper pair");
    assert_ne!(mf1.flag & 0x40, 0, "first in pair");
}

#[test]
fn unmapped_mate_inherits_coord() {
    let r1 = mapped_read(100, 100, false);
    let r2 = unmapped_read();

    let mf1 = mate_fields(&r1, &r2, true, false);
    let mf2 = mate_fields(&r2, &r1, false, false);

    // r1 is mapped; r2 is unmapped
    assert_eq!(mf1.flag & 0x8, 0x8, "r1 sees mate unmapped");
    assert_eq!(mf2.flag & 0x4, 0x4, "r2 is unmapped");

    // The unmapped mate (r2) should have placed_pos = Some(100) so the CLI
    // emits RNAME=ref and POS=101 instead of RNAME=* POS=0.
    assert_eq!(mf2.placed_pos, Some(100));
    assert_eq!(mf2.rnext, "=");
    assert_eq!(mf2.pnext, 101);
}

#[test]
fn both_unmapped_pair() {
    let r1 = unmapped_read();
    let r2 = unmapped_read();

    let mf1 = mate_fields(&r1, &r2, true, false);
    let mf2 = mate_fields(&r2, &r1, false, false);

    assert_ne!(mf1.flag & 0x4, 0, "r1 unmapped");
    assert_ne!(mf1.flag & 0x8, 0, "r1 mate unmapped");
    assert_eq!(mf1.rnext, "*");
    assert_eq!(mf1.pnext, 0);
    assert_eq!(mf1.tlen, 0);
    assert_eq!(mf1.placed_pos, None);

    assert_ne!(mf2.flag & 0x4, 0, "r2 unmapped");
    assert_ne!(mf2.flag & 0x8, 0, "r2 mate unmapped");
    assert_eq!(mf2.rnext, "*");
    assert_eq!(mf2.pnext, 0);
}

#[test]
fn write_paired_record_output_format() {
    let r1 = mapped_read(100, 100, false);
    let r2 = mapped_read(300, 100, true);
    let dist = InsertSizeDistribution::with_params(300.0, 50.0);
    let proper = is_proper_pair(&r1, &r2, &dist);

    let mf1 = mate_fields(&r1, &r2, true, proper);
    let mf2 = mate_fields(&r2, &r1, false, proper);

    let mut buf1: Vec<u8> = Vec::new();
    write_paired_record(&mut buf1, "read1", &r1, &mf1).unwrap();
    let line1 = String::from_utf8(buf1).unwrap();
    let f1 = parse_sam_fields(line1.trim_end_matches('\n'));

    assert_eq!(f1[0], "read1", "QNAME");
    let flag1: u16 = f1[1].parse().unwrap();
    assert_ne!(flag1 & 0x1, 0, "paired");
    assert_ne!(flag1 & 0x40, 0, "read1 bit");
    assert_eq!(f1[2], "ref", "RNAME for mapped read");
    assert_eq!(f1[3], "101", "POS = position+1");
    assert_eq!(f1[6], "=", "RNEXT = =");
    assert_eq!(f1[7], "301", "PNEXT = r2 position+1");
    let tlen1: i64 = f1[8].parse().unwrap();
    assert_ne!(tlen1, 0, "TLEN nonzero");
    assert_eq!(f1[9], "*", "SEQ stays *");
    assert_eq!(f1[10], "*", "QUAL stays *");

    let mut buf2: Vec<u8> = Vec::new();
    write_paired_record(&mut buf2, "read1", &r2, &mf2).unwrap();
    let line2 = String::from_utf8(buf2).unwrap();
    let f2 = parse_sam_fields(line2.trim_end_matches('\n'));

    let tlen2: i64 = f2[8].parse().unwrap();
    assert_eq!(tlen1, -tlen2, "TLEN antisymmetric");

    let flag2: u16 = f2[1].parse().unwrap();
    assert_ne!(flag2 & 0x80, 0, "read2 bit set");
    assert_eq!(f2[6], "=", "RNEXT = =");
}

#[test]
fn write_paired_record_unmapped_mate_convention() {
    // Unmapped r2 should get RNAME=ref, POS=r1+1 via placed_pos
    let r1 = mapped_read(99, 100, false);
    let r2 = unmapped_read();

    let mf2 = mate_fields(&r2, &r1, false, false);
    assert_eq!(mf2.placed_pos, Some(99), "placed_pos = r1's 0-based pos");

    let mut buf: Vec<u8> = Vec::new();
    write_paired_record(&mut buf, "read2", &r2, &mf2).unwrap();
    let line = String::from_utf8(buf).unwrap();
    let f = parse_sam_fields(line.trim_end_matches('\n'));

    assert_eq!(f[2], "ref", "unmapped-mate inherits RNAME=ref");
    assert_eq!(f[3], "100", "POS = placed_pos+1 = 100");
    assert_eq!(f[4], "0", "MAPQ 0 for unmapped");
    assert_eq!(f[5], "*", "CIGAR * for unmapped");

    let flag: u16 = f[1].parse().unwrap();
    assert_ne!(flag & 0x4, 0, "unmapped bit set");
    assert_ne!(flag & 0x80, 0, "second-in-pair bit set");
}
