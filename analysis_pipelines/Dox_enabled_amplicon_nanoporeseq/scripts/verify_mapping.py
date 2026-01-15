import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

GENOME_FASTA = "/Volumes/guttman/genomes/mm10/fasta/mm10.fa"
PRIMER_FILE = "ValidatedPrimers.fa"

def find_primer_genomic(primer_seq):
    start, end = 103460373 - 20000, 103483233 + 20000
    region = f"chrX:{start}-{end}"
    cmd = ["samtools", "faidx", GENOME_FASTA, region]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    seq = "".join(result.stdout.strip().split("\n")[1:]).upper()
    pos = seq.find(primer_seq)
    if pos != -1:
        return ("+", start + pos)
    pos_rc = seq.find(str(Seq(primer_seq).reverse_complement()))
    if pos_rc != -1:
        return ("-", start + pos_rc)
    return (None, None)

primers = list(SeqIO.parse(PRIMER_FILE, "fasta"))
pairs = [
    ("Exon", "AC_XistExAmp_5SNPs-F", "AC_XistExAmp_5SNPs-R"),
    ("Intron", "GG_XistInAmp1_4SNPs-F", "GG_XistInAmp1_4SNPs-R")
]

for name, f_id, r_id in pairs:
    f_seq = str(next(p for p in primers if p.id == f_id).seq).upper()
    r_seq = str(next(p for p in primers if p.id == r_id).seq).upper()

    s_f, p_f = find_primer_genomic(f_seq)
    s_r, p_r = find_primer_genomic(r_seq)

    print(f"\n--- {name} ---")
    print(f"F ({f_id}): {s_f} at {p_f}")
    print(f"R ({r_id}): {s_r} at {p_r}")
    if p_f and p_r:
        print(f"Span: chrX:{min(p_f, p_r)}-{max(p_f, p_r)}")
