
from Bio import SeqIO
import os

GB_FILE = "ref_seq/xist-ref-seq-with-snps-annot-1712-2295.gb"
FASTA_FILE = "ref_seq/xist_amplicon_ref.fa"

def convert_gb_to_fasta():
    if not os.path.exists(GB_FILE):
        print(f"Error: {GB_FILE} not found.")
        return

    records = list(SeqIO.parse(GB_FILE, "genbank"))
    if not records:
        print(f"Error: No records found in {GB_FILE}")
        return

    # We only have one sequence in the GB file
    record = records[0]
    # Use a clean ID for alignment
    record.id = "Xist_Amplicon"
    record.description = "Xist amplicon reference (1712-2295) with SNPs"

    SeqIO.write(record, FASTA_FILE, "fasta")
    print(f"Converted {GB_FILE} to {FASTA_FILE}")

if __name__ == "__main__":
    convert_gb_to_fasta()
