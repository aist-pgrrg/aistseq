import sys
from Bio import SeqIO

def calculate_kilobases(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            total_length += len(record.seq)
    
    kilobases = total_length / 1000
    return kilobases

def main():
    if len(sys.argv) != 2:
        print("Usage: python calculate_fasta.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    kb = calculate_kilobases(fasta_file)
    return kb

if __name__ == "__main__":
    size = main()
    print(size)
