from glob import glob

fasta_files = glob('./ncbi_dataset/data/*.1/*')
out_path = './phiX174.fasta'

from Bio import SeqIO
records = []
for path in fasta_files:
    for rec in SeqIO.parse(path, 'fasta'):
        records.append(rec)

SeqIO.write(records, out_path, 'fasta')