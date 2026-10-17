# filter GENCODE lncRNA to contain only sequences in EVL
from Bio import SeqIO

gencode_path = 'gencode.v50.lncRNA_transcripts.fa'
evl_path = 'EVL_human_lncRNA.txt'
out_path = 'human_lncRNA_filtered.fasta'

# get EVL accessions
accessions = set()
with open(evl_path, 'r') as f:
    for line in f.readlines():
        accessions.add(line[:-1])

# get GENCODE sequences matching accessions
records = []
for rec in SeqIO.parse(gencode_path, 'fasta'):
    desc = rec.description
    if desc.split('|', 1)[1].split('.', 1)[0] in accessions:
        records.append(rec)

# write
SeqIO.write(records, out_path, 'fasta')