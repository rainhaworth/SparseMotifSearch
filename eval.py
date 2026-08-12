# find all exact matches to set of motifs within target dataset
import utils
import argparse
import random

if __name__ == '__main__':
    print('made it here')
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='SILVA_SSU_tiny.fasta') # input FASTA file
    parser.add_argument('-m', type=str, default='out-motifs.csv') # input motif file
    parser.add_argument('-o', type=str, default='out-matches.csv') # output file
    parser.add_argument('-l', type=int, default=-1) # limit on seqs to load
    args = parser.parse_args()

    random.seed(42)

    # overwrite output file
    with open(args.o, 'w') as f:
        f.write('motif,hits,\n')

    # load seqs
    seqs = utils.load(args.f, args.l)

    print(len(seqs), 'seqs found')
    
    # load motifs
    motifs = []
    with open(args.m, 'r') as f:
        for line in f.readlines():
            if 'motif' in line: continue
            
            fields = line.split(',')
            if len(fields) == 0: break
            
            motifs.append(fields[0])

    print(len(motifs), 'motifs found')

    # find hits on entire dataset
    hits = utils.find_hits([(motif,None) for motif in motifs], seqs)
    with open(args.o, 'a')  as f:
        for i in range(len(motifs)):
            f.write(motifs[i] + ',' + str(len(hits[i])) + '\n')

    print('finished running')
