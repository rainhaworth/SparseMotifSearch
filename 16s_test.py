import numpy as np
from Bio import SeqIO
from train import mut_seqs
import polars as pl
import argparse
import pickle
from sklearn.metrics import confusion_matrix, balanced_accuracy_score

# polars GPU retrieval
def get_hits(seqs, motifs):
    exprs = [pl.col('seqs').str.contains(motifs[i].replace('-','[ACGU]')).alias(f'm_{i}') for i in range(len(motifs))]
    result = pl.LazyFrame({'seqs': seqs}).select(exprs).collect(engine='gpu')
    
    return result.to_numpy()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', type=str, default='./RNA/SILVA_138.2_SSURef_tax_silva.fasta')
    parser.add_argument('--train', type=str, default='./RNA/SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta') # for train ids
    parser.add_argument('-m', type=str, default='/RNA/16s.csv') # motifs
    parser.add_argument('-c', type=str, default='clf-r0.1.pickle') # classifier binary
    parser.add_argument('-r', type=float, default=0.1) # mutation rate
    args = parser.parse_args()

    # load motifs
    motifs = []
    with open(args.m, 'r') as f:
        for line in f.readlines():
            if line[0] == 'm': continue
            motifs.append(line.split(',', 1)[0])

    # load classifier
    with open(args.c, 'rb') as f:
        clf = pickle.load(f)

    # get train ids
    train_ids = set()
    for rec in SeqIO.parse(args.train, 'fasta'):
        train_ids.add(str(rec.description))

    # get ground truth test samples
    gt = []
    for rec in SeqIO.parse(args.test, 'fasta'):
        if rec.description in train_ids: continue
        gt.append(str(rec.seq))

    # mutate
    mut = mut_seqs(gt, args.r)

    # make dataset
    X = np.concat([get_hits(gt, motifs), get_hits(mut, motifs)], axis=0)
    y = [1] * len(gt) + [0] * len(gt)

    # run classifier and print metrics
    y_pred = clf.predict(X)
    print(confusion_matrix(y, y_pred))
    print('BA:', balanced_accuracy_score(y, y_pred))
        