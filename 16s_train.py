import numpy as np
from Bio import SeqIO
from train import mut_seqs
import polars as pl
import argparse
import pickle
from sklearn.linear_model import LogisticRegression
import pickle

# polars GPU retrieval
def get_hits(seqs, motifs):
    exprs = [pl.col('seqs').str.contains(motifs[i].replace('-','[ACGU]')).alias(f'm_{i}') for i in range(len(motifs))]
    result = pl.LazyFrame({'seqs': seqs}).select(exprs).collect(engine='gpu')
    
    return result.to_numpy()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--train', type=str, default='./RNA/SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta')
    parser.add_argument('-m', type=str, default='./RNA/16s.csv') # output file
    args = parser.parse_args()

    print('loading seqs and motifs')
    
    fn_mot = args.m
    fn_train = args.train    
    rs = [0.1, 0.05, 0.01]
    
    motifs = []
    with open(fn_mot, 'r') as f:
        for line in f.readlines():
            if line[0] == 'm': continue
            motifs.append(line.split(',', 1)[0])
    
    train_gt = []
    train_ids = []
    for rec in SeqIO.parse(fn_train, 'fasta'):
        train_gt.append(str(rec.seq))
        train_ids.append(str(rec.description))

    print('finding hits in ground truth set')

    train_gt_hits = get_hits(train_gt, motifs)

    print('finding hits in mutated sets and training classifiers')
    
    train_muts = [mut_seqs(train_gt, r) for r in rs]
    for i in range(len(rs)):
        train_mut = train_muts[i]
        train_mut_hits = get_hits(train_mut, motifs)
        
        X = np.concat([train_gt_hits, train_mut_hits], axis=0)
        y = [1] * len(train_gt) + [0] * len(train_gt)
        
        clf = LogisticRegression(max_iter=200)
        clf.fit(X, y)
        with open('clf-r' + str(rs[i]) + '.pickle', 'wb') as f:
            pickle.dump(clf, f)

    print('done')
        