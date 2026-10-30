
# full search algorithm
import utils
import argparse
import random
from tqdm import tqdm
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import IsolationForest
from sklearn.metrics import confusion_matrix, accuracy_score
from find import search, make_search_parser
import re

def mut_seqs(seqs, r=0.1):
    assert 0.0 <= r <= 1.0
    chars = list(set(''.join(seqs)))
    seqs_mut = []

    for seq in seqs:
        n_mut = int(len(seq) * r)
        seq_list = list(seq)

        mutation_indices = random.sample(range(len(seq)), n_mut)

        for index in mutation_indices:
            current_base = seq_list[index]
            valid_bases = list(chars)
            if current_base not in chars:
                continue
            else:
                valid_bases.remove(current_base) # Ensure the new base is different.
                new_base = random.choice(valid_bases)
                seq_list[index] = new_base

        seqs_mut.append("".join(seq_list))  # Convert back to a string
    return seqs_mut

def find_seq_hits(seqs, motifs):
    # make regex for motif
    motif_strs = motifs

    regexes = []
    for motif_str in motif_strs:
        re_str = ''
        wildcards = 0
        for i in range(len(motif_str)):
            if motif_str[i] == '-': wildcards += 1
            else:
                if wildcards != 0:
                    re_str += '.{' + str(wildcards) + '}'
                    wildcards = 0
                re_str += motif_str[i]
        regex = re.compile(re_str)
        regexes.append(regex)
    
    # search
    hits = []
    for seq in seqs:
        hits.append([r.search(seq) != None for r in regexes])
    
    return hits 

def train_test_split(seqs, split=0.7):
    train_idx = random.sample(range(len(seqs)), int(split * len(seqs)))
    seqs_train = []
    seqs_test = []
    for i, seq in enumerate(seqs):
        if i in train_idx: seqs_train.append(seq)
        else: seqs_test.append(seq)
    return seqs_train, seqs_test

if __name__ == '__main__':
    parser = make_search_parser()
    parser.add_argument('-i', type=str, default='SILVA_SSU_tiny.fasta') # input fasta path
    parser.add_argument('-l', type=int, default=-1) # limit on seqs to load
    args = parser.parse_args()

    random.seed(42)

    # load
    seqs = utils.load(args.i, args.l)

    print(len(seqs), 'seqs loaded')

    # train-test split -> make mutated sets
    seqs_train_gt, seqs_test_gt = train_test_split(seqs)
    seqs_train_mut = mut_seqs(seqs_train_gt, args.r)
    seqs_test_mut = mut_seqs(seqs_test_gt, args.r)

    print(len(seqs_train_gt), len(seqs_test_gt))

    print('finding motifs')

    motifs = search(seqs_train_gt, args)

    if len(motifs) > 0:
        print('training classifier')

        def gt_mut_to_X_y(gt, mut):
            X = find_seq_hits(gt, motifs) + find_seq_hits(mut, motifs)
            y = [1] * len(gt) + [0] * len(mut)
            return X, y
        
        X_train, y_train = gt_mut_to_X_y(seqs_train_gt, seqs_train_mut)
        X_test, y_test = gt_mut_to_X_y(seqs_test_gt, seqs_test_mut)

        clf = LogisticRegression(max_iter=200)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        print(confusion_matrix(y_test, y_pred))
        print('acc:', accuracy_score(y_test, y_pred))

        print('finished running')