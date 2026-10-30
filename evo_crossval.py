
# full search algorithm
import utils
import argparse
import random
from find import search, make_search_parser
from train import find_seq_hits, train_test_split

if __name__ == '__main__':
    parser = make_search_parser()
    parser.add_argument('--nat', type=str, default='./phage/ncbi_dataset_phiX174/ncbi_phiX174.fasta') # natural seqs fasta
    parser.add_argument('--evo', type=str, default='./phage/ai_phiX174.fasta') # evo seqs fasta
    parser.add_argument('--cv', type=int, default=5) # cross-val fold
    args = parser.parse_args()

    random.seed(42)

    # load
    seqs_nat = utils.load(args.nat, -1)
    seqs_evo = utils.load(args.evo, -1)

    print('found', len(seqs_nat), 'natural seqs and', len(seqs_evo), 'Evo seqs')

    for i in range(args.cv):
        print('trial', i)
        trial_prefix = './evo-trial' + str(i) + '-w' + str(args.w)

        # train-test split
        seqs_train, seqs_val = train_test_split(seqs_nat)

        print('finding motifs')

        args.o = trial_prefix + '-mot.csv'
        motifs = search(seqs_train, args)

        print('counting hits in each seq')

        nat_hits = [sum(x) for x in find_seq_hits(seqs_train, motifs)]
        val_hits = [sum(x) for x in find_seq_hits(seqs_val, motifs)]
        evo_hits = [sum(x) for x in find_seq_hits(seqs_evo, motifs)]

        with open(trial_prefix + '-seq.txt', 'w') as f:
            f.write(','.join(str(x) for x in nat_hits) + '\n')
            f.write(','.join(str(x) for x in val_hits) + '\n')
            f.write(','.join(str(x) for x in evo_hits) + '\n')