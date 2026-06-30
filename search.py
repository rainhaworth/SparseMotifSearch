# full search algorithm
import utils
import argparse
import random

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, default='SILVA_SSU_tiny.fasta') # input file
    parser.add_argument('-o', type=str, default='out-motifs.csv') # output file
    parser.add_argument('-l', type=int, default=-1) # limit on seqs to load
    parser.add_argument('-b', type=int, default=10) # max batch size
    parser.add_argument('-n', type=int, default=10) # number of batches to process
    parser.add_argument('-w', type=int, default=10) # min weight
    parser.add_argument('-c', type=int, default=8) # min count
    parser.add_argument('--mode', default='fft', choices=['direct','fft','gpu'])
    args = parser.parse_args()

    random.seed(42)

    # overwrite output file
    with open(args.o, 'w') as f:
        f.write('motif,weight,len,hits,\n')

    match args.mode:
        case 'direct': aligner = utils.ShiftAlign
        case 'fft': aligner = utils.ShiftAlignFFT
        case 'gpu': quit()

    # load
    seqs = utils.load(args.i, args.l)

    # if seqs fit in one batch, solve exactly
    if len(seqs) <= args.b:
        print('solving exactly')
        # TODO: we should probably pass in the sequences directly along with a class then run the constructor in solve_batch
        batch = aligner(seqs)
        motifs = utils.solve_batch(batch, args.w, args.c)
        utils.write_motifs(args.o, motifs)
        print(len(motifs), 'motifs found')
    # otherwise, run multiple times
    else:
        for i in range(args.n):
            print('batch', i)
            seqs_batch = random.sample(seqs, args.b)
            batch = utils.ShiftAlignFFT(seqs_batch)
            motifs = utils.solve_batch(batch, args.w, args.c)

            # find hits on entire dataset
            for i in range(len(motifs)):
                hits = utils.find_hits(motifs[i], seqs)
                motifs[i] = (motifs[i][0], hits)

            utils.write_motifs(args.o, motifs)
            print(len(motifs), 'motifs found')

    print('finished running')
        