# full search algorithm
import utils
import argparse
import random
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, default='SILVA_SSU_tiny.fasta') # input file
    parser.add_argument('-o', type=str, default='out-motifs.csv') # output file
    parser.add_argument('-l', type=int, default=-1) # limit on seqs to load
    parser.add_argument('-b', type=int, default=20) # max batch size
    parser.add_argument('-c', type=float, default=0.75) # min conservation frac
    parser.add_argument('-n', type=int, default=10) # number of batches to process
    parser.add_argument('-w', type=int, default=10) # min weight
    parser.add_argument('--mode', default='dense', choices=['exact', 'dense', 'sparse'])
    # TODO: choice between refine or w_mer
    # TODO: choice of "seed" type (seq vs pair)
    args = parser.parse_args()

    random.seed(42)

    # overwrite output file
    with open(args.o, 'w') as f:
        f.write('motif,weight,len,hits,\n')

    # load
    seqs = utils.load(args.i, args.l)
    seqs_OH, chars = utils.prepare(seqs)

    print(len(seqs), 'seqs loaded')
    
    iters = args.n
    bsz = args.b
    w = args.w

    if args.mode == 'exact':
        print('not yet implemented')

    elif args.mode == 'dense':
        # if seqs fit in one batch, solve exactly + override iters
        if len(seqs) <= bsz:
            print('solving exactly')
            bsz = len(seqs)
            iters = 1

        c_int = int(args.c * bsz)
        print(bsz, c_int)

        for iter in range(iters):
            print('batch', iter)
            seqs_batch = random.sample(seqs_OH, bsz)

            # find motifs
            motifs = []
            t = tqdm(total=(bsz**2 - bsz) // 2)
            for i in range(bsz-1):
                for j in range(i+1,bsz):
                    t.update()
                    # get initial candidate and set of all seqs containing an exact match
                    motif = utils.align(seqs_batch[i], seqs_batch[j])
                    targets = [seqs_batch[k] for k in range(bsz) if k not in {i,j}] # TODO: profiling, replace this if slow

                    # adjust c_int to account for 2 existing matches
                    motif = utils.refine(motif, targets, chars, w, c_int - 2)
                    if motif is not None: motifs.append(motif)
            t.close()

            # deduplicate
            motifs = list(set(motifs))

            # find hits on entire dataset
            hits = utils.find_hits(motifs, seqs)

            utils.write_motifs(args.o, zip(motifs, hits))
            print(len(motifs), 'motifs found')
    else:
        # sparse mode
        c_int = int(args.c * len(seqs))

        motifs = []
        for i in tqdm(range(iters)):
            q_idx = random.randint(0,len(seqs))
            query = seqs_OH[q_idx]
            targets = [seqs_OH[i] for i in range(len(seqs_OH)) if i != q_idx]

            motifs_new = utils.get_w_mers(query, targets, chars, w, c_int)
            if motifs_new is not None: motifs += motifs_new

        motifs = list(set(motifs))
        hits = utils.find_hits(motifs, seqs)

        utils.write_motifs(args.o, zip(motifs, hits))
        print(len(motifs), 'motifs found')

        
        

    print('finished running')
        
