# full search algorithm
import utils
import argparse
import random
from tqdm import tqdm

# inner loop helper function; output is always list
def find_motifs(seqlist, idxs, algo, chars, w, c_int):
    assert 1 <= len(idxs) <= 2
    if len(idxs) == 1: motif = seqlist[idxs[0]]
    else: motif = utils.align(seqlist[idxs[0]], seqlist[idxs[1]])

    targets = [seqlist[k] for k in range(len(seqlist)) if k not in set(idxs)] # TODO: profiling, replace this if slow

    # adjust c_int to account for matches at dxs
    return algo(motif, targets, chars, w, c_int - len(idxs))

# run full search; define here for portability
def search(seqs, args):
    # convert args to local vars + map algo choice to actual function
    out_path = args.o

    c = args.c
    w = args.w
    density = args.density
    ref_type = args.reference
    n_sparse = args.n

    bsz = args.b
    n_batch = args.m

    if args.algorithm == 'iterative':
        algo = utils.refine
    else:
        algo = utils.get_w_mers


    # check some loose bounds
    assert 0.0 <= c <= 1.0
    assert w > 0
    assert n_sparse > 0
    assert bsz > 1
    assert n_batch > 0
    assert len(seqs) > 1


    print('preprocessing')

    seqs_OH, chars = utils.prepare(seqs)

    # if not batching or seqs all fit in one batch, adjust bsz and n_batch
    if len(seqs) <= bsz or not args.batch:
        bsz = len(seqs)
        n_batch = 1

    # convert c to int occurrence count
    c_int = int(c * bsz)


    print('finding motifs')

    motifs = []
    for batch_i in range(n_batch):
        print('batch', batch_i)

        # get batch
        if bsz == len(seqs): seqs_batch = seqs_OH
        else: seqs_batch = random.sample(seqs_OH, bsz)

        # sparse search
        if density == 'sparse':
            for i in tqdm(range(n_sparse)):
                if ref_type == 'pair': idxs = random.sample(range(bsz), 2)
                else: idxs = [random.randint(0, bsz-1)]

                motifs += find_motifs(seqs_batch, idxs, algo, chars, w, c_int)

        # dense search
        else:
            # easiest to use different implementations for pairwise alignment vs seq refs
            if ref_type == 'pair':
                t = tqdm(total=(bsz**2 - bsz) // 2)
                for i in range(bsz-1):
                    for j in range(i+1,bsz):
                        t.update()
                        motifs += find_motifs(seqs_batch, [i,j], algo, chars, w, c_int)
                t.close()

            else:
                for i in tqdm(range(bsz)):
                    motifs += find_motifs(seqs_batch, [i], algo, chars, w, c_int)

        # get unique motifs
        motifs = list(set(motifs))
        print('motif subtotal:', len(motifs))

    print('retrieving', len(motifs), 'motifs and writing to', out_path)

    hits = utils.find_hits(motifs, seqs)

    # overwrite file, output header
    with open(out_path, 'w') as f:
        f.write('motif,weight,len,batch_hits,\n')
    utils.write_motifs(out_path, zip(motifs, hits))
    return motifs

# minimal argparse to provide all args needed for search
def make_search_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', type=str, default='out-motifs.csv') # output motif path

    parser.add_argument('-a', '--algorithm', default='iterative', choices=['iterative', 'column'])
    parser.add_argument('-c', type=float, default=0.75) # c hat (conservation target, as fraction)
    parser.add_argument('-w', type=int, default=10) # w hat (exact or min weight)

    parser.add_argument('-d', '--density', default='dense', choices=['dense', 'sparse'])
    parser.add_argument('-n', type=int, default=10) # for sparse search: number of sequences to search over

    parser.add_argument('--batch', action='store_true')
    parser.add_argument('-b', type=int, default=20) # max batch size
    parser.add_argument('-m', type=int, default=10) # number of batches to process

    parser.add_argument('-r', '--reference', default='pair', choices=['seq', 'pair'])
    return parser

if __name__ == '__main__':
    parser = make_search_parser()
    parser.add_argument('-i', type=str, default='SILVA_SSU_tiny.fasta') # input fasta path
    parser.add_argument('-l', type=int, default=-1) # limit on seqs to load
    args = parser.parse_args()

    random.seed(42)

    # load
    seqs = utils.load(args.i, args.l)
    print(len(seqs), 'seqs loaded')

    # run search
    search(seqs, args) 

    print('finished running')
        
