# from proposal outline (obsidian), components:
    # format seqs - convert raw data to some desired format for max shift
    # max shift (pairwise alignment) - input 2 formatted seqs, output motif string in same format
        # optionally, output matching indices in seq 1 and offset of seq 2
        # should be able to handle both full seqs and motif candidates
    # solve full batch - run max shift on all pairs, then refine all motifs
    # exact match (eval) - find all exact motif occurrences within set of seqs
# many possible formats; group together max shift and formatting in the same class to make implementation easy

from Bio import SeqIO
import numpy as np
from collections import Counter
import re
from scipy.fft import next_fast_len

# load from file with optional limit on # seqs, output list of strings with no metadata
def load(fn, lim=-1):
    seqs = []
    for i, rec in enumerate(SeqIO.parse(fn, 'fasta')):
        seqs.append(str(rec.seq))
        if lim > 0 and i+1 == lim: break
    return seqs

# convert seq strs to one hot arrays, get char list + fft size
def prepare(seqs):
    out_seqs = [np.array(list(s.encode()), dtype=np.uint8) for s in seqs]
    chars = np.unique(np.concat(out_seqs))
    chars = np.concat([np.array([ord('-')]), chars])

    # precompute optimal FFT size 
    maxlen = len(max(seqs, key=len))
    fft_sz = next_fast_len(2*maxlen-1)

    # convert all seqs to one-hots
    # TODO: figure out whether we can recover shared positions in pre-computed FFTs
    out_seqs = [(s[None,:] == chars[:,None]).astype(np.float32) for s in out_seqs]

    return out_seqs, chars, fft_sz

# convert one hot arrays back to strs
def to_str(s, chars, trim=True):
    # get indices
    char_idxs = s.argmax(axis=0)
    # trim
    if trim:
        nz_idxs = np.nonzero(s.sum(axis=0))[0]
        char_idxs = char_idxs[nz_idxs[0]:nz_idxs[-1]+1]
    # convert
    return ''.join(chr(chars[c]) for c in char_idxs)

# align two one hot seqs, return shared candidate motif
def align(s1, s2, fft_sz, return_shift=False):
    # cross correlate; output lags with padding are [0, 1, ..., n-2, n-1, 0, ..., 0, 1-m, 2-m, ..., -2, -1]
    s1_fft = np.fft.rfft(s1, fft_sz)
    s2_fft = np.fft.rfft(s2, fft_sz)
    xcorr = np.fft.irfft((s1_fft * s2_fft.conj()).sum(axis=0), fft_sz)

    # get best cross correlation
    off = np.argmax(xcorr).astype(int)
    if off >= s1.shape[1]:
        off -= fft_sz

    # get shifts and length
    t1 = max(0, off)
    t2 = max(0, -off)
    m_len = min(s1.shape[1] - t1, s2.shape[1] - t2)

    # get motif
    motif = np.logical_and(s1[:,t1:t1+m_len], s2[:,t2:t2+m_len]).astype(np.float32)

    out = (motif, np.sum(motif).astype(int))
    if return_shift: out += (t1,)
    return out

# given a query seq and a batch of targets, align one vs all, count column freqs, and output a w-mer set
def get_w_mers(query, targets, chars, fft_sz, min_weight=10, min_count=10):
    aligns = [align(query, t, fft_sz, True) for t in targets]

    col_freqs = np.zeros(query.shape[1], np.float32)
    for i in range(len(aligns)):
        t = aligns[i][-1]
        m_len = aligns[i][0].shape[1]
        col_freqs[t:t+m_len] += aligns[i][0].sum(0)

    # find cols exceeding threshold
    col_valid = col_freqs >= min_count
    if col_valid.sum() < min_weight: return None

    # make w-mers
    query[:, ~col_valid] = 0.0
    w = min_weight
    valid_idx = np.nonzero(col_valid)[0]

    seq_str = to_str(query, chars, trim=False)
    w_mers = [seq_str[valid_idx[i]:valid_idx[i+w]+1] for i in range(len(valid_idx) - w)]
    return w_mers


# given a query seq and a batch of targets, align one vs all repeatedly to refine into a conserved motif
def refine(query, targets, chars, fft_sz, min_weight=10, min_count=10):
    motif = query
    matches = set()
    # refine
    while True:
        if motif[1] < min_weight: break

        # search for best submotif and all exact matches of current motif simultaneously
        best_sub = (None,0)
        for i in range(len(targets)):
            if i in matches: continue
            submotif = align(motif[0], targets[i], fft_sz)

            if submotif[1] == motif[1]: matches.add(i)
            elif submotif[1] > best_sub[1]: best_sub = submotif
        
        # emit motif on meeting criteria
        if len(matches) >= min_count:
            return to_str(motif[0], chars)
        
        # update and continue refinement
        motif = best_sub

    return None

# find exact matches for given motif set
def find_hits(motifs, seqs):
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
    hits = [set() for _ in regexes]
    for i, seq in enumerate(seqs):
        for j, regex in enumerate(regexes):
            if regex.search(seq) != None:
                hits[j].add(i)
    
    return hits 

# given list of motifs and file, append to file
def write_motifs(fn, motifs):
    with open(fn, 'a') as f:
        for m in motifs:
            # get weight
            w = Counter(m[0])
            w['-'] = 0
            # set fields
            out = m[0] + ','
            out += str(w.total()) + ','
            out += str(len(m[0])) + ','
            out += str(len(m[1])) + ',' #','.join(str(e) for e in m[1])
            f.write(out + '\n')
