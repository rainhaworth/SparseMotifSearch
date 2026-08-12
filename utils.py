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
from tqdm import tqdm
from collections import Counter
import re
#import torch
from scipy.fft import next_fast_len

# load from file with optional limit on # seqs, output list of strings with no metadata
def load(fn, lim=-1):
    seqs = []
    for i, rec in enumerate(SeqIO.parse(fn, 'fasta')):
        seqs.append(str(rec.seq))
        if lim > 0 and i+1 == lim: break
    return seqs

# simple formatting + alignment
# also a base/parent class; child classes should reimplement __init__, align_seqs, and to_str
class ShiftAlign:
    # convert to numpy byte arrays
    def __init__(self, seqs):
        self.seqs = [np.array(list(s.encode()), dtype=np.uint8) for s in seqs]
    
    def __len__(self):
        return len(self.seqs)
    
    def __getitem__(self, idx):
        return self.seqs[idx]

    # align any arbitrary seqs
    def align_seqs(self, s1, s2):
        # explore all shifts; can only shift seqs "backward" with slicing
        best = (None,0) # motif, score

        # score seqs + update best shift
        def update_best(best, s1_, s2_, flag):
            match = (s1_ == s2_)
            score = np.sum(match)
            if score > best[1]:
                return (np.where(match, s1_, 0), score)
            return best

        # shift s2 backward
        for i in range(len(s2)):
            # make sure lengths match
            s2_ = s2[i:i+len(s1)]
            s1_ = s1[:len(s2_)]

            best = update_best(best, s1_, s2_, 's2')

        # shift s1 backward, effectively shifting s2 forward
        for i in range(len(s1)):
            s1_ = s1[i:i+len(s2)]
            s2_ = s2[:len(s1_)]

            best = update_best(best, s1_, s2_, 's1')
        
        return best
    
    # align by indices
    def align_at_idx(self, i1, i2):
        assert i1 >= 0 and i2 >= 0 and i1 < len(self.seqs) and i2 < len(self.seqs)
        s1 = self.seqs[i1]
        s2 = self.seqs[i2]
        return self.align_seqs(s1, s2)
    
    # convert to str (from byte array)
    def to_str(self, s):
        # trim
        idxs = np.nonzero(s)[0]
        s = s[idxs[0]:idxs[-1]+1]
        # replace 0 with -
        s = np.where(s == 0, ord('-'), s)
        # convert to string
        return ''.join([chr(x) for x in s])

# shift align on GPU: convert to one-hot for fast matmul comparisons
# TODO: port to tensorflow; efficient pytorch implementation is not feasible, too many missing features
"""
class ShiftAlignGPU(ShiftAlign):
    def __init__(self, seqs):
        # get device; no multi-GPU support for now
        if torch.cuda.is_available():
            dev = 'cuda:0'
        else:
            dev = 'cpu'
        self.device = torch.device(dev)
        # token mapping
        charset = set().union(*[set(s) for s in seqs]) # list of sets of unique chars -> expand -> union
        c2t = {k: v for k,v in zip(charset, range(1,len(charset)+1))}
        c2t['-'] = 0
        # tokenize -> one hot; torch doesn't have ragged tensors so store as a list of tensors of dim (L, |Sigma|+1)
        seq_toks = [[c2t[c] for c in s] for s in seqs]
        self.seqs = [torch.nn.functional.one_hot(torch.LongTensor(st)).half().to(self.device) for st in seq_toks]
        # decoder for to_str
        self.t2c = {k:v for v,k in c2t.items()}
    
    def align_seqs(self, s1, s2):
        # big matmul, shape (L_1, L_2), 1 for matches and 0 for mismatches
        algn_mat = s1 @ s2.T

        # we can now find the best shift as the diagonal with the largest sum. the hard part is computing the diagonal sums efficiently
        # torch can't retrieve all diagonals at once, so use list comprehension
        # some other ideas:
            # https://stackoverflow.com/a/35074207 -> replace np.indices with meshgrid
            # modify https://stackoverflow.com/a/74120776 to work with (n,m) matrices
        n, m = algn_mat.shape
        offsets = torch.arange(-n+1,m)
        diags = [torch.diagonal(algn_mat, offset=o) for o in offsets]
        diag_sums = torch.Tensor([t.sum() for t in diags]).to(self.device)

        # get best diagonal
        best = torch.argmax(diag_sums)
        best_offset = offsets[best]
        best_diag = diags[best]

        # retrieve motif candidate
        if best_offset >= 0:
            seq = s2[best_offset:best_offset+len(best_diag),:]
        else:
            seq = s1[-best_offset:-best_offset+len(best_diag),:]
        # convert mismatches to token 0 (-)
        zero_idx = torch.nonzero(1.0 - best_diag[:seq.shape[0]])
        tok_0 = torch.zeros(seq.shape[1], dtype=torch.half, device=self.device)
        tok_0[0] = 1.0
        seq[zero_idx] = tok_0
        return (seq, best_diag.sum())

    def to_str(self, s : torch.Tensor):
        s_np = s.numpy(force=True)
        s_toks = np.nonzero(s_np)[1]
        return ''.join(self.t2c[t] for t in s_toks)
"""

class ShiftAlignFFT(ShiftAlign):
    def __init__(self, seqs):
        # convert seqs to array, get unique chars
        seqs = [np.array(list(s.encode()), dtype=np.uint8) for s in seqs]
        self.chars = np.unique(np.concat(seqs))
        self.chars = np.concat([np.array([ord('-')]), self.chars])

        # precompute optimal FFT size 
        maxlen = len(max(seqs, key=len))
        self.fft_sz = next_fast_len(2*maxlen-1)

        # convert all seqs to one-hots; could pre-compute FFT but this way we can guarantee exact motif retrieval
        self.seqs = [(s[None,:] == self.chars[:,None]).astype(np.float32) for s in seqs]
    
    def align_seqs(self, s1, s2):
        # cross correlate; output lags with padding are [0, 1, ..., n-1, 0, ..., 0, 1-m, 2-m, ..., -1]
        s1_fft = np.fft.rfft(s1, self.fft_sz)
        s2_fft = np.fft.rfft(s2, self.fft_sz)
        xcorr = np.fft.irfft((s1_fft * s2_fft.conj()).sum(axis=0), self.fft_sz)

        # get best cross correlation
        off = np.argmax(xcorr).astype(int)
        if off >= s1.shape[1]:
            off -= self.fft_sz

        # get shift and length
        t1 = max(0, off)
        t2 = max(0, -off)
        m_len = min(s1.shape[1] - t1, s2.shape[1] - t2)

        # get motif
        motif = np.logical_and(s1[:,t1:t1+m_len],s2[:,t2:t2+m_len]).astype(np.float32)
        return (motif, np.sum(motif).astype(int))

    def to_str(self, s):
        # get indices
        char_idxs = s.argmax(axis=0)
        # trim
        nz_idxs = np.nonzero(s.sum(axis=0))[0]
        char_idxs = char_idxs[nz_idxs[0]:nz_idxs[-1]+1]
        # convert
        return ''.join(chr(self.chars[c]) for c in char_idxs)


# given an initialized aligner containing a batch of formatted seqs, a min size, and a min rate of occurrence, find all motifs   
def solve_batch(batch : ShiftAlign, min_weight=10, min_count=10, dedup=True):
    assert min_count <= len(batch)
    motifs = []
    
    # O(N^3L^2 * min(N,L)) search
    t = tqdm(total=(len(batch)**2 - len(batch)) // 2)
    for i in range(len(batch)-1):
        for j in range(i+1,len(batch)):
            t.update()
            # get initial candidate and set of all seqs containing an exact amtch
            motif = batch.align_at_idx(i,j)
            matches = set([i,j])
            
            # refine
            while True:
                if motif[1] < min_weight: break

                # search for best submotif and all exact matches of current motif simultaneously
                best_sub = (None,0)
                for k in range(len(batch)):
                    if k in matches: continue
                    submotif = batch.align_seqs(motif[0], batch[k])

                    if submotif[1] == motif[1]: matches.add(k)
                    elif submotif[1] > best_sub[1]: best_sub = submotif
                
                # emit motif on meeting criteria
                if len(matches) >= min_count:
                    motifs.append((batch.to_str(motif[0]), matches))
                    break
                
                # update and continue refinement
                motif = best_sub
    t.close()

    # deduplicate
    # TODO: if it's possible for the same motif to have different matches, aggregate all matches. otherwise this is fine
    if dedup:
        motifs = list(
            set(
                (m[0], frozenset(m[1])) for m in motifs
                )
            )

    return motifs

# find exact matches for given motif set
def find_hits(motifs, seqs):
    # make regex for motif
    motif_strs = [motif[0] for motif in motifs]

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
