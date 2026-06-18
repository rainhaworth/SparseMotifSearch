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


# load from file with optional limit on # seqs, output list of strings with no metadata
def load(fn, lim=-1):
    seqs = []
    for i, rec in enumerate(SeqIO.parse(fn, 'fasta')):
        seqs.append(str(rec.seq))
        if lim > 0 and i == lim: break
    return seqs

# simple formatting + alignment
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
        idxs = np.argwhere(s).squeeze()
        s = s[idxs[0]:idxs[-1]+1]
        # replace 0 with -
        s = np.where(s == 0, ord('-'), s)
        # convert to string
        return ''.join([chr(x) for x in s])


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
def find_hits(motif, seqs):
    # make regex for motif
    motif_str = motif[0]

    re_str = ''
    wildcards = 0
    for i in range(len(motif_str)):
        if motif_str[i] == '-': wildcards += 1
        else:
            if wildcards != 0:
                re_str += '.{' + str(wildcards) + '}'
                wildcards = 0
            re_str += motif_str[i]

    print(motif)
    print(re_str)
    regex = re.compile(re_str)
    
    # search
    hits = set()
    for i, s in enumerate(seqs):
        if regex.search(s) != None:
            hits.add(i)
    
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