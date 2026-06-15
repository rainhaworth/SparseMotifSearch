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


# load from file with optional limit on # seqs, output list of strings with no metadata
def load(fn, lim=-1):
    seqs = []
    for rec in SeqIO.parse(fn):
        seqs.append(rec.seq)
    return seqs

# simple formatting + alignment
class ShiftAlign:
    # convert to numpy byte arrays
    def __init__(self, seqs):
        self.seqs = [np.array(list(s.encode()), dtype=np.uint8) for s in seqs]

    # align seqs at indices
    def align(self, i1, i2):
        assert i1 >= 0 and i2 >= 0 and i1 < len(self.seqs) and i2 < len(self.seqs)

        # explore all shifts; can only shift seqs "backward" with slicing
        s1 = self.seqs[i1]
        s2 = self.seqs[i2]
        best = (None,0) # motif, score

        # score seqs + update best shift
        def update_best(best, s1_, s2_):
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

            best = update_best(best, s1_, s2_)

        # shift s1 backward, effectively shifting s2 forward
        for i in range(len(s1)):
            s1_ = s1[i:i+len(s2)]
            s2_ = s2[:len(s1_)]

            best = update_best(best, s1_, s2_)
        
        return best