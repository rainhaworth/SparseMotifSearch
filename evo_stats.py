# compute stats for all trials at single value of w
from glob import glob
import numpy as np
w = 100

files = glob('evo-trial*-w*-seq.txt')
trials_str = [open(f, 'r').readlines() for f in files if 'w' + str(w) in f]

# convert csvs to lists of ints, get size of each set
trials = [[[int(y) for y in x[:-1].split(',')] for x in l] for l in trials_str]
train_sz, val_sz, evo_sz = len(trials[0][0]), len(trials[0][1]), len(trials[0][2])
train_sz, val_sz, evo_sz

# get motif count, train hits, val hits, evo hits for each trial
trial_stats = np.zeros((len(trials), 4), int)

for i, trial in enumerate(trials):
    train, val, evo = trial
    a = max(train)
    tr = sum([x == a for x in train])
    va = sum([x == a for x in val])
    ev = sum([x == a for x in evo])

    trial_stats[i] = [a, tr, va, ev]

# convert to percentage, compute BA
trial_stats = trial_stats.astype(float)
trial_stats[:, 1:] /= [train_sz, val_sz, evo_sz]
trial_stats[:, 1:] *= 100
trial_stats = np.concat([trial_stats, ((trial_stats[:,2] + (100-trial_stats[:,3]))/2)[:,None]], axis=1)

# print stats for all trials + mean, std
print(trial_stats)
print('mean', trial_stats.mean(0))
print('stdev', trial_stats.std(0))