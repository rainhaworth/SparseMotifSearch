# protein motif "grammar" visualization
import re
import numpy as np
from Bio import SeqIO
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import random

# set to a value in 0-3 to generate different figures
fig_i = 2

fig_paths = [
    ('uniprotkb_LexA.fasta', 'LEXA_ECOLI', 'lexa.csv', 'P0A7C2.json'),
    ('uniprotkb_pgr.fasta', 'PRGR_HUMAN', 'pgr.csv', 'P06401.json'),
    ('uniprotkb_family_RRM_TET_family.fasta', 'FUS_HUMAN', 'tet.csv', 'P35637.json'),
    ('uniprotkb_family_Cu_Zn_superoxide_dismu.fasta', 'SODC_HUMAN', 'sod.csv', 'P00441.json')
]

# specify domain height, fig width, fig height, tick step
fig_mds = [
    (30, 5.4, 3, 50),
    (30, 9.6, 3, 200),
    (20, 9.6, 4, 100),
    (30, 3.8, 3, 50)
]

seq_path, seq_id, mot_path, md_path = fig_paths[fig_i]
md_whitelist = ['Domain', 'Zinc finger', 'Region', 'DNA binding', 'Active site', 'Site', 'Binding site']

# retrieve features matching whitelist, group by description, get all pos ranges, sort
md = json.load(open(md_path, 'r'))
feats = [x for x in md['features'] if x['type'] in md_whitelist]
feat_locs = defaultdict(list)
for feat in feats:
    loc = feat['location']
    s, e = loc['start']['value']-1, loc['end']['value']-1
    if feat['type'] == 'Binding site': fid = feat['ligand']['name']
    else: fid = feat['description']
    if 'AF' in fid or 'Mediates' in fid or 'Breakpoint' in fid or 'Modulating' in fid: continue
    if 'cleavage' in fid.lower(): fid = 'Cleavage site'
    feat_locs[fid].append(list(range(s, e+1)))
feat_locs = {a: b for a, b in sorted(feat_locs.items(), key = lambda x: x[1][0])}
assert len(feat_locs) > 0

# get ref seq + find n
n = 0
seq = None
for rec in SeqIO.parse(seq_path, 'fasta'):
    n+= 1
    if seq_id in rec.description and seq is None:
        seq = str(rec.seq)
assert seq is not None

# populate domain map
random.seed(42)
height = fig_mds[fig_i][0]
dommap = np.ones((len(feat_locs)*height, len(seq)))
for i, dom in enumerate(feat_locs.values()):
    c = random.random() * 0.9
    for j in dom:
        dommap[i*height+1:(i+1)*height-1, j] = c

# get motifs
motifs = []
with open(mot_path, 'r') as f:
    for line in f.readlines():
        if 'motif' in line: continue
        fields = line.split(',')
        motifs.append((fields[0], int(fields[3])))
assert len(motifs) > 0

# find hits
mot_res = [re.compile(mot[0].replace('-','.')) for mot in motifs]
hits = [m_r.search(seq) for m_r in mot_res]

# retrieve start pos, drop missing motifs, and sort
hits = [(x.start(), y) for x, y in zip(hits, motifs) if x is not None]
hits = sorted(hits, key=lambda x: x[0])
assert len(hits) > 0

# populate seq map
seqmap = np.ones((len(hits), len(seq)))
for i, (pos, mot) in enumerate(hits):
    cl = 1.0 - (mot[1] / n)
    for j, c in enumerate(mot[0]):
        if c != '-': seqmap[i, pos + j] = cl
fontsz = 15
fig, ax = plt.subplots(2, 1, 
    figsize=(fig_mds[fig_i][1], fig_mds[fig_i][2]),
    height_ratios=[dommap.shape[0]/seqmap.shape[0],1.0],
    sharex=True, layout='constrained'
)
fig.get_layout_engine().set(hspace=0) 
ax[0].imshow(dommap, cmap='gist_stern', vmin=0.0, vmax=1.0)
ax[0].set_yticks([x*height+height//2 for x in range(len(feat_locs))], [x.split('(')[0] for x in feat_locs.keys()])
ax[0].tick_params(left=False, bottom=False, labelsize=fontsz)
ax[1].imshow(seqmap, cmap='gray', vmin=0.0, vmax=1.0)
ax[1].set_yticks([])
ax[1].set_xticks(range(0,seqmap.shape[1],fig_mds[fig_i][3]))
ax[1].tick_params(bottom=False, labelsize=fontsz)
for i in [0,1]:
    ax[i].set_frame_on(False)

plt.xlabel('Position in ' + seq_id, fontsize=fontsz)
plt.grid()
plt.savefig('o_' + seq_id + '.png', dpi=300)