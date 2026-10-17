from Bio import SeqIO

in_path = 'lnc_clust2.fasta'
out_path_1 = 'lnc_largest_clust.fasta'
out_path_2 = 'lnc_2nd_largest_clust.fasta'

# get clusters; mmseqs2 cluster rep descriptions are duplicated, so first occurrence has no seq
clusters = dict()
c_id = None
cluster = []
for rec in SeqIO.parse(in_path, 'fasta'):
    if len(rec.seq) == 0:
        if c_id is not None: clusters[c_id] = cluster
        c_id = rec.description
        cluster = []
    else: cluster.append(rec.seq)

# sort by cluster length
sorted_keys = sorted(clusters, key=lambda x: len(clusters[x]), reverse=True)

# output top 2 clusters
for out_path, i in zip([out_path_1,out_path_2], [0,1]):
    with open(out_path, 'w') as f:
        for i, s in enumerate(clusters[sorted_keys[i]]):
            f.write('>' + str(i) + '\n' + str(s) + '\n')