# SparseMotifSearch
Methods for detecting sparse conserved motifs in groups of sequences

## Datasets

TODO

## Experiments

### Protein motif finding

```bash
python find.py -i uniprot_LexA.fasta -o LexA.csv --batch
python find.py -i uniprot_pgr.fasta -o pgr.csv -c 0.8 --batch -b 30
python find.py -i uniprot_Cu_Zn_SOD.fasta -o sod.csv -c 0.45 --batch
python find.py -i uniprot_TET.fasta -o tet.csv --batch
```

### LncRNA motif finding

```bash
python find.py -i human_lncRNA_filtered.fasta -o lnc_evl.csv -c 0.45 --batch -m 30
python find.py -i lncRNA_2nd_largest_clust.fasta -o lnc_2nd.csv -c 0.4 --batch
```

### 16S mutation detection

```bash
python find.py -i SILVA_SSU_Ref_NR_99_filtered.fasta -o 16s.csv -c 0.9 --batch -b 30 -m 20
```

### Evo-generated PhiX174 detection

TODO: convert notebook to script

search: always `-c 1.0`, weight is 20, 50, or 100, `-a column -r seq -d sparse`
