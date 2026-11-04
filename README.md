# SparseMotifSearch
Methods for detecting sparse conserved motifs in groups of sequences

## Datasets

TODO

## Experiments

### Protein motif finding

Run motif finding for each protein dataset with the following commands. These can be run without modification if the date at the end of the filename is removed.

```bash
python find.py -i proteins/uniprotkb_LexA.fasta -o proteins/lexa.csv --batch
python find.py -i proteins/uniprotkb_pgr.fasta -o proteins/pgr.csv -c 0.8 --batch -b 30
python find.py -i proteins/uniprotkb_family_Cu_Zn_superoxide_dismu.fasta -o proteins/sod.csv -c 0.45 --batch
python find.py -i proteins/uniprotkb_family_RRM_TET_family.fasta -o proteins/tet.csv --batch
```

We include the code used to generate our figures in `proteins/figures_proteins.py`.

### 16S mutation detection

First, filter the SSU Ref NR 99 database then find motifs with the following commands.

```bash
python filter_silva_seqs.py
python find.py -i SILVA_SSU_Ref_NR_99_tax_silva_filtered.fasta -o 16s.csv -c 0.9 --batch -b 30 -m 20
```

We then use two dedicated scripts implementing GPU-accelerated motif retrieval for training and testing. It is possible to run these scripts without a GPU by manually setting `engine=auto` in the `get_hits` function. The final output of the testing script is the balanced accuracy and confusion matrix. Note that the training script only needs to run once, while the testing script must run three times to cover all mutation rates.

```bash
python 16s_train.py
python 16s_test.py --train SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta --test SILVA_138.2_SSURef_tax_silva.fasta -m 16s.csv -c clf-r0.1.pickle -r 0.1
python 16s_test.py --train SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta --test SILVA_138.2_SSURef_tax_silva.fasta -m 16s.csv -c clf-r0.05.pickle -r 0.05
python 16s_test.py --train SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta --test SILVA_138.2_SSURef_tax_silva.fasta -m 16s.csv -c clf-r0.01.pickle -r 0.01
```

### Evo-generated PhiX174 detection

For this experiment, we use another dedicated script: `evo_crossval.py`. The commands below will run all 5 trials for each value of `w`. Each run generates one `.txt` file containing results and one `.csv` file containing motifs for each trial. The output `.txt` files contain the number of motifs found in each sequence in the training, natural test, and Evo test set on the first, second, and third lines, respectively. We label each sequence as predicted positive when the number of motifs contained matches the number of motifs found in the corresponding label. All other sequences are predicted negative. We use these labels to compute the TPR, FPR, and balanced accuracy reported in the paper.

```bash
python evo_crossval.py -a column -c 1.0 -w 20 -d sparse -r seq
python evo_crossval.py -a column -c 1.0 -w 50 -d sparse -r seq
python evo_crossval.py -a column -c 1.0 -w 100 -d sparse -r seq
```
