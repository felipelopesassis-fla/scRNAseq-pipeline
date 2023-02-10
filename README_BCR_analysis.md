## Scripts for BCR analyses in Assis, Hoehn et al. (submitted)
Kenneth B. Hoehn
kenneth.hoehn@yale.edu
1/28/23

### Processing 10X VDJ sequencing data

To align BCR sequences for all sequencing runs to the IMGT reference database, run
```
bash processVDJ.sh

```

### Analysis scripts

To combine BCR information with GEX information, identify clones, contruct clonal germlines, and identify convergent sequence clusters as well as cells with public heavy chains, run

```
Rscript identifyClones.R

```

To analyze clonal diversity, overlap, and cell type composition, run:

```
Rscript cloneAnalysis.R

```

To perform phylogenetic analysis and generate convergent cluster alignment figure, run

```
Rscript treeAnalysis.R

```
