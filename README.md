# seqdiff

This program implements Heng Li's formulas in
[A statistical framework for SNP calling, mutation discovery, association mapping and population
genetical parameter estimation from sequencing data]
(https://academic.oup.com/bioinformatics/article/27/21/2987/217423?login=true).

The goal of this program is to take as input a set of mapped reads and a reference diploid genome,
and estimate the fraction of nucleotides that mutated, where each nucleotide can take a continuous number
between 0 and 1, but likely 0, 0.5 or 1.

The important property of `seqdiff` is that it can be used for DNA methylation data, which allows Ts in
reads to map to either Cs or Ts in the reference genome. With this, one can use WGBS data to estimate
how much a dataset differs from the genome to which it was mapped.

Compile the program by running
```
make
make install
```
