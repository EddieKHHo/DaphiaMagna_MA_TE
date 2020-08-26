# *simMutAccumTE*

## Description
This program was used by ABC (Ho et al. 2020) to estimate the false positive and false negative rates of detecting tranposable element mutations when using TEFLoN (Adrion et al. 2017).
This program is a modified version of 'simpoolTE' from Adrion et al. (2018).
SimMutAccumTE simulates a mutation accumulation (MA) experiment by inserting and/or deleting a number of transposable elements (TEs) into a diploid genome of a focal MA line while leaving the ancestral (ANC) genomes intact. After mutations are simulated in the given genome, pIRS (Hu et al. 2012) is used to generate paired-end reads for all lines.

## Requirements
- Python v3.7
- pIRS (https://github.com/galaxy001/pirs)
- seqtk (https://github.com/lh3/seqtk)

## Usage

```
usage: simMutAccumTE.py <required> [optional] 
  -wd WD          <Full path to working directory (include prefix for new folders)>
  -pirs PIRSPATH  <Full path to pIRS (including /pirs at the end)>
  -c REFCHROM     <Chromosome to simulate in fasta format (must contain only one sequence)>
  -b BED          <TE annotation bed file>
  -nmut NMUT      <Number of TE mutations to simulate in MA line>
  -tmut TYPEMUT   <Type of mutation to simulate [1, 2, 3, 4]>
  -nhet NHET      <Number of shared heterozygous TE sites to simulate>
  -ncl NCL        <Number of non-focal MA lines to simulate (clones of Ancestor)>
  -mnlen MNLEN    [Minimum length of TEs to insert and delete]
  -mxlen MXLEN    [Maximum length of TEs to insert and delete]
  -r RANDSEED     [Seed for random number generator]
  -snp SNPRATE    [Rate of heterozygosity for haploid genome]
  -x COV          [Coverage for pIRS to simulate]
  -rlen RLEN      [Read length for pIRS to simulate]
  -insz INSIZE    [Insert size for pIRS to simulate]

```

SimMutAccumTE inserts or delete TEs from the TE annotation file (<em>b</em>) onto the genome fasta file (<em>c</em>).

First, pIRS (Hu et al. 2012) is used to generate a diploid ancestor (ANC) genome by duplicated the given genome fasta file and adding SNPs at a specified rate (<em>snp</em>); the fasta files of ANC is stored in /pirsAnc. pIRS requires that fasta file does not contain N's. 

Subsequently <em>nmut</em> TE mutations of a specified type (<em>tmut</em>) are simulated. Type 1 are novel TE insertions (0->1 gain), simulated by insertion a heterozygous TE (i.e., on one homolog of the genome) into the MA line. Type 2 are TE deletion from an ancestrally homozygous TE site (2->1 loss), simulated by inserting a homozygous TE (i.e. insertion of both homologs) on ANC and a heterozygous TE on the MA line. Type 3 are TE insertion onto an ancestrally heterozygous TE site (1->2 gain), simualted by inserting a heterozygous TE on the ANC and then a homozygous TE on the MA line. Type 4 are TE deletions from an ancestrally heterozygous TE site (1->0 loss), simualted by inserting a heterozygous TE on the ANC only. In addition a <em>nhet</em> number of shared TE heterozygous sites can be added onto both ANC and MA lines. 

Each TE insertion was chosen randomly from the annotation file and has a minimum and maximum length given by <em>mnlen</em> and <em>mxlen</em>, respectively. Insertions of TEs are separated by at least one read length (<em>rlen</em>) and flanked by a target site duplication with a mean length of 5 bp drawn from a Poisson distribution.



