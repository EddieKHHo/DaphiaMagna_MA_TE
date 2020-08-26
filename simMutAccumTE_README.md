# Description and usage of simMutAccumTE.py

## Description
This program was used by ABC (Ho et al. 2020) to estimate the false positive and false negative rates of detecting tranposable element mutations when using TEFLoN (Adrion et al. 2017).
This program is a modified version of 'simpoolTE' from Adrion et al. (2018).
SimMutAccumTE simulates a mutation accumulation (MA) experiment by inserting and/or deleting a number of transposable elements (TEs) into a diploid genome of a focal MA line while leaving the ancestral (ANC) and non-focal MA genomes intact. After mutations are simulated in the given genome, pIRS (Hu et al. 2012) is used to generate paired-end reads for all lines.

## Requirements
- Python v3.7
- pIRS (https://github.com/galaxy001/pirs)
- seqtk (https://github.com/lh3/seqtk)

## Usage
