#!/bin/sh

#### This is an example for SNP filtering using R script SNPfiltering.R
#### run `Rscript SNPfiltering.R --help` for more information
#### Required R packages: argparse, updog, XNomial, doParallel, foreach
#### This step is not mandatory. It might filter out a lot of SNPs especially when data coverage is low.

## sbatch --nodes=1 --ntasks=16 --mem=32000 --time=1-00:00:00 --partition=snowy,physical --job-name=z2 --export=ploidy=2 SNPfiltering.sh
## sbatch --nodes=1 --ntasks=16 --mem=32000 --time=1-00:00:00 --partition=snowy,physical --job-name=z4 --export=ploidy=4 SNPfiltering.sh

module load r/3.6.2   # load R module. Edit this as needed before run this script

## the input (-i) and output (-o) file should be unzipped raw text VCF file
Rscript SNPfiltering.R -i data/out2.vcf -o data/out2_filtered.vcf --parents P1:P2 --ploidy 2 -maxThreads 16

