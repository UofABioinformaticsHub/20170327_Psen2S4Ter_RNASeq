#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=10:00:00
#SBATCH --mem=8GB
#SBATCH -o /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.out
#SBATCH -e /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=4

## Modules
module load Jellyfish/2.2.6-foss-2016b

## Directories
PROJROOT=/data/biohub/20170327_Psen2S4Ter_RNASeq/data
JELLYDIR=${PROJROOT}/6_jellyfish

##----------------------------------------------------------------------------##
##                                Jellyfish                                   ##
##----------------------------------------------------------------------------##

## Count kmers in fasta file of coding sequences from genes classified as
## DE or not DE from differential expression testing using sample ribosomal
## RNA as predictor variable. Fasta file was generated in R. 

## Maximum hash size that fits in given memory was found using eg:
## jellyfish mem -m 10 --mem 8G

##############
## DE genes ##
##############

jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/dek10_counts.jf ${JELLYDIR}/deSeqs.fa

jellyfish dump -c ${JELLYDIR}/dek10_counts.jf > ${JELLYDIR}/dek10_dumps.txt

rm ${JELLYDIR}/dek10_counts.jf

##################
## Not DE genes ##
##################

jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/conk10_counts.jf ${JELLYDIR}/conSeqs.fa 

jellyfish dump -c ${JELLYDIR}/conk10_counts.jf > ${JELLYDIR}/conk10_dumps.txt

rm ${JELLYDIR}/conk10_counts.jf
