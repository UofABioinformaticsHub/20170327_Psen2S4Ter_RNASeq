#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH -o /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.out
#SBATCH -e /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=16

## Modules
module load Jellyfish/2.2.6-foss-2016b

## Directories
PROJROOT=/data/biohub/20170327_Psen2S4Ter_RNASeq/data
JELLYFQ=${PROJROOT}/7_jellyfishFq
RAWDIR=/data/biohub/170327MichaelLardelli_Nextseq/0_rawData

##----------------------------------------------------------------------------##
##                                Jellyfish                                   ##
##----------------------------------------------------------------------------##

## Count kmers in raw fastq files.

## Maximum hash size that fits in given memory was found using eg:
## jellyfish mem -m 10 --mem 8G

for R1 in ${RAWDIR}/fastq/*R1.fastq.gz
  do
    
    # R1=${RAWDIR}/fastq/10_Ps2Ex3M1_WT_6month_07_07_2016_F3_86_Fem_R1.fastq.gz
    echo -e "Counting kmers in ${R1}"

    ## Create output name
    out=${JELLYFQ}/$(basename ${R1%_R1.fastq.gz})
    echo -e "Output file will be ${out}_dumps.txt"

    jellyfish count -m 9 -s 8G -C -t ${CORES} -o ${out}_counts.jf <(zcat ${R1})
    jellyfish dump -c ${out}_counts.jf > ${out}_dumps.txt
    rm ${out}_counts.jf

  done
