#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=10:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20200310_rRNADepletion/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=16

## Modules
module load Salmon/1.1.0-foss-2016b

## Directories
PROJROOT=/data/biohub/20170327_Psen2S4Ter_RNASeq/data
TRIMDATA=${PROJROOT}/1_trimmedData
SALMON=${PROJROOT}/5_salmon

## Setup for salmon
mkdir -p ${SALMON}/quant

## Run salmon
for R1 in ${TRIMDATA}/fastq/*R1.fastq.gz
do

  BNAME=$(basename ${R1%_R1.fastq.gz})
  R2=${R1%_R1.fastq.gz}_R2.fastq.gz
  echo -e "file 1 is ${R1}"
  echo -e "file 2 is ${R2}"

  salmon quant \
    -i ${PROJROOT}/5_salmon/salmon_index_drerio98 \
    -l A \
    -1 ${R1} \
    -2 ${R2} \
    --validateMappings \
    --threads ${CORES} \
    -o ${SALMON}/quant/${BNAME}

done