#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH -o /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.out
#SBATCH -e /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

# Load modules
module load kallisto/0.43.1-foss-2016b

## Reference Files
REFS=/data/biorefs/reference_genomes/ensembl-release-98/danio_rerio/
IDX=/${REFS}/kallisto/Danio_rerio.GRCz11.cdna.primary_assembly.psen2S4Ter.idx

## Directories
PROJROOT=/data/biohub/20170327_Psen2S4Ter_RNASeq

## Setup for kallisto output
ALIGNDIR=${PROJROOT}/3_kallisto

## Now organise the input files
F1=$1
F2=${F1%_R1.fastq.gz}_R2.fastq.gz

## Organise the output files
OUTDIR=${ALIGNDIR}/$(basename ${F1%_R1.fastq.gz})
echo -e "Creating ${OUTDIR}"
mkdir -p ${OUTDIR}

echo -e "Currently aligning:\n\t${F1}\n\t${F2}"
echo -e "Output will be written to ${OUTDIR}"
kallisto quant \
	-b 50 \
	--rf-stranded \
	-t 1 \
	-i ${IDX} \
	-o ${OUTDIR} \
	${F1} ${F2} 
	