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
#SBATCH --mail-user=baerlachlan@gmail.com

## Cores
CORES=16

## Modules
module load FastQC/0.11.7
module load BWA/0.7.15-foss-2017a
module load SAMtools/1.9-foss-2016b

## Reference files
RRNA=/data/biorefs/rRNA/danio_rerio/bwa/danRer11

## Directories
PROJROOT=/data/biohub/20170327_Psen2S4Ter_RNASeq/data
TRIMDATA=${PROJROOT}/1_trimmedData

## Setup for BWA alignment
ALIGNDATABWA=${PROJROOT}/4_bwa
mkdir -p ${ALIGNDATABWA}/bam
mkdir -p ${ALIGNDATABWA}/fastq
mkdir -p ${ALIGNDATABWA}/log
mkdir -p ${ALIGNDATABWA}/FastQC

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference rRNA
##--------------------------------------------------------------------------------------------##

# ## Aligning and sorting
# for R1 in ${TRIMDATA}/fastq/*R1.fastq.gz                                                       
# do

#   out=${ALIGNDATABWA}/bam/$(basename ${R1%.fastq.gz})
#   echo -e "Output file will be ${out}"

#   ## Align and return reads as .bam
#   bwa mem -t ${CORES} ${RRNA} ${R1} \
#   | samtools view -u -h \
#   | samtools sort -o ${out}.sorted.bam

# done

# ## Run FastQC
# fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATABWA}/FastQC --noextract ${ALIGNDATABWA}/bam/*.bam

# ## Indexing, flagstat, and conversion of unmapped reads to fastq for further alignment
# for BAM in ${ALIGNDATABWA}/bam/*.bam                                                         
# do

#   outbam=${ALIGNDATABWA}/log/$(basename ${BAM%.sorted.bam})
#   outfastq=${ALIGNDATABWA}/fastq/$(basename ${BAM%.sorted.bam})
#   echo -e "Working on ${outbam}"

#   samtools index ${BAM}
#   samtools flagstat ${BAM} > ${outbam}.flagstat

#   echo -e "Now working on ${outfastq}"

#   ## Output only unmapped reads as fastq using -f 4
#   samtools fastq -f 4 -c 6 --threads ${CORES} ${BAM} > ${outfastq}.fastq

# done

gzip ${ALIGNDATABWA}/fastq/*.fastq