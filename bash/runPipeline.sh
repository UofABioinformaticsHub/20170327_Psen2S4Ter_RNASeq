#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time=36:00:00
#SBATCH --mem=32GB
#SBATCH -o /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.out
#SBATCH -e /data/biohub/20170327_Psen2S4Ter_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

## Clean run of the PSEN2 data.

## Cores
CORES=16

## Modules
module load FastQC/0.11.7
module load STAR/2.7.0d-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load cutadapt/1.14-foss-2016b-Python-2.7.13
module load Subread/1.5.2-foss-2016b

## Function for checking directories
checkAndMake () {
  echo "Checking if $1 exists"
  if [[ ! -d $1 ]]
  then 
    echo "Creating $1"
    mkdir -p $1
  fi
  
  if [[ -d $1 ]]
  then
    echo "Found $1"
  else
    echo "$1 could not be created or found"
    exit 1
  fi  
  
}

## Directories
PROJROOT=/data/biohub/20170327_Psen2S4Ter_RNASeq
REFS=/data/biorefs/reference_genomes/ensembl-release-98/danio_rerio/
if [[ ! -d ${REFS} ]]
then
  echo "Couldn't find ${REFS}"
  exit 1
if
GTF=${REFS}/Danio_rerio.GRCz11.98.chr.gtf.gz
if [[ ! -f ${GTF} ]]
then
  echo "Couldn't find ${GTF}"
  exit 1
if

# Raw Data
RAWDIR=${PROJROOT}/0_rawData
checkAndMake ${RAWDIR}
checkAndMake ${RAWDIR}/FastQC

## Trimmed 
TRIMDIR=${PROJROOT}/1_trimmedData
checkAndMake ${TRIMDIR}/fastq
checkAndMake ${TRIMDIR}/FastQC
checkAndMake ${TRIMDIR}/log

## Aligned
ALIGNDIR=${PROJROOT}/2_alignedData
checkAndMake ${ALIGNDIR}
checkAndMake ${ALIGNDIR}/bam
checkAndMake ${ALIGNDIR}/FastQC
checkAndMake ${ALIGNDIR}/log

echo "All directories checked and created"

##----------------------------------------------------------------------------##
##                              Initial FastQC                                ##
##----------------------------------------------------------------------------##

fastqc -t ${CORES} -o ${RAWDATA}/FastQC --noextract ${RAWDATA}/fastq/*fastq.gz

##----------------------------------------------------------------------------##
##                              Trimming                                      ##
##----------------------------------------------------------------------------##

for R1 in ${RAWDIR}/fastq/*R1.fastq.gz
  do
    R2=${R1%_R1.fastq.gz}_R2.fastq.gz
    echo -e "The R1 file should be ${R1}"
    echo -e "The R2 file should be ${R2}"

    ## Create output filenames
    out1=${TRIMDIR}/fastq/$(basename $R1)
    out2=${TRIMDIR}/fastq/$(basename $R2)
    BNAME=${TRIMDIR}/fastq/$(basename ${R1%_1.fq.gz})
    echo -e "Output file 1 will be ${out1}"
    echo -e "Output file 2 will be ${out2}"
    echo -e "Trimming:\t${BNAME}"
    
    LOG=${TRIMDIR}/log/${BNAME}.info
    
    cutadapt \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -o ${out1} \
      -p ${out2} \
      -m 35 \
      --trim-n \
      --max-n=1 \
      --nextseq-trim=30 \
      --info-file=${LOG} \
      ${R1} \
      ${R2} 
    
  done

fastqc -t ${CORES} -o ${TRIMDATA}/FastQC --noextract ${TRIMDATA}/fastq/*fastq.gz


##----------------------------------------------------------------------------##
##                                STAR Alignment                              ##                
##----------------------------------------------------------------------------##

## Aligning, filtering and sorting
for R1 in ${TRIMDATA}/fastq/*R1.fastq.gz
 do

 BNAME=$(basename ${R1%_R1.fastq.gz})
 R2=${R1%_R1.fastq.gz}_R2.fastq.gz
 echo -e "STAR will align:\t${R1}"
 echo -e "STAR will also align:\t${R2}"

  STAR \
    --runThreadN ${CORES} \
    --genomeDir ${REFS}/star \
    --readFilesIn ${R1} ${R2} \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix ${ALIGNDATA}/bam/${BNAME} \
    --outSAMtype BAM SortedByCoordinate

 done

## Move the log files into their own folder
mv ${ALIGNDATA}/bam/*out ${ALIGNDATA}/log
mv ${ALIGNDATA}/bam/*tab ${ALIGNDATA}/log

## Fastqc and indexing
for BAM in ${ALIGNDATA}/bam/*.bam
do
  fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATA}/FastQC --noextract ${BAM}
  samtools index ${BAM}
done


##----------------------------------------------------------------------------##
##                                featureCounts                               ##
##----------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
sampleList=`find ${ALIGNDATA}/bam -name "*out.bam" | tr '\n' ' '`

## Extract gtf for featureCounts
zcat ${GTF} > temp.gtf

## Running featureCounts on the sorted bam files
featureCounts -Q 10 \
  -s 2 \
  -T ${CORES} \
  -p \
  --fracOverlap 1 \
  -a temp.gtf \
  -o ${ALIGNDATA}/featureCounts/counts.out ${sampleList}

## Remove the temporary gtf
rm temp.gtf

 ## Storing the output in a single file
cut -f1,7- ${ALIGNDATA}/featureCounts/counts.out | \
  sed 1d > ${ALIGNDATA}/featureCounts/genes.out

##----------------------------------------------------------------------------##
##                                  kallisto                                  ##
##----------------------------------------------------------------------------##

## Aligning, filtering and sorting
for R1 in ${TRIMDATA}/fastq/*R1.fastq.gz
  do
    ${PROJROOT}/bash/singleKallisto.sh ${R1}
  done
  