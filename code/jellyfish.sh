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

## Maximum hash size that fits in given memory was found using eg:
## jellyfish mem -m 10 --mem 8G

##############
## DE genes ##
##############

jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/dek10_counts.jf ${JELLYDIR}/deSeqs.fa \

jellyfish dump -c ${JELLYDIR}/dek10_counts.jf > ${JELLYDIR}/dek10_dumps.txt

rm ${JELLYDIR}/dek10_counts.jf

##################
## Not DE genes ##
##################

jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/conk10_counts.jf ${JELLYDIR}/conSeqs.fa 

jellyfish dump -c ${JELLYDIR}/conk10_counts.jf > ${JELLYDIR}/conk10_dumps.txt

rm ${JELLYDIR}/conk10_counts.jf

################
## Test genes ##
################

# ## n = 500

# jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/test500k10_counts.jf ${JELLYDIR}/test500.fa

# jellyfish dump -c ${JELLYDIR}/test500k10_counts.jf > ${JELLYDIR}/test500k10_dumps.txt

# ## n = 1000

# jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/test1000k10_counts.jf ${JELLYDIR}/test1000.fa

# jellyfish dump -c ${JELLYDIR}/test1000k10_counts.jf > ${JELLYDIR}/test1000k10_dumps.txt

# ## n = 2000

# jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/test2000k10_counts.jf ${JELLYDIR}/test2000.fa

# jellyfish dump -c ${JELLYDIR}/test2000k10_counts.jf > ${JELLYDIR}/test2000k10_dumps.txt

###################
## Control genes ##
###################

# ## n = 500

# jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/control500k10_counts.jf ${JELLYDIR}/control500.fa

# jellyfish dump -c ${JELLYDIR}/control500k10_counts.jf > ${JELLYDIR}/control500k10_dumps.txt

# ## n = 1000

# jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/control1000k10_counts.jf ${JELLYDIR}/control1000.fa

# jellyfish dump -c ${JELLYDIR}/control1000k10_counts.jf > ${JELLYDIR}/control1000k10_dumps.txt

# ## n = 2000

# jellyfish count -m 10 -s 2G -t ${CORES} -C -o ${JELLYDIR}/control2000k10_counts.jf ${JELLYDIR}/control2000.fa

# jellyfish dump -c ${JELLYDIR}/control2000k10_counts.jf > ${JELLYDIR}/control2000k10_dumps.txt
