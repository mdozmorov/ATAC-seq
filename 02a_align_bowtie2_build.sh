#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov@vcu.edu
#PBS -N btbuild
#PBS -j oe
#PBS -q workq
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

source activate bowtie2

# Path to genome annotation files
## Human
DIRIN=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38ext/hg38ext.fa
## Mouse
# DIRIN=/home/mdozmorov/sequencing/data/ExtData/UCSC/mm10/mm10.fa

# Path to genome index
## Human
BOWTIE2_IDX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38ext/hg38ext
## Mouse
# BOWTIE2_IDX=/home/mdozmorov/sequencing/data/ExtData/UCSC/mm10/mm10

bowtie-build ${DIRIN} ${BOWTIE2_IDX}

conda deactivate
