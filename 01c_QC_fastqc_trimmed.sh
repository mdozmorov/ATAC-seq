#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M mdozmorov$vcu.edu
#PBS -m ae
#PBS -N fastqctrim
#PBS -j oe
#PBS -q workq
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

source activate trimgalore

# Number of CPUs, should match the ppn definition
THREADS=2

# Input folder
DIRIN=01_trimmed
# Output folder
DIROUT=01_trimmed_fastqc
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

# Only necessary if TrimGalore didn't output FastQC files
# for file in `find ${DIRIN} -type f -name "*.fq.gz" | sort`; do 
# 	fastqc -t ${THREADS} -o ${DIROUT} --noextract $file; 
# done

# Summarize quality metrics. I have to run it separately
# Output directory is a part of the output filename
FOUT=${DIROUT}.html
multiqc --filename ${FOUT} --outdir ${DIROUT} ${DIRIN}

conda deactivate
