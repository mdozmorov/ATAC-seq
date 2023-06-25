#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov$vcu.edu
#PBS -m ae
#PBS -N fastqcraw
#PBS -j oe
#PBS -q workq
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

source activate trimgalore

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=00_raw
# Output folder
DIROUT=00_raw_fastqc
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find ${DIRIN} -type f -name "*.fastq.gz" | sort`; do 
	fastqc -t ${THREADS} -o ${DIROUT} --noextract $file; 
done

# Summarize quality metrics. I have to run it separately
# Output directory is a part of the output filename
FOUT=${DIROUT}.html
multiqc --filename ${FOUT} --outdir ${DIROUT} ${DIROUT}

conda deactivate
