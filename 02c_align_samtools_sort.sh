#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N samsort
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
#PBS -q workq

cd $PBS_O_WORKDIR

source activate bowtie2

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=02_bowtie2
# Output folder
DIROUT=02_bowtie2
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find ${DIRIN} -name "*.bam" -type f | sort`; do
	FILEIN=$file;
	# Removes .bam and adds _sorted.bam to the file name
	FILEOUT=${DIROUT}/`basename $file .bam`_sorted.bam;
	samtools sort -o ${FILEOUT} -@ ${THREADS} ${FILEIN}
	samtools index ${FILEOUT}
	rm ${FILEIN}
done

conda deactivate
