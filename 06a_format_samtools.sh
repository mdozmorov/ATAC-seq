#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N samformat
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
#PBS -q workq

cd $PBS_O_WORKDIR

source activate samtools

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=03_bowtie2_sub
# Output folder
DIROUT=04_bedpe
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find ${DIRIN} -name "*_sorted_noDups_filt_noMT.bam" -type f | sort`; do
	# Removes _sorted_noDups_filt_noMT.bam and adds _namesorted_noDups_filt_noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename $file _sorted_noDups_filt_noMT.bam`_namesorted_noDups_filt_noMT.bam;
	# Sort by name
	samtools sort -n -o ${FILEOUT} -@ ${THREADS} ${file}
	# Fix mate
	FILEOUT1=${DIROUT}/`basename $file _sorted_noDups_filt_noMT.bam`_fixed.bam
	samtools fixmate ${FILEOUT} ${FILEOUT1}
	rm ${FILEOUT} ${FILEOUT}.bai
	# Convert to BEDPE
	FILEOUT2=${DIROUT}/`basename $file _sorted_noDups_filt_noMT.bam`_fixed.bedpe
	samtools view -bf 0x2 ${FILEOUT1} | bedtools bamtobed -i stdin -bedpe > ${FILEOUT2}
	rm ${FILEOUT1}
done

conda deactivate
