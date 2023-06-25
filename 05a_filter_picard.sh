#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mikhail.dozmorov@vcuhealth.org
#PBS -m ae
#PBS -N picard
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

source activate picard

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=03_bowtie2_sub
# Output folder
DIROUT=03_bowtie2_sub
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find $DIRIN -type f -name '*_sub_sorted_filt_noMT.bam' | sort`; do
	# Removes _sub_sorted_filt_noMT.bam and adds _noDups_filt_noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename $file _sub_sorted_filt_noMT.bam`_noDups_filt_noMT.bam;
	picard MarkDuplicates I=$file O=${FILEOUT} \
		M=${DIROUT}/`basename $file _sub_sorted_filt_noMT.bam`.txt REMOVE_DUPLICATES=true;
	# Sort and index
	FILEOUT1=${DIROUT}/`basename $file _sub_sorted_filt_noMT.bam`_sorted_noDups_filt_noMT.bam;
	samtools sort -o ${FILEOUT1} -@ ${THREADS} ${FILEOUT}
	samtools index ${FILEOUT1}
	rm ${FILEOUT} ${file} ${file}.bai
done

