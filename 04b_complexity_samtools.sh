#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N samsortnoMT
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
#PBS -q workq

cd $PBS_O_WORKDIR

source activate bowtie2

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=02_bowtie2
# Stats folder
DIRIN1=${DIRIN}_flagstat
# Output folder
DIROUT=03_bowtie2_sub
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find ${DIRIN} -name "*_sorted.bam_sorted_filt_noMT.bam" -type f | sort`; do
	# Get scaling factor
	FILEIN=${DIRIN1}/scaling_`basename $file _sorted.bam_sorted_filt_noMT.bam`.tsv
	FILESCALING=`cat ${FILEIN}`
	# Removes _noMT.bam and adds _sorted_noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename $file _sorted.bam_sorted_filt_noMT.bam`_sub_filt_noMT.bam;
	# If scaling factor is 1, simply copy the file
	if [ ${FILESCALING} -eq 1 ]
	then
		cp $file ${FILEOUT}
	else
		# Actual downsampling, see -s option at http://www.htslib.org/doc/samtools-view.html
		decimal_part=$(echo "$FILESCALING" | awk -F. '{print $2}')
		samtools view -h -b -s 1.${decimal_part} ${file} > ${FILEOUT}
	fi
	# Sort and index
	FILEOUT1=${DIROUT}/`basename $file _sorted.bam_sorted_filt_noMT.bam`_sub_sorted_filt_noMT.bam;
	samtools sort -o ${FILEOUT1} -@ ${THREADS} ${FILEOUT}
	samtools index ${FILEOUT1}
	rm ${FILEOUT}
done

conda deactivate
