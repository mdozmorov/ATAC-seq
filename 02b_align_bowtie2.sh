#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N bowtiealn
#PBS -j oe
#PBS -q workq
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

source activate bowtie2

# Path to genome index
## Human
BOWTIE2_IDX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38ext/hg38ext
## Mouse
# BOWTIE2_IDX=/home/mdozmorov/sequencing/data/ExtData/UCSC/mm10/mm10

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=01_trimmed
# Output folder
DIROUT=02_bowtie2
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

# Paired end
k=0
for file in `find ${DIRIN} -name "*.fq.gz" -type f | sort`; do
	if [ $k = 0 ]
	then
		s1=$file
		k+=1
	else
		s2=$file
		k=0
		# Adjust, removes the end of the filename from the second file in a pair
		FILEOUT=`basename ${s2} _R2_001_val_2.fq.gz`;
		bowtie2 --very-sensitive -X 1000 -x ${BOWTIE2_IDX} -1 ${s1} -2 ${s2} --threads ${THREADS} --mm \
			| samtools view -bS - > ${DIROUT}/${FILEOUT}.bam;
	fi
done

# # Single end
# for file in `find $DIRIN -name "*.fastq.gz" -type f | sort`; do
# 	FILEIN=$file;
# 	FILEOUT=`basename $file .fastq.gz`;
# 	bowtie2 --end-to-end --no-unal --no-mixed --no-discordant --dovetail --phred33 -p 12 -x ${BOWTIE2_IDX} -U $FILEIN -S ${DIROUT}/${FILEOUT}_bowtie2.sam > ${DIROUT}/${FILEOUT}_bowtie2.txt;
# done

conda deactivate
