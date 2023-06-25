#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N samsortfilt
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

for file in `find ${DIRIN} -name "*_sorted.bam" -type f | sort`; do

	# removes reads overlapping mitochondrial chromosome
	FILEIN=$file;
	# Removes _sorted.bam and adds _noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename $file _sorted.bam`_noMT.bam;
	samtools view -h ${FILEIN} | python removeChrom.py - - chrM | samtools view -bh - > ${FILEOUT}
	rm ${FILEIN}
	rm ${FILEIN}.bai

	# sort the chrM cleaned BAM files
	FILEIN=${FILEOUT};
	# Removes _noMT.bam and adds _sorted_noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename ${FILEIN} _noMT.bam`_sorted_noMT.bam;
	samtools sort -o ${FILEOUT} -@ ${THREADS} ${FILEIN}
	samtools index ${FILEOUT}
	rm ${FILEIN}

	# keep only properly paired reads (-f 3 option)
	FILEIN=${FILEOUT}
	# Removes _sorted_noMT.bam and adds _filt_noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename $file _sorted_noMT.bam`_filt_noMT.bam;
	samtools view -bh -f 3 ${FILEIN} > ${FILEOUT}
	rm ${FILEIN}
	rm ${FILEIN}.bai

	# sort filtered files
	FILEIN=${FILEOUT}
	# Removes _filt_noMT.bam and adds _sorted_filt_noMT.bam to the file name
	FILEOUT=${DIROUT}/`basename $file _filt_noMT.bam`_sorted_filt_noMT.bam;
	samtools sort -o ${FILEOUT} -@ ${THREADS} ${FILEIN}
	samtools index ${FILEOUT}
	rm ${FILEIN}

done

conda deactivate
