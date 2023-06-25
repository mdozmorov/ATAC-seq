#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N macs2
#PBS -j oe
#PBS -q workq
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

source activate macs2

# Number of CPUs, should match the ppn definition
THREADS=2

# Input folder
DIRIN=04_bedpe
# Output folder
DIROUT=05_macs2
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

# hg38 blacklisted: https://www.encodeproject.org/annotations/ENCSR636HFF/
# wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
FILEBLACKLISTED=ENCFF356LFX_hg38.bed

for file in `find $DIRIN -type f -name "*_minimal.bedpe"`; do
	FILEOUT=${DIROUT}/`basename $file _minimal.bedpe`
	macs2 callpeak -t $file -f BEDPE -n ${FILEOUT} -g hs --broad --broad-cutoff 0.05 --keep-dup all

	FILEOUT1=${DIROUT}/`basename $file _minimal.bedpe`_peaks_filt.broadPeak
	bedtools intersect -v -a ${FILEOUT}_peaks.broadPeak -b ${FILEBLACKLISTED} \
		| grep -P 'chr[\dXY]+[\t]' > ${FILEOUT1}

done

conda deactivate
