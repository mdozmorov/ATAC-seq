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
DIRIN1=05_macs2
# Output folder
DIROUT=05_macs2_pooled
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

# hg38 blacklisted: https://www.encodeproject.org/annotations/ENCSR636HFF/
# wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
FILEBLACKLISTED=ENCFF356LFX_hg38.bed

# Control, 2 replicates
FILEIN1=${DIRIN}/SAN_CTL2_S60_minimal.bedpe
FILEIN2=${DIRIN}/SAN_CTL3_S61_minimal.bedpe
FILEOUT=${DIROUT}/SAN_CTL_pool
macs2 callpeak -t ${FILEIN1} ${FILEIN2} -f BEDPE -n ${FILEOUT} -g hs --broad --broad-cutoff 0.05 --keep-dup all

FILEOUT1=${DIROUT}/SAN_CTL_pool_peaks_filt.broadPeak
bedtools intersect -v -a ${FILEOUT}_peaks.broadPeak -b ${FILEBLACKLISTED} \
	| grep -P 'chr[\dXY]+[\t]' > ${FILEOUT1}

FILEIN1=${DIRIN1}/SAN_CTL2_S60_peaks_filt.broadPeak
FILEIN2=${DIRIN1}/SAN_CTL3_S61_peaks_filt.broadPeak
FILEIN3=${DIROUT}/SAN_CTL_pool_peaks_filt.broadPeak
FILEOUT1=${DIROUT}/SAN_CTL_overlap_peaks_filt.broadPeak
bash naiveOverlapBroad.sh ${FILEIN1} ${FILEIN2} ${FILEIN3} > ${FILEOUT1}


# Treatment, 3 replicates
FILEIN1=${DIRIN}/SAN_DAC1_S62_minimal.bedpe
FILEIN2=${DIRIN}/SAN_DAC2_S63_minimal.bedpe
FILEIN3=${DIRIN}/SAN_DAC3_S64_minimal.bedpe
FILEOUT=${DIROUT}/SAN_DAC_pool
macs2 callpeak -t ${FILEIN1} ${FILEIN2} ${FILEIN3} -f BEDPE -n ${FILEOUT} -g hs --broad --broad-cutoff 0.05 --keep-dup all

FILEOUT1=${DIROUT}/SAN_DAC_pool_peaks_filt.broadPeak
bedtools intersect -v -a ${FILEOUT}_peaks.broadPeak -b ${FILEBLACKLISTED} \
	| grep -P 'chr[\dXY]+[\t]' > ${FILEOUT1}

FILEIN1=${DIRIN1}/SAN_DAC1_S62_peaks_filt.broadPeak
FILEIN2=${DIRIN1}/SAN_DAC2_S63_peaks_filt.broadPeak
FILEIN3=${DIRIN1}/SAN_DAC3_S64_peaks_filt.broadPeak
FILEIN4=${DIROUT}/SAN_DAC_pool_peaks_filt.broadPeak
FILEOUT1=${DIROUT}/SAN_DAC_overlap_peaks_filt.broadPeak
bash naiveOverlapBroad3.sh ${FILEIN1} ${FILEIN2} ${FILEIN3} ${FILEIN4} > ${FILEOUT1}


conda deactivate
