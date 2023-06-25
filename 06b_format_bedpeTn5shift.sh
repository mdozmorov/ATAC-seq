#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N bedpeshift
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
#PBS -q workq

cd $PBS_O_WORKDIR

source activate samtools

# Number of CPUs, should match the ppn definition
THREADS=1

# Input folder
DIRIN=04_bedpe
# Output folder
DIROUT=04_bedpe
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find ${DIRIN} -name "*_fixed.bedpe" -type f | sort`; do
	# Removes _fixed.bedpe and adds _tn5.bedpe to the file name
	FILEOUT=${DIROUT}/`basename $file _fixed.bedpe`_tn5.bedpe;
	# Shift BEDPE coordinates +4 and -5 bp, respectively, to compensate for Tn5 insertion in ATAC data
	bash bedpeTn5shift.sh ${file} > ${FILEOUT}
	# convert standard bedtools BEDPE format to "minimal" BEDPE format accepted by MACS2 for peak calling
	# FILEOUT1=${DIROUT}/`basename $file _fixed.bedpe`_minimal.bedpe;
	# bash bedpeMinimalConvert.sh ${FILEOUT} > ${FILEOUT1}
	# rm $file ${FILEOUT}
done

conda deactivate

# Must be run from the directory where the files are 
cd ${DIROUT}

for file in `ls *_tn5.bedpe`; do
	echo ${file}
	FILEOUT=`basename $file _tn5.bedpe`_minimal.bedpe
	bash ../bedpeMinimalConvert.sh ${file} > ${FILEOUT}
done

cd ..
