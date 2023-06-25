#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N trim
#PBS -q workq
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

source activate trimgalore

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=00_raw
# Output folder
DIROUT=01_trimmed
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

# Paired end
k=0
for file in `find ${DIRIN} -name "*.fastq.gz" -type f | sort`; do
	if [ $k = 0 ]
 	then
 		s1=$file
 		k+=1
 	else
 		s2=$file
 		k=0
        trim_galore --fastqc --paired -o ${DIROUT} --cores ${THREADS} ${s1} ${s2}
 	fi
done

# # Single end
# for file in `find $DIRIN -name "*.fastq.gz" -type f | sort`; do
# 	# java -jar trimmomatic-0.33.jar SE -threads 24 -phred33 $file $DIROUT"/"`basename $file .fastq.gz`.fastq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;
# 	trim_galore --fastqc -o $DIROUT --cores 4 $file
# done

conda deactivate
