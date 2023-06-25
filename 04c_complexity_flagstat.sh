#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N flagstat
#PBS -j oe
#PBS -q workq
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

source activate bowtie2

# Number of CPUs, should match the ppn definition
THREADS=24

# Input folder
DIRIN=03_bowtie2_sub
# Output folder
DIROUT=${DIRIN}_flagstat
# Make output folder if doesn't exist
mkdir -p ${DIROUT}

for file in `find ${DIRIN} -name "*_sub_sorted_filt_noMT.bam" -type f | sort`; do
    FOUT=${DIROUT}/flagstat_`basename ${file} _sub_sorted_filt_noMT.bam`.tsv
    samtools flagstat -@ ${THREADS} -O tsv $file > $FOUT; 
done

# Minimum depth
MINDEPTH=`head -q -n1 ${DIROUT}/flagstat*.tsv | cut -f1 | sort | head -n1`
echo "Minimum depth is" ${MINDEPTH}
# Scaling factor for each file
for file in `find ${DIROUT} -type f -name "flagstat*.tsv"`; do
    # Depth for the current file
    FILEDEPTH=`head -q -n1 ${file} | cut -f1`
    # Scaling factor for the current file
    FILESCALING=$( echo "scale=4; ${MINDEPTH} / ${FILEDEPTH}" | bc | awk '{printf "%f", $0}')
    # File to output scaling factor
    FOUT=${DIROUT}/`basename ${file} | sed 's/flagstat/scaling/'`
    echo ${FILESCALING} > ${FOUT}
done

conda deactivate
