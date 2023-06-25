# Scripts

High-performance computing scripts for the [ATAC-seq pipeline](https://github.com/reskejak/ATAC-seq) developed by [
Jake Reske](https://github.com/reskejak). See [below](#atac-seq-original-readme) how to cite this pipeline. Adapted for [PBS job scheduler](https://en.wikipedia.org/wiki/Portable_Batch_System), but easily modifiable for other systems. Each script implements steps and sub-steps outlined in the [workflow schema](https://github.com/reskejak/ATAC-seq#atacseq_workflowtxt).

Prerequisites: Tools must be installed within [Conda](https://docs.conda.io/en/latest/miniconda.html) environment. I prefer to install tools in separate environments. E.g., for [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml):

```{bash}
# Create an environment with the name (-n) bowtie2, answer yes (-y) to all prompts
conda create -n bowtie2 -y
# Activate this environment
conda activate bowtie2
# Install bowtie2 in this environment, answer yes (-y) to all prompts
# Google "conda install tool_name" for commands installing other tools
conda install -c bioconda bowtie2 -y
# Always good to have samtools installed
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y samtools==1.11 
# Deactivate the environment
conda deactivate
```

Check and adjust each script for the environment and tool names.

- [01a_QC_fastqc_raw.sh](01a_QC_fastqc_raw.sh) - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on raw FASTQ files, and [multiqc](https://multiqc.info/) summary. Requires `trimgalore` environment with `conda install -c bioconda trim-galore` and `conda install -c bioconda multiqc`.
  - Input: raw FASTQ files in the `00_raw` folder, `fastq.gz` extension. Paired-end (R1, R2 suffixes).
  - Output: `00_raw_fastqc` folder with HTML QC reports for individual files and the `00_raw_fastqc.html` file of multiqc summary.

- [01b_QC_trimgalore.sh](01b_QC_trimgalore.sh) - [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) adapter trimming
  - Input: raw FASTQ files in the `00_raw` folder. Paired-end. Requires `trimgalore` environment.
  - Output: `01_trimmed` folder with `.fq.gz` trimmed reads. Also contains `_trimming_report.txt`, `_fastqc.html`, and `_fastqc.zip` quality reports (TrimGalore! additionally runs FastQC).

- [01c_QC_fastqc_trimmed.sh](01c_QC_fastqc_trimmed.sh) - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on trimmed FASTQ files, and [multiqc](https://multiqc.info/) summary. Requires `trimgalore` environment.
  - Input: trimmed FASTQ files in the `01_trimmed` folder, `fq.gz` extension. Paired-end (val_1, val_2 suffixes).
  - Output: `01_trimmed_fastqc` folder with HTML QC reports for individual files and the `01_trimmed_fastqc.html` file of multiqc summary.

- [02a_align_bowtie2_build.sh](02a_align_bowtie2_build.sh) - [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) index building. Requires `bowtie2` environment with `conda install -c bioconda bowtie2` and samtools (see above).
  - Input: Path to genome FASTA file with all chromosomes, contigs, and decoy sequences.
  - Output: Path and basename of the bowtie2 index.

- [02b_align_bowtie2.sh](02b_align_bowtie2.sh) - [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment. Requires `bowtie2` environment.
  - Input: `01_trimmed` folder with `.fq.gz` trimmed reads.
  - Output: `02_bowtie2` folder with aligned `.bam` files.

- [02c_align_samtools_sort.sh](02c_align_samtools_sort.sh) - sort the bowtie2-aligned BAM files.
  - Input: `02_bowtie2` folder with aligned `.bam` files. Requires `bowtie2` environment.
  - Output: `02_bowtie2` folder with sorted and indexed `_sorted.bam` files. The original BAM files are deleted to save space.

- [03a_filter_samtools.sh](03a_filter_samtools.sh) - removes reads overlapping mitochondrial chromosome, keeps only properly paired reads, resorts and reindexes the files. Requires `bowtie2` environment.
  - Input: `02_bowtie2` folder with sorted and indexed `_sorted.bam` files.
  - Output: `02_bowtie2` folder with sorted properly paired reads, chrM cleaned `_sorted_filt_noTM.bam` files. The original sorted and indexed BAM files, and intermediate ones, are deleted to save space.

- [04a_complexity_flagstat.sh](04a_complexity_flagstat.sh) - calculate rescaling factors from flagstat stats. Requires `bowtie2` environment.
  - Input: `02_bowtie2` folder with sorted properly paired reads, chrM cleaned `_sorted_filt_noTM.bam` files.
  - Output: `02_bowtie2_flagstat` with `flagstat*` files with all statistics and `scaling*` files with scaling factors only.

- [04b_complexity_samtools.sh](04b_complexity_samtools.sh) - subsample, sort and index files.  Requires `bowtie2` environment.
  - Input: `02_bowtie2` folder with sorted properly paired reads, chrM cleaned `_sorted_filt_noTM.bam` files. Also, `02_bowtie2_flagstat` with `scaling*` files with scaling factors only.
  - Output: `03_bowtie2_sub` folder with `_sub_sorted_filt_noMT.bam;` files.

- [04c_complexity_flagstat.sh](04c_complexity_flagstat.sh) - double-check flagstat statistics and rescaling factors for subsampled files. They should show similar number of reads and scaling factors approximately equal to 1. Requires `bowtie2` environment.
  - Input: `03_bowtie2_sub` folder with sorted properly paired reads, chrM cleaned `_sorted_filt_noTM.bam` files.
  - Output: `03_bowtie2_sub_flagstat` with `flagstat*` files with all statistics and `scaling*` files with scaling factors only.

- [05a_filter_picard.sh](05a_filter_picard.sh) - removing duplicates with Picard, resort the files.
  - Input: `03_bowtie2_sub` folder with `_sub_sorted_filt_noMT.bam;` files.
  - Output: `03_bowtie2_sub` folder with `_sorted_noDups_filt_noMT.bam` files.

- [06a_format_samtools.sh](06a_format_samtools.sh) - name sort, fix mate and convert BAM files to BEDPE.
  - Input: `03_bowtie2_sub` folder with `_sorted_noDups_filt_noMT.bam` files.
  - Output: `04_bedpe` folder with `*_fixed.bedpe` files. Intermediate (name sorted, mate fixed) BAM files are deleted.

- [06b_format_bedpeTn5shift.sh](06b_format_bedpeTn5shift.sh) - shift BEDPE coordinates +4 and -5 bp, respectively, to compensate for Tn5 insertion in ATAC data as described by Buenrostro et al. 2013 (PMID: 24097267). Also, convert standard bedtools BEDPE format to "minimal" BEDPE format accepted by MACS2 for peak calling. Uses [bedpeMinimalConvert.sh](bedpeMinimalConvert.sh).
  - Input: `04_bedpe` folder with `*_fixed.bedpe` files. 
  - Output: `04_bedpe` folder with `*_tn5.bedpe` and `*_minimal.bedpe` files. `*_fixed.bedpe` files kept.

- [07a_peak_macs2.sh](07a_peak_macs2.sh) - call peaks from BEDPE files, remove those overlapping blacklist.
  - Input: `04_bedpe` folder with `*_minimal.bedpe` files. Human blacklist `ENCFF356LFX_hg38.bed`.
  - Output: `05_macs2` folder with MACS2 output and `*_filt.broadPeak` filtered broadpeak files

- [07b_peak_macs2_pooled.sh](07b_peak_macs2_pooled.sh) - calling pooled peaks and selecting those overlapping with individual replicates at least 50%. Uses [naiveOverlapBroad.sh](naiveOverlapBroad.sh) for two replicates and [naiveOverlapBroad3.sh](naiveOverlapBroad3.sh) for three replicates. Need to be manually edited.
  - Input: `04_bedpe` folder with `*_minimal.bedpe` files. `05_macs2` folder with MACS2 output and `*_filt.broadPeak` filtered broadpeak files.
  - Output: `05_macs2` folder with `*_overlap_peaks_filt.broadPeak` files.

# ATAC-seq (original README)

These scripts correspond to a (differential) ATAC-seq analysis workflow as described in [our recent report](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00342-y). It is largely based on the pipeline developed by [Anshul Kundaje's group (Stanford) and the ENCODE project](https://www.encodeproject.org/pipelines/ENCPL792NWO/).

If you use this methodology, please cite the following paper along with corresponding pipeline dependencies:

Jake J. Reske, Mike R. Wilson, and Ronald L. Chandler. 2020. ATAC-seq normalization method can significantly affect differential accessibility analysis and interpretation. *Epigenetics & Chromatin* **13**: 22.

Attempt to run each command individually or in blocks after editing to match your own data architecture.

### ATACseq_workflow.txt
**Generalized ATAC-seq data processing workflow intended for comparative analysis.** Stepwise bioinformatics process and example commands for analyzing ATAC-seq data from raw reads to calling peaks for downstream differential accessibility analysis. Consider “treat1” as an example mouse ATAC-seq Illumina paired-end library. Blue text denotes optional or conditional steps dependent on experimental design and desired output. Users seeking only to discover replicate-concordant accessible regions in a singular cell state may wish to call naïve overlapping peaks, though this step is not necessary for differential accessibility analysis.

![foo](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13072-020-00342-y/MediaObjects/13072_2020_342_Fig4_HTML.png)

### csaw_workflow.R
***csaw* workflow for multiple differential accessibility analyses in R.** Consider an experimental design with *n* = 2 biological replicates from two conditions: “treat” and “control”. Describes implementation of two possible normalization methods and use of either *MACS2* peaks or *de novo* locally enriched windows as query regions for output comparison; see [*csaw* manual](https://bioconductor.org/packages/release/bioc/html/csaw.html) for additional normalization frameworks. ***Note: updates to the behavior of certain csaw functions have required slight compatibility changes to the commands described graphically, so please reference the latest R script.***

![foo2](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13072-020-00342-y/MediaObjects/13072_2020_342_Fig6_HTML.png)
