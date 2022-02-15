#!/bin/bash

#SBATCH -c 4                # number of core to be used
#SBATCH -t 0-04:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=10000         # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB

# Get sample name
sample=${PWD##*/}

eautils=/path/to/eautils/ExpressionAnalysis-ea-utils-27a4809/clipper
echo $eautils
/path/to/FastQC/fastqc -version



# merge paired-end reads
date
cd ../../trim
echo "Merge paired-end reads..."
$eautils/fastq-join $sample.R1.fq $sample.R2.fq -o $sample.fq.
rm $sample.R1.fq $sample.R2.fq
echo "Done."


# zip trimmed sequences
date
echo "gzip fastq files..."
pigz -p4 $sample.fq.join
pigz -p4 $sample.fq.un1
pigz -p4 $sample.fq.un2
echo "Done."


# quality control via fastQC, plots
date
echo "Quality control for trimmed sequences of $sample..."
/path/to/FastQC/fastqc -t 3 $sample.fq.join.gz $sample.fq.un1.gz $sample.fq.un2.gz
# move to folder
mv $sample.*html ../fastqc/
rm $sample.*fastqc.zip
echo "Done."



date


