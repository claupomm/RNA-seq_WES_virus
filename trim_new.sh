#!/bin/bash

#SBATCH -c 2              # number of core to be used
#SBATCH -t 0-04:00        # estimated run-time in D-HH:MM
#SBATCH -p short          # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=5000        # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB

# Get sample name
sample=${PWD##*/}

eautils=/path/to/eautils/ExpressionAnalysis-ea-utils-27a4809/clipper
echo $eautils

# trim sequences
date
echo "Trim sequences for $sample..."
$eautils/fastq-mcf ../../illu_ad_sel.fa *_R1*.fastq.gz *_R2*.fastq.gz -o ../../trim/$sample.R1.fq -o ../../trim/$sample.R2.fq -q 20 -H
echo "Done."




