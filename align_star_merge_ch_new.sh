#!/bin/bash

#SBATCH -c 14                # number of core to be used
#SBATCH -t 0-04:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=40000         # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB

# Get sample name
sample=${PWD##*/}
star=/path/to/STAR-2.5.2a/bin/Linux_x86_64/STAR
genome=/path/to/star_index/Genomes/gatk/hg38/star_viruses
echo $star
echo $genome


# alignment via STAR, hg38+viruses hybrid genome
cd ../../bam
date
echo "Alignment via STAR, hg38+viruses hybrid genome, joined seqs..."
$star --outFileNamePrefix ../star/$sample.join. --genomeDir $genome --runThreadN 14 --readFilesIn ../trim/$sample.fq.join.gz --outStd BAM_SortedByCoordinate --chimSegmentMin 20 --outSAMtype BAM SortedByCoordinate --outSAMmode Full --outFilterMismatchNmax 10 --readFilesCommand zcat --twopassMode Basic > $sample.join.bam
echo "Done."
date
echo "Alignment via STAR, hg38+viruses hybrid genome, unmerged seqs..."
$star --outFileNamePrefix ../star/$sample.un. --genomeDir $genome --runThreadN 14 --readFilesIn ../trim/$sample.fq.un*.gz --outStd BAM_SortedByCoordinate --chimSegmentMin 20 --outSAMtype BAM SortedByCoordinate --outSAMmode Full --outFilterMismatchNmax 10 --readFilesCommand zcat --twopassMode Basic > $sample.un.bam
echo "Done."
date
echo "samtools merge for $sample..."
samtools merge -@ 14 $sample.j.s.bam $sample.join.bam $sample.un.bam
rm $sample.join.bam $sample.un.bam
samtools index -@ 14 $sample.j.s.bam
echo "Count all reads..."
samtools idxstats $sample.j.s.bam > $sample.j.idxstats
echo "Done."
date
echo "Count unique reads..."
samtools view -b -@ 14 -q 255 $sample.j.s.bam > $sample.u.bam
samtools index -@ 14 $sample.u.bam
samtools idxstats $sample.u.bam > $sample.u.idxstats
rm $sample.u.bam*
echo "Done."
date

rm ../trim/$sample.*gz
rm -r ../star/$sample.*._STARtmp ../star/$sample.*.Log.out ../star/$sample.*._STARgenome ../star/$sample.*._STARpass1
date

echo "Save chimeric reads + fasta files..."
mv ../star/$sample.join.Chimeric.out.sam .
mv ../star/$sample.un.Chimeric.out.sam .
samtools view -@ 14 -bS $sample.join.Chimeric.out.sam > $sample.join.Chimeric.out.bam
samtools sort -@ 14 $sample.join.Chimeric.out.bam -o $sample.join.Chimeric.bam
samtools view -@ 14 -bS $sample.un.Chimeric.out.sam > $sample.un.Chimeric.out.bam
samtools sort -@ 14 $sample.un.Chimeric.out.bam -o $sample.un.Chimeric.bam
rm $sample.join.Chimeric.out.sam $sample.un.Chimeric.out.sam $sample.join.Chimeric.out.bam $sample.un.Chimeric.out.bam
samtools merge -@ 14 $sample.Chimeric.bam $sample.join.Chimeric.bam $sample.un.Chimeric.bam
samtools index -@ 14 $sample.Chimeric.bam
rm $sample.join.Chimeric.bam $sample.un.Chimeric.bam
~/Programme/samtools-1.3/samtools fasta $sample.Chimeric.bam > $sample.Chimeric.fa
echo "Done."
date

