# RNA-seq_WES_virus

Commonly viral sequences are captured in NGS data extracted from human samples. In order to screen NGS data of human cell lines, this workflow has been developed. It was conducted on publicly available human CCLE (https://sites.broadinstitute.org/ccle) RNA-seq and whole-exome-sequencing (WES) data sets. In particular DLBCL, LCLL, multiple myeloma and some virus negative cell lines were screened for viral sequences. By evaluating with other wet lab methods, RNA-seq turned out to reliably detect existing viral existence in the samples. The analysis results have been published in Plos One (PMID: 30629668 PMCID: PMC6328144 DOI: 10.1371/journal.pone.0210404). CCLE sequencing data is available via https://portal.gdc.cancer.gov/legacy-archive/search/f - since NGS data is provided as bam files, these have to be translated back to fastq files before starting this workflow.



## Creating STAR aligner index files

Before starting the alignment via STAR, index files need to be set up. In order to obtain unique mapping hits, short read sequences were mapped to the human genome and selected viruses simultaneously. For this, fasta files for Homo sapiens were retrieved from GATK (hg38 including ebv) and for viruses from NCBI (hbv:NC_003977.2.fa, hcv:NC_004102.1.fa, hhv8:NC_009333.1.fa, hiv1:NC_001802.1.fa, hiv2:NC_001722.1.fa, hpv:NC_004500.1.fa, htlv1:NC_001436.1.fa, htlv2:NC_001488.1.fa, mlv.AF221065.1.fa, smrv:NC_001514.1.fa, xmrv:FN692043.2.fa) and merged to one fasta file, which served for producing the STAR index files.
```
mkdir /path/to/star_index/Genomes/gatk/hg38/star_viruses
cd /path/to/star_index/Genomes/gatk/hg38/star_viruses
cat  /path/to/Genomes/gatk/hg38/Homo_sapiens_assembly38.fasta /path/to/star_index/Genomes/ncbi/viruses/*fa > hg38.viruses.fa
sbatch -c 8 -p short --mem=300000 -J star -o star.out -e star.err <<EOF
#!/bin/sh
$star --runMode genomeGenerate --genomeDir . --genomeFastaFiles hg38.viruses.fa --runThreadN 8 --limitGenomeGenerateRAM 300000000000 &> build.out
EOF
```


## Preparing project folder and subfolders for RNA-seq
```
DIR=/path/to/project_folder
mkdir $DIR
cd $DIR
mkdir raw
raw=/path/to/fastq_files/cghub/rna
# link to fastq files; DLBCL, LCLL, multiple myeloma and some virus negative cell lines were selected
cells=(BL-41 BL-70 CA-46 CI-1 DAUDI DB DEL DOHH-2 EB-1 GRANTA-519 HD-MY-Z HH HT JEKO-1 KARPAS-299 KARPAS-422 KI-JK KM-H2 L-1236 L-428 L-540 MC-116 NAMALWA NU-DHL-1 NU-DUL-1 OCI-LY19 RAJI REC-1 RI-1 RL SR-786 SU-DHL-10 SU-DHL-1 SU-DHL-4 SU-DHL-5 SU-DHL-6 SU-DHL-8 SUP-M2 SUP-T1 U-937 WSU-DLCL2 KYSE-520 KYSE-70 697 ALL-SIL AML-193 BV-173 CMK CML-T1 EHEB EM-2 EOL-1 F-36P GDM-1 HEL HL-60 HPB-ALL JK-1 JURKAT JURL-MK1 JVM-2 JVM-3 K-562 KASUMI-1 KASUMI-2 KASUMI-6 KCL-22 KE-37 KG-1 KOPN-8 KU-812 KYO-1 LAMA-84 LOUCY M-07e ME-1 MEC-1 MEG-01 MHH-CALL-2 MHH-CALL-3 MHH-CALL-4 MOLM-13 MOLM-16 MOLM-6 MOLT-13 MOLT-16 MOLT-3 MONO-MAC-1 MONO-MAC-6 MUTZ-3 MUTZ-5 MV4-11 NALM-19 NALM-1 NALM-6 NB-4 NOMO-1 OCI-AML2 OCI-AML3 OCI-AML5 OCI-M1 P12-ICHIKAWA PEER PF-382 PL-21 Reh RPMI-8402 RS4-11 SEM SET-2 SIG-M5 SKM-1 SUP-B15 SUP-T11 TALL-1 TF-1 THP-1 KELLY SIMA HEP-3B HEP-G2 LXF-289 MSTO-211H SCLC-21H SK-MES-1 22RV1 DU-145  SK-MEL-1 SK-MEL-30 ML-1 S-117 BCP-1 AMO-1 EJM KARPAS-620 KMS-12-BM L-363 LP-1 MOLP-2 MOLP-8 NCI-H929 OPM-2 RPMI-8226 SK-MM-2)
ids=(G27374 G27333 G27359 G27371 G27352 G20488 G27318 G27375 G41717 G27335 G41694 G28867 G28835 G28842 G28063 G28013 G28046 G28039 G28032 G28022 G28040 G28058 G28578 G28542 G28540 G28543 G27519 G27498 G27521 G41675 G30588 G30631 G30627 G30609 G30641 G30573 G30556 G41688 G30561 G41727 G30568 G28054 G28067 G26223 G27288 G26259 G27347 G27390 G27344 G27378 G27324 G27380 G26244 G27388 G28844 G26200 G28828 G28000 G28080 G28010 G28008 G28068 G25242 G30554 G28002 G28043 G28049 G28009 G28028 G28086 G28055 G28037 G26177 G28006 G28007 G28048 G41748 G41728 G28551 G28609 G28532 G41692 G28555 G28544 G28576 G28565 G28623 G26187 G26207 G28621 G28530 G41677 G28531 G28567 G26180 G26239 G26205 G26238 G27236 G27258 G41747 G28600 G27500 G27485 G27530 G27250 G27474 G27290 G27232 G41659 G27458 G27510 G30569 G30618 G30610 G30593 G27241 G28082 G27482 G28888 G41665 G28047 G25204 G27237 G27245 G20502 G20506 G27463 G27243 G28541 G27483 G27358 G27213 G26178 G30567 G26182 G26221 G26193 G28550 G28560 G28571 G27263 G27229 G27201 )
for (( i=0; i<${#cells[*]}; i++ ))
do
    mkdir $DIR/raw/${cells[i]}
    cd $DIR/raw/${cells[i]}
    ln -s $raw/${ids[i]}.* .
done
```


## Run trimming, merging, and aligning on SLURM for RNA-seq
```
# Get sample names (string)
SAMPLES=$(find $DIR/raw/* -type d)
DIR=/path/to/project_folder
cd $DIR

# copy adapter sequences and needed shell scripts into this folder, too
# samples to be iterated need to be in single folder
cd $DIR
mkdir trim fastqc star bam
source ~/.bashrc
for SAMPLE in $SAMPLES
do
	cd $SAMPLE
	sample=${PWD##*/}
	# trimming
	RES=$(sbatch -J $sample.1 -o $sample.1.out -e $sample.1.err ../../trim_new.sh)
	# joining reads
	RES2=$(sbatch --dependency=afterok:${RES##* } -J $sample.2  -o $sample.2.out -e $sample.2.err ../../merge_reads_new.sh)
	# alignment via star for merged sequences for human and human+virus genomes, >50min
	sbatch --dependency=afterok:${RES2##* } -J $sample.3  -o $sample.3.out -e $sample.3.err ../../align_star_merge_ch_new.sh
done
```


## Surveying jobs and get statistics
```
ls -lh $DIR/raw/*/*.1.err
ls -lh $DIR/raw/*/*.2.err
ls -lh $DIR/raw/*/*.3.err
# remove files, if not needed
# rm $DIR/raw/*/*.1.err
# rm $DIR/raw/*/*.2.err
# rm $DIR/raw/*/*.3.err

# trimming statistics
grep "by an average of" $DIR/raw/L*/A.*.o*
# merging
grep "Total joined: " $DIR/raw/L*/*.o*
grep "Total reads: " $DIR/raw/L*/*.o*

# alignment statistics
grep "mapped reads %" $DIR/star/L*final*
# grep "mapped reads n" $DIR/star/L*final*
grep "Number of input reads" $DIR/star/L*final*
grep "multiple loci" $DIR/star/L*final*
```


## Preparing project folder and subfolders for WES
```
DIR=/path/to/project_folder
cd $DIR
# all DLBCL + LCLL, and at least one virus negative cell line in group to virus positive cell line
cells=(BL-41 BL-70 CI-1 GRANTA-519 HT L-1236 MINO NAMALWA NU-DHL-1 OCI-LY19 OCI-LY3 SU-DHL-10 SU-DHL-4 SU-DHL-5 SU-DHL-6 SU-DHL-8 SUP-HD1 SUP-M2 U-937 WSU-DLCL2 ALL-SIL AML-193 CMK CML-T1 DND-41 F-36P GDM-1 HPB-ALL JK-1 JURKAT JURL-MK1 K-562 KASUMI-2 KASUMI-6 KCL-22 KOPN-8 KYO-1 M-07e ME-1 MEC-1 MHH-CALL-2 MHH-CALL-3 MHH-CALL-4 MOLM-13 MOLM-16 MOLM-6 MOLT-3 MUTZ-5 NALM-19 NALM-1 NB-4 OCI-AML5 OCI-M1 PEER RPMI-8402 SEM SIG-M5 TF-1 HEP-3B DV-90 22RV1 LNCaP_clone_FGC AMO-1 EJM JJN-3 KARPAS-620 KMS-12-BM LP-1 MOLP-2 MOLP-8)
ids=(BL-41 BL-70 CI-1 GRANTA-519 HT L-1236 Mino NAMALWA NU-DHL-1 OCI-LY-19 OCI-LY3 SU-DHL-10 SU-DHL-4 SU-DHL-5 SU-DHL-6 SU-DHL-8 SUP-HD1 SUP-M2 U-937 WSU-DLCL2 ALL-SIL AML-193 CMK CML-T1 DND-41 F-36P GDM-1 HPB-ALL JK-1 JURKAT JURL-MK1 K-562 KASUMI-2 Kasumi-6 KCL-22 KOPN-8 KYO-1 M-07e ME-1 MEC-1 MHH-CALL-2 MHH-CALL-3 MHH-CALL-4 MOLM-13 MOLM-16 MOLM-6 MOLT-3 MUTZ-5 NALM-19 NALM-1 NB-4 OCI-AML5 OCI-M1 PEER RPMI-8402 SEM SIG-M5 TF-1 Hep_3B2.1-7 DV-90 22Rv1 LNCaP_clone_FGC AMO-1 EJM JJN-3 KARPAS-620 KMS-12-BM LP-1 MOLP-2 MOLP-8)

mkdir raw_wes fastqc_wes bam_wes star_wes
raw=/path/to/fastq_files/cghub/wes
for (( i=0; i<${#cells[*]}; i++ ))
do
mkdir $DIR/raw_wes/${cells[i]}
cd $DIR/raw_wes/${cells[i]}
if [ ${cells[i]} = "K-562" ]; then
    ln -s $raw/C835.${ids[i]}.* .
else
    ln -s $raw/C836.${ids[i]}.* .
fi
done
```


## Run trimming, merging, and alignment on SLURM for WES
```
# Get sample names (string)
SAMPLES=$(find $DIR/raw_wes/* -maxdepth 1 -type d)
DIR=/path/to/project_folder
cd $DIR

# copy adapter sequences and needed shell scripts into this folder, too
# samples to be iterated need to be in single folder
for SAMPLE in $SAMPLES
do
	cd $SAMPLE
	sample=${PWD##*/}
	# trimming
	RES=$(sbatch -J $sample.1 -o $sample.w1.out -e $sample.w1.err ../../trim_wes.sh)
	# joining reads
	RES2=$(sbatch --dependency=afterok:${RES##* } -J $sample.w2  -o $sample.w2.out -e $sample.w2.err ../../merge_reads_wes.sh)
	# alignment via star for merged sequences for human and human+virus genomes, >50min
	sbatch --dependency=afterok:${RES2##* } -J $sample.w3  -o $sample.w3.out -e $sample.w3.err ../../align_star_merge_ch_wes.sh
done
```


## Summarise and visualise data

```
# analysis via R, summarise data
DIR=/path/to/project_folder
cd $DIR

# info on sequencing + mapping, star statistics
mkdir tables plots
R --file="stats_seq_rna.R" --save > R.out # rna-seq
# unique mapping hits
R --file="virus_hits.R" >> R.out # rna-seq
R --file="virus_hits_wes.R" # wes

```
