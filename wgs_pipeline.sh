#!/usr/bin/env bash
# =============================================================================
# Bacterial Whole-Genome Sequencing (WGS) Workflow — Illumina + Oxford Nanopore
# -----------------------------------------------------------------------------
# A hybrid bacterial-genomics training workflow covering:
#   * Short-read (Illumina) QC + de novo assembly  : FastQC, fastp, BBMap, SPAdes, Bandage
#   * Long-read (Oxford Nanopore) QC + assembly     : NanoPlot, Flye, dnaapler
#   * Reference-based variant calling               : Bowtie2, SAMtools, GATK, BCFtools, FreeBayes
#   * Structural-variant calling                    : Manta, DELLY, CNVkit
#   * Annotation / typing                           : snpEff, VEP, Prokka, Bakta, AMRFinderPlus, abricate
#
# Data: public NCBI SRA accessions (Agrobacterium, Shigella, E. coli).
# Completed as a training assignment — see README.md for attribution.
#
# NOTE: This file documents the commands run during the training. It is a
# reference / notebook, not a turn-key script — paths, sample names and
# accessions should be adapted to your working directory before running.
# Server addresses, logins and personal paths have been redacted (<REMOTE_HOST>,
# <user>).
# =============================================================================

# Bari Bacterial WGS Workflow_Oxford Nanopore_Data

# Download anaconda 
wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh

# Install anaconda
chmod +x Anaconda3-2024.06-1-Linux-x86_64.sh

# Check Anaconda Version
conda --version

# Check Python Version (Installed via Anaconda)
python --version

# Initialize Anaconda by running:
conda init
source ~/.bashrc

# Create channel and conda environment
conda config --show channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n wgs_env -y

# Install and activate software
conda install bioconda::kingfisher
conda install bioconda::bbmap 
conda install bioconda::fastqc 
conda install bioconda::fastp 
conda install bioconda::spades 
conda install bioconda::bandage

# Downloading Raw Data from NCBI
SRR24759617
SRR24759618
SRR10506616
SRR28839724
SRR32486128
SRR32486129

# To check the folder files
ls nanopore/
ls illumina/

mv SRR24759617_1.fastq.gz Agro_i_R1.fastq.gz
mv SRR24759617_2.fastq.gz Agro_i_R2.fastq.gz
mv SRR10506616_1.fastq.gz Shig1_i_R1.fastq.gz
mv SRR10506616_2.fastq.gz Shig1_i_R2.fastq.gz
mv SRR28839724_1.fastq.gz Shig2_i_R1.fastq.gz
mv SRR28839724_2.fastq.gz Shig2_i_R2.fastq.gz

# Software Installation
conda install bioconda::kingfisher

# Using Kingfisher to download reads from SRA (or any other NCBI database)
kingfisher get -r $biosample -m ena-ascp ena-ftp aws-http prefetch

# make list of the target data, save it in name sra_list.tab
sra_list.tab
SRR24759617	Agro_illumina
SRR24759618	Agro_nanopore
SRR10506616	Shig1_illumina
SRR32486128	Shig1_nanopore
SRR28839724	Shig2_illumina
SRR32486129	Shig2_nanopore

# make a script individual_kf.sh for downloading these genomes
individual_kf.sh
#!/bin/bash

biosample=$1
strain=$2

kingfisher get -r $biosample -m ena-ascp ena-ftp aws-http prefetch

if [ -f ${biosample}_1.fastq ]; then
        #mv ${biosample}_1.fastq ${strain}_R1.fastq
        #mv ${biosample}_2.fastq ${strain}_R2.fastq
        gzip ${biosample}_R1.fastq
        gzip ${biosample}_R2.fastq
fi

mv ${biosample}_1.fastq.gz ${strain}_R1.fastq.gz
mv ${biosample}_2.fastq.gz ${strain}_R2.fastq.gz

# Make sure to make it executable:
chmod +x individual_kf.sh

# make another bash script that will handle the individual_kf.sh for all strains. The name of this script is run_kf_handle.sh:

run_kf_handle.sh
​​#!/bin/bash

while read line; do
        biosample=`echo -e "$line" | cut -f 1 -d '      '`
        strain=`echo -e "$line" | cut -f 2 -d ' '`
        echo "bash individual_kf.sh $biosample $strain"
        bash individual_kf.sh $biosample $strain
done < sra_list.tab

# Make sure to make it executable:
chmod +x run_kf_handle.sh

# To run this script, just do:

./run_kf_handle.sh

# organize the downloaded file in dedicated directory. So let’s organize our current project directory:

mkdir raw
mv *gz raw

cd raw

mkdir illumina nanopore
mv *illumina.fastq.gz illumina
mv *nanopore.fastq.gz nanopore

ls

## De Novo Genome Assembly from Short-Read Data
# Software Installation
conda install bioconda::bbmap -y
conda install bioconda::fastqc -y
conda install bioconda::fastp -y
conda install bioconda::spades -y 
conda install bioconda::bandage -y

conda install bioconda::bakta -y
conda install bioconda::prokka -y
conda install bioconda::kraken2 -y
conda install bioconda::antismash -y
conda install bioconda::ncbi-amrfinderplus -y
conda install bioconda::abricate -y

# check the quality of the Illumina reads first
# If you are just trying to test out the tools without taking long time to run them, you can resample the reads of your fastq files using reformat.sh from bbmap. The following commands will randomly sample 10% reads from both of your fastq files (forward and reverse):
reformat.sh in1=strain1_R1.fq.gz in2=strain1_R2.fq.gz out1=strain1_sampled_R1.fq.gz out2=strain1_sampled_R2.fq.gz samplerate=0.1

# determining taxonomy: identifying and characterizing the samples based on their k-mer profiles.
#!/bin/bash
#unset PYTHONPATH

for infile in `ls -1 ./raw/*.fastq.gz`; do
        dataset=`basename "$infile" | sed 's/.fastq.gz//g'`
        echo -e "$infile $dataset"
        sendsketch.sh -da in=$infile out=./sketches/${dataset}.sketch
done

# Once the .sketch files have been generated, we can summarize those using the following script.
#!/bin/bash

# Create the header for summary file
echo 'strainID  ANI     gSize   taxName' > sketch_summary.tsv

# Loop through sketch report and extract necessary column from the first hit
for sketch in `ls -d *.sketch`;
do
        strainid=`basename "$sketch" | sed 's/.sketch//g'`
        ani=`tail -n +4 $sketch | head -n 1 | cut -f 3`
        gSize=`tail -n +4 $sketch | head -n 1 | cut -f 10`
        taxa=`tail -n +4 $sketch | head -n 1 | cut -f 12`
        echo "${strainid}       ${ani}  ${gSize}        ${taxa}" >> sketch_summary.tsv
done

# Quality-Assessment of FastQ: Normally, we want to visualize the quality of the raw-reads before we start. For that, we can use FastQC:
fastqc strain1_R1.fastq.gz
fastqc strain1_R2.fastq.gz

fastqc Agro_i_R1.fastq.gz
fastqc Agro_i_R1.fastq.gz

#Or, using bash wildcard:
fastqc *fastq.gz

# use another program names FastP to pre-process the reads further for assembly. We can remove any reads which has overall quality below some threshold:
# fastp outputs has fastp.fastq.gz suffix.
fastp -i strain1_R1.fastq.gz -I strain1_R2.fastq.gz -o strain1_R1.fastp.fastq.gz -O strain1_R2.fastp.fastq.gz -e 28
fastp -i Agro_i_R1.fastq.gz -I Agro_i_R2.fastq.gz -o Agro_i_R1.fastp.fastq.gz -O Agro_i_R1.fastp.fastq.gz -e 28

# Now, we can repeat FastQC again:Once we are happy with the quality of the reads, we can move to assembly!
fastqc *fastp.fastq.gz

# Sequence Assembly: 
spades.py -o output_directory --phred-offset 33 --careful -1 strain1_R1.fastp.fastq.gz -2 strain1_R2.fastp.fastq.gz -k 21,33,55,77,99
spades.py -o Agro_i_spades_output --phred-offset 33 --careful -1 Agro_i_R1.fastp.fastq.gz -2 Agro_i_R2.fastp.fastq.gz -k 21,33,55,77,99

# But if we have limited resources, we can specify the run with fewer ram, fewer k-mer estimation, etc:
# # There are two main files that we need later: contigs.fasta which contains all the assembled contigs, and assembly_graph.fastg, which contains the graph data that we can visualize.
spades.py -o output_directory --phred-offset 33 -1 strain1_R1.fastp.fastq.gz -2 strain1_R2.fastp.fastq.gz -k 21,33,55 -t 4 -m 10 --only-assembler

# Getting assembly statistics
# There are multiple tools for doing that. A quick and dirty approach is to use an old perl program called assemblathon.
wget https://raw.githubusercontent.com/KorfLab/Assemblathon/refs/heads/master/assemblathon_stats.pl
wget https://raw.githubusercontent.com/ucdavis-bioinformatics/assemblathon2-analysis/refs/heads/master/FAlite.pm

# Now run assemblathon:
# Quast is also good tool, specially if there are any particular reference you want to compare the assembly with. 
perl assemblathon_stats.pl spades/output_directory/contigs.fasta
perl -I. assemblathon_stats.pl Agro_i_spades_output/contigs.fasta > Agro_i_assembly_stats.txt

# Visualizing assembly graph
bandage image assembly_graph.fastg output.png
Bandage image assembly_graph.fastg Agro_graph.png
ls Agro_graph.png

# again use the `sendsketch.sh` to check the assembled contigs individually, and look for any potential contamination. However, we will be using the mode=sequence to search the each contigs individually, and determine its taxonomy:
sendsketch.sh mode=sequence in=spades/output_directory/contigs.fasta out=sketches/contigs.sketc
sendsketch.sh mode=sequence in=contigs.fasta out=contigs.sketch

## long read assembly

# Installation
conda install bioconda::flye -y
conda install bioconda::nanoplot -y 
conda install matplotlib seaborn -y
conda install bioconda::dnaapler -y 

# Consolidating FastQ files
#!/bin/bash
for folder in `ls -d barcode*`; do
        newfile="${folder}.fastq.gz"
        cat $folder/* > $newfile
        echo "$folder $newfile"
done

# Checking taxonomy
sendsketch.sh -da in=strain1.fastq.gz out=./sketches/strain1.sketch

# Quality of the Nanopore Reads
NanoPlot -t 2 -o dataset/strain1 -p strain1.hex --loglength -f png --plots hex --title strain1 --fastq strain1.fastq.gz

# Our command
NanoPlot -t 12 -o dataset/SRR32486128 -p SRR32486128.hex --loglength -f png --plots hex --title SRR32486128 --fastq nanopore/SRR32486128.fastq.gz --legacy hex
# Warning
WARNING: hex as part of --plots has been deprecated and will be ignored. To get the hex output, rerun with --legacy hex.
# image
imgcat SRR24759618.hexWeightedLogTransformed_HistogramReadlength.png

# generated NanoPlot output 
#!/bin/bash
echo -e "dataset        number of reads total bases     mean length     median length   mean quality read N50   longest read"
for infile in `ls -1 ../raw/*gz`; do
     	   dataset=`echo -e "$infile" | sed 's/.fastq.gz//g;s/^.*\///g'`
     	   n50=`grep "N50" ./$dataset/${dataset}.hexNanoStats.txt | sed 's/^.*\s//g'`
     	   meanlength=`grep "Mean read length:" ./$dataset/${dataset}.hexNanoStats.txt | sed 's/^.*\s//g'`
     	   meanquality=`grep "Mean read quality:" ./$dataset/${dataset}.hexNanoStats.txt | sed 's/^.*\s//g'`
     	   medianlength=`grep "Median read length:" ./$dataset/${dataset}.hexNanoStats.txt | sed 's/^.*\s//g'`
     	   numreads=`grep "Number of reads:" ./$dataset/${dataset}.hexNanoStats.txt | sed 's/^.*\s//g'`
     	   totalbases=`grep "Total bases:" ./$dataset/${dataset}.hexNanoStats.txt | sed 's/^.*\s//g'`
     	   longestread=`grep -A 1 "Top 5 longest reads" ./$dataset/${dataset}.hexNanoStats.txt | tail -n 1 | sed 's/^1:\s//g'`
     	   echo -e "$dataset       $numreads       $totalbases     $meanlength     $medianlength   $meanquality    $n50    $longestread"
Done

# Assembly and Rotation
flye --iterations 5 -t $threads -o flye/flyehq_strain1 --scaffold --read-error 0.03 --nano-hq strain1.fastq.gz

# Explantion
You run 5 passes to carefully piece it together (--iterations 5)
You work with 12 assistants to speed it up (-t 12)
You save the reconstructed book in a folder (-o flye/flyehq_SRR32486128)
You also try to order the chapters and infer missing words (--scaffold)
You assume each page scan has a 3% chance of error (--read-error 0.03)
Your input is a ZIP file of page scans (--nano-hq SRR32486128.fastq.gz)

# use dnaapler tool to rotate circular contigs.
#find non-circular contigs and linear chromosomes and put in ignore list
#rotate the rest

echo "Make ignore list" 
touch flyehq_${dataset}/rotate_ignore_list

awk '$4 == "N"' flyehq_strain1/assembly_info.txt | cut -f 1 > flyehq_strain1/rotate_ignore_list

echo "Ignore list created" 

dnaapler all --ignore flyehq_strain1/rotate_ignore_list --input flyehq_strain1/assembly.fasta --output flyehq_strain1/rotated --prefix strain1 --threads 8 -f

# Assembly Statistics
assemblysize=`grep 'Total length:' flye/strain1/flye.log | sed 's/\tTotal length:\t//g'`
coverage=`grep 'Mean coverage:' flye/strain1/flye.log | sed 's/\tMean coverage:\t//g'`
fragments=`grep 'Fragments:' flye/strain1/flye.log | sed 's/\tFragments:\t//g'`
N50=`grep 'Fragments N50:' flye/strain1/flye.log | sed 's/\tFragments N50:\t//g'`

echo "$strain   $assemblysize   $coverage       $fragments      $N50"

# Visualize the Contigs
bandage image assembly_graph.fastg output.png

# Install and activate software
conda install -c bioconda sra-tools -y
conda install -c bioconda fastqc -y
conda install -c bioconda multiqc -y
conda install -c bioconda trimmomatic -y
conda install -c bioconda cutadapt -y
conda install -c bioconda seqkit -y
conda install -c bioconda samtools -y
conda install -c bioconda gatk -y
conda install -c bioconda freebayes -y
conda install -c bioconda -c conda-forge delly -y

conda create -n manta_env python=2.7
conda activate manta_env
conda install -c bioconda manta -y
conda install -c bioconda samtools bedtools
wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

conda install -c bioconda cnvkit -y

conda install -c bioconda annovar -y
conda install -c bioconda ensembl-vep -y

conda install -c bioconda r-clusterprofiler

# Create fresh conda evironemnt for GATk
conda create -n gatk_env
conda activate gatk_env
conda install -c bioconda -c conda-forge gatk4

# optional
conda install -c bioconda gatk4 -y # If you want latest version, sometimes not provide latest version 
conda install -c bioconda gatk=4.3.0.0 -y # If you want version specific

# To add read group in sorte.bam file need this software
conda install -c conda-forge openjdk=17 -y
conda install -c bioconda picard -y

# Or in one line code
conda create -n gatk_env -c bioconda -c conda-forge gatk4 openjdk=8 -y
conda activate gatk_env
gatk --version

# Download and extract the SRA files by SRAToolkit
prefetch SRR5384356

ls *.sra
ls -lh *.sra
mv SRR5384356 SRR
mv SRR5384356.sra ..
rm -v *.sra
rm -r */
tar xjvf # for unzip

# Split the pair-end files by fatsq-dump
fastq-dump --split-files SRR5384356.sra

# Quality check of fastq files by FastQC
fastqc -t 14 -o fastqc_results SRR5384356_1.fastq SRR5384356_2.fastq 

# Download fastqc and multiqc reports                               
scp <user>@<REMOTE_HOST>:~/wgs/multiqc_report.html C:/Users/<user>/Downloads/

# Trimming with Trimmomatic
# Download the TruSeq3-PE-2.fa adapter 
wget https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE-2.fa

ls -lh TruSeq3-PE-2.fa

trimmomatic PE -phred33 -threads 12 \
  SRR5384356_1.fastq SRR5384356_2.fastq \
  SRR5384356_1_paired.fastq SRR5384356_1_unpaired.fastq \
  SRR5384356_2_paired.fastq SRR5384356_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Check the trimmed files
# Check size
ls -lh SRR5384356_1_paired.fastq
ls -lh SRR5384356_1_unpaired.fastq
ls -lh SRR5384356_2_paired.fastq
ls -lh SRR5384356_2_unpaired.fastq
# Check output files
head SRR5384356_1_paired.fastq
head SRR5384356_1_unpaired.fastq
head SRR5384356_2_paired.fastq
head SRR5384356_2_unpaired.fastq

# Check the trimmed fastq or combined multiqc files again
fastqc SRR5384356_1_paired.fastq SRR5384356_1_unpaired.fastq SRR5384356_2_paired.fastq SRR5384356_2_unpaired.fastq -o fastqc_output/
multiqc fastqc_output/

# View the FastQC reports for individual files (you can open them one by one in a browser)
firefox fastqc_output/SRR5384356_1_paired_fastqc.html
firefox fastqc_output/SRR5384356_2_paired_fastqc.htm
firefox multiqc_report.html

# Download Fastqc or MultiQC html file in local Server
pscp <user>@<REMOTE_HOST>:/omyk/ 5_PRJNA1000995/trimmed_reads/ SRR25482873_2_val_2_fastqc.html C:\Users\a\Downloads\

# Bowtie2 download
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.4/
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.4/bowtie2-2.5.4-linux-x86_64.zip

### Alignment with Reference Genome
# Download reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz

# Unzip the fna and gff files 
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.gff.gz

# Rename the fna and gff files
mv GCF_000005845.2_ASM584v2_genomic.fna ref.ecoli.fna
mv GCF_000005845.2_ASM584v2_genomic.gff ref.ecoli.gff

# Build the Reference Genome Index (Bowtie2 index building)
bowtie2-build ref.ecoli.fna reference_index

# This will create following files 
reference_index.1.bt2
reference_index.2.bt2
reference_index.3.bt2
reference_index.4.bt2
reference_index.rev.1.bt2
reference_index.rev.2.bt2

# Check them
ls -lh reference_index.*

# Alignment to Reference genome by Bowtie2
bowtie2 -x reference_index -1 SRR5384356_1_paired.fastq -2 SRR5384356_2_paired.fastq -S SRR5384356.sam

# If want to download sam file to check in genomw browser
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356_aligned.sam C:/Users/<user>/Downloads/
# SAM to BAM Conversion and Sorting by Samtool
samtools view -Sb SRR5384356.sam > SRR5384356.bam
# Sort bam file
samtools sort SRR5384356.bam -o SRR5384356_sorted.bam
# Index bam file
samtools index SRR5384356_sorted.bam
# Download bam file to visualize in IGV
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356.bam C:\Users\a\Downloads\

# Remove PCR duplicates
samtools rmdup SRR5384356_sorted.bam SRR5384356_sorted_rmdup.bam
# Idenx dmdup.bam files
samtools index SRR5384356_sorted_rmdup.bam

# Index genome file with samtool
samtools faidx ref.ecoli.fasta

# Index genome file for GATK
gatk CreateSequenceDictionary \
  -R ref.ecoli.fna

# To add read group in sorte.bam file need this software
picard AddOrReplaceReadGroups \
  I=SRR5384356_sorted_rmdup.bam \
  O=SRR5384356_sorted_rmdup_with_rg.bam \
  RGID=1 \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU=unit1 \
  RGSM=SRR5384356
Or
# Add read groups
gatk AddOrReplaceReadGroups \
  -I SRR5384356_sorted.bam \
  -O SRR5384356_rg.bam \
  -RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM sample

# Index the rg.bam file 
samtools index SRR5384356_sorted_rmdup_with_rg.bam

# Variant calling with GATK
gatk HaplotypeCaller \
  -R ref.ecoli.fna \
  -I SRR5384356_sorted_rmdup_with_rg.bam \
  -O SRR5384356.vcf.gz

# Index the bam
samtools index SRR5384356_dedup.bam
# Mark duplicates
gatk MarkDuplicates \
  -I SRR5384356_rg.bam \
  -O SRR5384356_dedup.bam \
  -M dedup_metrics.txt

# Variant calling with bcftools/ Call SNPs
# Need samtoll indexing of fna and sorted bam file
samtools index SRR5384356_sorted_rmdup_with_rg.bam # this will generate a .bai
samtools faidx ref.ecoli.fna # this will create ref.ecoli.fna.fai

bcftools mpileup -f ref.ecoli.fna SRR5384356_sorted_rmdup_with_rg.bam | bcftools call -mv -Oz -o SRR5384356_variants.vcf.gz
bcftools mpileup -Ov -f my.fasta my-sorted.bam | bcftools call -mv -Ov -v variants.vcf
bcftools stats 

# Filter vcf file/Variant filtering
# Filtration of SNP
bcftools view -v snps SRR5384356.vcf >snp.vcf
bcftools filter -i 'QUAL>30 && DP>10' SRR5384356_variants.vcf.gz -o SRR5384356_filtered.vcf.gz
# Filter quality socre
bcftools filter -i 'QUAL>=20' snp.vcf -o snp_qual.vcf
# Filter based on Depth Coverage DP
bcftools filter -i 'DP>=20' snp_qual.vcf -o snp_qual_dp.vcf

# Freebayes for variant calling
/home/<user>/Desktop/SNP_Samtools/freebayes/build/freebayes -f Ecoli.fasta -p1 Ecoli_sorted.bam > variants_freebayes.vcf

# Annontate vcf files
conda install -c bioconda snpeff -y # For annoatation install snpeff
snpEff ann -v ecoli SRR5384356_filtered.vcf.gz > SRR5384356_filtered_annotated.vcf

# Snpff
sudo apt install openjdk-21-jdk
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
java -jar ~/wgs/snpEff/snpEff.jar
java -jar ~/wgs/snpEff/SnpSift.jar

mkdir -p ./Softwares/snpEff/data/K12

find ~ -name "snpEff.jar"
sudo find / -name "snpEff.jar"
sudo find / -type f -name "snpEff.jar" 2>/dev/null
snpEff -h
snpsift -h

# Download SnpEff Databases

# Download cufflink
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
sudo apt install cufflinks
gffread -F -y protein.fa -g sequences.fa genes.gff
gffread --version

K12.genome : Ecoli genome (K12)

java -jar snpEff.jar build -gff3 -v K12
java -jar ./Softwares/snpEff/snpEff.jar -v K12 snp_qual_dp.vcf >annotated.vcf

java -jar ./snpEff/snpEff -jar -v K12 snp_qual_dp.vcf >annotated.vcf
scp <user>@<REMOTE_HOST>:/home/<user>/wgs/snpEff/snpEff_summary.html C:\Users\a\Downloads\

scp <user>@<REMOTE_HOST>:/home/<user>/wgs/snpEff_summary.html C:\Users\a\Downloads\

bcftools stats snp_qual_dp.vcf >stat.txt

bcftools mpileup -Ov -f Ecoli.fasta Ecoli_sorted.bam | bcftools call -mv -Ov -o variants.vcf

# Manta Install new env
# Ensure the BAM is indexed
samtools index SRR5384356_sorted.bam

# Run
configManta.py \
--bam SRR5384356_sorted.bam \
--referenceFasta ref.ecoli.fna \
--runDir manta_sortedbam_run

# Then
cd manta_sortedbam_run
./runWorkflow.py -m local -j 4

# output
diploidSV.vcf.gz

# Check output file
ls results/variants/
zgrep -v "^#" diploidSV.vcf.gz | less

# Download
scp <user>@<REMOTE_HOST>:~/wgs/manta_sortedbam_run/results/variants/candidateSV.vcf.gz C:\Users\a\Downloads\

# Index sorted BAM file
samtools index SRR5384356_sorted.bam

# Run DELLY for structural variant calling
# Deletion file check
delly call -t DEL -g ref.ecoli.fna -o SRR5384356.del.delly.bcf SRR5384356_sorted.bam

# Duplication file check
delly call -t DUP -g ref.ecoli.fna -o SRR5384356.dup.delly.bcf SRR5384356_sorted.bam

# Convert BCF to VCF
bcftools view SRR5384356.del.delly.bcf > SRR5384356.del.delly.vcf
bcftools view SRR5384356.dup.delly.bcf > SRR5384356.dup.delly.vcf

delly filter -f germline -o SRR5384356.delly.vcf -s sample.txt SRR5384356.delly.vcf
# Merge your existing variant calls (e.g., DEL and DUP) into one common site list:
delly merge -o sites.cnv.bcf SRR5384356.del.delly.bcf SRR5384356.dup.delly.bcf

# Then, you can use this merged site file to genotype SVs in the same or other samples
delly call -g ref.ecoli.fna -v sites.cnv.bcf -o SRR5384356.cnv.bcf SRR5384356_sorted.bam

# Download files for IGV
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356_sorted.bam "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356_sorted.bam.bai "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/ref.ecoli.fna "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/ref.ecoli.fna.fai "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356.cnv.vcf "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356.del.delly.vcf "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356.dup.delly.vcf "C:\Users\a\Downloads\"
scp <user>@<REMOTE_HOST>:~/wgs/SRR5384356_freebayes.vcf "C:\Users\a\Downloads\"

# Additional reference commands (compilation / troubleshooting notes)
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.4/bowtie2-2.5.4-linux-x86_64.zip/download
https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2
mkdir build && cd build
# Ensure that CC and CXX are updated to target compiler if needed, e.g.:
#     export CC=/path/to/cc
#     export CXX=/path/to/c++
../manta-${MANTA_VERSION}.release_src/configure --jobs=4 --prefix=/path/to/install
make -j4 install
apt-get update -qq

apt-get install -qq bzip2 gcc g++ make python python-software-properties zlib1g-dev

apt-get update -qq && apt-get install -qq bzip2 gcc g++ make python3 python3-pip software-properties-common zlib1g-dev

apt-get install -qq python3 python3-pip
apt-get install -qq software-properties-common

E: Could not open lock file /var/lib/dpkg/lock-frontend - open (13: Permission denied)
E: Unable to acquire the dpkg frontend lock (/var/lib/dpkg/lock-frontend), are you root?

sudo rm /var/lib/dpkg/lock-frontend
sudo rm /var/cache/apt/archives/lock

sudo apt update please

sudo dpkg --configure -a
sudo apt-get update

sudo rm /var/lib/dpkg/lock-frontend
sudo rm /var/cache/apt/archives/lock

wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.release_src.tar.bz2
tar -xjf

mkdir build && cd build

export CC=/opt/rh/devtoolset-3/root/usr/bin/gcc
export CXX=/opt/rh/devtoolset-3/root/usr/bin/g++

bowtie2-build Ecoli.fasta index

bowtie2-align-s -x index -1 SRR._R1.fastq -2 SRR_R2.fastq -S Ecoli.sam

/home/<user>/WGS/softwares/sratoolkit.3.2.1-ubuntu64/bin/fastq-dump --split-files SRR5384356

./softwares/bowtie2-2.5.4-linux-x86_64/bowtie2-build ecoli.fasta index

./softwares/bowtie2-2.5.4-linux-x86_64/bowtie2-align-s -x index -1 SRR5384356_1_paired.fastq -2 SRR5384356_2_paired.fastq -S ecoli.sam

samtools view -S -b ecoli.sam > ecoli.bam
samtools sort ecoli.bam -o  ecoli_sorted.bam

samtools faidx ecoli.fasta 

bcftools mpileup -Ov -f ecoli.fasta ecoli_sorted.bam | bcftools call -mv -Ov -o

sudo apt-get install samtools

bcftools

https://github.com/dellytools/delly.git
git clone --recursive

sudo apt-get install libhts-dev libcurl4-openssl-dev
sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev

sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev

libcurl4-openssl-dev : Conflicts: libcurl4-gnutls-dev but 8.5.0-2ubuntu10.6 is to be installed

sudo apt-get remove libcurl4-gnutls-dev
sudo apt-get install libcurl4-openssl-dev

if [ -r src/htslib/Makefile ]; then cd src/htslib && autoreconf -i && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && make && make lib-static && cd ../../ && touch .htslib; fi
/bin/sh: 1: autoreconf: not found
make: *** [Makefile:57: .htslib] Error 127