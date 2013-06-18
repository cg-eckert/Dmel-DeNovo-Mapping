#!/bin/bash
cd /
#mkdir /data
#mount -o ro /dev/xvdf /data

#The de novo assembly is in the file /data/dmel_trinity/Trinity.fasta

#The read fastq files are in /data/OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq
#OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq
#OREf_SAMm_w_GTCCGC_L006_R1_001.fastq
#OREf_SAMm_w_GTCCGC_L006_R2_001.fastq

#Start with some read QC
#mkdir /mnt/QC
#/usr/local/share/FastQC/fastqc /data/*.fastq --outdir=/mnt/QC

#All other files will be written to /mnt/map
mkdir /mnt/map

#Interleave the paired end reads in preparation for fastx trimming
cd /data
python /usr/local/share/khmer/sandbox/interleave.py OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq > /mnt/map/vg.combined.fastq
python /usr/local/share/khmer/sandbox/interleave.py OREf_SAMm_w_GTCCGC_L006_R1_001.fastq OREf_SAMm_w_GTCCGC_L006_R2_001.fastq > /mnt/map/w.combined.fastq

cd /mnt/map
#Use the FASTX toolkit to trim off bases over 70 (replace with new number)
fastx_trimmer -Q33 -l 70 -i vg.combined.fastq > vg.trimmed.fq
fastx_trimmer -Q33 -l 70 -i w.combined.fastq > w.trimmed.fq

#Split the paired reads
python /usr/local/share/khmer/sandbox/split-pe.py vg.trimmed.fq
python /usr/local/share/khmer/sandbox/split-pe.py w.trimmed.fq

echo "**********QC DONE************"

# Index the de novo transcriptome
bowtie-build /data/dmel_trinity/Trinity.fasta /mnt/map/denovo_bowtie

# Map with bowtie
echo "******** MAPPING WITH BOWTIE*********"
cd /mnt/map
bowtie denovo_bowtie vg.trimmed.fq.1 vg.1.sam
bowtie denovo_bowtie vg.trimmed.fq.2 vg.2.sam
bowtie denovo_bowtie w.trimmed.fq.1 w.1.sam
bowtie denovo_bowtie w.trimmed.fq.2 w.2.sam

# Convert SAM files to indexed BAM files
echo "******converting SAM to BAM********"
samtools view -Sb vg.1.sam > vg.1.temp.bam
samtools view -Sb vg.2.sam > vg.2.temp.bam
samtools view -Sb w.1.sam > w.1.temp.bam
samtools view -Sb w.2.sam > w.2.temp.bam
samtools sort -f vg.1.temp.bam vg.1.bam
samtools sort -f vg.2.temp.bam vg.2.bam
samtools sort -f w.1.temp.bam w.1.bam
samtools sort -f w.2.temp.bam w.2.bam
samtools index vg.1.bam
samtools index vg.2.bam
samtools index w.1.bam
samtools index w.2.bam
rm *.temp.bam
echo "********* SAMs have been BAMed *************"


# Make sure bedtools is installed
# Using bedtools to calculate read counts
echo "********bedtools analysis starting*********"
cd /mnt/map
coverageBed -s -abam vg.1.bam -b /data/dmel_trinity/final_dmel.bed > vg1.counts.txt
coverageBed -s -abam vg.2.bam -b /data/dmel_trinity/final_dmel.bed > vg2.counts.txt
coverageBed -s -abam w.1.bam -b /data/dmel_trinity/final_dmel.bed > w1.counts.txt
coverageBed -s -abam w.2.bam -b /data/dmel_trinity/final_dmel.bed > w2.counts.txt

echo "********bedtools analysis FINISHED********"
















