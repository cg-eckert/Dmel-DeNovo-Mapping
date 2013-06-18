#!/bin/bash
cd /
#mkdir /data
#mount -o ro /dev/xvdf /data

#The de novo assembly is in the file /data/dmel_trinity/Trinity.fasta

#The read fastq files are in /data/OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq
#OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq
#OREf_SAMm_w_GTCCGC_L006_R1_001.fastq
#OREf_SAMm_w_GTCCGC_L006_R2_001.fastq

#mkdir /mnt/map

#Make little read files with 10,000 reads each and write to /mnt/map. Note that we use the same seed for the paired reads so that the same reads are sampled
cd /data
seqtk sample -s 99 OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq 10000 > /mnt/map/vg1_10k_1.fastq
seqtk sample -s 99 OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq 10000 > /mnt/map/vg1_10k_2.fastq
seqtk sample -s 53 OREf_SAMm_w_GTCCGC_L006_R1_001.fastq 10000 > /mnt/map/w_10k_1.fastq
seqtk sample -s 53 OREf_SAMm_w_GTCCGC_L006_R2_001.fastq 10000 > /mnt/map/w_10k_2.fastq

cd /mnt
#Start with some read QC
mkdir QC
#/usr/local/share/FastQC/fastqc /data/*.fastq --outdir=/mnt/QC
/usr/local/share/FastQC/fastqc /mnt/map/*.fastq --outdir=/mnt/QC

#All other files will be written to /mnt/map

#Interleave the paired end reads in preparation for fastx trimming
#cd /data
#python /usr/local/share/khmer/sandbox/interleave.py OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq > /mnt/map/vg.combined.fastq
#python /usr/local/share/khmer/sandbox/interleave.py OREf_SAMm_w_GTCCGC_L006_R1_001.fastq OREf_SAMm_w_GTCCGC_L006_R2_001.fastq > /mnt/map/w.combined.fastq

#Little read files: Interleave the paired end reads in preparation for fastx trimming
cd /mnt/map
python /usr/local/share/khmer/sandbox/interleave.py vg1_10k_1.fastq vg1_10k_2.fastq > vg.combined.fastq
python /usr/local/share/khmer/sandbox/interleave.py w_10k_1.fastq w_10k_1.fastq > w.combined.fastq


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
bowtie -S -p 2 denovo_bowtie -1 vg.trimmed.fq.1 -2 vg.trimmed.fq.2 vg.paired.sam
bowtie -S -p 2 denovo_bowtie -1 w.trimmed.fq.1 -2 w.trimmed.fq.2 w.paired.sam

# Convert SAM files to indexed BAM files
echo "******converting SAM to BAM********"
samtools view -bS vg.paired.sam > vg.temp.bam
samtools view -bS w.paired.sam > w.temp.bam

samtools sort -f vg.temp.bam vg.paired.bam
samtools sort -f w.temp.bam w.paired.bam

samtools index vg.paired.bam
samtools index w.paired.bam

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
















