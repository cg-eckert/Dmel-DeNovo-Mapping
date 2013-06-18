#!/bin/bash
cd /
#mkdir /data
#mount -o ro /dev/xvdf /data

#The de novo assembly is in the file /data/dmel_trinity/Trinity.fasta

#The read fastq files are in /data
DATA1=/data/OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq
DATA2=/data/OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq
DATA3=/data/OREf_SAMm_w_GTCCGC_L006_R1_001.fastq
DATA4=/data/OREf_SAMm_w_GTCCGC_L006_R2_001.fastq

#Start with some read QC
#mkdir /mnt/QC
#/usr/local/share/FastQC/fastqc /data/*.fastq --outdir=/mnt/QC

#All other files will be written to /mnt/map
mkdir /mnt/map

#Interleave the paired end reads in preparation for fastx trimming
cd /data
python /usr/local/share/khmer/sandbox/interleave.py $DATA1 $DATA2 > /mnt/map/vg.combined.fastq
python /usr/local/share/khmer/sandbox/interleave.py $DATA3 $DATA4 > /mnt/map/w.combined.fastq

cd /mnt/map

for i in vg w ; do
  #Use the FASTX toolkit to trim off bases over 70 (replace with new number)
  fastx_trimmer -Q33 -l 70 -i $i.combined.fastq > $i.trimmed.fq
  #Split the paired reads 
  python /usr/local/share/khmer/sandbox/split-pe.py $i.trimmed.fq
done

echo "**********QC DONE************"

# Index the de novo transcriptome
bowtie-build /data/dmel_trinity/Trinity.fasta /mnt/map/denovo_bowtie

# Map with bowtie
echo "******** MAPPING WITH BOWTIE*********"
cd /mnt/map
bowtie denovo_bowtie -1 vg.trimmed.fq.1 -2 vg.trimmed.fq.2 vg.sam
bowtie denovo_bowtie -1 w.trimmed.fq.1 -2 w.trimmed.fq.2 w1.sam

# Convert SAM files to indexed BAM files
cp /data/dmel_trinity/Trinity.fasta /mnt/map/Trinity.fasta
echo "******converting SAM to BAM********"
samtools faidx Trinity.fasta
for i in vg1 vg2 w1 w2 ; do
  samtools view -Sb $i.sam > $i.temp.bam
  samtools sort -f $i.temp.bam $i.bam
  samtools index $i.bam
done
rm *.temp.bam
echo "********* SAMs have been BAMed *************"

# Make sure bedtools is installed
# Using bedtools to calculate read counts
echo "********bedtools analysis starting*********"
cd /mnt/map
for i in vg1 vg2 w1 w2 ; do
  coverageBed -s -abam $i.bam -b /data/dmel_trinity/final_dmel.bed > $i.counts.txt
done

echo "********bedtools analysis FINISHED********"
















