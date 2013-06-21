#!/bin/bash
cd /
#mkdir /data
#mount -o ro /dev/xvdf /data

mkdir /mnt/map

#The de novo assembly is in the file /data/dmel_trinity/Trinity.fasta

#The read fastq.gz files are in /data2. We need to copy them to a directory where we can unzip them
cd /data2
cp OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq.gz /mnt/map/vg1_1.fastq.gz
cp OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq.gz /mnt/map/vg1_2.fastq.gz
cp OREf_SAMm_w_GTCCGC_L006_R1_001.fastq.gz /mnt/map/w1_1.fastq.gz
cp OREf_SAMm_w_GTCCGC_L006_R2_001.fastq.gz /mnt/map/w1_2.fastq.gz
cp SAMf_OREm_vg1_ACTGAT_L004_R1_001.fastq.gz /mnt/map/vg2_1.fastq.gz
cp SAMf_OREm_vg1_ACTGAT_L004_R2_001.fastq.gz /mnt/map/vg2_2.fastq.gz
cp SAMf_OREm_w_CAGATC_L005_R1_001.fastq.gz /mnt/map/w2_1.fastq.gz
cp SAMf_OREm_w_CAGATC_L005_R2_001.fastq.gz /mnt/map/w2_2.fastq.gz

cd /mnt/map
gunzip *.gz

#for i in *.fastq; do seqtk sample -s 11 $i 10000 > ${i/.fastq/S.fastq}; done

VG1_1=/mnt/map/vg1_1.fastq
VG1_2=/mnt/map/vg1_2.fastq
VG2_1=/mnt/map/vg2_1.fastq
VG2_2=/mnt/map/vg2_2.fastq
W1_1=/mnt/map/w1_1.fastq
W1_2=/mnt/map/w1_2.fastq
W2_1=/mnt/map/w2_1.fastq
W2_2=/mnt/map/w2_2.fastq

#VG1_1=vg1_1S.fastq
#VG1_2=vg1_2S.fastq
#VG2_1=vg2_1S.fastq
#VG2_2=vg2_2S.fastq
#W1_1=w1_1S.fastq
#W1_2=w1_2S.fastq
#W2_1=w2_1S.fastq
#W2_2=w2_2S.fastq


#Start with some read QC
cd /mnt
rm -r QC
mkdir /mnt/QC
/usr/local/share/FastQC/fastqc /mnt/map/*.fastq --outdir=/mnt/QC
#/usr/local/share/FastQC/fastqc /mnt/map/*S.fastq --outdir=/mnt/QC

cd /mnt/map
#Interleave the paired end reads in preparation for fastx trimming

python /usr/local/share/khmer/sandbox/interleave.py $VG1_1 $VG1_2 > /mnt/map/vg1.combined.fastq
python /usr/local/share/khmer/sandbox/interleave.py $VG2_1 $VG2_2 > /mnt/map/vg2.combined.fastq
python /usr/local/share/khmer/sandbox/interleave.py $W1_1 $W1_2 > /mnt/map/w1.combined.fastq
python /usr/local/share/khmer/sandbox/interleave.py $W2_1 $W2_2 > /mnt/map/w2.combined.fastq


cd /mnt/map

for i in vg1 vg2 w1 w2 ; do
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
bowtie -S -p 2 denovo_bowtie -1 vg1.trimmed.fq.1 -2 vg1.trimmed.fq.2 vg1.sam
bowtie -S -p 2 denovo_bowtie -1 vg2.trimmed.fq.1 -2 vg2.trimmed.fq.2 vg2.sam
bowtie -S -p 2 denovo_bowtie -1 w1.trimmed.fq.1 -2 w1.trimmed.fq.2 w1.sam
bowtie -S -p 2 denovo_bowtie -1 w2.trimmed.fq.1 -2 w2.trimmed.fq.2 w2.sam

# Convert SAM files to indexed BAM files
cp /data/dmel_trinity/Trinity.fasta /mnt/map/Trinity.fasta
echo "******converting SAM to BAM********"
samtools faidx Trinity.fasta
for i in vg1 vg2 w1 w2
do
  samtools view -Sb $i.sam > $i.temp.bam
  samtools sort $i.temp.bam $i.sorted
  samtools index $i.bam
done

echo "********* SAMs have been BAMed *************"

# Make sure bedtools is installed
# Using bedtools to calculate read counts
echo "********bedtools analysis starting*********"
cd /mnt/map
cp /data/dmel_trinity/final_dmel.bed temp_dmel.bed
sed '$d' temp_dmel.bed > final_dmel.bed
bedtools multicov -q 30 -p -bams vg1.sorted.bam vg2.sorted.bam w1.sorted.bam w2.sorted.bam -bed final_dmel.bed > transcriptome_counts.txt

echo "********bedtools analysis FINISHED********"

#The end product of this pipeline is the file transcriptome_counts.txt
#For thorough data exploration and plotting, it's probably best to analyze the file using the Rscript.R script on your local computer.

















