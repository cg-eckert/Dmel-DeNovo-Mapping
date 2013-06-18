#!/bin/bash
cd /
#mkdir /data
#mount -o ro /dev/xvdf /data

mkdir /mnt/map

#The de novo assembly is in the file /data/dmel_trinity/Trinity.fasta

#The read fastq files are in /data
#DATA1=/data/OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq
#DATA2=/data/OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq
#DATA3=/data/OREf_SAMm_w_GTCCGC_L006_R1_001.fastq
#DATA4=/data/OREf_SAMm_w_GTCCGC_L006_R2_001.fastq

cd /data
seqtk sample -s 11 OREf_SAMm_vg1_CTTGTA_L005_R1_001.fastq 10000 > /mnt/map/vg1.fastq
seqtk sample -s 11 OREf_SAMm_vg1_CTTGTA_L005_R2_001.fastq 10000 > /mnt/map/vg2.fastq
seqtk sample -s 11 OREf_SAMm_w_GTCCGC_L006_R1_001.fastq 10000 > /mnt/map/w1.fastq
seqtk sample -s 11 OREf_SAMm_w_GTCCGC_L006_R2_001.fastq 10000 > /mnt/map/w2.fastq


DATA1=/mnt/map/vg1.fastq
DATA2=/mnt/map/vg2.fastq
DATA3=/mnt/map/w1.fastq
DATA4=/mnt/map/w2.fastq

#Start with some read QC
cd /mnt
rm -r QC
mkdir /mnt/QC
/usr/local/share/FastQC/fastqc /data/*.fastq --outdir=/mnt/QC

#All other files will be written to /mnt/map

mkdir /mnt/map

#Interleave the paired end reads in preparation for fastx trimming

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
bowtie -S -p 2 denovo_bowtie -1 vg.trimmed.fq.1 -2 vg.trimmed.fq.2 vg.sam
bowtie -S -p 2 denovo_bowtie -1 w.trimmed.fq.1 -2 w.trimmed.fq.2 w.sam

# Convert SAM files to indexed BAM files
cp /data/dmel_trinity/Trinity.fasta /mnt/map/Trinity.fasta
echo "******converting SAM to BAM********"
samtools faidx Trinity.fasta
for i in vg w ; do
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
cp /data/dmel_trinity/final_dmel.bed temp_dmel.bed
sed '$d' temp_dmel.bed > final_dmel.bed
bedtools multicov -q 30 -p -bams vg.bam w.bam -bed final_dmel.bed > transcriptome_counts.txt

echo "********bedtools analysis FINISHED********"

R --no-save < Rscript.R

















