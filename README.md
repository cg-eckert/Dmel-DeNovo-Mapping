This repository contains the scripts necessary to compare transcript expression between lines of Drosophila melanogaster (vg = something with initials vg vs. w = wild type).
There are two biological replicates of each line (vg1, vg2 and w1, w2).
The reads are paired-end illumina reads.
The reads are mapped back to a transcriptome assembled de novo from these reads.

The first script (xxx) performs quality control on the reads, maps them back to the de novo transcriptome, then counts the reads mapping to each contig.
The second script (Rscript.R) analyses the count data in the file transcriptome_counts.txt. This is designed to work on a local computer for full data exploration.

Good luck!

Chris, Sara, David & Michael

MSU NGS summer course 2013
