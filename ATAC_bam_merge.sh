## this script is to filter bam files based on consensus peaks 
## script was run on HPC and made by Huifang Yuan
## this output file will be used for tracks visualization in Gviz and parallel with 5 of 8 consensus peaks

## step 1: filter bam files based on consensus peaks
## step 2:  merge alll filtered bam files into one bam file using smatools merge


#!/bin/bash

for i in /QRISdata/Q2677/07.ATAC-2025/00.my.analysis/03.mapping/03.last.bam/01.nat_light/02.constant.last/*.bam
do 
  echo "Processing $i"
  outfile="${i%.bam}.filtered.bam"
  samtools view -@ 4 -L 01.light.merge.sort.withID.bed -b "$i" > "$outfile"
done

################
nohup samtools merge 002.constant_all.filtered.bam /QRISdata/Q2677/07.ATAC-2025/00.my.analysis/03.mapping/03.last.bam/01.nat_light/02.constant.last/*.filtered.bam & 


