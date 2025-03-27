#!/bin/bash 
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma_possorted_genome_bam.bam
wget https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma_possorted_genome_bam.bam.bai

samtools view -h Parent_SC5v1_Human_Glioblastoma_possorted_genome_bam.bam chr8:102238559-102238874 | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b > RRM2B_Isoforms.bam
samtools index RRM2B_Isoforms.bam

#Remove reads with >7000 bp or <2000 bp distance between the two exons
samtools view RRM2B_Isoforms.bam |grep -Eo "[0-9]+N" |sed "s/N//g" |awk '$1 >2000 && $1 < 7000'|sed 's/$/N/g' >Patterns_to_keep.txt
samtools view -H RRM2B_Isoforms.bam >RRM2B_Isoform1+2.sam
samtools view RRM2B_Isoforms.bam |grep -f Patterns_to_keep.txt >>RRM2B_Isoform1+2.sam
samtools view -bh RRM2B_Isoform1+2.sam > RRM2B_Isoform1+2.bam
samtools index RRM2B_Isoform1+2.bam

#Separate into isoform 1 and 2
samtools view -bh RRM2B_Isoform1+2.bam chr8:102232149-102232304 >RRM2B_Isoform1.bam
samtools index RRM2B_Isoform1.bam

samtools view -bh RRM2B_Isoform1+2.bam chr8:102238559-102238822 >RRM2B_Isoform2.bam
samtools index RRM2B_Isoform2.bam

#Get the cell barcodes from the bam
samtools view Isoform1.bam |grep -oE "CB:Z:.+-1" |sort |uniq >Isoform1.txt

samtools view Isoform2.bam |grep -oE "CB:Z:.+-1" |sort |uniq >Isoform2.txt
