#!/bin/bash

#Create and activate ambient
#mamba create -y -n SeqsAnalyses cutadapt fastqc multiqc bwa samtools igv qualimap flash2 parallel fastx_toolkit
#mamba activate SeqsAnalyses

#echo "getting current directory"
WD=$(pwd)
#File in which I'm going to work
file_name=$WD/data/30-1-02-EL3_R1_all.fastq.gz

#Quality report
mkdir -p $WD/reports
mkdir -p $WD/logs/fastqc
mkdir -p $WD/logs/multiqc
mkdir -p $WD/reports/multiqc

find $WD/data -name '*R1*' | sed 's/.fastq.gz$//' | parallel -j $(nproc) \
	"fastqc $WD/data/{/.}.fastq.gz -o \
	$WD/reports 2>&1 | tee \
	$WD/logs/fastqc/{/.}.log"

find $WD/data -name '*R2*' | sed 's/.fastq.gz$//' | parallel -j $(nproc) \
	"fastqc $WD/data/{/.}.fastq.gz -o \
	$WD/reports 2>&1 | tee \
	$WD/logs/fastqc/{/.}.log"
	
multiqc "$WD/reports/" -f -m fastqc -o \
	"$WD/reports/multiqc/" 2>&1 | tee \
	"$WD/logs/multiqc/Report.log"

#Remove barcodes, primers and low quality bases
#Auxiliar directories
mkdir -p $WD/aux/cutadapt/R1_results $WD/aux/cutadapt/too_short/Reg1 \
       $WD/aux/cutadapt/untrimmed/Reg1  $WD/logs/cutadapt/Reg1

PR1=GAATGTTGGTGACATACTTGCT
PR2=TTGCCATGATCAAATTGACC
LB=10


#This line uses cutadapt to remove primers, barcodes and low quality ends. 
#Just works for R1 reads, need to integrate the R2 reads.
#Requiers file_name
cutadapt -q 25 -g ^$PR1 -O 13 -o $WD/aux/cutadapt/R1_results/$file_name.fastq.gz --too-short-output \
	$WD/aux/cutadapt/too_short/Reg1/sample_less_m50.fastq --untrimmed-output $WD/aux/cutadapt/untrimmed/Reg1/sample_untrimmed.fastq \
	$file_name
	
