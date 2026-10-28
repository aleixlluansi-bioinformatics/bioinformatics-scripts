#!/bin/bash
#SBATCH --job-name=fastqc_test
#SBATCH --output=logs/fastqc_test_%j.out
#SBATCH --error=logs/fastqc_test_%j.err
#SBATCH --qos=class_a
#SBATCH --account=res23021
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8G

#3 muestras ~15min

module load java fastqc/0.11.9

cd /gpfs/projects/res23021/metagenomica/neuroblastoma_heces/

mkdir -p output_fastqc

fastqc /gpfs/home/sjd/sjd685439/raw_data/*.fastq.gz -o output_fastqc -t 4 #ajustar segun numero muestras

