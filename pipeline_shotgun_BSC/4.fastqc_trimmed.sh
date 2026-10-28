#!/bin/bash
#SBATCH --job-name=fastqc_test
#SBATCH --output=logs/fastqc_test_%j.out
#SBATCH --error=logs/fastqc_test_%j.err
#SBATCH --qos=class_a
#SBATCH --account=res23021
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8G

module load java fastqc/0.11.9

cd /gpfs/archive/resdata/res23021/metagenomica/prueba_aleix/output_fastqc/trimmed

mkdir -p output_fastqc/trimmed

fastqc clean_data/*_val_*.fq.gz -o output_fastq/trimmed -t 4 #ajustar segun numero muestras
