
#################################
# Pipeline Shotgun Metagenomica #
#################################



##### PREPROCESSING ####



#!/bin/bash


#1. Install required packages: fastqc, muiltiqc, bowtie2, samtools, trim_galore

conda create -n preprocess_env -c bioconda trim-galore

conda activate preprocess_env

conda install -c bioconda bowtie2

 pip install multiqc

conda install bioconda::fastqc

conda install bioconda::samtools



#2. Create directories for project data

mkdir -p project/{raw_data,clean_data,host_removed,qc_reports,bowtie2_index,logs}


#3. Run fastqc and multiqc for quality control check

fastqc raw_data/*.fastq -o qc_reports/ -t 8 
multiqc qc_reports/ -o qc_reports/  #Outputs will be stored at qc_reports folder
                                               
                                               
#4. Cut with trimGalore
 
mkdir -p clean_data logs

for fq1 in raw_data/*_1.fastq; do   sample=$(basename "$fq1" _1.fastq);   fq2="raw_data/${sample}_2.fastq";   if [[ -f "$fq2" ]]; then     trim_galore       --paired "$fq1" "$fq2"       --quality 20       --length 75       --max_n 2       -o clean_data/       > logs/${sample}_trim_galore.log 2>&1;   else     echo "No paired file found for sample $sample";   fi; done  #Outputs will be stored at qc_reports folder

less clean_data/ERR11594211_1.fastq_trimming_report.txt  #Open report



#5. Fastqc and multiqc on trimmed data

mkdir -p qc_reports/trimmed
fastqc clean_data/*_val_*.fq -o qc_reports/trimmed/ -t 8
multiqc qc_reports/trimmed/ -o qc_reports/trimmed/ #Outputs will be stored at qc_reports/trimmed folder


#6. Download and build Bowtie2 index if already hasn't

mkdir -p bowtie2_index
cd bowtie2_index

# Descarga índice preconstruido (Bowtie2)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# Construir índice
bowtie2-build hg19.fa hg19


#6. Contaminant and host DNA was identification with Bowtie2 (v2.3.4.3)28 using the–sensitive–local parameter, allowing confident removal of the Instrument's spike-ins and human-associated reads (hg19 human genome release)

for fq1 in clean_data/*_1_val_1.fq; do
  sample=$(basename "$fq1" _1_val_1.fq)
  fq2="clean_data/${sample}_2_val_2.fq"

  if [[ -f "$fq2" ]]; then
    echo "✅ Procesando muestra: $sample"

    bowtie2 -x bowtie2_index/hg19 \
      --sensitive-local \
      -1 "$fq1" -2 "$fq2" \
      --un-conc "host_removed/${sample}_nohuman.fq" \
      -S /dev/null \
      2> logs/${sample}_bowtie2.log

  else
    echo "⚠️  No se encontró el archivo pareado para: $sample"
  fi
done




#7. Read statistics of all samples (number of reads, number of bases, minimum, median, and maximum read length per sample)
# Use read_stats.sh stored at project folder

chmod +x read_stats.sh
./read_stats.sh


#8. Filter Metagenomes with ≥5 M reads were included in the analysis
# Use filter_samples.sh stored at "project" folder

chmod +x filter_samples.sh
./filter_samples.sh




####### Species-level profiling ######

#We will use MetaPhlAn 4.0.329 with default parameters and the latest database version from http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/. Current latest version: mpa_vJan25_CHOCOPhlAnSGB_202503.tar

#Use Specieslevel_profiling.sh stored at "project" folder 

chmod +x Specieslevel_profiling.sh
./Specieslevel_profiling.sh


#merge profile.txt of each sample (output from metaphlan)
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt


###### End of preprocessing steps #####


#metaphlan has some extra steps for visualizing taxonomic profiles. See github metaphlan website.



                                                                                  
