#!/bin/bash
#SBATCH --job-name=trimgalore_all
#SBATCH --output=logs/trimgalore_all_%j.out
#SBATCH --error=logs/trimgalore_all_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=res23021
#SBATCH --qos=class_a

#3 muestras ~1h 

# Cargar módulo y entorno
module load miniforge
source activate trim_galore

# Definir directorios
RAW_DIR="/gpfs/archive/resdata/res23021/metagenomica/prueba_aleix/raw_data"
OUT_DIR="/gpfs/archive/resdata/res23021/metagenomica/prueba_aleix/clean_data"
LOG_DIR="/gpfs/archive/resdata/res23021/metagenomica/prueba_aleix/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Loop para procesar archivos PE
for fq1 in ${RAW_DIR}/*_1.fastq.gz; do
    sample=$(basename "$fq1" _1.fastq.gz)
    fq2="${RAW_DIR}/${sample}_2.fastq.gz"

    if [[ -f "$fq2" ]]; then
        echo "Procesando muestra: $sample"
        trim_galore \
            --paired "$fq1" "$fq2" \
            --quality 20 \
            --length 75 \
            --max_n 2 \
            -j 8 \
            -o "$OUT_DIR" \
            > "${LOG_DIR}/${sample}_trim_galore.log" 2>&1
    else
        echo "No paired file found for sample $sample" >> "${LOG_DIR}/missing_pairs.log"
    fi
done


