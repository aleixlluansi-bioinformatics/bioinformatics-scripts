#!/bin/bash
#SBATCH --job-name=bowtie2_nohuman
#SBATCH --output=logs/bowtie2_nohuman_%j.out
#SBATCH --error=logs/bowtie2_nohuman_%j.err
#SBATCH --qos=class_a
#SBATCH --account=res23021
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G

# Cargar módulo
module load bowtie2/2.3.2

# Crear carpetas si no existen
mkdir -p host_removed logs

# Ruta al índice de Bowtie2
INDEX="bowtie2_index/hg19"

# Loop sobre muestras
for fq1 in clean_data/*_1_val_1.fq.gz; do
  sample=$(basename "$fq1" _1_val_1.fq.gz)
  fq2="clean_data/${sample}_2_val_2.fq.gz"

  if [[ -f "$fq2" ]]; then
    echo "✅ Procesando muestra: $sample"

    bowtie2 -x "$INDEX" \
      --sensitive-local \
      -1 "$fq1" -2 "$fq2" \
      --un-conc "host_removed/${sample}_nohuman.fq" \
      -p 8 \
      -S /dev/null \
      2> "logs/${sample}_bowtie2.log"
  else
    echo "⚠️  No se encontró el archivo pareado para: $sample" >> logs/missing_pairs.log
  fi
done

