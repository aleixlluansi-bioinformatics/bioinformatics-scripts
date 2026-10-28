#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --output=logs/metaphlan_%j.out
#SBATCH --error=logs/metaphlan_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=24
#SBATCH --account=res23021
#SBATCH --qos=class_a


#~20 minut per sample with 24 cpus. Alternative could be using more tasks for many samples.


# Cargar entorno
module load miniforge
eval "$(conda shell.bash hook)"
conda activate metagen-env

cd /gpfs/archive/resdata/res23021/metagenomica/prueba_aleix


mkdir -p metaphlan_output

# Ruta a la base de datos de MetaPhlAn
db_dir="mpa_vJan25_CHOCOPhlAnSGB_202503"
index_name="mpa_vJan25_CHOCOPhlAnSGB_202503"

# Crear directorio de salida si no existe
mkdir -p metaphlan_output

# Leer los nombres de muestras de alta calidad
tail -n +2 high_quality_samples.tsv | cut -f1 | while read sample; do
  fq1="host_removed/${sample}_nohuman.1.fq"
  fq2="host_removed/${sample}_nohuman.2.fq"

  if [[ -f "$fq1" && -f "$fq2" ]]; then
    echo "🔄 Procesando muestra: $sample"

    metaphlan \
      "${fq1},${fq2}" \
      --input_type fastq \
      --db_dir "$db_dir" \
      --index "$index_name" \
      --mapout metaphlan_output/${sample}_bowtie2.bz2 \
      --nproc 24 \
      --output_file metaphlan_output/${sample}_profile.txt

  else
    echo "⚠️ Archivo pareado faltante para muestra: $sample"
  fi
done
