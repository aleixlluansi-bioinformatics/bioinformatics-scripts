#!/bin/bash
#SBATCH --job-name=filter_samples
#SBATCH --output=logs/filter_samples_%j.out
#SBATCH --error=logs/filter_samples_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=res23021
#SBATCH --qos=class_a

# Cambiar al directorio donde están los archivos
cd /gpfs/archive/resdata/res23021/metagenomica/prueba_aleix

input="read_stats_summary.tsv"
output="high_quality_samples.tsv"
min_reads=5000000

echo -e "#sample\ttotal_reads" > "$output"

# Obtener lista de muestras únicas
samples=$(awk 'NR>1 {print $1}' "$input" | sort | uniq)

for sample in $samples; do
  r1_reads=$(awk -v s="$sample" '$1==s && $2=="R1" {print $4}' "$input")
  r2_reads=$(awk -v s="$sample" '$1==s && $2=="R2" {print $4}' "$input")

  # Validar que sean números enteros
  if [[ "$r1_reads" =~ ^[0-9]+$ && "$r2_reads" =~ ^[0-9]+$ ]]; then
    total=$((r1_reads + r2_reads))

    # Comparación: mínimo de lecturas por par, no por R1 o R2 por separado
    if (( total >= min_reads * 2 )); then
      echo -e "$sample\t$total" >> "$output"
    fi
  fi
done

echo "✅ Muestras filtradas guardadas en: $output"
