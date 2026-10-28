#!/bin/bash
#SBATCH --job-name=multiqc_test
#SBATCH --output=logs/multiqc_test_%j.out
#SBATCH --error=logs/multiqc_test_%j.err
#SBATCH --account=res23021
#SBATCH --qos=class_a
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G

# Cargar entorno Conda según instrucciones support BSC
module load miniforge
source activate metagen-env

# Entrar al directorio del proyecto
cd /gpfs/archive/resdata/res23021/metagenomica/prueba_aleix/

# Crear carpeta para reporte si no existe
mkdir -p multiqc_reports/trimmed

# Ejecutar MultiQC sobre los resultados de FastQC
multiqc output_fastqc/ -o multiqc_reports/

multiqc multiqc_reports/trimmed/ -o multiqc_reports/trimmed/ 
