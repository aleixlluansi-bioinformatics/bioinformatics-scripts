#!/bin/bash
#SBATCH --job-name=read_stats
#SBATCH --output=logs/read_stats_%j.out
#SBATCH --error=logs/read_stats_%j.err
#SBATCH --account=res23021
#SBATCH --qos=class_a
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

cd /gpfs/archive/resdata/res23021/metagenomica/prueba_aleix

output_file="read_stats_summary.tsv"
echo -e "#sample\tfile\tbases\treads\tmin_len\tmedian_len\tmean_len\tmax_len" > "$output_file"

calculate_stats() {
    local fq=$1
    local sample=$2
    local file_label=$3

    local lengths_file=$(mktemp)
    local tmp_stats=$(mktemp)

    # Inicializar min_length para gawk
    local min_length=-1

    gawk -v min_length="$min_length" -v lengths_file="$lengths_file" '
    BEGIN {
        total_bases=0; total_reads=0; max_length=0;
    }
NR%4==2 {
        len=length($0);
        total_bases+=len;
        total_reads++;
        if (min_length == -1 || len < min_length) min_length=len;
        if (len > max_length) max_length=len;
        print len >> lengths_file
    }
    END {
        if (total_reads > 0) {
            mean_len=total_bases/total_reads;
        } else {
            mean_len=0; min_length=0; max_length=0;
        }
        print total_bases, total_reads, min_length, mean_len, max_length
    }' "$fq" > "$tmp_stats"

    read bases reads min mean max < "$tmp_stats"

    median=$(sort -n "$lengths_file" | awk '{
        a[NR]=$1
    }
    END {
        if (NR % 2 == 0) {
            print (a[NR/2] + a[NR/2+1]) / 2
        } else {
            print a[(NR+1)/2]
        }
    }')

    echo -e "${sample}\t${file_label}\t${bases}\t${reads}\t${min}\t${median}\t${mean}\t${max}" >> "$output_file"

    rm "$lengths_file" "$tmp_stats"
}

# Procesar todos los pares en host_removed/
for fq1 in host_removed/*_nohuman.1.fq; do
    fq2="${fq1%.1.fq}.2.fq"
    sample=$(basename "$fq1" _nohuman.1.fq)

    if [[ -f "$fq2" ]]; then
        echo "📊 Procesando $sample"
        calculate_stats "$fq1" "$sample" "R1"
        calculate_stats "$fq2" "$sample" "R2"
    else
        echo "⚠️  Pares faltantes para $sample"
    fi
done

echo "✅ Estadísticas generadas en: $output_file"

