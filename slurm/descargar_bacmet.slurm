#!/bin/bash
#SBATCH --job-name=bacmet_db
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=bacmet_db_%j.out
#SBATCH --error=bacmet_db_%j.err

set -euo pipefail

DB="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre/analysis/bacmet_db"
URL="http://bacmet.biomedicine.gu.se/download"

mkdir -p "$DB"
cd "$DB"

# Descargar las proteínas y sus anotaciones
for FILE in \
    BacMet2_EXP_database.fasta \
    BacMet2_EXP.753.mapping.txt
do
    curl --fail --location \
        --connect-timeout 30 \
        --max-time 300 \
        --retry 3 \
        "$URL/$FILE" \
        -o "${FILE}.part"

    test -s "${FILE}.part"
    mv "${FILE}.part" "$FILE"
done

# Verificar que el archivo descargado tiene formato FASTA
if ! head -n 1 BacMet2_EXP_database.fasta | grep -q '^>'; then
    echo "ERROR: la descarga no parece un FASTA." >&2
    exit 1
fi

echo "Número de secuencias:"
grep -c '^>' BacMet2_EXP_database.fasta

echo "Primeras líneas de las anotaciones:"
head -n 3 BacMet2_EXP.753.mapping.txt

sha256sum \
    BacMet2_EXP_database.fasta \
    BacMet2_EXP.753.mapping.txt > checksums.sha256

date -u > fecha_descarga.txt

echo "Descarga completada en $DB"

