#!/bin/bash
#SBATCH --job-name=diamond_bacmet
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/bacmet_%A_%a.out
#SBATCH --error=logs/bacmet_%A_%a.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"

IMAGE="$BASE/analysis/containers/diamond_2.1.9.sif"
DB="$BASE/analysis/bacmet_db/BacMet2_EXP.dmnd"
INPUT="$BASE/secuencias_extraidas/faa"
OUTPUT="$BASE/analysis/resultados_diamond_bacmet"

LIST="${1:?Debes indicar la ruta de lista.txt}"

[[ -s "$LIST" ]] || {
    echo "ERROR: lista inexistente o vacía: $LIST" >&2
    exit 1
}

mapfile -t SAMPLES < <(
    sed 's/\r$//' "$LIST" | awk 'NF {print $1}'
)

if (( ${#SAMPLES[@]} != 100 )); then
    echo "ERROR: se esperaban 100 IDs; hay ${#SAMPLES[@]}" >&2
    exit 1
fi

SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
FAA="$INPUT/${SAMPLE}.faa"

for FILE in "$IMAGE" "$DB" "$FAA"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: archivo inexistente o vacío: $FILE" >&2
        exit 1
    }
done

RUN="$OUTPUT/$SLURM_ARRAY_JOB_ID"
mkdir -p "$OUTPUT"

echo "Muestra: $SAMPLE"
echo "Inicio: $(date)"

singularity exec --cleanenv \
    "$IMAGE" diamond version \
    > "$OUTPUT/${SAMPLE}.version.txt"

# Búsqueda exploratoria: conservar todas las referencias
# que cumplan el E-value para revisar y filtrar después.
singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    diamond blastp \
        --query "$FAA" \
        --db "$DB" \
        --out "$OUTPUT/${SAMPLE}.bacmet.raw.tsv.tmp" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --very-sensitive \
        --evalue 1e-5 \
        --max-target-seqs 0 \
        --outfmt 6 \
            qseqid sseqid pident length \
            qstart qend sstart send \
            evalue bitscore qlen slen \
    2> "$OUTPUT/${SAMPLE}.diamond.log"

mv "$OUTPUT/${SAMPLE}.bacmet.raw.tsv.tmp" \
   "$OUTPUT/${SAMPLE}.bacmet.raw.tsv"

echo "Coincidencias: $(wc -l < "$OUTPUT/${SAMPLE}.bacmet.raw.tsv")"
echo "Finalizado: $(date)"

