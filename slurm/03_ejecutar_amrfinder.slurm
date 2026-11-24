#!/bin/bash
#SBATCH --job-name=amrfinder
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=amrfinder_%j.out
#SBATCH --error=amrfinder_%j.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"

IMAGE="$BASE/analysis/containers/amrfinder_4.2.7_db2026-08-07.1.sif"
DB="$BASE/analysis/amrfinder_db/2026-08-07.1"
INPUT="$BASE/secuencias_extraidas/fna"
OUTPUT="$BASE/analysis/resultados_amrfinder"

[[ -f "$IMAGE" ]] || { echo "No existe: $IMAGE" >&2; exit 1; }
[[ -d "$DB" ]] || { echo "No existe: $DB" >&2; exit 1; }

shopt -s nullglob
FILES=("$INPUT"/*.fna)

if (( ${#FILES[@]} == 0 )); then
    echo "No se encontraron archivos .fna en $INPUT" >&2
    exit 1
fi

mkdir -p "$OUTPUT"

# Registrar las versiones utilizadas
singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    amrfinder -d "$DB" -V \
    > "$OUTPUT/versiones_${SLURM_JOB_ID}.txt" 2>&1

for FNA in "${FILES[@]}"; do
    SAMPLE=$(basename "$FNA" .fna)

    echo "Iniciando: $SAMPLE — $(date)"

    singularity exec --cleanenv \
        --bind "$BASE:$BASE" \
        "$IMAGE" \
        amrfinder \
            -n "$FNA" \
            -d "$DB" \
            --threads "$SLURM_CPUS_PER_TASK" \
            --name "$SAMPLE" \
            -o "$OUTPUT/${SAMPLE}.amrfinder.tsv" \
        2> "$OUTPUT/${SAMPLE}.amrfinder.log"

    echo "Finalizado: $SAMPLE — $(date)"
done

echo "Análisis completado: ${#FILES[@]} genomas"


