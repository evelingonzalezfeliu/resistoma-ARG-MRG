#!/bin/bash
#SBATCH --job-name=integronfinder
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/integronfinder_%A_%a.out
#SBATCH --error=logs/integronfinder_%A_%a.err

set -euo pipefail
#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"
IMAGE="$BASE/analysis/containers/integronfinder_2.0.6.sif"
LIST="${1:?Indica la ruta de lista.txt}"
TASK="${SLURM_ARRAY_TASK_ID:?Envía el trabajo con sbatch --array}"

[[ -s "$LIST" ]] || { echo "Lista inexistente o vacía: $LIST" >&2; exit 1; }
mapfile -t SAMPLES < <(awk '{sub(/\r$/, ""); if (NF) print $1}' "$LIST")
SAMPLE="${SAMPLES[$((TASK - 1))]:?No hay un aislado para esta tarea}"

FNA="$BASE/secuencias_extraidas/fna/$SAMPLE.fna"
OUT="$BASE/analysis/resultados_integronfinder/$SAMPLE"

for FILE in "$IMAGE" "$FNA"; do
    [[ -s "$FILE" ]] || { echo "Falta el archivo o está vacío: $FILE" >&2; exit 1; }
done

mkdir -p "$OUT"
cd "$OUT"

singularity exec --cleanenv --bind "$BASE:$BASE" "$IMAGE" \
    integron_finder \
        --cpu "$SLURM_CPUS_PER_TASK" \
        --linear \
        --outdir "$OUT" \
        "$FNA"

