#!/bin/bash
#SBATCH --job-name=mefinder
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/mefinder_%A_%a.out
#SBATCH --error=logs/mefinder_%A_%a.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"
IMAGE="$BASE/analysis/containers/mobileelementfinder.sif"
LIST="${1:?Indica la ruta de lista.txt}"
TASK="${SLURM_ARRAY_TASK_ID:?Envía el trabajo con sbatch --array}"

[[ -s "$LIST" ]] || { echo "Lista inexistente o vacía: $LIST" >&2; exit 1; }

# Seleccionar el ID correspondiente a esta tarea; ignorar líneas vacías.
mapfile -t SAMPLES < <(awk '{sub(/\r$/, ""); if (NF) print $1}' "$LIST")
INDEX=$((TASK - 1))
SAMPLE="${SAMPLES[$INDEX]:?No hay un aislado para esta tarea}"

FNA="$BASE/secuencias_extraidas/fna/${SAMPLE}.fna"
OUT="$BASE/analysis/resultados_mobileelementfinder/$SAMPLE"

for FILE in "$IMAGE" "$FNA"; do
    [[ -s "$FILE" ]] || { echo "Falta el archivo o está vacío: $FILE" >&2; exit 1; }
done

mkdir -p "$OUT/tmp"

echo "Aislado: $SAMPLE | Inicio: $(date)"

singularity exec --cleanenv --bind "$BASE:$BASE" "$IMAGE" \
    mefinder find \
    --contig "$FNA" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --temp-dir "$OUT/tmp" \
    --gff \
    "$OUT/$SAMPLE"

echo "Aislado: $SAMPLE | Fin: $(date)"
