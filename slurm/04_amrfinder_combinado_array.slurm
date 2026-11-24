#!/bin/bash
#SBATCH --job-name=amrfinder_plus
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/amrfinder_%A_%a.out
#SBATCH --error=logs/amrfinder_%A_%a.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"

IMAGE="$BASE/analysis/containers/amrfinder_4.2.7_db2026-08-07.1.sif"
DB="$BASE/analysis/amrfinder_db/2026-08-07.1"
INPUT="$BASE/secuencias_extraidas"
OUTPUT="$BASE/analysis/resultados_amrfinder_combinado_plus"

# Lista de muestras recibida al enviar el trabajo
LIST="${1:?Falta la ruta de la lista de muestras}"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[[ -n "$SAMPLE" ]] || {
    echo "ERROR: muestra vacía para tarea $SLURM_ARRAY_TASK_ID" >&2
    exit 1
}

FNA="$INPUT/fna/${SAMPLE}.fna"
FAA="$INPUT/faa/${SAMPLE}.faa"
GFF="$INPUT/gff/${SAMPLE}.gff"

for FILE in "$IMAGE" "$FNA" "$FAA" "$GFF"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: archivo inexistente o vacío: $FILE" >&2
        exit 1
    }
done

[[ -d "$DB" ]] || {
    echo "ERROR: no existe la base: $DB" >&2
    exit 1
}

mkdir -p "$OUTPUT"

# Separar cada ejecución para conservar los resultados anteriores
RUN="$OUTPUT/$SLURM_ARRAY_JOB_ID"
mkdir -p "$RUN"

echo "Muestra: $SAMPLE"
echo "Nodo: $(hostname)"
echo "Inicio: $(date)"

singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    amrfinder \
        -n "$FNA" \
        -p "$FAA" \
        -g "$GFF" \
        -d "$DB" \
        --plus \
        --threads "$SLURM_CPUS_PER_TASK" \
        --name "$SAMPLE" \
        -o "$RUN/${SAMPLE}.amrfinder.tsv.tmp" \
    2> "$RUN/${SAMPLE}.amrfinder.log"

# Publicar la tabla solo si el programa terminó correctamente
mv "$RUN/${SAMPLE}.amrfinder.tsv.tmp" \
   "$RUN/${SAMPLE}.amrfinder.tsv"

echo "Finalizado: $(date)"
