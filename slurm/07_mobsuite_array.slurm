#!/bin/bash
#SBATCH --job-name=mobsuite
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mobsuite_%A_%a.out
#SBATCH --error=logs/mobsuite_%A_%a.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"

IMAGE="$BASE/analysis/containers/mobsuite_3.1.9.sif"
DB="$BASE/analysis/mobsuite_db"
INPUT="$BASE/secuencias_extraidas/fna"
OUTPUT="$BASE/analysis/resultados_mobsuite"

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

DUPLICATES=$(printf '%s\n' "${SAMPLES[@]}" | sort | uniq -d)
if [[ -n "$DUPLICATES" ]]; then
    echo "ERROR: IDs duplicados:" >&2
    printf '%s\n' "$DUPLICATES" >&2
    exit 1
fi

SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
FNA="$INPUT/${SAMPLE}.fna"
SAMPLE_OUTPUT="$OUTPUT/$SAMPLE"

for FILE in "$IMAGE" "$FNA"; do
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

# Evitar sobrescribir o mezclar ejecuciones
if [[ -e "$SAMPLE_OUTPUT" ]]; then
    echo "ERROR: ya existe la carpeta: $SAMPLE_OUTPUT" >&2
    exit 1
fi

echo "Muestra: $SAMPLE"
echo "Nodo: $(hostname)"
echo "Inicio: $(date)"

singularity exec --cleanenv \
    "$IMAGE" \
    mob_recon --version

singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    mob_recon \
        --infile "$FNA" \
        --outdir "$SAMPLE_OUTPUT" \
        --database_directory "$DB" \
        --num_threads "$SLURM_CPUS_PER_TASK" \
        --sample_id "$SAMPLE"

echo "Finalizado correctamente: $(date)"
