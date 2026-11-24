#!/bin/bash
#SBATCH --job-name=isescan
#SBATCH --array=1-100%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/isescan_%A_%a.out
#SBATCH --error=logs/isescan_%A_%a.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"

IMAGE="$BASE/analysis/containers/isescan_1.7.3.sif"
INPUT="$BASE/secuencias_extraidas/fna"
OUTPUT="$BASE/analysis/resultados_isescan"

LIST="${1:?Debes indicar la ruta de lista.txt}"

[[ -s "$LIST" ]] || {
    echo "ERROR: lista inexistente o vacía: $LIST" >&2
    exit 1
}

# Leer IDs e ignorar líneas vacías y retornos de Windows
mapfile -t SAMPLES < <(
    sed 's/\r$//' "$LIST" | awk 'NF {print $1}'
)

if (( ${#SAMPLES[@]} != 100 )); then
    echo "ERROR: se esperaban 100 IDs; hay ${#SAMPLES[@]}" >&2
    exit 1
fi

# Evitar que dos tareas escriban en la misma carpeta
DUPLICATES=$(printf '%s\n' "${SAMPLES[@]}" | sort | uniq -d)
if [[ -n "$DUPLICATES" ]]; then
    echo "ERROR: hay IDs duplicados:" >&2
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

mkdir -p "$SAMPLE_OUTPUT"

# Directorio de trabajo independiente para cada aislado
cd "$SAMPLE_OUTPUT"

echo "Muestra: $SAMPLE"
echo "Nodo: $(hostname)"
echo "Inicio: $(date)"

singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    isescan.py --version

singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    isescan.py \
        --seqfile "$FNA" \
        --output "$SAMPLE_OUTPUT" \
        --nthread "$SLURM_CPUS_PER_TASK"

echo "Finalizado correctamente: $(date)"
