#!/bin/bash
#SBATCH --job-name=mobsuite_db
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/mobsuite_db_%j.out
#SBATCH --error=logs/mobsuite_db_%j.err

set -euo pipefail

#module load singularity

BASE="/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre"

IMAGE="$BASE/analysis/containers/mobsuite_3.1.9.sif"
DB="$BASE/analysis/mobsuite_db"

[[ -s "$IMAGE" ]] || {
    echo "ERROR: no existe el contenedor: $IMAGE" >&2
    exit 1
}

mkdir -p "$DB"

echo "Inicio: $(date)"

singularity exec --cleanenv \
    --bind "$BASE:$BASE" \
    "$IMAGE" \
    mob_init -d "$DB"

echo "Inicialización finalizada: $(date)"
du -sh "$DB"

