#!/bin/bash
#SBATCH --job-name=amrfinder_db
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=amrfinder_db_%j.out
#SBATCH --error=amrfinder_db_%j.err

set -euo pipefail

# Ajusta el nombre del módulo si corresponde
#module load singularity

cd "$SLURM_SUBMIT_DIR"

IMAGE="$PWD/containers/amrfinder_4.2.7_db2026-08-07.1.sif"

if [[ ! -f "$IMAGE" ]]; then
    echo "ERROR: no se encuentra el contenedor: $IMAGE" >&2
    exit 1
fi

mkdir -p amrfinder_db

singularity exec --cleanenv \
    --bind "$PWD:/work" \
    "$IMAGE" \
    amrfinder_update -d /work/amrfinder_db

# Verificar las versiones del programa y la base descargada
singularity exec --cleanenv \
    --bind "$PWD:/work" \
    "$IMAGE" \
    amrfinder -d /work/amrfinder_db/latest -V
