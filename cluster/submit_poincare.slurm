#!/bin/bash
#SBATCH --job-name=poincare
#SBATCH --partition=academic
#SBATCH --qos=academic
#SBATCH --array=113,114,115
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/poincare_%A_%a.out
#SBATCH --error=logs/poincare_%A_%a.err

# Go to the project directory
cd /home/raventos/SdS/cluster || exit 1

# Activate the Python environment
source .venv/bin/activate

# Make the project directory importable
export PYTHONPATH="$PWD:$PYTHONPATH"

echo "Starting Poincare plot for configuration ${SLURM_ARRAY_TASK_ID}"

python -u compute_poincare.py --config "${SLURM_ARRAY_TASK_ID}"

status=$?

echo "Finished configuration ${SLURM_ARRAY_TASK_ID} with exit code ${status}"

exit "$status"