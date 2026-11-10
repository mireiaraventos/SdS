#!/bin/bash
#SBATCH --job-name=manifold
#SBATCH --partition=academic
#SBATCH --qos=academic
#SBATCH --array=23
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/manifold_%A_%a.out
#SBATCH --error=logs/manifold_%A_%a.err

cd /home/raventos/SdS/cluster || exit 1

source .venv/bin/activate

# Allow configuration scripts to import test_manifold_config.py
export PYTHONPATH="$PWD:$PYTHONPATH"

echo "Starting configuration ${SLURM_ARRAY_TASK_ID}"

python -u "configurations/two_pert${SLURM_ARRAY_TASK_ID}.py"
status=$?

echo "Finished configuration ${SLURM_ARRAY_TASK_ID} with exit code ${status}"
exit "$status"