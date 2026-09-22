#!/bin/bash
#SBATCH --job-name=tpiv_sweep
#SBATCH --partition=academic
#SBATCH --qos=academic
#SBATCH --array=1-40%40
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/tpiv_%A_%a.out
#SBATCH --error=logs/tpiv_%A_%a.err


module load gcc/13.2.0 julia/1.10.4

cd ~/TPIV2
mkdir -p logs num_results

DELTAS=(0.007 0.015 0.004 0.00001 0.00005) #(0.005 0.001 0.1 0.05 1.0 0.2 0.0001 0.0002 0.0003 0.0005 0.002 0.003 0.01 0.02 0.03 0.3 0.5)
KPARAS=(1.0e5 1.0e6 1.0e7 1.0e8 1.0e9 1.0e10 1.0e11 1.0e12)

NX=201
NUM_K=${#KPARAS[@]}

TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

DELTA_INDEX=$((TASK_ID / NUM_K))
K_INDEX=$((TASK_ID % NUM_K))

DELTA=${DELTAS[$DELTA_INDEX]}
KPARA=${KPARAS[$K_INDEX]}

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "delta = $DELTA"
echo "k_para = $KPARA"
echo "nx = $NX"

julia --project=. SI_simulate_unif_grid_cluster.jl $DELTA $KPARA $NX
