#!/bin/bash
#SBATCH --job-name=tpiv_2I_sweep
#SBATCH --partition=academic
#SBATCH --qos=academic
#SBATCH --array=1-100%60
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/tpiv_%A_%a.out
#SBATCH --error=logs/tpiv_%A_%a.err

module load gcc/13.2.0 julia/1.10.4

cd ~/SdS/cluster
mkdir -p logs num_results

# Symmetric two-island cases: delta1 = delta2
DELTAS=(0.7 0.6 0.2 0.1)
# DELTAS=(0.5)

# k_para = 10^k_exp
KEXPS=(4.0 6.0 8.0 10.0 12.0)

# Iota values
# IOTAS=(2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8)
IOTAS=(3.0)

NX=201

NUM_DELTAS=${#DELTAS[@]}
NUM_K=${#KEXPS[@]}
NUM_IOTAS=${#IOTAS[@]}

# Convert SLURM task ID (1-based) to 0-based
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

# Determine parameter indices
IOTA_INDEX=$((TASK_ID % NUM_IOTAS))
K_INDEX=$(((TASK_ID / NUM_IOTAS) % NUM_K))
DELTA_INDEX=$((TASK_ID / (NUM_IOTAS * NUM_K)))

DELTA=${DELTAS[$DELTA_INDEX]}
KEXP=${KEXPS[$K_INDEX]}
IOTA=${IOTAS[$IOTA_INDEX]}

# Compute k_para = 10^KEXP
KPARA=$(python3 -c "print(10**float('$KEXP'))")

echo "=============================="
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "delta1 = $DELTA"
echo "delta2 = $DELTA"
echo "k_exp = $KEXP"
echo "k_para = $KPARA"
echo "iota = $IOTA"
echo "nx = $NX"
echo "=============================="

# Run Julia simulation:
# ARGS = delta1 delta2 k_para nx iota
julia --project=. 2I_simulate_unif_grid_cluster_Neumann.jl $DELTA $DELTA $KPARA $NX $IOTA
