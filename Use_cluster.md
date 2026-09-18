# Run simulations using FaADE in EPFL clusters

This file contains the information needed to run simulations using `FaADE` with EPFL's cluster `jed`.

## 1. Contents/files
### 1.1 Julia project
The files `Manifest.toml`and `Project.toml` contain the Julia environment to use `FaADE` to run the simulations.
### 1.2 One simulation
The file `SI_simulate_unif_grid_cluster.jl` solves the anisotropic diffusion equation for a single island using an uniform grid.

- The perpendicular diffusion coefficient is fixed to $\kappa_\perp=1.0$
- Time integration goes from $t_0=0$ to $t_f=0.01$, with a timestep $\Delta t=10^{-4}$.
- The amplitude of the island $\delta$, the parallel diffusion coefficient $\kappa_\parallel$ and the resolution of the grid $n_x=n_y$ are introduced as arguments.

The results are saved as `num_results/SI_k$(exp_k)_tf$(t_f)_d$(δ)_res$(nx)_dT.jld2` (in the folder `num_results/`, with `exp_k`=$log_{10}(\kappa_\parallel/\kappa_\perp)$, `t_f`=0.01, `δ` is the amplitude and `nx` the grid resolution). What is saved:
- `soln`: the whole solution
- `x_mid`: a vector of the x coordinates in the middle of the island ($y=\pi$)
- `u_mid`: the temperature profile with respect to x in the middle of the island ($y=\pi$)
- `δ`
- `t_f`


To execute it, run 
```
julia --project=. SI_simulate_unif_grid_cluster.jl <delta> <k_para> <nx>
```
For further details on how to execute it, see section 2.2.


### 1.3 Multiple simulations
The file `run_sweep.slurm` contains the code to run multiple simulations in parallel in the cluster. It executes `SI_simulate_unif_grid_cluster.jl` for all the amplitudes specified in `DELTAS` and the parallel diffusion coefficients in `KPARAS`, always with a 201x201 grid. To change grid resolution for all sweep simulations, modify `NX` directly inside `run_sweep.slurm`.

The first lines of the file define how the cluster should schedule and execute the batch of simulations:
- `--array=1-40%40`: Defines the total number of simulations and how many run simultaneously. For large number of simulations, recommended to execute them "in two runs" (ex: for 88 simulations, `--array=1-88%44`, so the first 44 simulations will run in parallel and then the other 44). Update 40 to: `(number of DELTAS) × (number of KPARAS)`
- `--time=05:00:00`: Maximum allowed runtime per simulation. Increase if jobs are too long.
- `--cpus-per-task=2`: CPU cores per simulation. Increase only if needed.
- `--mem=8G`: Memory per simulation. Increase if jobs fail due to memory limits.
- Output and errors are saved in:
    - `logs/*.out` → normal output
    - `logs/*.err` → errors/debugging







## 2. Instructions on how to run simulations using the cluster
### 2.1 Connecting EPFL JED Cluster to VS Code (Remote Development Setup)

#### 2.1.1. Prerequisites

Before starting, make sure you have:
- EPFL HPC account with JED access
- Working SSH login from terminal (PowerShell for Windows):
    ```
    ssh username@jed.hpc.epfl.ch
    ```

If this works in PowerShell/Terminal, VS Code should also work.


#### 2.1.2. Remote SSH in VS Code
1. Install the extension `Remote - SSH` in VS Code.
2. Open the command palette by pressing `Ctrl + Shift + P`
3. Select `Remote-SSH: Connect to Host...`
4. Enter `username@jed.hpc.epfl.ch`
5. If it asks “select SSH configuration file”, choose usual Windows default `C:\Users\YOURNAME\.ssh\config`.

    Your config file should contain:
    ``` 
    Host jed.hpc.epfl.ch
        HostName jed.hpc.epfl.ch
        User username 
    ```
6. When VS Code asks “what platform is the remote host?”, choose Linux
7. Enter your EPFL password

Only need to do all of this the first time. The next times, you can directly access your cluster profile by going through points 2., 3., choosing `jed.hpc.epfl.ch` and then 7.

Once inside your cluster profile, go to your desired working directory on the cluster. You are now editing directly on the cluster, any file you edit/save in VS Code is already on JED.

#### 2.1.3. Working with Julia
Each time we want to use Julia we need to load it: 
```
module load gcc/13.2.0 julia/1.10.4
```
(choose desired version).

### 2.2 Running one simulation
See section 1.2 for further details on the file `SI_simulate_unif_grid_cluster.jl`.

Users must have these files inside the main project directory (inside the cluster):
```
WorkingDirectory/
├── Project.toml
├── Manifest.toml
├── SI_simulate_unif_grid_cluster.jl
├── run_sweep.slurm
└── num_results/
```

1. Activate the project environment - This tells Julia to use the exact environment defined by these files. 
    ```
    julia --project=.
    ```
    Important: Always run Julia commands with `--project=.` to ensure the correct environment is used.
2. Install exact package versions: inside julia run
    ```
    using Pkg
    Pkg.instantiate()
    ```
    This step is usually only required the first time (or after changing `Project.toml` / `Manifest.toml`).
3. Verify environment (optional but recommended):
    ```
    Pkg.status()
    Base.active_project()
    ```
    We'd expect `FaADE` to point to the patched GitHub fork `https://github.com/mireiaraventos/FaADE.jl.git`.

4. Exit Julia and run a test simulation (with low resolution, so it doesn't take too long):
    ```
    julia --project=. SI_simulate_unif_grid_cluster.jl 0.01 1e5 21
    ```
    where the arguments are 
    ```
    δ       = 0.01
    k_para  = 1e5
    nx = ny = 21
    ```
    and the results are saved in `num_results/`.

### 2.3 Running multiple simulations with Slurm
For parameter sweeps, use `run_sweep.slurm`. See section 1.3 for further details on this file.

Submit:
```bash
sbatch run_sweep.slurm
```
Useful commands:
- Check running jobs: `squeue -u $USER`
- Check specific job: `squeue -j JOBID`
- Cancel all your jobs: `scancel -u $USER`


### 2.4 Downloading results to local computer

From local PowerShell (Windows):
```bash
scp "username@jed.hpc.epfl.ch:/home/username/WorkingDirectory/num_results/*.jld2" "C:\Users\YOURNAME\DesiredFolder\"
```
This downloads all `.jld2` result files from the cluster.