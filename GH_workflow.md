Workflow of the GitHub repository SdS:

```
SdS/
├── cluster/
│   ├── Project.toml          ← tracked (Julia env)
│   ├── Manifest.toml         ← tracked (Julia env)
│   ├── simulation code       ← tracked
│   ├── Slurm scripts         ← tracked
│   └── num_results/          ← ignored
│
├── local/
|   ├── .venv/                ← ignored (Python env)
|   ├── requirements.txt      ← tracked (Python env)
│   ├── Project.toml          ← tracked (Julia env)
│   ├── Manifest.toml         ← tracked (Julia env)
│   ├── analysis/code         ← tracked
│   └── num_results_loc/      ← ignored
│
├── literature/               ← ignored
├── .gitignore                ← tracked
└── ...
```

### Work locally 
Should only edit `local/`. The results generated go to `local/num_results_loc`, but notice that this folder isn't loaded into GitHub.

### Work from the cluster
Should only edit `cluster/`. The results generated go to `cluster/num_results`, but notice that this folder isn't loaded into GitHub. Also, the results generated in the cluster should be loaded into `local/num_results_loc`. 