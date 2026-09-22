Workflow of the GitHub repository SdS:

```
SdS/
├── cluster/
│   ├── Project.toml          ← tracked
│   ├── Manifest.toml         ← tracked
│   ├── simulation code       ← tracked
│   ├── Slurm scripts         ← tracked
│   └── num_results/          ← ignored
│
├── local/
│   ├── Project.toml          ← tracked
│   ├── Manifest.toml         ← tracked
│   ├── analysis/code         ← tracked
│   └── num_results_loc/      ← ignored
│
├── literature/               ← ignored
├── .gitignore                ← tracked
└── ...
```