Workflow of the GitHub repository SdS:

```
SdS/
├── cluster/
|   ├── .venv/                ← ignored (Python env)
|   ├── requirements.txt      ← tracked (Python env)
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
The modified PyOculus library is in a separate repository:
```
pyoculus/
└── ...
```

### Work locally 
Should only edit `local/`. The results generated go to `local/num_results_loc`, but notice that this folder isn't loaded into GitHub. Should always pull at the beginning of the session and push at the end.

### Work from the cluster
Should only edit `cluster/`. The results generated go to `cluster/num_results`, but notice that this folder isn't loaded into GitHub. Also, the results generated in the cluster should be loaded into `local/num_results_loc`. Should always start the session by doing a pull and push all the changes at the end.

### Edit forked PyOculus package
Both the local and the cluster `.venv` work with the forked `pyoculus` library. 

Any modifications should be done locally, and then commited and pushed into GH in the `slab-model` branch (important to make sure any changes are made in the `slab-model`branch!!!). 

Then, to use the latest version of the package in the cluster, we need to run:
```
cd ~/pyoculus
git pull
```
