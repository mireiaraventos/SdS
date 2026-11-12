import json
import numpy as np
from pathlib import Path

from pyoculus.fields.slab_bfield_two_pert import SlabBfieldTwoPert
from pyoculus.maps import ToroidalBfieldSection
from pyoculus.solvers import PoincarePlot

# Paths on the cluster
cluster_dir = Path("/home/raventos/SdS/cluster")
results_dir = cluster_dir / "poincare_plots"
results_dir.mkdir(parents=True, exist_ok=True)

poincare_plots = {}

def compute_poincare(config_id, BC=1, npts_initial=50, npts_crossings=150):

    # Load configuration parameters from JSON
    params_path = (
        cluster_dir / "manifold_results"
        / f"config_{config_id}" / "parameters.json"
    )

    with open(params_path, "r", encoding="utf-8") as f:
        params = json.load(f)

    # Reconstruct field and map
    field_params = params.get("field", params)

    field = SlabBfieldTwoPert(**field_params)

    mymap = ToroidalBfieldSection(
        field,
        domain=[(0, 2), (0, 2*np.pi)]
    )

    # Define initial conditions
    if BC == 1:
        pplot = PoincarePlot.with_linspace(
            mymap,
            np.array([0.05, np.pi]),
            np.array([1.95, np.pi]),
            npts=npts_initial
        )

    elif BC == 2:
        pplot = PoincarePlot.with_segments(
            mymap,
            np.array([
                [0.05, np.pi/2],
                [1.95, np.pi/2],
                [0.05, 3*np.pi/2],
                [1.95, 3*np.pi/2]
            ]),
            neps=[npts_initial, npts_initial],
            connected=False
        )

    else:
        raise ValueError("BC must be 1 or 2.")

    # Compute Poincare trajectories
    pplot.compute(npts=npts_crossings)

    # Store in memory
    poincare_plots[config_id] = pplot

    # Save computed coordinates
    np.savez_compressed(
        results_dir / f"pplot_{config_id}.npz",
        hits=pplot._hits,
        BC=BC,
        npts_initial=npts_initial,
        npts_crossings=npts_crossings
    )

    print(f"Saved Poincare plot for configuration {config_id}")

    return pplot

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=int, required=True)
    args = parser.parse_args()

    config_id = args.config

    # Set the configurations requiring BC=2
    bc2_configs = {21, 22, 23}  # Adjust to your actual configurations

    BC = 2 if config_id in bc2_configs else 1

    compute_poincare(
        config_id=config_id,
        BC=BC,
        npts_initial=50,
        npts_crossings=150
    )