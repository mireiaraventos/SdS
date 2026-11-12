import numpy as np
from pyoculus.fields.slab_bfield_two_pert import SlabBfieldTwoPert
from pyoculus.maps import ToroidalBfieldSection
from pyoculus.solvers import PoincarePlot, FixedPoint, Manifold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from time import perf_counter

import json
from pathlib import Path

# Configuration number
config_id = 23

# Results directory, independent of the working directory
cluster_dir = Path(__file__).resolve().parent.parent
results_dir = cluster_dir / "manifold_results" / f"config_{config_id}"
results_dir.mkdir(parents=True, exist_ok=True)

field_params = {
    "m1": 2,
    "n1": 1,
    "delta1": 0.5,
    "m2": 1,
    "n2": -1,
    "delta2": 0.5,
    "iota_prime": 4.0,
    "x0": 1.0
}

with open(results_dir / "parameters.json", "w") as f:
    json.dump(field_params, f, indent=4)

print(f"Parameters saved to {results_dir / 'parameters.json'}")

field = SlabBfieldTwoPert(**field_params)

xs1 = field.x0 + field.n1 / (field.iota_prime * field.m1)
xs2 = field.x0 + field.n2 / (field.iota_prime * field.m2)

map = ToroidalBfieldSection(field, domain=[(0,2), (0, 2*np.pi)], type="dopri5", rtol=1e-10, nsteps=100000)
pplot = PoincarePlot.with_linspace(
    map,
    np.array([0.05, np.pi]),
    np.array([1.95, np.pi]),
    npts=50
)

# Follow each initial condition for 200 crossings
pplot.compute(npts=150)

# Fixed points:
# X-points island 1:
fp1 = FixedPoint(map)
fp1.find(t = 2, guess = np.array([1.25, 0.]), method="scipy.root")  

fp2 = FixedPoint(map)
fp2.find(t = 2, guess = np.array(fp1.coords[1]), method="scipy.root")
# fp2.find(t = 2, guess = np.array([1.21, np.pi]), method="scipy.root")

fp3 = FixedPoint(map)
fp3.find(t = 2, guess = np.array([1.287, 6.28]), method="scipy.root")

# X-points island 2:
fp4 = FixedPoint(map)
fp4.find(t = 1, guess = np.array([0.7, 0.]), method="scipy.root")

fp5 = FixedPoint(map)
fp5.find(t = 1, guess = np.array(fp4.coords[1]), method="scipy.root")  

# O-points island 1:
ofp1 = FixedPoint(map)
ofp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([1.3, np.pi]), method="scipy.root") 

ofp2 = FixedPoint(map)
ofp2.find_with_iota(n = field.n1, m = field.m1, guess = np.array(ofp1.coords[1]), method="scipy.root")

# O-points island 2:
ofp3 = FixedPoint(map)
ofp3.find_with_iota(n = field.n2, m = field.m2, guess = np.array([xs2, 3]), method="scipy.root")  


# Manifolds:
# Island 1 bottom:
manifold1 = Manifold(
    map, fp2, fp1,
    dir1=ofp1.coords[0] - fp2.coords[0],
    dir2=ofp1.coords[0] - fp1.coords[0]
)

manifold2 = Manifold(
    map, fp2, fp1,
    dir1=ofp1.coords[0] - fp2.coords[0],
    dir2=ofp1.coords[0] - fp1.coords[0],
    first_stable=False
)

# Island 1 top:
manifold3 = Manifold(
    map, fp2, fp3,
    dir1=ofp2.coords[0] - fp2.coords[0],
    dir2=ofp2.coords[0] - fp3.coords[0]
)

manifold4 = Manifold(
    map, fp2, fp3,
    dir1=ofp2.coords[0] - fp2.coords[0],
    dir2=ofp2.coords[0] - fp3.coords[0],
    first_stable=False
)


# Island 2:
manifold5 = Manifold(
    map, fp4, fp5,
    dir1=ofp3.coords[0] - fp4.coords[0],
    dir2=ofp3.coords[0] - fp5.coords[0]
)

manifold6 = Manifold(
    map, fp4, fp5,
    dir1=ofp3.coords[0] - fp4.coords[0],
    dir2=ofp3.coords[0] - fp5.coords[0],
    first_stable=False
)

for i, manifold in enumerate([manifold1, manifold2, manifold3, manifold4, manifold5, manifold6], start=1):
    manifold.vector_s = np.real_if_close(manifold.vector_s)
    manifold.vector_u = np.real_if_close(manifold.vector_u)


# Test manifold configurations:# Directory for this configuration's results

manifolds = [
    manifold1, manifold2, manifold3,
    manifold4, manifold5, manifold6
]

errors = [4e-3, 4e-3, 4e-3, 4e-3, 4e-3, 4e-3]

manifold_arrays = []
unwrapped_arrays = []
save_data = {}

for i, (manifold, err) in enumerate(zip(manifolds, errors), start=1):
    t0 = perf_counter()

    print(f"Starting manifold {i}", flush=True)

    manifold.find_clinics_exhaustive(
        first_guess_eps_s=1e-6,
        first_guess_eps_u=1e-6,
        ERR=err
    )

    print(f"Clinic search: {(perf_counter() - t0) / 60:.1f} min", flush=True)
    t1 = perf_counter()

    stable, unstable = manifold.compute_adaptive2(tol=4e-3, max_periods=20, max_points=1000, n_u_extra=2, n_s_extra=2)

    print(f"Adaptive computation: {(perf_counter() - t1) / 60:.1f} min", flush=True)

    stable_u = stable.copy()
    unstable_u = unstable.copy()

    stable_u[:, 1] = np.unwrap(stable[:, 1], period=2*np.pi)
    unstable_u[:, 1] = np.unwrap(unstable[:, 1], period=2*np.pi)

    manifold_arrays.append((stable, unstable))
    unwrapped_arrays.append((stable_u, unstable_u))

    save_data[f"stable{i}"] = stable
    save_data[f"unstable{i}"] = unstable
    save_data[f"stable_unwrapped{i}"] = stable_u
    save_data[f"unstable_unwrapped{i}"] = unstable_u

    # Save after every completed manifold
    np.savez_compressed(
        results_dir / "manifolds.npz",
        **save_data
    )

    print(f"Manifold {i} completed and saved", flush=True)


# Plot Poincare section and manifolds
fig, ax = pplot.plot(
    marker=".",
    s=0.5,
    xlim=[0.0, 2.0],
    plottype=None
)
  
for i, (stable, unstable) in enumerate(unwrapped_arrays, start=1):

    ax.plot(
        stable[:, 0], stable[:, 1],
        color="blue",
        label="Stable" if i == 1 else None,
        marker=".",
        markersize=3,
        linestyle=""
    )

    ax.plot(
        unstable[:, 0], unstable[:, 1],
        color="red",
        label="Unstable" if i == 1 else None,
        marker=".",
        markersize=3,
        linestyle=""
    )

ax.legend()

fig.savefig(results_dir / "manifolds.png", dpi=300, bbox_inches="tight")

plt.close(fig)