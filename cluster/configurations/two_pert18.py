import numpy as np
from pyoculus.fields.slab_bfield_two_pert import SlabBfieldTwoPert
from pyoculus.maps import ToroidalBfieldSection
from pyoculus.solvers import PoincarePlot, FixedPoint, Manifold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import json
from pathlib import Path

# Configuration number
config_id = 18

# Results directory, independent of the working directory
cluster_dir = Path(__file__).resolve().parent.parent
results_dir = cluster_dir / "manifold_results" / f"config_{config_id}"
results_dir.mkdir(parents=True, exist_ok=True)

field_params = {
    "m1": 1,
    "n1": 1,
    "delta1": 0.5,
    "m2": 1,
    "n2": -1,
    "delta2": 0.5,
    "iota_prime": 3.0,
    "x0": 1.0
}

with open(results_dir / "parameters.json", "w") as f:
    json.dump(field_params, f, indent=4)

print(f"Parameters saved to {results_dir / 'parameters.json'}")

field = SlabBfieldTwoPert(**field_params)

xs1 = field.x0 + field.n1 / (field.iota_prime * field.m1)
xs2 = field.x0 + field.n2 / (field.iota_prime * field.m2)

map = ToroidalBfieldSection(field, domain=[(0, 2), (0, 2*np.pi)])

pplot = PoincarePlot.with_linspace(
    map,
    np.array([0.05, np.pi]),
    np.array([1.95, np.pi]),
    npts=50
)

# Follow each initial condition for 200 crossings
pplot.compute(npts=150)

fp2 = FixedPoint(map)
fp2.find_with_iota(n = field.n1, m = field.m1, guess = np.array([xs1, 0.01]), method="scipy.root")   

fp1 = FixedPoint(map)
fp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array(fp2.coords[1]), method="scipy.root")  

fp3 = FixedPoint(map)
fp3.find_with_iota(n = field.n2, m = field.m2, guess = np.array([xs2, 0.1]), method="scipy.root")  

fp4 = FixedPoint(map)
fp4.find(t=1, guess = np.array([0.684, 6.15]), method="scipy.root")   

ofp1 = FixedPoint(map)
ofp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([xs1, 3]), method="scipy.root")

ofp2 = FixedPoint(map)
ofp2.find_with_iota(n = field.n2, m = field.m2, guess = np.array([xs2, 3]), method="scipy.root")

# Manifolds:

manifold1 = Manifold(
    map, fp1, fp2,
    dir1=ofp1.coords[0] - fp1.coords[0],
    dir2=ofp1.coords[0] - fp2.coords[0]
)
manifold2 = Manifold(
    map, fp1, fp2,
    dir1=ofp1.coords[0] - fp1.coords[0],
    dir2=ofp1.coords[0] - fp2.coords[0],
    first_stable=False
)
manifold3 = Manifold(
    map, fp3, fp4,
    dir1=ofp2.coords[0] - fp3.coords[0],
    dir2=ofp2.coords[0] - fp4.coords[0]
)
manifold4 = Manifold(
    map, fp3, fp4,
    dir1=ofp2.coords[0] - fp3.coords[0],
    dir2=ofp2.coords[0] - fp4.coords[0],
    first_stable=False
)

for i, manifold in enumerate([manifold1, manifold2, manifold3, manifold4], start=1):
    manifold.vector_s = np.real_if_close(manifold.vector_s)
    manifold.vector_u = np.real_if_close(manifold.vector_u)


# Test manifold configurations:# Directory for this configuration's results

manifold1.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.3e-3)
stable1, unstable1 = manifold1.compute_adaptive2(tol=3e-3, n_u_extra=0, n_s_extra=0)
    
manifold2.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.3e-3)
stable2, unstable2 = manifold2.compute_adaptive2(tol=3e-3, n_u_extra=0, n_s_extra=0)

manifold3.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.3e-3)
stable3, unstable3 = manifold3.compute_adaptive2(tol=3e-3, n_u_extra=0, n_s_extra=0)

manifold4.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.3e-3)
stable4, unstable4 = manifold4.compute_adaptive2(tol=3e-3, n_u_extra=0, n_s_extra=0)

stable_unwrapped1 = stable1.copy()
unstable_unwrapped1 = unstable1.copy()
stable_unwrapped1[:, 1] = np.unwrap(stable1[:, 1], period=2*np.pi)
unstable_unwrapped1[:, 1] = np.unwrap(unstable1[:, 1], period=2*np.pi)

stable_unwrapped2 = stable2.copy()
unstable_unwrapped2 = unstable2.copy()
stable_unwrapped2[:, 1] = np.unwrap(stable2[:, 1], period=2*np.pi)
unstable_unwrapped2[:, 1] = np.unwrap(unstable2[:, 1], period=2*np.pi)

stable_unwrapped3 = stable3.copy()
unstable_unwrapped3 = unstable3.copy()
stable_unwrapped3[:, 1] = np.unwrap(stable3[:, 1], period=2*np.pi)
unstable_unwrapped3[:, 1] = np.unwrap(unstable3[:, 1], period=2*np.pi)

stable_unwrapped4 = stable4.copy()
unstable_unwrapped4 = unstable4.copy()
stable_unwrapped4[:, 1] = np.unwrap(stable4[:, 1], period=2*np.pi)
unstable_unwrapped4[:, 1] = np.unwrap(unstable4[:, 1], period=2*np.pi)

# Save manifold arrays
np.savez_compressed(
    results_dir / "manifolds.npz",
    stable1=stable1,
    unstable1=unstable1,
    stable2=stable2,
    unstable2=unstable2,
    stable3=stable3,
    unstable3=unstable3,
    stable4=stable4,
    unstable4=unstable4,
    stable_unwrapped1=stable_unwrapped1,
    unstable_unwrapped1=unstable_unwrapped1,
    stable_unwrapped2=stable_unwrapped2,
    unstable_unwrapped2=unstable_unwrapped2,
    stable_unwrapped3=stable_unwrapped3,
    unstable_unwrapped3=unstable_unwrapped3,
    stable_unwrapped4=stable_unwrapped4,
    unstable_unwrapped4=unstable_unwrapped4
)

print(f"Manifold results saved to {results_dir / 'manifolds.npz'}")


fig, ax = pplot.plot(marker=".", s=0.5, xlim=[.0, 2.], plottype=None)

ax.plot(stable_unwrapped1[:, 0], stable_unwrapped1[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped1[:, 0], unstable_unwrapped1[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped2[:, 0], stable_unwrapped2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped2[:, 0], unstable_unwrapped2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped4[:, 0], stable_unwrapped4[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped4[:, 0], unstable_unwrapped4[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped3[:, 0], stable_unwrapped3[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped3[:, 0], unstable_unwrapped3[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.legend()

# Save figure
fig.savefig(results_dir / "manifolds.png", dpi=300, bbox_inches="tight")

print(f"Figure saved to {results_dir / 'manifolds.png'}")

plt.close(fig)