import numpy as np
from pyoculus.fields.slab_bfield_two_pert import SlabBfieldTwoPert
from pyoculus.maps import ToroidalBfieldSection
from pyoculus.solvers import PoincarePlot, FixedPoint, Manifold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
import traceback
from test_manifold_config import test_manifold_configuration

field = SlabBfieldTwoPert(
    m1=1,
    n1=1,
    delta1=0.5,
    m2=1,
    n2=-1,
    delta2=0.5,
    iota_prime=2.0,
    x0=1.0,
)

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

# Fixed points:
fp1 = FixedPoint(map)
fp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([1.4, 0.2]), method="scipy.root")   

fp2 = FixedPoint(map)
fp2.find_with_iota(n = field.n1, m = field.m1, guess = np.array([1.5, 2*np.pi]), method="scipy.root")  

fp3 = FixedPoint(map)
fp3.find_with_iota(n = field.n2, m = field.m2, guess = np.array([0.5, 0.2]), method="scipy.root")  

fp4 = FixedPoint(map)
fp4.find_with_iota(n = field.n2, m = field.m2, guess = np.array([0.48, 2*np.pi]), method="scipy.root")  

ofp1 = FixedPoint(map)
ofp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([1.4, 3]), method="scipy.root")

ofp2 = FixedPoint(map)
ofp2.find_with_iota(n = field.n2, m = field.m2, guess = np.array([0.4, 3]), method="scipy.root")


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

# Configuration number
config_id = 12  # Change to 2, 3, etc. in the other scripts

manifold1.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.5e-3)
stable1_2, unstable1_2 = manifold1.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

stable_unwrapped1_2 = stable1_2.copy()
unstable_unwrapped1_2 = unstable1_2.copy()

stable_unwrapped1_2[:, 1] = np.unwrap(stable1_2[:, 1], period=2*np.pi)
unstable_unwrapped1_2[:, 1] = np.unwrap(unstable1_2[:, 1], period=2*np.pi)
    
manifold2.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.5e-3)
stable2_2, unstable2_2 = manifold2.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

manifold3.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.5e-3)
stable3_2, unstable3_2 = manifold3.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

manifold4.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.5e-3, eps_res=20)
stable4_2, unstable4_2 = manifold4.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

stable_unwrapped4_2 = stable4_2.copy()
unstable_unwrapped4_2 = unstable4_2.copy()

stable_unwrapped4_2[:, 1] = np.unwrap(stable4_2[:, 1], period=2*np.pi)
unstable_unwrapped4_2[:, 1] = np.unwrap(unstable4_2[:, 1], period=2*np.pi)


# Results directory, independent of the working directory
cluster_dir = Path(__file__).resolve().parent.parent
results_dir = cluster_dir / "manifold_results" / f"config_{config_id}"
results_dir.mkdir(parents=True, exist_ok=True)

# Save manifold arrays
np.savez_compressed(
    results_dir / "manifolds_opt.npz",
    stable1=stable1_2,
    unstable1=unstable1_2,
    stable2=stable2_2,
    unstable2=unstable2_2,
    stable3=stable3_2,
    unstable3=unstable3_2,
    stable4=stable4_2,
    unstable4=unstable4_2,
    stable_unwrapped1=stable_unwrapped1_2,
    unstable_unwrapped1=unstable_unwrapped1_2,
    stable_unwrapped4=stable_unwrapped4_2,
    unstable_unwrapped4=unstable_unwrapped4_2
)

print(f"Manifold results saved to {results_dir / 'manifolds_opt.npz'}")

fig, ax = pplot.plot(marker=".", s=0.5, xlim=[.0, 2.], plottype=None)

ax.plot(stable_unwrapped1_2[:, 0], stable_unwrapped1_2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped1_2[:, 0], unstable_unwrapped1_2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable2_2[:, 0], stable2_2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable2_2[:, 0], unstable2_2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped4_2[:, 0], stable_unwrapped4_2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped4_2[:, 0], unstable_unwrapped4_2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable3_2[:, 0], stable3_2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable3_2[:, 0], unstable3_2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

# Save figure
fig.savefig(results_dir / "manifolds_opt.png", dpi=300, bbox_inches="tight")

print(f"Figure saved to {results_dir / 'manifolds_opt.png'}")

plt.close(fig)