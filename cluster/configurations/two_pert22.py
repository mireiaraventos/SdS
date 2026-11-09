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
    m1=2,
    n1=1,
    delta1=0.3,
    m2=1,
    n2=-1,
    delta2=0.3,
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
# X-points island 1:
fp1 = FixedPoint(map)
fp1.find(t = 2, guess = np.array([1.25, 0.]), method="scipy.root")  

fp2 = FixedPoint(map)
fp2.find(t = 2, guess = np.array(fp1.coords[1]), method="scipy.root")

fp3 = FixedPoint(map)
fp3.find(t = 2, guess = np.array(fp1.coords[2]), method="scipy.root")

# X-points island 2:
fp4 = FixedPoint(map)
fp4.find_with_iota(n = field.n2, m = field.m2, guess = np.array([xs2, 0.0]), method="scipy.root")  

fp5 = FixedPoint(map)
fp5.find(t = 1, guess = np.array([xs2, 0.0]), method="scipy.root")  

# O-points island 1:
ofp2 = FixedPoint(map)
ofp2.find_with_iota(n = field.n1, m = field.m1, guess = np.array([1.3, 4.7]), method="scipy.root")

ofp1 = FixedPoint(map)
ofp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array(ofp2.coords[1]), method="scipy.root")

# O-points island 2:
ofp3 = FixedPoint(map)
ofp3.find_with_iota(n = field.n2, m = field.m2, guess = np.array([xs2, np.pi]), method="scipy.root")


# Manifolds:
# Island 1 top:
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

# Island 1 bottom:
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

# Configuration number
config_id = 22  # Change to 2, 3, etc. in the other scripts

manifold1.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=2.2e-3)
stable1_2, unstable1_2 = manifold1.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

manifold2.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=2.2e-3)
stable2_2, unstable2_2 = manifold2.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

manifold3.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=2.2e-3)
stable3_2, unstable3_2 = manifold3.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

manifold4.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=2.2e-3)
stable4_2, unstable4_2 = manifold4.compute_adaptive2(tol=3e-3, n_u_extra=1, n_s_extra=1)

manifold5.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=2.2e-3)
stable5_2, unstable5_2 = manifold5.compute_adaptive2(tol=3e-3, n_u_extra=0, n_s_extra=0)

manifold6.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=2.2e-3)
stable6_2, unstable6_2 = manifold6.compute_adaptive2(tol=3e-3, n_u_extra=0, n_s_extra=0)

stable_unwrapped1_2 = stable1_2.copy()
unstable_unwrapped1_2 = unstable1_2.copy()
stable_unwrapped1_2[:, 1] = np.unwrap(stable1_2[:, 1], period=2*np.pi)
unstable_unwrapped1_2[:, 1] = np.unwrap(unstable1_2[:, 1], period=2*np.pi)

stable_unwrapped2_2 = stable2_2.copy()
unstable_unwrapped2_2 = unstable2_2.copy()
stable_unwrapped2_2[:, 1] = np.unwrap(stable2_2[:, 1], period=2*np.pi)
unstable_unwrapped2_2[:, 1] = np.unwrap(unstable2_2[:, 1], period=2*np.pi)

stable_unwrapped3_2 = stable3_2.copy()
unstable_unwrapped3_2 = unstable3_2.copy()
stable_unwrapped3_2[:, 1] = np.unwrap(stable3_2[:, 1], period=2*np.pi)
unstable_unwrapped3_2[:, 1] = np.unwrap(unstable3_2[:, 1], period=2*np.pi)

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
    stable5=stable5_2,
    unstable5=unstable5_2,
    stable6=stable6_2,
    unstable6=unstable6_2,
    stable_unwrapped1=stable_unwrapped1_2,
    unstable_unwrapped1=unstable_unwrapped1_2,
    stable_unwrapped2=stable_unwrapped2_2,
    unstable_unwrapped2=unstable_unwrapped2_2,
    stable_unwrapped3=stable_unwrapped3_2,
    unstable_unwrapped3=unstable_unwrapped3_2,
    stable_unwrapped4=stable_unwrapped4_2,
    unstable_unwrapped4=unstable_unwrapped4_2
)

print(f"Manifold results saved to {results_dir / 'manifolds_opt.npz'}")

fig, ax = pplot.plot(marker=".", s=0.5, xlim=[.0, 2.], plottype=None)

ax.plot(stable_unwrapped1_2[:400, 0], stable_unwrapped1_2[:400, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped1_2[:2000, 0], unstable_unwrapped1_2[:2000, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped2_2[:2000, 0], stable_unwrapped2_2[:2000, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped2_2[:500, 0], unstable_unwrapped2_2[:500, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped3_2[:270, 0], stable_unwrapped3_2[:270, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped3_2[:230, 0], unstable_unwrapped3_2[:230, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable_unwrapped4_2[:500, 0], stable_unwrapped4_2[:500, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable_unwrapped4_2[:320, 0], unstable_unwrapped4_2[:320, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable5_2[:1000, 0], stable5_2[:1000, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable5_2[:, 0], unstable5_2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable6_2[:, 0], stable6_2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable6_2[:, 0], unstable6_2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

# Save figure
fig.savefig(results_dir / "manifolds_opt.png", dpi=300, bbox_inches="tight")

print(f"Figure saved to {results_dir / 'manifolds_opt.png'}")

plt.close(fig)