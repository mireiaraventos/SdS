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
    delta1=0.05,
    m2=1,
    n2=-1,
    delta2=0.05,
    iota_prime=2.0,
    x0=1.0,
)

xs1 = field.x0 + field.n1 / (field.iota_prime * field.m1)
xs2 = field.x0 + field.n2 / (field.iota_prime * field.m2)

map = ToroidalBfieldSection(field, domain=[(0, 2), (0, 2*np.pi)])

pplot = PoincarePlot.with_segments(
    map,
    np.array([
        [0.05, np.pi/2],
        [1.95, np.pi/2],
        [0.05, 3*np.pi/2],
        [1.95, 3*np.pi/2]
    ]),
    neps=[50, 50],
    connected=False,
)

# Follow each initial condition for 200 crossings
pplot.compute(npts=150)

# Fixed points:
fp1 = FixedPoint(map)
fp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([xs1, 0.0]), method="scipy.root")   

fp2 = FixedPoint(map)
fp2.find_with_iota(n = field.n1, m = field.m1, guess = np.array(fp1.coords[1]), method="scipy.root")  

fp3 = FixedPoint(map)
fp3.find_with_iota(n = field.n1, m = field.m1, guess = np.array(fp1.coords[2]), method="scipy.root")  

fp4 = FixedPoint(map)
fp4.find_with_iota(n = field.n2, m = field.m2, guess = np.array([0.5, 0.0]), method="scipy.root")  

fp5 = FixedPoint(map)
fp5.find(t = 1, guess = np.array([0.5, 0.0]), method="scipy.root")  

ofp1 = FixedPoint(map)
ofp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([xs1, 4]), method="scipy.root")

ofp2 = FixedPoint(map)
ofp2.find_with_iota(n = field.n1, m = field.m1, guess = np.array([xs1, 2]), method="scipy.root")

ofp3 = FixedPoint(map)
ofp3.find_with_iota(n = field.n2, m = field.m2, guess = np.array([xs2, 3]), method="scipy.root")


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
config_id = 21  # Change to 2, 3, etc. in the other scripts

manifold1.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.5e-3)
stable1, unstable1 = manifold1.compute_adaptive2(tol=2e-3, n_u_extra=1, n_s_extra=1)
    
manifold2.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6, ERR=1.5e-3)
stable2, unstable2 = manifold2.compute_adaptive2(tol=2e-3,  n_u_extra=1, n_s_extra=1)

manifold3.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6)
stable3, unstable3 = manifold3.compute_adaptive2(tol=2e-3, n_u_extra=1, n_s_extra=1)

manifold4.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6)
stable4, unstable4 = manifold4.compute_adaptive2(tol=2e-3, n_u_extra=1, n_s_extra=1)

manifold5.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6)
stable5, unstable5 = manifold5.compute_adaptive2(tol=2e-3, n_u_extra=1, n_s_extra=1)

manifold6.find_clinics_exhaustive(first_guess_eps_s=1e-6, first_guess_eps_u=1e-6)
stable6, unstable6 = manifold6.compute_adaptive2(tol=2e-3, n_u_extra=1, n_s_extra=1)

# Results directory, independent of the working directory
cluster_dir = Path(__file__).resolve().parent.parent
results_dir = cluster_dir / "manifold_results" / f"config_{config_id}"
results_dir.mkdir(parents=True, exist_ok=True)

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
    stable5=stable5,
    unstable5=unstable5,
    stable6=stable6,
    unstable6=unstable6
)

print(f"Manifold results saved to {results_dir / 'manifolds.npz'}")

fig, ax = pplot.plot(marker=".", s=0.5, xlim=[.0, 2.], plottype=None)

ax.plot(stable1[:, 0], stable1[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable1[:, 0], unstable1[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable2[:, 0], stable2[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable2[:, 0], unstable2[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable3[:, 0], stable3[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable3[:, 0], unstable3[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable4[:, 0], stable4[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable4[:, 0], unstable4[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable5[:, 0], stable5[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable5[:, 0], unstable5[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

ax.plot(stable6[:, 0], stable6[:, 1], color="blue", label="Stable", marker=".", markersize=3, linestyle="")
ax.plot(unstable6[:, 0], unstable6[:, 1], color="red", label="Unstable", marker=".", markersize=3, linestyle="")

# Save figure
fig.savefig(results_dir / "manifolds.png", dpi=300, bbox_inches="tight")

print(f"Figure saved to {results_dir / 'manifolds.png'}")

plt.close(fig)