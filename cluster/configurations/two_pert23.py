import numpy as np
from pyoculus.fields.slab_bfield_two_pert import SlabBfieldTwoPert
from pyoculus.maps import ToroidalBfieldSection
from pyoculus.solvers import PoincarePlot, FixedPoint, Manifold
import matplotlib.pyplot as plt

from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
import traceback
from test_manifold_config import test_manifold_configuration

field = SlabBfieldTwoPert(
    m1=2,
    n1=1,
    delta1=0.5,
    m2=1,
    n2=-1,
    delta2=0.5,
    iota_prime=4.0,
    x0=1.0,
)

xs1 = field.x0 + field.n1 / (field.iota_prime * field.m1)
xs2 = field.x0 + field.n2 / (field.iota_prime * field.m2)

map = ToroidalBfieldSection(field, domain=[(0, 2), (0, 2*np.pi)])


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
    # dir1=ofp1.coords[0] - fp2.coords[0],
    # dir2=ofp1.coords[0] - fp1.coords[0],
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

# Configuration number
config_id = 23  # Change to 2, 3, etc. in the other scripts

# Results directory, independent of the working directory
cluster_dir = Path(__file__).resolve().parent.parent
results_dir = cluster_dir / "manifold_test_results" / f"config_{config_id}"
results_dir.mkdir(parents=True, exist_ok=True)

manifolds = [manifold1, manifold2, manifold3, manifold4]

for i, manifold in enumerate(manifolds, start=1):
    output_file = results_dir / f"manifold_{i}.txt"

    with open(output_file, "w", encoding="utf-8", buffering=1) as f:
        with redirect_stdout(f), redirect_stderr(f):
            print(f"Configuration {config_id} - Manifold {i}", flush=True)

            try:
                test_manifold_configuration(
                    map_obj=map,
                    manifold=manifold,
                    t=2,
                    run_clinics=True,
                    clinic_kwargs={
                        "first_guess_eps_s": 1e-6,
                        "first_guess_eps_u": 1e-6,
                        "ERR": 2.2e-3
                    }
                )
            except Exception:
                traceback.print_exc()

    print(f"Finished manifold {i}: {output_file}", flush=True)