import numpy as np
from pyoculus.fields.slab_bfield_two_pert import SlabBfieldTwoPert
from pyoculus.maps import ToroidalBfieldSection
from pyoculus.solvers import PoincarePlot, FixedPoint, Manifold
import matplotlib.pyplot as plt

from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
import traceback
from test_manifold_config import test_manifold_configuration
import json

config_id = 212

cluster_dir = Path(__file__).resolve().parent.parent

params_path = (
    cluster_dir / "configurations" / "parameters"
    / f"parameters_{config_id}.json"
)

with open(params_path, "r", encoding="utf-8") as f:
    field_params = json.load(f)

field = SlabBfieldTwoPert(**field_params)

xs1 = field.x0 + field.n1 / (field.iota_prime * field.m1)
xs2 = field.x0 + field.n2 / (field.iota_prime * field.m2)

map = ToroidalBfieldSection(field, domain=[(0,2), (0, 2*np.pi)], type="dopri5", rtol=1e-10, nsteps=200000)


# Fixed points:
# X-points island 1:
fp1 = FixedPoint(map)
fp1.find(t = 2, guess = np.array([1.1, 0.]), method="scipy.root")  

fp2 = FixedPoint(map)
fp2.find(t = 2, guess = np.array(fp1.coords[1]), method="scipy.root")
# fp2.find(t = 2, guess = np.array([1.21, np.pi]), method="scipy.root")

fp3 = FixedPoint(map)
fp3.find(t = 2, guess = np.array([1.09934463, 6.28318478]), method="scipy.root")

# X-points island 2:
fp4 = FixedPoint(map)
fp4.find(t = 1, guess = np.array([xs2, 0.]), method="scipy.root")

fp5 = FixedPoint(map)
fp5.find(t = 1, guess = np.array([0.733, 6.23]), method="scipy.root")  

# O-points island 1:
ofp1 = FixedPoint(map)
ofp1.find_with_iota(n = field.n1, m = field.m1, guess = np.array([xs1, np.pi/2]), method="scipy.root") 

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


# Test manifold configurations

results_dir = (
    cluster_dir / "manifold_test_results" / f"config_{config_id}"
)
results_dir.mkdir(parents=True, exist_ok=True)

output_file = results_dir / "manifold_4_eps_test.txt"

with open(output_file, "w", encoding="utf-8", buffering=1) as f:
    with redirect_stdout(f), redirect_stderr(f):

        print("=" * 60)
        print(f"Configuration: {config_id}")
        print("Manifold: 4")
        print(f"Field parameters: {field_params}")
        print("nsteps: 200000")
        print("ERR: 4e-3")
        print("=" * 60, flush=True)

        eps1 = [2e-8, 5e-8, 1e-8, 5e-8]
        eps2 = [2e-8, 5e-8, 5e-8, 1e-8]
        for i in range(4):

            print("\n" + "=" * 60)
            print(f"Testing epsilons = {eps1[i]:.1e}, {eps2[i]:.1e}")
            print("=" * 60, flush=True)

            try:
                # Create a fresh manifold for each test
                manifold4 = Manifold(
                    map, fp2, fp3,
                    dir1=ofp2.coords[0] - fp2.coords[0],
                    dir2=ofp2.coords[0] - fp3.coords[0],
                    first_stable=False
                )

                manifold4.vector_s = np.real_if_close(
                    manifold4.vector_s
                )
                manifold4.vector_u = np.real_if_close(
                    manifold4.vector_u
                )

                manifold4.find_clinics_exhaustive(
                    first_guess_eps_s=eps1[i],
                    first_guess_eps_u=eps2[i],
                    ERR=4e-3
                )

                print(
                    f"ClinicSet empty: {manifold4.clinics.is_empty}",
                    flush=True
                )

                if not manifold4.clinics.is_empty:
                    print(
                        f"nint_pair: {manifold4.clinics.nint_pair}",
                        flush=True
                    )

            except Exception:
                traceback.print_exc()

print(f"Finished manifold 4 tests: {output_file}", flush=True)