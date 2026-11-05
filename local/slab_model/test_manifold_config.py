
import numpy as np
import traceback

def numerical_jacobian(map_obj, x, t=2, h=1e-7):
    x = np.asarray(x, dtype=float)
    J = np.zeros((2, 2))

    for j in range(2):
        dx = np.zeros(2)
        dx[j] = h

        fp = np.asarray(map_obj.f(t, x + dx), dtype=float)
        fm = np.asarray(map_obj.f(t, x - dx), dtype=float)

        diff = fp - fm
        diff[1] = (diff[1] + np.pi) % (2*np.pi) - np.pi

        J[:, j] = diff / (2*h)

    return J


def test_manifold_configuration(
    map_obj,
    manifold,
    t=2,
    eps_values=None,
    run_clinics=False,
    clinic_kwargs=None
):
    if eps_values is None:
        eps_values = np.logspace(-4, -10, 7)

    print("=" * 65)
    print("MANIFOLD CONFIGURATION DIAGNOSTIC")
    print("=" * 65)

    # --------------------------------------------------
    # 1. Fixed-point Jacobians
    # --------------------------------------------------

    for label, fp in [
        ("Stable reference FP", manifold.rfp_s),
        ("Unstable reference FP", manifold.rfp_u)
    ]:
        x = np.asarray(fp, dtype=float)

        J_num = numerical_jacobian(map_obj, x, t=t)
        J_df = np.asarray(map_obj.df(t, x), dtype=float)

        rel_error = (
            np.linalg.norm(J_df - J_num)
            / np.linalg.norm(J_num)
        )

        print(f"\n--- {label} ---")
        print("Coordinates:", x)
        print("Jacobian relative error:", rel_error)
        print("det(J_df):", np.linalg.det(J_df))
        print("Eigenvalues:", np.linalg.eigvals(J_df))

    # --------------------------------------------------
    # 2. Manifold eigenvalues and eigenvectors
    # --------------------------------------------------

    print("\n--- Manifold eigendata ---")
    print("lambda_s:", manifold.lambda_s)
    print("lambda_u:", manifold.lambda_u)
    print("vector_s:", manifold.vector_s)
    print("vector_u:", manifold.vector_u)

    # --------------------------------------------------
    # 3. Linear-regime checks
    # --------------------------------------------------

    print("\n--- Linear-regime errors ---")
    print(
        f"{'epsilon':>12} "
        f"{'stable error':>18} "
        f"{'unstable error':>18}"
    )

    for eps in eps_values:
        try:
            err_s = manifold.error_linear_regime(
                eps,
                manifold.rfp_s,
                manifold.vector_s,
                direction=-1
            )
        except Exception:
            err_s = np.nan

        try:
            err_u = manifold.error_linear_regime(
                eps,
                manifold.rfp_u,
                manifold.vector_u,
                direction=1
            )
        except Exception:
            err_u = np.nan

        print(f"{eps:12.1e} {err_s:18.6e} {err_u:18.6e}")

    # --------------------------------------------------
    # 4. Clinic search
    # --------------------------------------------------

    if run_clinics:
        print("\n--- Clinic search ---")

        if clinic_kwargs is None:
            clinic_kwargs = {
                "first_guess_eps_s": 1e-6,
                "first_guess_eps_u": 1e-6,
                "ERR": 2.2e-3
            }

        try:
            manifold.find_clinics_exhaustive(**clinic_kwargs)
            print("Clinic search completed successfully.")

            if hasattr(manifold, "clinics"):
                print("Clinics:", manifold.clinics)

        except Exception as e:
            print("Clinic search failed.")
            print("Error:", str(e))
            traceback.print_exc()

    print("\n" + "=" * 65)
