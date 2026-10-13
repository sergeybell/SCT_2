"""Fit slope = kappa * (-2*pi*xi/C) + b on the sextupole scans."""

import numpy as np

# pred = -2*pi*xi/C, then slopeT, slopeC. Magnetic from the summary,
# plus the large-J phase average where it was repeated.
mag = np.array(
    [
        [6.852035e-03, 8.373388e-03, 8.212726e-03],
        [-7.699560e-03, -7.278178e-03, -7.281440e-03],
        [-2.221941e-02, -2.291348e-02, -2.307313e-02],
    ]
)
# electrostatic: summary points and the S=0.0025 phase average
es = np.array(
    [
        [1.730483e-02, 1.158699e-02, -7.649001e-03],
        [3.193883e-03, -5.759835e-03, -1.816941e-02],
        [-4.120843e-04, -1.010524e-02, -2.081090e-02],
        [-1.099608e-02, -2.315520e-02, -2.875255e-02],
        [-3.932617e-02, -5.791095e-02, -4.989642e-02],
        [-6.766645e-02, -9.268831e-02, -7.107790e-02],
    ]
)


def fit(rows, name):
    pred, st, sc = rows.T
    a = np.column_stack([pred, np.ones_like(pred)])
    kt, bt = np.linalg.lstsq(a, st, rcond=None)[0]
    kc, bc = np.linalg.lstsq(a, sc, rcond=None)[0]
    print(f"== {name} ==")
    print(f"time  kappa={kt:.4f} offset={bt:.6e}  rms={np.std(st-kt*pred-bt):.3e}")
    print(f"path  kappa={kc:.4f} offset={bc:.6e}  rms={np.std(sc-kc*pred-bc):.3e}")
    print("residuals T", np.array2string(st - kt * pred - bt, precision=3))
    # value at xi=0
    print(f"at xi=0: slopeT={bt:.6e} slopeC={bc:.6e}")


fit(mag, "magnetic, summary S=0,0.002,0.004")
fit(es, "electrostatic")
