"""Plots for the OV 3 3 0 magnetic_2 mapping test."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DAT = Path(__file__).resolve().parents[1] / "dat" / "magnetic_2"
OUT = DAT

COLS_INT = [
    "SGx1", "SGx2", "SGy1", "SGy2",
    "chrom_x", "chrom_y", "chrom_ell_x", "chrom_ell_y",
    "eta0", "eta1",
    "MU_ref", "quad_Kx", "quad_Ky", "quad_Kz",
    "alpha_0", "alpha_1", "Cs", "LTx", "LTy",
]
COLS_S0 = ["SGx1", "SGx2", "SGy1", "SGy2", "spin_ref", "spin_x", "spin_y", "spin_z"]
COLS_S1 = ["SGx1", "SGx2", "SGy1", "SGy2", "spin_ref_1", "spin_x_1", "spin_y_1", "spin_z_1"]


def _read(path: Path, names: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", names=names, engine="python", skiprows=1)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col].astype(str).str.replace("D", "E"), errors="coerce")
    return df.dropna().reset_index(drop=True)


def load() -> pd.DataFrame:
    a = _read(DAT / "integrals_ov330.dat", COLS_INT)
    s0 = _read(DAT / "particle_spin_tune_ov330.dat", COLS_S0)
    s1 = _read(DAT / "particle_spin_tune_1_ov330.dat", COLS_S1)
    keys = ["SGx1", "SGx2", "SGy1", "SGy2"]
    return a.merge(s0, on=keys).merge(s1, on=keys)


def tss_ok(df: pd.DataFrame) -> pd.DataFrame:
    """Drop points where TSS returned MU=0 (spin tune not extracted)."""
    return df[df.MU_ref != 0].reset_index(drop=True)


def _fit_spin(df: pd.DataFrame, ycol: str, xcols: tuple[str, ...]) -> tuple[np.ndarray, float]:
    x = np.column_stack([np.ones(len(df)), df[list(xcols)].to_numpy()])
    y = df[ycol].to_numpy()
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    return coef, r2


def plot_dnu_vs_xi(df: pd.DataFrame) -> Path:
    xvars = [
        ("chrom_x", r"$\xi_{J,x}$ from MAP(5)"),
        ("chrom_y", r"$\xi_{J,y}$ from MAP(5)"),
        ("eta1", r"$\eta_1$ from MAP(5)"),
    ]
    yvars = [
        ("spin_x", r"$\Delta\nu_{s,x}$ (P0, $x=1$ mm)"),
        ("spin_y", r"$\Delta\nu_{s,y}$ (P0, $y=1$ mm)"),
        ("spin_z", r"$\Delta\nu_{s,z}$ (P0, $\delta=10^{-4}$)"),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(11.5, 10.5), sharex="col", sharey="row")
    fig.suptitle(
        r"magnetic_2, OV 3 3 0: $\Delta\nu_s$ of offset particles vs $\xi_J$, not vs Berz $\xi_B$",
        fontsize=12,
    )
    nat = df[(df.SGx1 == 0) & (df.SGx2 == 0) & (df.SGy1 == 0)]
    for i, (ycol, ylab) in enumerate(yvars):
        for j, (xcol, xlab) in enumerate(xvars):
            ax = axes[i, j]
            ax.scatter(df[xcol], df[ycol], s=22, c="#1f4e79", alpha=0.75, zorder=2)
            if len(nat):
                ax.scatter(
                    nat[xcol], nat[ycol], s=70, c="#c0392b", zorder=3,
                    label="sextupoles off",
                )
            m, b = np.polyfit(df[xcol], df[ycol], 1)
            xx = np.linspace(df[xcol].min(), df[xcol].max(), 50)
            ax.plot(xx, m * xx + b, color="#e67e22", lw=1)
            if i == 2:
                ax.set_xlabel(xlab)
            if j == 0:
                ax.set_ylabel(ylab)
            ax.grid(True, alpha=0.3)
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    fig.tight_layout()
    out = OUT / "fig_ov330_dnu_vs_xi.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_p0_p1(df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    fig.suptitle(r"P0 vs P1 $\Delta\nu_s$ (same $\xi_J$)", fontsize=12)
    pairs = [
        ("spin_x", "spin_x_1", r"$x$-probe"),
        ("spin_y", "spin_y_1", r"$y$-probe"),
        ("spin_z", "spin_z_1", r"$\delta$-probe"),
    ]
    for ax, (a, b, title) in zip(axes, pairs):
        ax.scatter(df[a], df[b], s=22, c="#1f4e79", alpha=0.75)
        lims = [
            min(df[a].min(), df[b].min()),
            max(df[a].max(), df[b].max()),
        ]
        ax.plot(lims, lims, "k--", lw=0.8)
        ax.set_title(title)
        ax.set_xlabel("P0")
        ax.set_ylabel("P1")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = OUT / "fig_ov330_p0_p1.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_link(df: pd.DataFrame) -> Path:
    xcols = ("chrom_x", "chrom_y", "eta1")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    fig.suptitle(r"Linear map $\Delta\nu_s \approx a\,\xi_{J,x}+b\,\xi_{J,y}+c\,\eta_1+d$", fontsize=12)
    for ax, plane in zip(axes, ("x", "y", "z")):
        col = f"spin_{plane}"
        coef, r2 = _fit_spin(df, col, xcols)
        x = np.column_stack([np.ones(len(df)), df[list(xcols)].to_numpy()])
        pred = x @ coef
        y = df[col].to_numpy()
        ax.scatter(y, pred, s=22, c="#6c3483", alpha=0.75)
        lims = [min(y.min(), pred.min()), max(y.max(), pred.max())]
        ax.plot(lims, lims, "r--", lw=0.8)
        ax.set_title(fr"plane {plane}, $R^2={r2:.4f}$")
        ax.set_xlabel(r"measured $\Delta\nu_s$")
        ax.set_ylabel(r"predicted from $\xi_J$")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = OUT / "fig_ov330_link.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def summary(df: pd.DataFrame) -> str:
    n_all = len(df)
    n_fail = int((df.MU_ref == 0).sum())
    df = tss_ok(df)
    nat = df[(df.SGx1 == 0) & (df.SGx2 == 0) & (df.SGy1 == 0)].iloc[0]
    lines = [
        f"points: {len(df)}  (TSS fail MU=0 dropped: {n_fail} of {n_all})",
        f"natural xi_J = ({nat.chrom_x:.6f}, {nat.chrom_y:.6f})",
        f"natural xi_ell = ({nat.chrom_ell_x:.6f}, {nat.chrom_ell_y:.6f})",
        f"natural eta0, eta1 = {nat.eta0:.6f}, {nat.eta1:.6f}",
        f"xi_J range x: [{df.chrom_x.min():.3f}, {df.chrom_x.max():.3f}]"
        f"  y: [{df.chrom_y.min():.3f}, {df.chrom_y.max():.3f}]",
        f"natural Delta nu_s P0 (x,y,z) = {nat.spin_x:.6e}, {nat.spin_y:.6e}, {nat.spin_z:.6e}",
    ]
    for plane in ("x", "y", "z"):
        coef, r2 = _fit_spin(df, f"spin_{plane}", ("chrom_x", "chrom_y", "eta1"))
        lines.append(
            f"spin_{plane}: R^2={r2:.6f}  "
            f"a(xi_x)={coef[1]:.4e}  a(xi_y)={coef[2]:.4e}  a(eta1)={coef[3]:.4e}  d={coef[0]:.4e}"
        )
    return "\n".join(lines)


def main() -> None:
    raw = load()
    print(summary(raw))
    df = tss_ok(raw)
    for fn in (plot_dnu_vs_xi, plot_p0_p1, plot_link):
        print("wrote", fn(df))


if __name__ == "__main__":
    main()
