#!/usr/bin/env python3
"""Make the Bagpipes SED-fit example used on the SAGUI website.

This is intentionally a lightweight diagnostic fit for one representative
regional SED, not a replacement for a full production SED-fitting pipeline.
It reads the flux-preserving SAGUI regional photometry table, fits a simple
Bagpipes model to Sagui-8 region 8, and writes both the plot and the
model-vs-data table.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy.optimize import differential_evolution, least_squares

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")

import bagpipes  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
PHOTOMETRY_PATH = Path(
    "/Users/rd23aag/Documents/GitHub/crp8_segmentation/"
    "results/flux_per_region/SED_flux_wide_sagui8.csv"
)
THROUGHPUT_PATH = ROOT / "inst/extdata/throughput_nircam.csv"
OUTPUT_DIR = ROOT / "outputs/website_sedfit"
FILTER_DIR = OUTPUT_DIR / "bagpipes_filters"
IMAGE_PATH = ROOT / "images/examples/sagui8_region8_bagpipes_sedfit.png"
TABLE_PATH = OUTPUT_DIR / "sagui8_region8_bagpipes_model_photometry.csv"

REGION_ID = 8
REDSHIFT = 0.6222
SYSTEMATIC_FRAC = 0.06

BANDS = [
    "F090W",
    "F115W",
    "F150W",
    "F182M",
    "F200W",
    "F210M",
    "F277W",
    "F335M",
    "F356W",
    "F410M",
    "F430M",
    "F444W",
    "F460M",
    "F480M",
]

CENWAVE_UM = {
    "F090W": 0.902,
    "F115W": 1.154,
    "F150W": 1.501,
    "F182M": 1.845,
    "F200W": 1.990,
    "F210M": 2.093,
    "F277W": 2.786,
    "F335M": 3.365,
    "F356W": 3.563,
    "F410M": 4.092,
    "F430M": 4.280,
    "F444W": 4.421,
    "F460M": 4.624,
    "F480M": 4.834,
}

# Approximate wavelength-ordered colors. These are not meant to reproduce
# human vision in the NIR; they simply make the filter ordering intuitive.
FILTER_COLORS = {
    "F090W": "#3F4EA1",
    "F115W": "#315FB9",
    "F150W": "#2577D9",
    "F182M": "#1FA2D6",
    "F200W": "#28B5B5",
    "F210M": "#39BF75",
    "F277W": "#90D94A",
    "F335M": "#E4D642",
    "F356W": "#F0C145",
    "F410M": "#EE9A3A",
    "F430M": "#E97832",
    "F444W": "#D85A2A",
    "F460M": "#C9432B",
    "F480M": "#A9362D",
}

SAGUI_BLUE = "#213E60"
SAGUI_ORANGE = "#E68C3A"


def write_bagpipes_filters(throughput: pd.DataFrame) -> list[str]:
    FILTER_DIR.mkdir(parents=True, exist_ok=True)
    filt_paths = []

    for band in BANDS:
        curve = throughput.loc[throughput["filter"] == band]
        path = FILTER_DIR / f"{band}.dat"
        np.savetxt(
            path,
            np.column_stack(
                [
                    curve["wavelength"].to_numpy(dtype=float) * 1e4,
                    curve["throughput"].to_numpy(dtype=float),
                ]
            ),
        )
        filt_paths.append(str(path))

    return filt_paths


def read_region_photometry() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = pd.read_csv(PHOTOMETRY_PATH)
    row = table.loc[table["region"] == REGION_ID]
    if row.empty:
        raise ValueError(f"Region {REGION_ID} not found in {PHOTOMETRY_PATH}")

    row = row.iloc[0]
    obs_ujy = np.array([row[band] for band in BANDS], dtype=float) * 1e6
    err_ujy = np.array([row[f"{band}_err"] for band in BANDS], dtype=float) * 1e6
    err_ujy = np.sqrt(err_ujy**2 + (SYSTEMATIC_FRAC * np.maximum(obs_ujy, 1e-9)) ** 2)
    wave_um = np.array([CENWAVE_UM[band] for band in BANDS], dtype=float)

    return wave_um, obs_ujy, err_ujy


def make_components(params: np.ndarray) -> dict:
    log_mass, age_gyr, log_tau_gyr, log_metallicity, av = params

    return {
        "redshift": REDSHIFT,
        "exponential": {
            "age": age_gyr,
            "tau": 10**log_tau_gyr,
            "massformed": log_mass,
            "metallicity": 10**log_metallicity,
        },
        "dust": {"type": "Calzetti", "Av": av},
    }


def fit_bagpipes_model(
    filt_paths: list[str],
    obs_ujy: np.ndarray,
    err_ujy: np.ndarray,
) -> tuple[bagpipes.model_galaxy, np.ndarray, float]:
    spec_wavs = np.linspace(7600, 56000, 1400)
    start = np.array([10.2, 3.0, 0.0, np.log10(0.5), 0.4])

    model = bagpipes.model_galaxy(
        make_components(start),
        filt_list=filt_paths,
        spec_wavs=spec_wavs,
        spec_units="mujy",
        phot_units="mujy",
    )

    def model_photometry(params: np.ndarray) -> np.ndarray:
        try:
            model.update(make_components(params))
            phot = np.asarray(model.photometry, dtype=float)
            if not np.all(np.isfinite(phot)):
                return np.full_like(obs_ujy, 1e30)
            return phot
        except Exception:
            return np.full_like(obs_ujy, 1e30)

    def residuals(params: np.ndarray) -> np.ndarray:
        return (model_photometry(params) - obs_ujy) / err_ujy

    bounds = (
        [7.5, 0.08, -1.3, np.log10(0.005), 0.0],
        [12.5, 7.2, 1.0, np.log10(2.5), 3.5],
    )

    starts = [
        start,
        [10.0, 1.0, -0.3, np.log10(1.0), 0.2],
        [10.5, 5.0, 0.3, np.log10(0.2), 0.8],
        [9.8, 2.0, -0.6, np.log10(0.4), 0.0],
        [10.8, 6.5, 0.6, np.log10(1.0), 1.2],
    ]

    fits = [
        least_squares(
            residuals,
            trial,
            bounds=bounds,
            max_nfev=400,
            xtol=1e-5,
            ftol=1e-5,
            gtol=1e-5,
        )
        for trial in starts
    ]
    best = min(fits, key=lambda fit: np.sum(fit.fun**2))

    # Short global pass to avoid a visually poor local minimum.
    global_fit = differential_evolution(
        lambda x: np.sum(residuals(x) ** 2),
        list(zip(*bounds)),
        maxiter=20,
        popsize=8,
        polish=False,
        seed=13,
        workers=1,
        updating="immediate",
    )
    polished = least_squares(
        residuals,
        global_fit.x,
        bounds=bounds,
        max_nfev=500,
        xtol=1e-5,
        ftol=1e-5,
        gtol=1e-5,
    )
    if np.sum(polished.fun**2) < np.sum(best.fun**2):
        best = polished

    model.update(make_components(best.x))
    rms = float(np.sqrt(np.mean(residuals(best.x) ** 2)))
    return model, best.x, rms


def write_model_table(
    wave_um: np.ndarray,
    obs_ujy: np.ndarray,
    err_ujy: np.ndarray,
    model_ujy: np.ndarray,
    rms: float,
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(
        {
            "band": BANDS,
            "wave_um": wave_um,
            "obs_ujy": obs_ujy,
            "err_ujy": err_ujy,
            "model_ujy": model_ujy,
            "pull": (obs_ujy - model_ujy) / err_ujy,
            "fit_rms_sigma": rms,
        }
    )
    table.to_csv(TABLE_PATH, index=False)


def plot_sed(
    throughput: pd.DataFrame,
    wave_um: np.ndarray,
    obs_ujy: np.ndarray,
    err_ujy: np.ndarray,
    model: bagpipes.model_galaxy,
    rms: float,
) -> None:
    IMAGE_PATH.parent.mkdir(parents=True, exist_ok=True)

    spectrum = np.asarray(model.spectrum, dtype=float)
    model_ujy = np.asarray(model.photometry, dtype=float)
    y_values = np.r_[obs_ujy, model_ujy]
    y_min = float(np.nanmin(y_values) * 0.55)
    y_max = float(np.nanmax(y_values) * 1.28)
    filter_base = y_min * 1.08
    filter_height = (y_max - y_min) * 0.20

    fig, ax = plt.subplots(figsize=(8.2, 4.7), dpi=220)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for band in BANDS:
        curve = throughput.loc[throughput["filter"] == band]
        x = curve["wavelength"].to_numpy(dtype=float)
        y = curve["throughput"].to_numpy(dtype=float)
        if np.nanmax(y) > 0:
            y = y / np.nanmax(y)
        ax.fill_between(
            x,
            filter_base,
            filter_base + filter_height * y,
            color=FILTER_COLORS[band],
            alpha=0.28,
            linewidth=0,
            zorder=0,
        )

    visible = (
        (spectrum[:, 0] >= 7500)
        & (spectrum[:, 0] <= 56000)
        & np.isfinite(spectrum[:, 1])
    )
    ax.plot(
        spectrum[visible, 0] / 1e4,
        spectrum[visible, 1],
        color=SAGUI_BLUE,
        linewidth=2.35,
        label="Bagpipes model",
        zorder=2,
    )
    ax.plot(wave_um, model_ujy, color=SAGUI_BLUE, linewidth=1.1, alpha=0.25, zorder=3)
    ax.scatter(wave_um, model_ujy, s=22, color=SAGUI_BLUE, zorder=4)
    ax.errorbar(
        wave_um,
        obs_ujy,
        yerr=err_ujy,
        fmt="o",
        ms=7.4,
        mfc="white",
        mec=SAGUI_ORANGE,
        mew=1.9,
        ecolor="#424242",
        elinewidth=1.25,
        capsize=2.8,
        label="Observed photometry",
        zorder=5,
    )

    ax.text(
        0.03,
        0.94,
        "Sagui-8, region 8",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        color=SAGUI_BLUE,
    )
    ax.text(
        0.03,
        0.875,
        f"RMS = {rms:.2f}\u03c3",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        color="#555555",
    )

    ax.set_xlim(0.78, 5.15)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("Observed wavelength [\u00b5m]", fontsize=12)
    ax.set_ylabel("Flux density [\u00b5Jy]", fontsize=12)
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.tick_params(colors="#1b1b1b", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#222222")
        spine.set_linewidth(0.9)

    fig.tight_layout(pad=0.55)
    fig.savefig(IMAGE_PATH, facecolor="white")
    plt.close(fig)


def main() -> None:
    if not PHOTOMETRY_PATH.exists():
        raise FileNotFoundError(PHOTOMETRY_PATH)
    if not THROUGHPUT_PATH.exists():
        raise FileNotFoundError(THROUGHPUT_PATH)

    throughput = pd.read_csv(THROUGHPUT_PATH)
    filt_paths = write_bagpipes_filters(throughput)
    wave_um, obs_ujy, err_ujy = read_region_photometry()
    model, params, rms = fit_bagpipes_model(filt_paths, obs_ujy, err_ujy)
    write_model_table(wave_um, obs_ujy, err_ujy, np.asarray(model.photometry), rms)
    plot_sed(throughput, wave_um, obs_ujy, err_ujy, model, rms)

    print(f"Wrote {IMAGE_PATH}")
    print(f"Wrote {TABLE_PATH}")
    print("Best-fit parameters:", params)
    print(f"RMS = {rms:.3f} sigma")


if __name__ == "__main__":
    main()
