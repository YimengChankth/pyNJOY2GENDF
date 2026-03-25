from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


def interp1(x: np.ndarray, y: np.ndarray, xq: np.ndarray | float) -> np.ndarray:
    """1D linear interpolation with clamping to endpoints.

    Like the Fortran MATH_interp1 used in MONTE_extractmacro:
    - inside bounds: linear interpolation
    - out of bounds: clamp to nearest endpoint
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xq_arr = np.asarray(xq, dtype=float)

    # numpy.interp requires increasing x.
    if x.size < 2:
        return np.full_like(xq_arr, y[0] if y.size else 0.0, dtype=float)

    if x[0] <= x[-1]:
        return np.interp(xq_arr, x, y, left=y[0], right=y[-1])

    # Decreasing x: reverse
    xr = x[::-1]
    yr = y[::-1]
    return np.interp(xq_arr, xr, yr, left=yr[0], right=yr[-1])


@dataclass
class IsoTempInterpolated:
    isoname: str
    ng: int
    nsig0: int
    sig0_grid: np.ndarray  # (nsig0,)

    # Arrays after temperature interpolation: (ng, nsig0)
    sigt: np.ndarray
    sigf: np.ndarray
    kerma: np.ndarray
    flux: np.ndarray

    # Scalars/vectors
    nubar: np.ndarray  # (ng,)
    chi: np.ndarray  # (ng,)

    # Elastic scatter (first temp only, but kept here after temp interpolation)
    nl: int
    f_e: np.ndarray
    t_e: np.ndarray
    sige: np.ndarray  # (nl, nonze, nsig0)

    # Inelastic scatter (typically stored without temp axis)
    f_i: np.ndarray
    t_i: np.ndarray
    sigi: np.ndarray  # (nl, nonzi, nsig0)


@dataclass
class IsoSig0Interpolated:
    isoname: str
    ng: int

    # Arrays after sig0 interpolation: (ng,)
    sigt: np.ndarray
    sigf: np.ndarray
    kerma: np.ndarray

    # Scalars/vectors
    nubar: np.ndarray
    chi: np.ndarray

    nl: int
    f_e: np.ndarray
    t_e: np.ndarray
    sige: np.ndarray  # (nl, nonze)

    f_i: np.ndarray
    t_i: np.ndarray
    sigi: np.ndarray  # (nl, nonzi)

    flux: Optional[np.ndarray] = None  # (ng,) optional


def interpolate_temperature(
    *,
    isoname: str,
    temps: np.ndarray,
    sig0_grid: np.ndarray,
    sigt: np.ndarray,
    sigf: np.ndarray,
    kerma: np.ndarray,
    flux: np.ndarray,
    nubar: np.ndarray,
    chi: np.ndarray,
    temp_query: float,
    # scattering
    f_e: np.ndarray,
    t_e: np.ndarray,
    sige: np.ndarray,
    f_i: np.ndarray,
    t_i: np.ndarray,
    sigi: np.ndarray,
) -> IsoTempInterpolated:
    """Interpolate temperature-dependent arrays to a single temperature."""

    temps = np.asarray(temps, dtype=float)
    sig0_grid = np.asarray(sig0_grid, dtype=float)

    ntemp, ng, nsig0 = sigt.shape

    sigt_out = np.zeros((ng, nsig0), dtype=float)
    sigf_out = np.zeros((ng, nsig0), dtype=float)
    flux_out = np.zeros((ng, nsig0), dtype=float)

    for g in range(ng):
        for k in range(nsig0):
            sigt_out[g, k] = float(interp1(temps, sigt[:, g, k], temp_query))
            sigf_out[g, k] = float(interp1(temps, sigf[:, g, k], temp_query))
            flux_out[g, k] = float(interp1(temps, flux[:, g, k], temp_query))

    # Kerma is stored without temperature axis in the original Fortran
    kerma_out = np.asarray(kerma, dtype=float)

    nl = int(sige.shape[1]) if sige.ndim == 4 else int(sige.shape[0])
    # Expected input sige shape: (ntemp, nl, nonze, nsig0)
    if sige.ndim == 4:
        _, nl, nonze, _ = sige.shape
        sige_out = np.zeros((nl, nonze, nsig0), dtype=float)
        for ell in range(nl):
            for j in range(nonze):
                for k in range(nsig0):
                    sige_out[ell, j, k] = float(interp1(temps, sige[:, ell, j, k], temp_query))
    else:
        sige_out = np.asarray(sige, dtype=float)
        nl = sige_out.shape[0]

    # Inelastic scatter is temperature-independent, just pass through
    sigi_out = np.asarray(sigi, dtype=float)

    return IsoTempInterpolated(
        isoname=isoname,
        ng=ng,
        nsig0=nsig0,
        sig0_grid=sig0_grid,
        sigt=sigt_out,
        sigf=sigf_out,
        kerma=kerma_out,
        flux=flux_out,
        nubar=np.asarray(nubar, dtype=float),
        chi=np.asarray(chi, dtype=float),
        nl=nl,
        f_e=np.asarray(f_e, dtype=int),
        t_e=np.asarray(t_e, dtype=int),
        sige=sige_out,
        f_i=np.asarray(f_i, dtype=int),
        t_i=np.asarray(t_i, dtype=int),
        sigi=sigi_out,
    )


def interpolate_sig0(iso_t: IsoTempInterpolated, sig0_by_group: np.ndarray) -> IsoSig0Interpolated:
    """Interpolate all temperature-corrected data to a single sigma0 per group.

    Uses log10(sig0) interpolation, clamped to endpoints.
    """

    sig0_by_group = np.asarray(sig0_by_group, dtype=float)
    if sig0_by_group.shape != (iso_t.ng,):
        raise ValueError(f"sig0_by_group must be shape ({iso_t.ng},), got {sig0_by_group.shape}")

    x = np.log10(np.asarray(iso_t.sig0_grid, dtype=float))

    sigt = np.zeros((iso_t.ng,), dtype=float)
    sigf = np.zeros((iso_t.ng,), dtype=float)
    kerma = np.zeros((iso_t.ng,), dtype=float)
    flux = np.zeros((iso_t.ng,), dtype=float)

    for g in range(iso_t.ng):
        xq = float(np.log10(max(sig0_by_group[g], 1.0)))
        sigt[g] = float(interp1(x, iso_t.sigt[g, :], xq))
        sigf[g] = float(interp1(x, iso_t.sigf[g, :], xq))
        kerma[g] = float(interp1(x, iso_t.kerma[g, :], xq))
        flux[g] = float(interp1(x, iso_t.flux[g, :], xq))

    # Elastic scattering: sige(ell, j, k) -> interpolate over k using group-dependent sig0
    nl, nonze, nsig0 = iso_t.sige.shape
    sige = np.zeros((nl, nonze), dtype=float)
    for ell in range(nl):
        for j in range(nonze):
            g = int(iso_t.f_e[j]) - 1  # 1-based -> 0-based
            xq = float(np.log10(max(sig0_by_group[g], 1.0)))
            sige[ell, j] = float(interp1(x, iso_t.sige[ell, j, :], xq))

    # Inelastic: sigi(ell, j, k) -> use from-group index to pick sig0
    nl_i, nonzi, nsig0_i = iso_t.sigi.shape
    sigi = np.zeros((nl_i, nonzi), dtype=float)
    for ell in range(nl_i):
        for j in range(nonzi):
            g = int(iso_t.f_i[j]) - 1
            xq = float(np.log10(max(sig0_by_group[g], 1.0)))
            sigi[ell, j] = float(interp1(x, iso_t.sigi[ell, j, :], xq))

    return IsoSig0Interpolated(
        isoname=iso_t.isoname,
        ng=iso_t.ng,
        sigt=sigt,
        sigf=sigf,
        kerma=kerma,
        nubar=iso_t.nubar,
        chi=iso_t.chi,
        nl=int(iso_t.nl),
        f_e=iso_t.f_e,
        t_e=iso_t.t_e,
        sige=sige,
        f_i=iso_t.f_i,
        t_i=iso_t.t_i,
        sigi=sigi,
        flux=flux,
    )
