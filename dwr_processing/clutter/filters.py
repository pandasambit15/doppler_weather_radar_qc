"""
dwr_processing.clutter.filters
================================

Clutter identification and removal for IMD DWR data.

Two strategies are supported:

* **Single-polarisation:** Gabella statistical filter
  (Gabella & Notarpietro, 2002).
* **Dual-polarisation:** Fuzzy-logic polarimetric echo classifier
  (Schuur et al., 2003 / wradlib implementation).

Both routines return a NumPy masked array where non-meteorological
echoes are masked and set to NaN.
"""

import numpy as np
import wradlib as wrl


def gabella_filter(
    rawdbz: np.ndarray,
    wsize: int = 5,
    thrsnorain: float = 0.0,
    tr1: float = 10.0,
    n_p: int = 8,
    tr2: float = 1.3,
) -> np.ndarray:
    """
    Apply the Gabella clutter filter to a reflectivity sweep.

    Suitable for **single-polarisation** radars where only the base
    moments (Z, V, W) are available.

    Parameters
    ----------
    rawdbz : np.ndarray
        Reflectivity array of shape ``(num_azimuth, num_gates)`` in dBZ.
    wsize : int, optional
        Side length (in gates) of the sliding analysis window. Default 5.
    thrsnorain : float, optional
        Reflectivity threshold below which a pixel is considered no-rain.
        Default 0.0 dBZ.
    tr1 : float, optional
        First texture-ratio threshold. Default 10.0.
    n_p : int, optional
        Minimum number of non-rain neighbours required for clutter
        classification. Default 8.
    tr2 : float, optional
        Second texture-ratio threshold. Default 1.3.

    Returns
    -------
    np.ndarray
        Clutter-corrected reflectivity array (same shape as input).
        Flagged gates are set to NaN.
    """
    clutter = wrl.clutter.filter_gabella(
        rawdbz,
        wsize=wsize,
        thrsnorain=thrsnorain,
        tr1=tr1,
        n_p=n_p,
        tr2=tr2,
    )
    mdata = np.ma.array(rawdbz.copy(), mask=clutter)
    mdata[mdata.mask] = np.nan
    return mdata


def polarimetric_filter(
    rawdbz: np.ndarray,
    rhohv: np.ndarray,
    phidp: np.ndarray,
    velocity: np.ndarray,
    zdr: np.ndarray,
    clutter: np.ndarray = None,
    weights: dict = None,
    thresh: float = 0.5,
    despeckle_nrange: int = 5,
) -> np.ndarray:
    """
    Apply a fuzzy-logic polarimetric clutter filter to a reflectivity sweep.

    Uses co-polar correlation coefficient (ρ_hv), differential phase (Φ_DP),
    differential reflectivity (Z_DR), and Doppler velocity as discriminating
    features. If no pre-computed clutter map is supplied, the Gabella filter
    is run first to provide an initial clutter estimate.

    Parameters
    ----------
    rawdbz : np.ndarray
        Reflectivity array (num_azimuth × num_gates), dBZ.
    rhohv : np.ndarray
        Co-polar correlation coefficient ρ_hv (same shape).
    phidp : np.ndarray
        Differential phase Φ_DP in degrees (same shape).
    velocity : np.ndarray
        Radial velocity in m/s (same shape).
    zdr : np.ndarray
        Differential reflectivity Z_DR in dB (same shape).
    clutter : np.ndarray, optional
        Pre-computed boolean clutter map. If ``None``, the Gabella
        filter is applied first.
    weights : dict, optional
        Feature weights for the fuzzy classifier. Defaults to
        ``{"zdr": 0.5, "rho": 0.5, "rho2": 0.3, "phi": 0.3,
           "dop": 0.3, "map": 0.5}``.
    thresh : float, optional
        Fuzzy-score threshold above which a gate is classified as
        clutter. Default 0.5.
    despeckle_nrange : int, optional
        Window size for linear despeckle after classification. Default 5.

    Returns
    -------
    np.ndarray
        Clutter-corrected reflectivity array (same shape as input).
        Clutter gates are set to NaN.
    """
    if weights is None:
        weights = {
            "zdr": 0.5,
            "rho": 0.5,
            "rho2": 0.3,
            "phi": 0.3,
            "dop": 0.3,
            "map": 0.5,
        }

    if clutter is None:
        clutter = wrl.clutter.filter_gabella(
            rawdbz, wsize=5, thrsnorain=0.0, tr1=10.0, n_p=8, tr2=1.3
        )

    # Mask flagged gates
    mdata = np.ma.array(rawdbz.copy(), mask=clutter)
    mdata[mdata.mask] = np.nan

    dat = {
        "rho": rhohv,
        "phi": phidp,
        "ref": mdata,
        "dop": velocity,
        "zdr": zdr,
        "map": clutter,
    }

    cmap, _ = wrl.clutter.classify_echo_fuzzy(dat, weights=weights, thresh=thresh)
    rdata = np.ma.array(mdata.copy(), mask=cmap)
    rdata[rdata.mask] = np.nan
    rdata = wrl.dp.linear_despeckle(rdata, despeckle_nrange)

    return rdata
