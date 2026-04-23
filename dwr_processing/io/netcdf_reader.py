"""
dwr_processing.io.netcdf_reader
================================

NetCDF ingestion utilities for IMD DWR Level-2 data
distributed via MOSDAC (https://mosdac.gov.in/).
"""

import wradlib as wrl


def read_imd_netcdf(filepath: str) -> dict:
    """
    Read an IMD DWR Level-2 NetCDF file using wradlib's generic reader.

    Parameters
    ----------
    filepath : str
        Absolute or relative path to the .nc file.

    Returns
    -------
    dict
        wradlib generic NetCDF dictionary with keys ``'dimensions'``,
        ``'variables'``, and ``'attributes'``.
    """
    return wrl.io.read_generic_netcdf(filepath)


def detect_polarisation_type(raw: dict) -> str:
    """
    Determine whether a radar volume is single- or dual-polarisation.

    Detection is based on the presence *and* non-emptiness of the three
    mandatory dual-pol variables: ZDR, RHOHV, and PHIDP.

    Parameters
    ----------
    raw : dict
        Output of :func:`read_imd_netcdf`.

    Returns
    -------
    str
        ``'dual_pol'`` or ``'single_pol'``.
    """
    dual_pol_vars = ("ZDR", "RHOHV", "PHIDP")
    variables = raw.get("variables", {})

    if all(v in variables for v in dual_pol_vars):
        if all(variables[v]["data"].size != 0 for v in dual_pol_vars):
            return "dual_pol"

    return "single_pol"


def extract_scan_metadata(raw: dict) -> dict:
    """
    Extract key scan geometry and site metadata from a raw NetCDF dict.

    Parameters
    ----------
    raw : dict
        Output of :func:`read_imd_netcdf`.

    Returns
    -------
    dict
        Dictionary with fields:
        ``site_lat``, ``site_lon``, ``site_alt``,
        ``num_elev``, ``num_bins``, ``num_azim``,
        ``range_resol``, ``nyquist``, ``elevation``, ``sweep``.
    """
    v = raw["variables"]
    d = raw["dimensions"]

    return {
        "site_lat": v["siteLat"]["data"].item(),
        "site_lon": v["siteLon"]["data"].item(),
        "site_alt": v["siteAlt"]["data"].item(),
        "num_elev": d["sweep"]["size"],
        "num_bins": d["bin"]["size"],
        "num_azim": d["radial"]["size"],
        "range_resol": v["gateSize"]["data"],
        "azimuth": v["radialAzim"]["data"],
        "elevation": v["radialElev"]["data"],
        "elevation_angle": v["elevationAngle"]["data"],
        "elevation_list": v["elevationList"]["data"],
        "sweep": v["elevationNumber"]["data"],
        "nyquist": v["nyquist"]["data"],
        "max_range": v["unambigRange"]["data"],
        "first_gate": v["firstGateRange"]["data"],
    }


def get_scan_mode(num_elev: int) -> str:
    """
    Determine the scan strategy from the number of elevation sweeps.

    Parameters
    ----------
    num_elev : int
        Number of elevation angles in the volume scan.

    Returns
    -------
    str
        ``'long_range'`` (≤3 sweeps) or ``'short_range'`` (10 sweeps).

    Raises
    ------
    ValueError
        If ``num_elev`` does not match a known scan mode.
    """
    if num_elev <= 3:
        return "long_range"
    if num_elev == 10:
        return "short_range"
    raise ValueError(
        f"Unrecognised number of elevation sweeps: {num_elev}. "
        "Expected ≤3 (long-range) or 10 (short-range)."
    )
