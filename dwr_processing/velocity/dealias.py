"""
dwr_processing.velocity.dealias
================================

Custom bidirectional velocity dealiasing algorithm for IMD DWR data.

Algorithm overview
------------------
1. **Reference radial detection** — Identify a radial near the zero-velocity
   isodop by searching for a sign change in azimuthally-averaged mean
   Doppler velocity across consecutive reliable radials.

2. **First-pass dealiasing** — Propagate phase corrections clockwise and
   anti-clockwise from the reference radial, in both azimuthal and radial
   directions, using the Nyquist folding correction:
       v_corrected = v_observed + k × 2 × v_n
   where k = round((v_ref − v_obs) / (2 × v_n)).

3. **Second-pass dealiasing** — Retry uncorrected gates with an expanded
   search window (M = 3 additional radials).

Key tuning parameters
---------------------
alpha : float
    Velocity shear threshold expressed as a fraction of the Nyquist
    velocity v_n. Gates with shear > alpha × v_n are considered
    unreliable. Default 0.8.
beta : float
    Low-velocity fraction of v_n used for reference radial selection.
    Radials with mean |V| < beta × v_n are candidate reference radials.
    Default 0.3.
"""

import numpy as np


def _find_initial_radial(
    velocity: np.ndarray,
    v_n: float,
    num_azim: int,
    num_bins: int,
    alpha: float,
    beta: float,
) -> tuple:
    """
    Identify the reference (initial) radial for velocity dealiasing.

    The reference radial lies near the zero-velocity isodop — the boundary
    where the mean Doppler velocity transitions between positive and negative.

    Parameters
    ----------
    velocity : np.ndarray
        Raw radial velocity array (num_azimuth × num_gates), m/s.
    v_n : float
        Nyquist velocity, m/s.
    num_azim : int
        Number of azimuths (radials).
    num_bins : int
        Number of range gates.
    alpha : float
        Shear threshold fraction of v_n.
    beta : float
        Low-velocity gate fraction of v_n.

    Returns
    -------
    initial_radial : int or -999
        Index of the reference radial, or -999 if not found.
    v_ref : float or nan
        Mean velocity at the reference radial.
    radial_flag : np.ndarray
        Boolean reliability mask (1 = reliable radial).
    """
    radial_flag = np.ones(num_azim)
    N1 = np.zeros(num_azim)
    v_m1 = np.empty(num_azim)
    v_m1[:] = np.nan
    initial_radial = -999
    v_ref = np.nan

    # --- Pass 1: flag unreliable radials by shear ---
    for a in range(num_azim):
        gate_value = []
        for r in range(num_bins):
            if np.isfinite(velocity[a, r]):
                gate_value.append(velocity[a, r])
            else:
                for m in range(len(gate_value) - 1):
                    if abs(gate_value[m + 1] - gate_value[m]) > alpha * v_n:
                        radial_flag[a] = 0.0
                gate_value = []

    # --- Compute mean low-velocity per reliable radial ---
    ng = 0
    for index in range(num_azim):
        if radial_flag[index] == 1.0:
            good_radial_vel = []
            for x in range(num_bins):
                if np.isfinite(velocity[index, x]) and abs(velocity[index, x]) < beta * v_n:
                    good_radial_vel.append(velocity[index, x])
                    N1[index] += 1
            if good_radial_vel:
                v_m1[index] = np.mean(good_radial_vel)
                ng += 1

    # --- Detect sign-change transition over 4 consecutive reliable azimuths ---
    if ng > 1:
        for l in range(num_azim - 3):
            if radial_flag[l] == 1.0:
                # Positive → Negative transition
                if v_m1[l] >= 0 and v_m1[l + 1] >= 0:
                    if v_m1[l + 2] < 0 and v_m1[l + 3] < 0:
                        if N1[l + 1] > N1[l + 2]:
                            initial_radial = l + 1
                            v_ref = v_m1[l + 1]
                        else:
                            initial_radial = l + 2
                            v_ref = v_m1[l + 2]
                # Negative → Positive transition
                elif v_m1[l] < 0 and v_m1[l + 1] < 0:
                    if v_m1[l + 2] >= 0 and v_m1[l + 3] >= 0:
                        if N1[l + 1] > N1[l + 2]:
                            initial_radial = l + 1
                            v_ref = v_m1[l + 1]
                        else:
                            initial_radial = l + 2
                            v_ref = v_m1[l + 2]

    # --- Fallback: relaxed gate-count search ---
    if initial_radial == -999:
        lim = 40
        while lim >= 5:
            N0 = np.zeros(num_azim)
            v_m0 = np.empty(num_azim)
            v_m0[:] = np.nan
            rad_flag_0 = np.ones(num_azim)
            ng_0 = 0

            for a in range(num_azim):
                nn0 = 0
                gate_value_0 = []
                for r in range(num_bins):
                    if np.isfinite(velocity[a, r]):
                        nn0 += 1
                        gate_value_0.append(velocity[a, r])
                    else:
                        if nn0 >= lim:
                            for m in range(len(gate_value_0) - 1):
                                if abs(gate_value_0[m + 1] - gate_value_0[m]) > alpha * v_n:
                                    rad_flag_0[a] = 0.0
                        gate_value_0 = []

            good_rad_vel0 = []
            for index in range(num_azim):
                if rad_flag_0[index] == 1.0:
                    for x in range(num_bins):
                        if np.isfinite(velocity[index, x]):
                            good_rad_vel0.append(velocity[index, x])
                            N0[index] += 1
                    if good_rad_vel0:
                        v_m0[index] = np.mean(good_rad_vel0)
                    if np.isfinite(v_m0[index]) and v_m0[index] < beta * v_n:
                        ng_0 += 1
                        initial_radial = index
                        v_ref = v_m0[index]

            if ng_0 > 1:
                min_vm0 = v_m0[initial_radial]
                for pos in range(num_azim):
                    if rad_flag_0[pos] == 1.0 and np.isfinite(v_m0[pos]):
                        if v_m0[pos] < min_vm0:
                            min_vm0 = v_m0[pos]
                            initial_radial = pos
                            v_ref = v_m0[pos]

            if initial_radial != -999:
                break
            lim -= 5

    return initial_radial, v_ref, radial_flag


def _correct_initial_radial(
    velocity: np.ndarray,
    corr_vel: np.ndarray,
    gate_flag: np.ndarray,
    loc: int,
    loc1: int,
    loc2: int,
    v_ref: float,
    v_n: float,
    num_bins: int,
) -> None:
    """
    Apply Nyquist correction to the reference radial and its immediate
    azimuthal neighbours (in-place).
    """
    # Reference radial
    for r in range(num_bins):
        if np.isfinite(velocity[loc, r]):
            k0 = np.rint((v_ref - velocity[loc, r]) / (2 * v_n))
            corr_vel[loc, r] = velocity[loc, r] + k0 * 2 * v_n
            gate_flag[loc, r] = 1

    # Immediate neighbours ±1 in azimuth
    for r in range(num_bins):
        if np.isfinite(corr_vel[loc, r]) and np.isfinite(velocity[loc, r]):
            for neighbour in (loc1, loc2):
                if np.isfinite(velocity[neighbour, r]):
                    k = np.rint((corr_vel[loc, r] - velocity[neighbour, r]) / (2 * v_n))
                    if np.isfinite(k):
                        corr_vel[neighbour, r] = velocity[neighbour, r] + k * 2 * v_n
                        gate_flag[neighbour, r] = 1
                else:
                    gate_flag[neighbour, r] = 0


def dealias_velocity(
    velocity: np.ndarray,
    v_n: float,
    range_resolution: float,
    alpha: float = 0.8,
    beta: float = 0.3,
    fill_value: float = None,
) -> np.ndarray:
    """
    Bidirectional radial velocity dealiasing for a single PPI sweep.

    The algorithm identifies a reference radial near the zero-velocity
    isodop and propagates Nyquist-interval corrections clockwise and
    anti-clockwise in both azimuthal and radial directions. A second pass
    resolves residual aliasing using an expanded search window.

    Parameters
    ----------
    velocity : np.ndarray
        Raw radial velocity array of shape ``(num_azimuth, num_gates)``
        in m/s. Missing values should be NaN (or indicated by
        ``fill_value``).
    v_n : float
        Nyquist velocity in m/s.
    range_resolution : float
        Gate spacing in metres (used to compute the 20-km radial
        anchor window).
    alpha : float, optional
        Shear threshold as a fraction of v_n. Default 0.8.
    beta : float, optional
        Low-velocity fraction of v_n for reference radial selection.
        Default 0.3.
    fill_value : float, optional
        If provided, gates equal to ``fill_value`` are replaced with NaN
        before dealiasing.

    Returns
    -------
    np.ndarray
        Dealiased velocity array of shape ``(num_azimuth, num_gates)``.
        Uncorrectable gates remain as NaN.
    """
    velocity = velocity.astype(float).copy()
    if fill_value is not None:
        velocity[velocity == fill_value] = np.nan

    num_azim, num_bins = velocity.shape

    initial_radial, v_ref, _ = _find_initial_radial(
        velocity, v_n, num_azim, num_bins, alpha, beta
    )

    if initial_radial == -999:
        print("  [dealias] No reference radial found — returning original velocities.")
        return velocity

    print(f"  [dealias] Reference radial: azimuth index {initial_radial}")

    # --- Initialise ---
    gate_flag = np.zeros((num_azim, num_bins), dtype=int)
    corr_vel = velocity.copy()

    if initial_radial == 359:
        loc, loc1, loc2 = initial_radial, initial_radial - 1, 0
    else:
        loc, loc1, loc2 = initial_radial, initial_radial - 1, initial_radial + 1

    _correct_initial_radial(velocity, corr_vel, gate_flag, loc, loc1, loc2, v_ref, v_n, num_bins)

    jump = int(np.rint(10_000 / range_resolution))  # ~20-km span in gates

    # -----------------------------------------------------------------------
    # FIRST PASS — azimuthal sweep, both directions
    # -----------------------------------------------------------------------
    for direction in ("cw", "acw"):
        if direction == "cw":
            azimuths = range(loc2 + 1, loc2 + 1 + 180)
        else:
            azimuths = range(loc1 - 1, loc1 - 1 - 180, -1)

        # --- Azimuthal direction ---
        for theta_raw in azimuths:
            theta = theta_raw % num_azim
            good_vel = []
            for r in range(num_bins):
                if np.isfinite(velocity[theta, r]):
                    prev = [(theta - p) % num_azim for p in range(1, 4)]
                    candidates = [corr_vel[p, r] for p in prev if np.isfinite(corr_vel[p, r]) and gate_flag[p, r] == 1]
                    if len(candidates) == 3:
                        if (abs(candidates[1] - candidates[0]) < alpha * v_n
                                or abs(candidates[2] - candidates[1]) < alpha * v_n):
                            v_r = np.mean(candidates)
                            if np.isfinite(v_r) and abs(velocity[theta, r] - v_r) < alpha * v_n:
                                k = np.rint((v_r - velocity[theta, r]) / (2 * v_n))
                                if np.isfinite(k):
                                    corr_vel[theta, r] = velocity[theta, r] + k * 2 * v_n
                                    gate_flag[theta, r] = 1

        # --- Radial direction ---
        initial_gate = -999
        for theta_raw in azimuths:
            theta = theta_raw % num_azim

            for j in range(1, int(num_bins / jump)):
                bin_min = (j - 1) * jump
                bin_max = (j + 1) * jump
                seg_vel = []
                seg_bin = []
                contiguous = True

                for x in range(bin_min, min(bin_max, num_bins)):
                    if np.isfinite(corr_vel[theta, x]):
                        seg_vel.append(corr_vel[theta, x])
                        seg_bin.append(x)
                    else:
                        contiguous = False

                if contiguous and len(seg_vel) >= jump:
                    for vi in range(2, len(seg_vel) - 2):
                        for gi in range(vi - 2, min(vi + 2, len(seg_vel) - 1)):
                            if gate_flag[theta, vi] == 1:
                                if (gate_flag[theta, gi] == 1 and gate_flag[theta, gi - 1] == 1
                                        and abs(seg_vel[gi] - seg_vel[gi - 1]) < alpha * v_n):
                                    prevs = [(theta - p) % num_azim for p in range(1, 4)]
                                    az_consistent = all(
                                        gate_flag[p, vi] == 1 for p in prevs
                                    ) and all(
                                        abs(corr_vel[prevs[k_], vi] - corr_vel[prevs[k_ + 1], vi]) < alpha * v_n
                                        for k_ in range(2)
                                    )
                                    if az_consistent:
                                        initial_gate = seg_bin[vi]

            if initial_gate != -999:
                for pos in range(initial_gate + 1, num_bins):
                    if gate_flag[theta, pos] == 0 and np.isfinite(corr_vel[theta, pos]):
                        gv = [corr_vel[theta, rb] for rb in range(pos - 1, max(initial_gate - 2, -1), -1)
                              if gate_flag[theta, rb] == 1][:3]
                        if gv:
                            k_un = np.rint((np.mean(gv) - corr_vel[theta, pos]) / (2 * v_n))
                            if np.isfinite(k_un):
                                corr_vel[theta, pos] = corr_vel[theta, pos] + k_un * 2 * v_n
                                gate_flag[theta, pos] = 1

                for pos in range(initial_gate - 1, 2, -1):
                    if gate_flag[theta, pos] == 0 and np.isfinite(corr_vel[theta, pos]):
                        gv = [corr_vel[theta, rb] for rb in range(pos + 1, min(initial_gate + 3, num_bins))
                              if gate_flag[theta, rb] == 1][:3]
                        if gv:
                            k_un = np.rint((np.mean(gv) - corr_vel[theta, pos]) / (2 * v_n))
                            if np.isfinite(k_un):
                                corr_vel[theta, pos] = corr_vel[theta, pos] + k_un * 2 * v_n
                                gate_flag[theta, pos] = 1

    # -----------------------------------------------------------------------
    # SECOND PASS — expanded azimuthal search for residual aliasing
    # -----------------------------------------------------------------------
    for direction in ("cw", "acw"):
        if direction == "cw":
            azimuths = range(loc2 + 1, loc2 + 1 + 180)
        else:
            azimuths = range(loc1 - 1, loc1 - 1 - 180, -1)

        for theta_raw in azimuths:
            theta = theta_raw % num_azim
            for r in range(num_bins):
                if np.isfinite(velocity[theta, r]) and gate_flag[theta, r] == 0:
                    M = 3
                    sa = theta - 1
                    ea = sa - 3
                    rounds = 0
                    while rounds <= M:
                        good_vel = [
                            corr_vel[p % num_azim, r]
                            for p in range(sa, ea, -1)
                            if np.isfinite(corr_vel[p % num_azim, r]) and gate_flag[p % num_azim, r] == 1
                        ]
                        v_r = np.mean(good_vel) if good_vel else np.nan
                        if len(good_vel) == 3 and np.isfinite(v_r):
                            if (abs(good_vel[1] - good_vel[0]) < alpha * v_n
                                    and abs(good_vel[2] - good_vel[1]) < alpha * v_n):
                                if abs(velocity[theta, r] - v_r) < alpha * v_n:
                                    k = np.rint((v_r - velocity[theta, r]) / (2 * v_n))
                                    if np.isfinite(k):
                                        corr_vel[theta, r] = velocity[theta, r] + k * 2 * v_n
                                        gate_flag[theta, r] = 1
                                    rounds = M + 1  # exit inner loop
                                    continue
                            sa -= 1
                            ea = sa - 3
                        rounds += 1

    return corr_vel
