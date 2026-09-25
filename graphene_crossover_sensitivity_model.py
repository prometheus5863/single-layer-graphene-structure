"""
graphene_crossover_sensitivity_model.py

HOW FAR CAN THE p/n CROSSOVER MOVE BEFORE CHAPTER 6'S CONCLUSIONS CHANGE?

Both papers this repo takes `W_CROSS_CHEM = 5.4 eV` from write it with a
tilde:

  Giovannetti et al., Phys. Rev. Lett. 101, 026803 (2008)  [arXiv:0802.2267]
  Khomyakov et al.,  Phys. Rev. B   79, 195425 (2009)      [arXiv:0902.1203]

    "the crossover from p-type to n-type doping occurs for a metal work
     function of ~5.4 eV"

Since 2026-09-18 this repo has read that tilde as three significant figures.
This module puts an unknown scalar offset `delta` on the crossover and asks,
for all 21 asymmetric metal pairs, how hard each one amplifies it.

WHY A SCALAR OFFSET, GIVEN THAT 2026-09-21 BUILT A PER-METAL CROSSOVER
---------------------------------------------------------------------
Because the per-metal form cannot reach the pairs that matter.  It needs the
metal-graphene separation d_eq, tabulated for only Cu, Au and Pt of this
repo's seven metals (Ti/Ni/Pd are chemisorbed, where Section 6.12 Result 2
showed the anchored exponential diverges; Cr has no tabulated d_eq).  Every
one of the five headline pairs of the 2026-09-20 retraction table -- Au/Pd,
Ti/Cr, Ni/Au, Cr/Cu, Ni/Pd -- contains at least one metal outside that set.

A scalar offset has the wrong SHAPE (it is not per-metal) but the right
REACH (it needs no d_eq).  The two are complementary; neither supersedes the
other.  See notes/2026-09-22-how-far-can-the-crossover-move.md Section 1.

The offset is also algebraically identical to a common shift of every entry
of METAL_WORK_FUNCTIONS, since dW_i = W_i - w_cross.  So the same sweep
bounds a systematic error in the work-function table -- a live worry, since
Ni's entry already carries a 4.9-5.35 eV spread in its own comment.

PRE-REGISTRATION
----------------
Four predictions (P1-P4) and one derivation (D5) were committed to git in
notes/2026-09-22-how-far-can-the-crossover-move.md BEFORE this file existed
(note in 235b420).  check_predictions() scores them mechanically.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS, W_GRAPHENE
from graphene_photodetector_signed_carrier_model import (
    net_response, signed_offset, W_CROSS_CHEM, N_POINTS,
)
from graphene_photodetector_nonuniform_illumination_model import (
    collection_kernel, g_uniform, g_shadow_mask,
    response_under, response_per_incident,
)

# ---------------------------------------------------------------------
# The perturbation
# ---------------------------------------------------------------------
DELTA_NOMINAL = 0.2      # eV -- band adopted in Section 2 of the note
DELTA_OUTER = 0.5        # eV -- outer sweep
H_DIFF = 1e-3            # eV -- central-difference step for the derivative
#
# H_DIFF IS 3.5 DECADES TOO COARSE.  Measured 2026-09-25 by
# graphene_default_scale_audit.py (RESULT 3): the h -> 0 plateau of
# d ln|N| / d delta runs from ~1e-10 to ~3e-7 eV, and 1e-3 sits well above its
# top.  The value is LEFT UNCHANGED so that every number this module has
# already published stays reproducible and can be compared against rather than
# silently replaced -- see the repo's standing rule about superseded numbers.
# Use `sensitivity_converged()` below for the converged derivative, and
# `H_DIFF_CONVERGED` for the step it uses.
#
# The consequences, measured: per-pair |S| moves by up to 44% (Ti-Ni), and
# prediction P2 ("the largest |S| at delta = 0.1 belongs to Au-Pd"), scored
# PASS on 2026-09-21, is FALSIFIED -- converged, the pair is Ti-Cu.  P1's
# substance (the straddling/same-side separation) survives: 30.4 shipped,
# 26.7 converged, against a claimed factor of 10.
H_DIFF_CONVERGED = 1e-8  # eV -- inside the measured plateau, 2026-09-25
#
# 1e-8 AND NOT 1e-7, AND THE REASON IS A MEASUREMENT.  A first attempt used
# 1e-7 with a 1e-6 convergence requirement and RAISED on two pairs: Cr-Ni at
# 3.98e-5 and Cu-Au at 5.07e-5, against 1e-9..1e-7 for the other nineteen.
# The ratio is exactly second-order -- both drop by 100x when h drops by 10x --
# so it is ordinary truncation and those two pairs simply have a third
# derivative ~100x the rest.  At 1e-8 the worst estimate over all 21 pairs is
# 8.8e-7.  THE CONSEQUENCE IS A DESIGN POINT: no single global step is right
# for every pair, and the only reason this is visible at all is that
# sensitivity_converged() returns its convergence estimate instead of
# asserting one.

# Metals, low work function first, so pair ordering is deterministic.
METALS = sorted(METAL_WORK_FUNCTIONS, key=lambda m: METAL_WORK_FUNCTIONS[m])


def all_asymmetric_pairs():
    """The 21 unordered pairs (A, B) with W_A < W_B."""
    out = []
    for i, a in enumerate(METALS):
        for b in METALS[i + 1:]:
            out.append((a, b))
    return out


def w_cross_of(delta):
    return W_CROSS_CHEM + delta


def straddles(pair, delta=0.0):
    """True when the two contacts sit on opposite sides of the crossover."""
    wc = w_cross_of(delta)
    a, b = pair
    return (METAL_WORK_FUNCTIONS[a] - wc) * (METAL_WORK_FUNCTIONS[b] - wc) < 0.0


def N_of(pair, delta=0.0, e_applied=0.0):
    """Signed net charge to contact A per absorbed photon, uniform light."""
    a, b = pair
    return net_response(METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b],
                        e_applied=e_applied, w_cross=w_cross_of(delta))


def mean_abs_dW(pair, delta=0.0):
    """<|dW|>, the input magnitude the amplification ratio is taken against."""
    wc = w_cross_of(delta)
    a, b = pair
    return 0.5 * (abs(signed_offset(METAL_WORK_FUNCTIONS[a], wc))
                  + abs(signed_offset(METAL_WORK_FUNCTIONS[b], wc)))


def sensitivity(pair, h=H_DIFF):
    """
    THE DEFAULT STEP IS KNOWN TOO COARSE -- see H_DIFF above and
    `sensitivity_converged()`.  Kept at 1e-3 deliberately: this function's
    outputs are quoted in Chapter 6 and in the 2026-09-21/09-23 prediction
    records, and changing the default would rewrite them in place.

    S = d ln|N| / d delta  [1/eV]   and
    A = (d|N|/|N|) / (d<|dW|>/<|dW|>)   [dimensionless, 2026-09-21's ratio]

    Central difference at delta = 0.  Returns (S, A, N0).
    """
    Np, Nm = N_of(pair, +h), N_of(pair, -h)
    N0 = N_of(pair, 0.0)
    dlnN = (abs(Np) - abs(Nm)) / (2.0 * h) / abs(N0)
    wp, wm, w0 = mean_abs_dW(pair, +h), mean_abs_dW(pair, -h), mean_abs_dW(pair)
    dlnW = (wp - wm) / (2.0 * h) / w0
    A = abs(dlnN / dlnW) if dlnW != 0.0 else np.inf
    return dlnN, A, N0


def sensitivity_converged(pair, h=H_DIFF_CONVERGED, rtol=1e-3, decades=2):
    """
    The same quantities as `sensitivity()`, at a step inside the measured
    plateau, with TWO convergence diagnostics returned alongside -- and only
    one of them is allowed to be the pass criterion.

    Returns (S, A, N0, conv_local, conv_plateau).

        conv_local   = |S(h) - S(2h)| / |S(h)|            REPORTED ONLY
        conv_plateau = max over k = 1..decades of
                       |S(h) - S(h/10**k)| / |S(h)|       THE PASS CRITERION

    WHY THE OBVIOUS CHECK IS NOT THE CRITERION, MEASURED 2026-09-25.  The
    first version of this function used step-doubling alone, the pattern
    `graphene_sensitivity_audit.log_sensitivity` has used since 2026-09-22.
    At the shipped H_DIFF = 1e-3 that estimate is BLIND:

        pair    true error vs converged S     step-doubling estimate
        Ti-Ni              44.38 %                   2.33e-06
        Cu-Ni              10.16 %                   5.48e-05
        Cr-Pd              31.99 %                   3.96e-04

    Ti-Ni under-reports its own error by a factor of 190 000.  The mechanism
    is not mysterious: step-doubling measures dE/d(log h) of the error curve
    E(h) = S(h) - S(0), so it reads zero wherever that curve is STATIONARY,
    and S(h) for Ti-Ni is flat to four digits from h = 1e-3 to 1e-2 while
    sitting 44% away from its limit.  A shipped default has no reason to avoid
    a stationary point of its own error curve, and if it lands on one, every
    local self-consistency test certifies it.

    So the criterion walks DOWNWARD instead: a step is accepted only if the
    derivative is unchanged at h/10 and h/100.  That is an anchored
    comparison, not a self-consistency one, and it is the thing that fires.
    `rtol = 1e-3` IS PLACED IN A MEASURED GAP, not chosen for roundness --
    which matters in a function written on the day this repo learned that a
    threshold is a default and a default is a claim about scale.  The two
    populations it has to separate were measured over all 21 pairs:

        worst conv_plateau at H_DIFF_CONVERGED (must pass)   9.60e-05
        cheapest conv_plateau at H_DIFF among the 11 pairs
        whose |S| moves by more than 5%      (must fire)     6.36e-02

    -- a clean separation of 2.82 decades, and 1e-3 sits inside it with an
    order of margin below the firing population and an order above the passing
    one.  At that value: 11 of 11 moved pairs fire at H_DIFF, zero escape, and
    all 21 pass at H_DIFF_CONVERGED.  The residual 9.6e-05 belongs to Pd-Pt,
    a straddling pair whose |S| is 0.0077, where walking two decades down to
    h = 1e-10 reaches the cancellation floor rather than any truncation.

    Raises rather than returning a flag, so an unconverged derivative cannot
    be received silently.
    """
    def _S(hh):
        Np, Nm = N_of(pair, +hh), N_of(pair, -hh)
        return (abs(Np) - abs(Nm)) / (2.0 * hh) / abs(N_of(pair, 0.0))

    S = _S(h)
    denom = abs(S) if S != 0.0 else 1.0
    conv_local = abs(S - _S(2.0 * h)) / denom
    conv_plateau = max(abs(S - _S(h / 10.0 ** k)) / denom
                       for k in range(1, decades + 1))
    if conv_plateau > rtol:
        raise ValueError(
            f"d ln|N|/d delta for {pair} is not on a plateau at h={h:g}: "
            f"worst deviation over {decades} decades below is "
            f"{conv_plateau:.3e} > rtol={rtol:g} "
            f"(the step-doubling estimate here is {conv_local:.3e}, which is "
            f"why step-doubling is not the criterion)")
    wp, wm, w0 = mean_abs_dW(pair, +h), mean_abs_dW(pair, -h), mean_abs_dW(pair)
    dlnW = (wp - wm) / (2.0 * h) / w0
    A = abs(S / dlnW) if dlnW != 0.0 else np.inf
    return S, A, N_of(pair, 0.0), conv_local, conv_plateau


# =====================================================================
# Validations -- four, all against EXACTLY known values.
# =====================================================================
def validate_shift_invariance(deltas=None, verbose=True):
    """
    EXACT #1 (up to float rounding, and the rounding is PREDICTED).

    Moving the crossover by delta and moving every work function by -delta
    are the same perturbation, because only dW = W - w_cross enters:

        N(W_A, W_B ; 5.4 + delta)  ==  N(W_A - delta, W_B - delta ; 5.4)

    Exact in real arithmetic.  NOT expected bitwise in floating point:
    (4.33) - (5.4 + 0.1) and (4.33 - 0.1) - 5.4 round differently.  The note
    predicted <= 1e-15 absolute in advance; claiming bitwise here would have
    been wrong, and saying so before running is the point.
    """
    if deltas is None:
        deltas = np.linspace(-DELTA_OUTER, DELTA_OUTER, 11)
    worst, worst_case, bitwise = 0.0, None, 0
    total = 0
    for pair in all_asymmetric_pairs():
        a, b = pair
        WA, WB = METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b]
        for d in deltas:
            lhs = net_response(WA, WB, e_applied=0.0, w_cross=W_CROSS_CHEM + d)
            rhs = net_response(WA - d, WB - d, e_applied=0.0,
                               w_cross=W_CROSS_CHEM)
            dev = abs(lhs - rhs)
            total += 1
            if dev == 0.0:
                bitwise += 1
            if dev > worst:
                worst, worst_case = dev, (a, b, float(d), lhs, rhs)
    if verbose:
        print("\nValidation 1 (EXACT in exact arithmetic): crossover shift == "
              "common work-function shift")
        print(f"  {total} (pair, delta) combinations; bitwise identical in "
              f"{bitwise}/{total}")
        if worst_case:
            a, b, d, l, r = worst_case
            print(f"  worst: {a}/{b} at delta={d:+.3f}  "
                  f"N={l:+.12f} vs {r:+.12f}")
        print(f"  worst |lhs - rhs|: {worst:.3e}   "
              f"(predicted <= 1e-15 in the note, BEFORE running)")
    return worst, bitwise, total


def validate_delta_zero_reduction(verbose=True):
    """
    EXACT #2.  At delta = 0 this module must reproduce the 2026-09-18 signed
    model BITWISE -- it is the identical call, so anything but 0.0 means the
    wrapper is not neutral.
    """
    worst = 0.0
    for pair in all_asymmetric_pairs():
        a, b = pair
        ref = net_response(METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b],
                           e_applied=0.0, w_cross=W_CROSS_CHEM)
        worst = max(worst, abs(N_of(pair, 0.0) - ref))
    if verbose:
        print("\nValidation 2 (EXACT, bitwise): delta = 0 reproduces the "
              "2026-09-18 signed model")
        print(f"  21 asymmetric pairs; worst |difference|: {worst:.3e} "
              f"({'BITWISE' if worst == 0.0 else 'NOT BITWISE'})")
    return worst


def validate_symmetric_at_every_delta(deltas=None, verbose=True):
    """
    EXACT #3.  Identical contacts at zero bias give an antisymmetric E(x) for
    ANY offset whatever, so N == 0 must hold across the whole sweep, not just
    at delta = 0 where 2026-09-18 checked it.
    """
    if deltas is None:
        deltas = np.linspace(-DELTA_OUTER, DELTA_OUTER, 11)
    worst, worst_case = 0.0, None
    for m in METALS:
        W = METAL_WORK_FUNCTIONS[m]
        for d in deltas:
            N = net_response(W, W, e_applied=0.0, w_cross=W_CROSS_CHEM + d)
            if abs(N) > worst:
                worst, worst_case = abs(N), (m, float(d))
    if verbose:
        print("\nValidation 3 (EXACT): symmetric pair -> N == 0 at EVERY delta")
        print(f"  {len(METALS)}x{len(deltas)} cases; worst |N|: {worst:.3e}"
              + (f"  ({worst_case[0]} at delta={worst_case[1]:+.2f})"
                 if worst_case else ""))
    return worst


def validate_charge_conjugation_at_every_delta(deltas=None, verbose=True):
    """
    EXACT #4.  Reflecting every work function through the SHIFTED crossover,
    W -> 2*(5.4 + delta) - W, negates both offsets, which flips E -> -E and
    exchanges the two carriers exactly:  N(-dW) == -N(dW).

    2026-09-18 got exactly 0.000e+00 on 27 cases at delta = 0.  It must stay
    exact as the crossover moves; a wrapper that leaked a magnitude or mixed
    the carriers would fail here while passing Validations 2 and 3.
    """
    if deltas is None:
        deltas = np.linspace(-DELTA_OUTER, DELTA_OUTER, 11)
    worst, worst_case, exact, total = 0.0, None, 0, 0
    for d in deltas:
        wc = W_CROSS_CHEM + d
        for a in METALS:
            for b in METALS:
                WA, WB = METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b]
                N = net_response(WA, WB, e_applied=0.0, w_cross=wc)
                Nf = net_response(2 * wc - WA, 2 * wc - WB,
                                  e_applied=0.0, w_cross=wc)
                dev = abs(N + Nf)
                total += 1
                if dev == 0.0:
                    exact += 1
                if dev > worst:
                    worst, worst_case = dev, (a, b, float(d))
    if verbose:
        print("\nValidation 4 (EXACT): charge conjugation at EVERY delta")
        print(f"  {total} cases; exactly 0.000e+00 in {exact}/{total}")
        print(f"  worst |N + N_flip|: {worst:.3e}"
              + (f"  ({worst_case[0]}/{worst_case[1]} at "
                 f"delta={worst_case[2]:+.2f})" if worst_case else ""))
    return worst, exact, total


# =====================================================================
# Results
# =====================================================================
def sensitivity_table(verbose=True):
    """
    RESULT 1.  S and A for all 21 asymmetric pairs, split by sign structure.
    This is the 21-pair test of the 65-fold gap 2026-09-21 measured on three.
    """
    rows = []
    for pair in all_asymmetric_pairs():
        S, A, N0 = sensitivity(pair)
        rows.append((pair, straddles(pair), S, A, N0))
    rows.sort(key=lambda r: -r[3])
    strad = [r for r in rows if r[1]]
    same = [r for r in rows if not r[1]]
    if verbose:
        print("\nRESULT 1: crossover sensitivity of all 21 asymmetric pairs "
              "(zero bias, uniform light)")
        print(f"{'pair':<10}{'sign':>10}{'N(0)':>12}"
              f"{'S=dln|N|/dd':>14}{'A (ratio)':>12}")
        for (a, b), st, S, A, N0 in rows:
            print(f"{a+'/'+b:<10}{'straddle' if st else 'same-sign':>10}"
                  f"{N0:>+12.4f}{S:>+14.4f}{A:>12.3f}")
        mx_s = max(r[3] for r in strad)
        mn_m = min(r[3] for r in same)
        print(f"\n  straddling pairs ({len(strad)}): A in "
              f"[{min(r[3] for r in strad):.3f}, {mx_s:.3f}]")
        print(f"  same-sign  pairs ({len(same)}): A in "
              f"[{mn_m:.3f}, {max(r[3] for r in same):.3f}]")
        print(f"  separation ratio min(same-sign)/max(straddling) = "
              f"{mn_m / mx_s:.2f}")
        print("  2026-09-21 measured 4.55 (Cu/Au) vs 0.07-0.11 (Cu/Pt, Au/Pt)")
        print("  from a sample of THREE pairs under a PER-METAL perturbation.")
    return rows


def sign_invariance_scan(lo=-3.0, hi=3.0, n=601, verbose=True):
    """
    RESULT 2, and the FALSIFICATION of derivation D5.

    D5 (pre-registered, note 235b420) asserted that N changes sign exactly
    when w_cross moves BETWEEN a pair's two work functions, because that is
    the only way sign(dW_A * dW_B) can change, and derived a table of flip
    windows from arithmetic on METAL_WORK_FUNCTIONS alone.

    Both halves are wrong, and the first draft of this function hid it.

    THE BUG, RECORDED BECAUSE IT IS THE INSTRUCTIVE PART.  That draft
    bisected between delta = 0 and the midpoint of D5's window WITHOUT ever
    checking that the root was bracketed.  When the sign does not change,
    such a loop walks its lower bound up to the upper bound and returns the
    ENDPOINT.  It printed 21 roots.  Every one of them equalled
    (W_A + W_B)/2 - 5.4 -- the bracket endpoint -- and every one of them
    looked like a physical result, including a "nearest sign flip at
    -0.015 eV" that would have gone straight into the thesis.  It was caught
    by hand-checking one row against what the field actually does at
    dW_A = -dW_B, not by any check in the code.  A bracket test is now the
    first thing this function does.

    WHAT IS ACTUALLY TRUE.  Write dW_A = m - s, dW_B = m + s with
    s = (W_B - W_A)/2 and m = (W_A + W_B)/2 - w_cross.  A uniform crossover
    offset moves m and leaves s EXACTLY fixed, so

        E(x) = m [f(L-x) - f(x)]  +  s [f(x) + f(L-x)],   f(x) = lambda^-1
                                                          (1 + x/lambda)^-2

    splits into an antisymmetric part carrying m -- which contributes zero to
    N by the same argument as Validation 3 -- and a symmetric part carrying
    s, which does not change at all.  The sign of N is therefore set by
    sign(s), i.e. by which contact has the higher work function, and no
    scalar offset can touch it.  m only rescales |N|, going to zero as
    |m| -> infinity: that limit IS the near-cancellation of Section 6.11.2.

    So the scan below is a test of an exactly known answer: ZERO sign changes.
    """
    deltas = np.linspace(lo, hi, n)
    changes, total, smallest = [], 0, (np.inf, None)
    for pair in all_asymmetric_pairs():
        vals = np.array([N_of(pair, float(d)) for d in deltas])
        total += len(vals)
        sgn = np.sign(vals)
        flips = int(np.count_nonzero(sgn[1:] != sgn[:-1]))
        if flips:
            changes.append((pair, flips))
        mn = float(np.min(np.abs(vals)))
        if mn < smallest[0]:
            smallest = (mn, (pair, float(deltas[int(np.argmin(np.abs(vals)))])))
        # every value must carry the sign the s-term predicts
        assert np.all(vals < 0.0), f"{pair}: sign(s) argument violated"
    if verbose:
        print("\nRESULT 2 (falsifies derivation D5): the response sign is "
              "INVARIANT under any scalar crossover offset")
        print(f"  scanned delta in [{lo}, {hi}] eV, {n} points x 21 pairs "
              f"= {total} evaluations")
        print(f"  pairs showing a sign change: {len(changes)}  "
              f"(D5 predicted 21, with the nearest at -0.29 eV)")
        print(f"  smallest |N| anywhere in the scan: {smallest[0]:.4e} "
              f"({smallest[1][0][0]}/{smallest[1][0][1]} at "
              f"delta = {smallest[1][1]:+.3f} eV) -- approached, never crossed")
        print("  Exact reason: dW_A - dW_B = W_A - W_B is invariant under a "
              "uniform offset, and")
        print("  the symmetric part of E(x) that it multiplies is what sets "
              "sign(N).  See docstring.")
    return changes, total, smallest


def validate_straddling_input_invariance(verbose=True):
    """
    EXACT #5, discovered while scoring P1 rather than planned.

    For a pair that STRADDLES the crossover, the mean offset magnitude is

        <|dW|> = ((w_cross - W_A) + (W_B - w_cross)) / 2 = (W_B - W_A) / 2

    which is EXACTLY independent of w_cross.  So the 2026-09-21 amplification
    ratio A = (d|N|/|N|) / (d<|dW|>/<|dW|>) divides by an exactly zero input
    change for every straddling pair: A is not merely large, it is UNDEFINED.

    This is why P1 could not be scored as written, and it is the session's
    main methodological result: A is a property of the (pair, perturbation)
    couple, not of the pair.  Checked here to machine precision over the
    outer sweep.
    """
    worst, worst_case = 0.0, None
    checks = 0
    for pair in all_asymmetric_pairs():
        if not straddles(pair):
            continue
        a, b = pair
        ref = 0.5 * (METAL_WORK_FUNCTIONS[b] - METAL_WORK_FUNCTIONS[a])
        for d in np.linspace(-0.2, 0.2, 21):
            if not straddles(pair, float(d)):
                continue
            dev = abs(mean_abs_dW(pair, float(d)) - ref)
            checks += 1
            if dev > worst:
                worst, worst_case = dev, (a, b, float(d))
    if verbose:
        print("\nValidation 5 (EXACT, found mid-session): a straddling pair's "
              "<|dW|> is independent of w_cross")
        print(f"  {checks} (pair, delta) checks; worst deviation from "
              f"(W_B - W_A)/2: {worst:.3e}")
        if worst_case:
            print(f"  worst: {worst_case[0]}/{worst_case[1]} at "
                  f"delta={worst_case[2]:+.2f}")
        print("  => the 2026-09-21 amplification ratio is UNDEFINED for every "
              "straddling pair under this perturbation.")
    return worst, checks


def mask_gain_under_delta(deltas=(-0.2, -0.1, 0.0, 0.1, 0.2), verbose=True):
    """
    RESULT 3.  Section 6.11.2's table -- best perfect shadow mask scored PER
    INCIDENT photon, divided by |N_uniform| -- recomputed at each crossover
    offset.  The claim under test is "five pairs gain more than 2x".

    The Section 6.9.2 ceiling identity |N[g]| <= max|k| is re-checked here at
    every delta rather than assumed to survive the perturbation.
    """
    pairs = all_asymmetric_pairs()
    table = {}
    ceiling_violations = 0
    ceiling_checks = 0
    for d in deltas:
        wc = W_CROSS_CHEM + d
        row = {}
        for pair in pairs:
            a, b = pair
            WA, WB = METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b]
            x, k, _ = collection_kernel(WA, WB, e_applied=0.0, w_cross=wc)
            un = response_under(k, x, g_uniform(x))
            best = max(abs(response_per_incident(
                k, x, g_shadow_mask(x, T=0.0, masked=m)))
                for m in ("left", "right"))
            ceiling = float(np.max(np.abs(k)))
            ceiling_checks += 2
            for g in (g_uniform(x), g_shadow_mask(x, T=0.0, masked="left")):
                if abs(response_under(k, x, g)) > ceiling * (1 + 1e-12):
                    ceiling_violations += 1
            row[pair] = (un, best,
                         best / abs(un) if abs(un) > 1e-12 else np.inf)
        table[d] = row
    if verbose:
        head = [p for p in pairs
                if table[0.0][p][2] > 2.0 or any(table[d][p][2] > 2.0
                                                 for d in deltas)]
        head.sort(key=lambda p: -table[0.0][p][2])
        print("\nRESULT 3: Section 6.11.2's 'five pairs gain more than 2x', "
              "re-scored as the crossover moves")
        print("  (best perfect mask, PER INCIDENT photon, / |N_uniform|)")
        print(f"{'pair':<10}" + "".join(f"{('d='+format(d,'+.2f')):>12}"
                                        for d in deltas))
        for p in head:
            cells = []
            for d in deltas:
                g = table[d][p][2]
                cells.append("         inf" if not np.isfinite(g)
                             else f"{g:>12.3f}")
            print(f"{p[0]+'/'+p[1]:<10}" + "".join(cells))
        counts = {d: sum(1 for p in pairs if table[d][p][2] > 2.0)
                  for d in deltas}
        print("\n  pairs above 2x: " +
              ", ".join(f"delta={d:+.2f} -> {counts[d]}" for d in deltas))
        print(f"  Section 6.9.2 ceiling |N[g]| <= max|k|: "
              f"{ceiling_violations} violations in {ceiling_checks} checks")
    return table, ceiling_violations, ceiling_checks


def headline_pair_band(pair=("Ti", "Pt"), verbose=True):
    """RESULT 4.  How much does Chapter 6's headline pair move over the band?"""
    N0 = abs(N_of(pair, 0.0))
    ds = np.linspace(-DELTA_NOMINAL, DELTA_NOMINAL, 41)
    vals = np.array([abs(N_of(pair, float(d))) for d in ds])
    rel = (vals - N0) / N0
    if verbose:
        print(f"\nRESULT 4: {pair[0]}/{pair[1]} over the nominal band "
              f"|delta| <= {DELTA_NOMINAL} eV")
        print(f"  |N| at delta=0 : {N0:.4f}")
        print(f"  range over band: {vals.min():.4f} to {vals.max():.4f}")
        print(f"  worst relative change: {np.abs(rel).max()*100:.2f}%  "
              f"(at delta = {ds[np.argmax(np.abs(rel))]:+.3f} eV)")
    return ds, vals, float(np.abs(rel).max())


# =====================================================================
# Mechanical scoring of the pre-registered predictions
# =====================================================================
def check_predictions(verbose=True):
    rows = sensitivity_table(verbose=False)
    strad = [r for r in rows if r[1]]
    same = [r for r in rows if not r[1]]
    # P1 as WRITTEN, in the ratio A
    max_sA = max(r[3] for r in strad)
    min_mA = min(r[3] for r in same)
    p1_written = (max_sA < 1.0) and (min_mA > 2.0) and (min_mA / max_sA > 10.0)
    # P1's SUBSTANCE, in the well-defined statistic S
    max_sS = max(abs(r[2]) for r in strad)
    min_mS = min(abs(r[2]) for r in same)
    p1_substance = min_mS / max_sS > 10.0

    h = H_DIFF

    def S_at(pair, d0):
        Np, Nm = N_of(pair, d0 + h), N_of(pair, d0 - h)
        N0 = N_of(pair, d0)
        return (abs(Np) - abs(Nm)) / (2 * h) / abs(N0)

    at01 = sorted(((abs(S_at(p, 0.1)), p) for p in all_asymmetric_pairs()),
                  reverse=True)
    p2 = at01[0][1] == ("Au", "Pd")

    _, _, worst_rel = headline_pair_band(("Ti", "Pt"), verbose=False)
    p3 = worst_rel < 0.05

    table, _, _ = mask_gain_under_delta(verbose=False)
    base = {p for p in all_asymmetric_pairs() if table[0.0][p][2] > 2.0}
    p4 = any({p for p in all_asymmetric_pairs() if table[d][p][2] > 2.0}
             != base for d in (-0.2, -0.1, 0.1, 0.2))

    changes, _, _ = sign_invariance_scan(verbose=False)
    d5 = len(changes) > 0          # D5 claimed 21 pairs flip

    if verbose:
        print("\n" + "=" * 70)
        print("SCORING THE PRE-REGISTERED PREDICTIONS (note 235b420)")
        print("=" * 70)
        print(f"P1 as written (ratio A), clean, >10x  : "
              f"{'HELD' if p1_written else 'FALSIFIED'}")
        print(f"     max A(straddling) = {max_sA:.4f} (need < 1.0) "
              f"-- UNDEFINED, see Validation 5")
        print(f"     min A(same-sign)  = {min_mA:.4f} (need > 2.0)")
        print(f"P1 substance, in S = dln|N|/d(delta)  : "
              f"{'HELD' if p1_substance else 'FALSIFIED'}")
        print(f"     max |S|(straddling) = {max_sS:.4f} /eV")
        print(f"     min |S|(same-sign)  = {min_mS:.4f} /eV")
        print(f"     separation ratio    = {min_mS / max_sS:.1f}x "
              f"(2026-09-21 measured 65x on three pairs)")
        print(f"P2 largest sensitivity at d=+0.1 is Au/Pd : "
              f"{'HELD' if p2 else 'FALSIFIED'}  (scored in S, since A is "
              f"undefined)")
        print(f"     winner {at01[0][1][0]}/{at01[0][1][1]} "
              f"(|S| = {at01[0][0]:.3f}); runner-up "
              f"{at01[1][1][0]}/{at01[1][1][1]} (|S| = {at01[1][0]:.3f})")
        print(f"P3 Ti/Pt moves < 5% over |delta|<=0.2  : "
              f"{'HELD' if p3 else 'FALSIFIED'}")
        print(f"     worst relative change = {worst_rel*100:.2f}%")
        print(f"P4 the >2x membership changes          : "
              f"{'HELD' if p4 else 'FALSIFIED'}")
        print(f"     at delta=0: {sorted(a+'/'+b for a, b in base)}")
        for d in (-0.2, -0.1, 0.1, 0.2):
            s = {p for p in all_asymmetric_pairs() if table[d][p][2] > 2.0}
            if s != base:
                print(f"     at delta={d:+.2f}: "
                      f"{sorted(a+'/'+b for a, b in s)}")
        print(f"D5 sign flips exist (21 predicted)     : "
              f"{'CONFIRMED' if d5 else 'FALSIFIED'}  "
              f"({len(changes)} pairs flip anywhere in +-3 eV)")
    return {"P1_written": p1_written, "P1_substance": p1_substance,
            "P2": p2, "P3": p3, "P4": p4, "D5": d5}


def make_plots(fname="crossover_sensitivity.png"):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))

    # (a) amplification, straddling vs same-sign
    rows = sensitivity_table(verbose=False)
    labels = [f"{a}/{b}" for (a, b), *_ in rows]
    A = [r[3] for r in rows]
    cols = ["#1f77b4" if r[1] else "#d62728" for r in rows]
    axes[0].bar(range(len(rows)), A, color=cols)
    axes[0].set_yscale("log")
    axes[0].set_xticks(range(len(rows)))
    axes[0].set_xticklabels(labels, rotation=90, fontsize=7)
    axes[0].set_ylabel(r"amplification $A$  (dimensionless)")
    axes[0].set_title("(a) all 21 asymmetric pairs\n"
                      "blue = straddling, red = same-sign")
    axes[0].axhline(1.0, color="k", ls=":", lw=0.8)

    # (b) |N| vs delta for a representative few
    ds = np.linspace(-DELTA_OUTER, DELTA_OUTER, 201)
    for pair, style in ((("Ti", "Pt"), "-"), (("Cu", "Pt"), "-"),
                        (("Au", "Pd"), "--"), (("Ni", "Pd"), "--"),
                        (("Cr", "Cu"), "--")):
        y = [N_of(pair, float(d)) for d in ds]
        axes[1].plot(ds, y, style, lw=1.6, label=f"{pair[0]}/{pair[1]}")
    axes[1].axvspan(-DELTA_NOMINAL, DELTA_NOMINAL, color="0.85", zorder=0)
    axes[1].axhline(0.0, color="k", lw=0.8)
    axes[1].set_xlabel(r"crossover offset $\delta$  (eV)")
    axes[1].set_ylabel(r"$N$ (charge to contact A / absorbed photon)")
    axes[1].set_title("(b) response vs crossover offset\n"
                      "shaded = nominal $|\\delta|\\leq0.2$ eV")
    axes[1].legend(fontsize=8)

    # (c) mask gain of the five headline pairs vs delta
    dlist = np.linspace(-0.25, 0.25, 11)
    table, _, _ = mask_gain_under_delta(deltas=tuple(float(d) for d in dlist),
                                        verbose=False)
    for pair in (("Au", "Pd"), ("Ti", "Cr"), ("Ni", "Au"),
                 ("Cr", "Cu"), ("Ni", "Pd"), ("Ti", "Pt")):
        y = [table[float(d)][pair][2] for d in dlist]
        axes[2].plot(dlist, y, "o-", ms=3, lw=1.4,
                     label=f"{pair[0]}/{pair[1]}")
    axes[2].axhline(2.0, color="k", ls=":", lw=1.0)
    axes[2].set_yscale("log")
    axes[2].set_xlabel(r"crossover offset $\delta$  (eV)")
    axes[2].set_ylabel("best mask / uniform (per incident photon)")
    axes[2].set_title("(c) Section 6.11.2's '>2x' claim\n"
                      "dotted line = the 2x threshold")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"\nwrote {fname}")
    return fname


if __name__ == "__main__":
    print("=" * 70)
    print("CROSSOVER SENSITIVITY: how far can ~5.4 eV move?")
    print("=" * 70)
    validate_shift_invariance()
    validate_delta_zero_reduction()
    validate_symmetric_at_every_delta()
    validate_charge_conjugation_at_every_delta()
    validate_straddling_input_invariance()
    sensitivity_table()
    sign_invariance_scan()
    mask_gain_under_delta()
    headline_pair_band()
    check_predictions()
    make_plots()
