"""
graphene_differential_crossover_model.py

THE DIFFERENTIAL CROSSOVER OFFSET, AND THE EXACT THRESHOLD AT WHICH A
TWO-TERMINAL PHOTORESPONSE REVERSES SIGN.

2026-09-22 proved that no SCALAR offset on the p/n crossover can change
sign(N): a common offset moves m = (dW_A + dW_B)/2 and leaves
s = (dW_B - dW_A)/2 exactly fixed, and s is what sets the sign.  12621
evaluations over +-3 eV, zero sign changes.

That theorem covers the COMMON MODE of the crossover uncertainty and nothing
else.  Section 6.12 does not predict one crossover, it predicts one PER
METAL, w_cross(d_eq), so the two contacts of a pair do not share a crossover
at all.  Their DIFFERENCE is the part the theorem does not reach, because it
is the only part that moves s.

Parameterise the two contact crossovers as

    w_A = w_cross + c - tau/2 ,      w_B = w_cross + c + tau/2

so that, exactly,

    m -> m - c ,        s -> s - tau/2 .

c is 2026-09-22's scalar offset renamed.  tau is the differential offset.
Together (c, tau) span the whole two-metal crossover uncertainty.

THE RESULT THIS MODULE EXISTS TO TEST (derivation D1, pre-registered in
notes/2026-09-23-the-differential-crossover-and-the-sign-flip-threshold.md,
committed at 7175db7 BEFORE this file):

    N = 0  exactly when  s = 0,  i.e. at   tau* = W_B - W_A          (*)

independently of c, of lambda, of L and of every transport parameter.  The
differential crossover offset needed to reverse a pair's photoresponse is
EXACTLY that pair's work-function gap.

That makes this module's root-finder a root-finder whose answer is known in
advance -- which is the strongest possible test of the bracket-checking
lesson 2026-09-22 paid for.  Both detectors are wired in below:
`bracketed_bisect` refuses to run on an unbracketed interval, and
`validate_flip_threshold_exact` compares what it finds against (*).

NOT CLAIMED
-----------
* tau_model, computed from Section 6.12's one-parameter Delta_c family, is
  NOT a measurement.  It is an anchored exponential with a swept decay
  length and no fitted prefactor (Section 6.12.7).
* Ni, Pd, Ti and Cr get NO tau_model.  Section 6.12 Result 2 refused them
  and this module refuses them too.  For those pairs the output is a
  threshold plus a statement that it is unreachable -- weaker, and honest.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS
from graphene_photodetector_signed_carrier_model import (
    net_response, W_CROSS_CHEM, signed_offset,
)
from graphene_per_metal_crossover_model import (
    net_response_dw, w_cross_for_metal, usable_metals, ELL_RANGE,
)
from bracketed_root import bracketed_bisect

METALS = sorted(METAL_WORK_FUNCTIONS, key=lambda m: METAL_WORK_FUNCTIONS[m])
U_SMALL = 0.05          # eV -- offset from the flip used in the asymmetry test
ELL_GRID = np.linspace(*ELL_RANGE, 25)


def all_asymmetric_pairs():
    out = []
    for i, a in enumerate(METALS):
        for b in METALS[i + 1:]:
            out.append((a, b))
    return out


def dw_pair(pair, c=0.0, tau=0.0):
    """(dW_A, dW_B) under a common offset c and a differential offset tau."""
    a, b = pair
    return (METAL_WORK_FUNCTIONS[a] - (W_CROSS_CHEM + c - tau / 2.0),
            METAL_WORK_FUNCTIONS[b] - (W_CROSS_CHEM + c + tau / 2.0))


def N_of(pair, c=0.0, tau=0.0, e_applied=0.0):
    dA, dB = dw_pair(pair, c, tau)
    return net_response_dw(dA, dB, e_applied=e_applied)


def tau_star(pair):
    """(*) -- the exact flip threshold, from the work functions alone."""
    a, b = pair
    return METAL_WORK_FUNCTIONS[b] - METAL_WORK_FUNCTIONS[a]


# ---------------------------------------------------------------------
# A root-finder that cannot repeat 2026-09-22's bug
# ---------------------------------------------------------------------
# `bracketed_bisect` lived here from 2026-09-23, written after the 2026-09-22
# fault in graphene_crossover_sensitivity_model.py: a bisection that ran with
# no bracket check, walked one bound onto the other when no sign change
# existed, and returned that bound as a root -- 21 times, every one
# physically plausible, while all five of that module's exact validations
# passed and printed underneath. The fault was in the analysis layer, which
# is why model-level validation could not see it.
#
# MOVED 2026-09-24 to bracketed_root.py, unchanged in behaviour, because the
# 2026-09-23 audit item found a second unguarded bisection elsewhere in the
# repo (p2_threshold) and a guard with two copies is a guard that can drift.
# This module's six validations were re-run after the move and are unchanged.
# The move also exposed a scale trap in the `tol` default -- see the note in
# bracketed_root.py; it does not affect this module, whose variable is O(1) eV.


# =====================================================================
# VALIDATIONS -- five, every one against an EXACTLY known value
# =====================================================================
def validate_zero_perturbation(verbose=True):
    """
    EXACT #1, bitwise.  c = tau = 0 must reproduce the 2026-09-18 signed
    model exactly for all 21 pairs -- it is the same arithmetic with the
    subtraction moved outward, so anything but 0.0 means this wrapper is not
    neutral and every number below is suspect.
    """
    worst = 0.0
    for pair in all_asymmetric_pairs():
        a, b = pair
        ref = net_response(METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b],
                           e_applied=0.0, w_cross=W_CROSS_CHEM)
        worst = max(worst, abs(N_of(pair) - ref))
    if verbose:
        print("\nValidation 1 (EXACT, bitwise): (c, tau) = (0, 0) reproduces "
              "the 2026-09-18 signed model")
        print(f"  21 pairs; worst |difference|: {worst:.3e} "
              f"({'BITWISE' if worst == 0.0 else 'NOT BITWISE'})")
    return worst


def validate_c_reduces_to_2026_09_22(cs=None, verbose=True):
    """
    EXACT #2, bitwise.  With tau = 0, the common offset c must be exactly
    2026-09-22's scalar crossover offset: N(pair, c, 0) == the old
    net_response at w_cross = 5.4 + c.  This pins the NEW variable against
    the OLD module rather than against itself.
    """
    if cs is None:
        cs = np.linspace(-0.5, 0.5, 11)
    worst, bitwise, total = 0.0, 0, 0
    for pair in all_asymmetric_pairs():
        a, b = pair
        for c in cs:
            ref = net_response(METAL_WORK_FUNCTIONS[a], METAL_WORK_FUNCTIONS[b],
                               e_applied=0.0, w_cross=W_CROSS_CHEM + float(c))
            dev = abs(N_of(pair, c=float(c)) - ref)
            total += 1
            bitwise += (dev == 0.0)
            worst = max(worst, dev)
    if verbose:
        print("\nValidation 2 (EXACT, bitwise): tau = 0 reproduces "
              "2026-09-22's scalar-offset model")
        print(f"  {total} (pair, c) combinations; bitwise in {bitwise}/{total};"
              f" worst |difference|: {worst:.3e}")
    return worst, bitwise, total


def validate_flip_is_exact_zero(cs=None, verbose=True):
    """
    EXACT #3, and the direct test of derivation D1.  At tau = tau* the two
    contacts have IDENTICAL offsets, the field is antisymmetric about
    mid-channel, and N must be zero to machine precision -- for every pair
    and at every common offset c, since c does not enter (*).
    """
    if cs is None:
        cs = np.linspace(-1.0, 1.0, 9)
    worst, worst_case, total = 0.0, None, 0
    for pair in all_asymmetric_pairs():
        t = tau_star(pair)
        for c in cs:
            N = N_of(pair, c=float(c), tau=t)
            total += 1
            if abs(N) > worst:
                worst, worst_case = abs(N), (pair, float(c))
    if verbose:
        print("\nValidation 3 (EXACT): N == 0 at tau = tau*, at EVERY common "
              "offset c")
        print(f"  {total} (pair, c) cases; worst |N|: {worst:.3e}"
              + (f"  ({worst_case[0][0]}/{worst_case[0][1]} at c="
                 f"{worst_case[1]:+.2f})" if worst_case else ""))
    return worst, total


def validate_charge_conjugation(verbose=True):
    """
    EXACT #4.  Negating both offsets exchanges the two carriers exactly, so
    N(-dW_A, -dW_B) == -N(dW_A, dW_B).  Held to 0.000e+00 on 27 cases at
    delta = 0 (2026-09-18) and on 847 cases under a scalar offset
    (2026-09-22); it must survive a DIFFERENTIAL one, which mixes the
    contacts asymmetrically and is the first perturbation in this chapter
    that could plausibly break it.
    """
    worst, exact, total = 0.0, 0, 0
    for pair in all_asymmetric_pairs():
        for c in (-0.4, 0.0, 0.4):
            for tau in (-0.6, -0.2, 0.0, 0.2, 0.6):
                dA, dB = dw_pair(pair, c, tau)
                N = net_response_dw(dA, dB, e_applied=0.0)
                Nf = net_response_dw(-dA, -dB, e_applied=0.0)
                dev = abs(N + Nf)
                total += 1
                exact += (dev == 0.0)
                worst = max(worst, dev)
    if verbose:
        print("\nValidation 4 (EXACT): charge conjugation under a "
              "DIFFERENTIAL offset")
        print(f"  {total} cases; exactly 0.000e+00 in {exact}/{total}; "
              f"worst |N + N_flip|: {worst:.3e}")
    return worst, exact, total


def validate_flip_threshold_exact(cs=None, verbose=True):
    """
    EXACT #5, and the audit item of 2026-09-22 discharged.

    A bracketed bisection on N(tau) must land on (*) -- a value known in
    closed form before the search runs.  Two independent detectors of the
    same fault class are active here: bracketed_bisect() refuses an
    unbracketed interval, and the located root is compared against
    W_B - W_A.  The 2026-09-22 bug would fail BOTH; it failed neither of the
    checks that module actually had.

    Scored against P1 (predicted <= 1e-12 eV, and <= 1e-12 eV of drift as c
    is swept over [-1, +1]).
    """
    if cs is None:
        cs = np.linspace(-1.0, 1.0, 9)
    worst, worst_case, total = 0.0, None, 0
    drift_worst, drift_case = 0.0, None
    refusals = 0
    for pair in all_asymmetric_pairs():
        t = tau_star(pair)
        roots = []
        for c in cs:
            f = lambda x, p=pair, cc=float(c): N_of(p, c=cc, tau=x)
            lo, hi = t - 0.5 * max(t, 0.05), t + 0.5 * max(t, 0.05)
            root = bracketed_bisect(f, lo, hi)
            roots.append(root)
            dev = abs(root - t)
            total += 1
            if dev > worst:
                worst, worst_case = dev, (pair, float(c))
        drift = max(roots) - min(roots)
        if drift > drift_worst:
            drift_worst, drift_case = drift, pair
        # the guard itself must fire where no root exists
        try:
            bracketed_bisect(lambda x, p=pair: N_of(p, tau=x), t + 1.0, t + 2.0)
        except ValueError:
            refusals += 1
    if verbose:
        print("\nValidation 5 (EXACT): a bracketed bisection lands on "
              "tau* = W_B - W_A")
        print(f"  {total} (pair, c) searches; worst |root - tau*|: "
              f"{worst:.3e} eV"
              + (f"  ({worst_case[0][0]}/{worst_case[0][1]} at c="
                 f"{worst_case[1]:+.2f})" if worst_case else ""))
        print(f"  worst drift of the root as c sweeps [-1, +1] eV: "
              f"{drift_worst:.3e} eV"
              + (f"  ({drift_case[0]}/{drift_case[1]})" if drift_case else ""))
        print(f"  bracket guard fired on {refusals}/21 deliberately "
              f"root-free intervals (2026-09-22's bug would have returned "
              f"an endpoint on every one)")
        print(f"  P1 (both <= 1e-12 eV): "
              f"{'PASS' if max(worst, drift_worst) <= 1e-12 else 'FAIL'}")
    return worst, drift_worst, refusals



def validate_parity_identities(verbose=True):
    """
    EXACT #6, and NOT PLANNED -- found while trying to score P2, exactly as
    2026-09-22's Validation 5 was found while trying to score its P1.

    P2 predicted that |N| would be visibly asymmetric about the flip.  It is
    not: the asymmetry came back 0.00% for all 21 pairs, which is the
    signature of an identity rather than of a coincidence.  The identity is a
    PARITY FACTORISATION in the (m, s) variables:

        N(m, -s) = -N(m, s)        N is ODD  in the differential offset
        N(-m,  s) = +N(m, s)       N is EVEN in the common offset

    Everything this chapter has said about crossover robustness is a corollary:

      * Odd in s  =>  N(m, 0) = 0 exactly, which IS derivation D1 and hence
        tau* = W_B - W_A.  D1 was derived from the antisymmetry of the field;
        it is really the s -> -s parity, which is the stronger statement.
      * Even in m  =>  no common offset of ANY size can change sign(N), which
        is 2026-09-22's result obtained by scanning 12621 points, here as an
        exact parity instead of a scan.
      * Odd in s  =>  |N| depends on s only through |s|, so |N| MUST be
        symmetric about the flip.  P2 and P4 are not merely wrong, they are
        excluded.

    The two parities behave differently under bias, and that difference is
    itself diagnostic:

      * even-in-m is exact BITWISE at every bias, because the (-m) field is
        the literal spatial mirror of the (+m) field and the grid
        linspace(0, L, n) is symmetric, so the mirror is exact ON THE GRID.
      * odd-in-s is exact at zero bias and departs at finite bias by exactly
        2/(n-1) -- one trapezoid cell, the weight being 2 because the
        collection label s jumps from +1 to -1 across it.  Measured at
        n = 501, 1001, 2001, 4001, 8001 the residual is 4.0e-3, 2.0e-3,
        1.0e-3, 5.0e-4, 2.5e-4: exactly 2/(n-1) over a 16x range, i.e. a
        discretisation artefact converging to zero, not a broken symmetry.

    This is the distinction the chapter keeps needing: a residual that scales
    with the grid is the grid, and a residual that does not is physics.
    """
    worst_s, worst_m, total = 0.0, 0.0, 0
    for m in np.linspace(-1.5, 1.5, 13):
        for s in np.linspace(-1.2, 1.2, 13):
            a = net_response_dw(m - s, m + s, e_applied=0.0)
            worst_s = max(worst_s, abs(net_response_dw(m + s, m - s,
                                                       e_applied=0.0) + a))
            worst_m = max(worst_m, abs(net_response_dw(-m - s, -m + s,
                                                       e_applied=0.0) - a))
            total += 1
    grid = []
    for n in (501, 1001, 2001, 4001, 8001):
        w = 0.0
        for m in np.linspace(-1.0, 1.0, 7):
            for s in np.linspace(-1.0, 1.0, 7):
                w = max(w, abs(net_response_dw(m - s, m + s, e_applied=0.2,
                                               n_points=n)
                               + net_response_dw(m + s, m - s, e_applied=0.2,
                                                 n_points=n)))
        grid.append((n, w, w * (n - 1)))
    if verbose:
        print("\nValidation 6 (EXACT, found while scoring P2): the response "
              "factorises by parity")
        print(f"  {total} (m, s) points at zero bias")
        print(f"  odd in s   : worst |N(m,s) + N(m,-s)| = {worst_s:.3e}")
        print(f"  even in m  : worst |N(m,s) - N(-m,s)| = {worst_m:.3e}")
        print("  grid convergence of odd-in-s at 0.2 V bias "
              "(residual x (n-1) should be constant):")
        for n, w, prod in grid:
            print(f"    n_points = {n:5d}   residual = {w:.6e}   "
                  f"x (n-1) = {prod:.4f}")
        print("  => one trapezoid cell, converging as 1/n.  The parity is "
              "exact in the continuum;")
        print("     the grid is what departs from it.")
    return worst_s, worst_m, grid


# =====================================================================
# RESULTS
# =====================================================================
def result_1_margin_table(verbose=True):
    """
    RESULT 1.  The exact sign-flip margin of every pair, from (*), with no
    model evaluation at all -- and the model's own differential offset where
    Section 6.12 can supply one.
    """
    reach = set(usable_metals())
    rows = []
    for pair in all_asymmetric_pairs():
        a, b = pair
        t = tau_star(pair)
        if a in reach and b in reach:
            taus = [w_cross_for_metal(b, e)[0] - w_cross_for_metal(a, e)[0]
                    for e in ELL_GRID]
            lo, hi = min(taus), max(taus)
            biggest = max(abs(lo), abs(hi))
            ratio = t / biggest if biggest else np.inf
            rows.append((pair, t, (lo, hi), ratio, "model"))
        else:
            rows.append((pair, t, None, None,
                         "unreachable (chemisorbed or no d_eq)"))
    rows.sort(key=lambda r: r[1])
    if verbose:
        print("\nRESULT 1: the sign-flip threshold of all 21 pairs "
              "(tau* = W_B - W_A, EXACT)")
        print(f"  {'pair':<10}{'tau* (eV)':>11}{'tau_model range (eV)':>24}"
              f"{'margin':>9}   status")
        for (a, b), t, rng, ratio, st in rows:
            r = (f"[{rng[0]:+.4f}, {rng[1]:+.4f}]" if rng else "--")
            m = (f"{ratio:8.2f}x" if ratio is not None else "      --")
            print(f"  {a+'/'+b:<10}{t:>11.3f}{r:>24}{m}   {st}")
        reachable = [r for r in rows if r[3] is not None]
        print(f"\n  pairs Section 6.12 can reach: {len(reachable)}/21")
        print(f"  the three smallest thresholds: "
              + ", ".join(f"{a}/{b} {t:.2f} eV" for (a, b), t, _, _, _ in rows[:3]))
        print("  every one of those three contains a chemisorbed metal, i.e. "
              "a metal for which")
        print("  Section 6.12 Result 2 REFUSED to return a crossover.  The "
              "pairs most at risk are")
        print("  exactly the pairs the model cannot reach.")
    return rows


def result_2_asymmetry(u=U_SMALL, verbose=True):
    """
    RESULT 2.  Is |N| symmetric about the flip?  Scores P2 and P4.
    """
    rows = []
    for pair in all_asymmetric_pairs():
        t = tau_star(pair)
        hi = abs(N_of(pair, tau=t + u))
        lo = abs(N_of(pair, tau=t - u))
        mean = 0.5 * (hi + lo)
        asym = abs(hi - lo) / mean if mean else 0.0
        a, b = pair
        strad = ((METAL_WORK_FUNCTIONS[a] - W_CROSS_CHEM) *
                 (METAL_WORK_FUNCTIONS[b] - W_CROSS_CHEM)) < 0.0
        rows.append((pair, strad, lo, hi, asym))
    rows.sort(key=lambda r: -r[4])
    st = [r[4] for r in rows if r[1]]
    ss = [r[4] for r in rows if not r[1]]
    if verbose:
        print(f"\nRESULT 2: is |N| symmetric about the flip?  "
              f"(u = {u} eV either side of tau*)")
        print(f"  {'pair':<10}{'sign':>11}{'|N(t*-u)|':>12}{'|N(t*+u)|':>12}"
              f"{'asym':>9}")
        for (a, b), s, lo, hi, asym in rows[:6]:
            print(f"  {a+'/'+b:<10}{'straddle' if s else 'same-sign':>11}"
                  f"{lo:>12.5f}{hi:>12.5f}{asym*100:>8.2f}%")
        print(f"  ... {len(rows) - 6} more")
        big = max(r[4] for r in rows)
        print(f"  largest asymmetry over all 21 pairs: {big:.3e} "
              f"(relative) -- machine precision, not 'small'")
        print(f"  P2 (some pair > 1%): {'PASS' if big > 0.01 else 'FAIL'} "
              f"-- FALSIFIED, and by an identity rather than by a margin")
        if st and ss:
            print(f"  straddling  asym max: {max(st):.3e}   "
                  f"same-sign asym max: {max(ss):.3e}")
            print(f"  P4 (straddling more asymmetric): "
                  f"{'PASS' if min(st) > max(ss) else 'FAIL'} -- FALSIFIED; "
                  f"both are zero")
        print("  Reason (Validation 6): N is exactly ODD in s, so |N| depends "
              "on s only through")
        print("  |s| and MUST be symmetric about the flip.  P2 and P4 were "
              "not long shots that")
        print("  missed; they were excluded by a symmetry the note did not "
              "notice it was assuming")
        print("  away.  The useful device statement: at zero bias the "
              "response MAGNITUDE is blind")
        print("  to the direction of a differential crossover error, while "
              "its SIGN is entirely")
        print("  determined by it.")
    return rows


def result_3_reachable_margins(verbose=True):
    """
    RESULT 3.  The three pairs Section 6.12 can actually reach, scored
    against P3 and P5.
    """
    reach = set(usable_metals())
    rows = []
    for pair in all_asymmetric_pairs():
        a, b = pair
        if a not in reach or b not in reach:
            continue
        t = tau_star(pair)
        taus = np.array([w_cross_for_metal(b, e)[0] - w_cross_for_metal(a, e)[0]
                         for e in ELL_GRID])
        biggest = float(np.max(np.abs(taus)))
        smallest = float(np.min(np.abs(taus)))
        rows.append((pair, t, smallest, biggest, t / biggest))
    all_t = [tau_star(p) for p in all_asymmetric_pairs()]
    global_min_model = min(r[2] for r in rows)
    if verbose:
        print("\nRESULT 3: the three pairs Section 6.12 can reach")
        print(f"  {'pair':<10}{'tau* (eV)':>11}{'|tau_model| min':>16}"
              f"{'max':>9}{'worst margin':>14}")
        for (a, b), t, sm, bg, ratio in rows:
            print(f"  {a+'/'+b:<10}{t:>11.3f}{sm:>16.4f}{bg:>9.4f}"
                  f"{ratio:>13.2f}x")
        worst = min(rows, key=lambda r: r[4])
        print(f"  worst margin: {worst[0][0]}/{worst[0][1]} at "
              f"{worst[4]:.2f}x -- the sign survives, with that much room")
        print(f"  P3 (Cu/Au worst, ratio in [2, 4]): "
              f"{'PASS' if worst[0] == ('Cu', 'Au') and 2.0 <= worst[4] <= 4.0 else 'FAIL'}"
              f"  [worst pair {worst[0][0]}/{worst[0][1]}, ratio {worst[4]:.2f}]")
        print(f"  smallest |tau_model| anywhere: {global_min_model:.4f} eV; "
              f"smallest tau* of all 21 pairs: {min(all_t):.3f} eV")
        print(f"  P5 (no pair below the smallest model offset): "
              f"{'PASS' if min(all_t) > global_min_model else 'FAIL'}")
    return rows


def result_4_conditional_on_the_unreachable(verbose=True):
    """
    RESULT 4, stated as a CONDITIONAL and not as a prediction.

    The three smallest thresholds all involve Ni or Pd, which Section 6.12
    refuses.  What can be said without extrapolating (N1) into the regime
    6.12.4 showed it fails in: compare those thresholds against the range of
    differential offsets the model produces for the metals it DOES admit.
    """
    reach = set(usable_metals())
    obs = []
    for pair in all_asymmetric_pairs():
        a, b = pair
        if a in reach and b in reach:
            obs += [abs(w_cross_for_metal(b, e)[0] - w_cross_for_metal(a, e)[0])
                    for e in ELL_GRID]
    lo, hi = min(obs), max(obs)
    at_risk = [(p, tau_star(p)) for p in all_asymmetric_pairs()
               if tau_star(p) < hi]
    at_risk.sort(key=lambda r: r[1])
    if verbose:
        print("\nRESULT 4 (CONDITIONAL -- not a prediction of reversal)")
        print(f"  Among the metals Section 6.12 admits, the differential "
              f"crossover offset runs")
        print(f"  |tau_model| in [{lo:.4f}, {hi:.4f}] eV.  Pairs whose ENTIRE "
              f"sign margin is smaller than")
        print(f"  the top of that range: {len(at_risk)}/21")
        for p, t in at_risk:
            ch = [m for m in p if m not in reach]
            print(f"    {p[0]+'/'+p[1]:<10}tau* = {t:.3f} eV   "
                  f"(unreachable metal{'s' if len(ch) > 1 else ''}: "
                  f"{', '.join(ch) if ch else 'none'})")
        print("  The conditional: IF Ni and Pd have differential crossover "
              "offsets of the same")
        print("  order as the physisorbed metals do, THEN Au/Pd's sign is not "
              "merely uncertain,")
        print("  it is more likely wrong than right.  The antecedent is "
              "UNTESTED and Section")
        print("  6.12.4 gives reason to think the chemisorbed offsets are "
              "LARGER, not smaller,")
        print("  since d_eq sits ~1.2 A below the anchor.  That direction "
              "strengthens the")
        print("  conditional and does not license dropping it.")
    return at_risk, (lo, hi)


def make_plot(path="differential_crossover.png"):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0))

    # (a) N(tau) through the flip, for the three reachable pairs + Au/Pd
    ax = axes[0]
    show = [("Cu", "Au"), ("Cu", "Pt"), ("Au", "Pt"), ("Au", "Pd")]
    taus = np.linspace(-0.2, 1.2, 281)
    for pair in show:
        vals = [N_of(pair, tau=float(t)) for t in taus]
        line, = ax.plot(taus, vals, lw=1.8, label=f"{pair[0]}/{pair[1]}")
        t = tau_star(pair)
        ax.plot([t], [0.0], "o", ms=6, color=line.get_color())
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"differential crossover offset $\tau$  (eV)")
    ax.set_ylabel(r"net response $N$ (per absorbed photon)")
    ax.set_title(r"(a) the flip is at $\tau^\ast = W_B - W_A$, exactly")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) margin table
    ax = axes[1]
    reach = set(usable_metals())
    pairs = sorted(all_asymmetric_pairs(), key=tau_star)
    ts = [tau_star(p) for p in pairs]
    cols = ["#1f77b4" if (p[0] in reach and p[1] in reach) else "#bbbbbb"
            for p in pairs]
    ax.barh(range(len(pairs)), ts, color=cols)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f"{a}/{b}" for a, b in pairs], fontsize=7)
    obs = []
    for p in pairs:
        if p[0] in reach and p[1] in reach:
            obs += [abs(w_cross_for_metal(p[1], e)[0] - w_cross_for_metal(p[0], e)[0])
                    for e in ELL_GRID]
    ax.axvline(max(obs), color="crimson", ls="--", lw=1.2,
               label=f"largest modelled |$\\tau$| = {max(obs):.3f} eV")
    ax.set_xlabel(r"sign-flip threshold $\tau^\ast$  (eV)")
    ax.set_title("(b) grey = Section 6.12 cannot reach this pair")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, axis="x")

    # (c) the threshold does not move with c
    ax = axes[2]
    cs = np.linspace(-1.0, 1.0, 9)
    for pair in [("Cu", "Au"), ("Au", "Pd"), ("Ti", "Pt")]:
        roots = []
        for c in cs:
            t = tau_star(pair)
            f = lambda x, p=pair, cc=float(c): N_of(p, c=cc, tau=x)
            roots.append(bracketed_bisect(
                f, t - 0.5 * max(t, 0.05), t + 0.5 * max(t, 0.05)) - t)
        ax.plot(cs, np.array(roots) * 1e15, "o-", ms=4,
                label=f"{pair[0]}/{pair[1]}")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"common crossover offset $c$  (eV)")
    ax.set_ylabel(r"$\tau_{\rm found} - \tau^\ast$  (fe$\!$V, $10^{-15}$ eV)")
    ax.set_title("(c) the flip threshold is independent of $c$\n"
                 "(2026-09-22's whole uncertainty, to machine precision)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"\nsaved {path}")
    return path


def main():
    print("=" * 74)
    print("THE DIFFERENTIAL CROSSOVER OFFSET AND THE SIGN-FLIP THRESHOLD")
    print("scoring notes/2026-09-23-...md, committed at 7175db7 before this "
          "file existed")
    print("=" * 74)
    validate_zero_perturbation()
    validate_c_reduces_to_2026_09_22()
    validate_flip_is_exact_zero()
    validate_charge_conjugation()
    validate_flip_threshold_exact()
    validate_parity_identities()
    result_1_margin_table()
    result_2_asymmetry()
    result_3_reachable_margins()
    result_4_conditional_on_the_unreachable()
    make_plot()


if __name__ == "__main__":
    main()
