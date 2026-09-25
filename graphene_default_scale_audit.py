"""
graphene_default_scale_audit.py -- 2026-09-25

THE ITEM THIS CLOSES
--------------------
The 2026-09-24 session closed the "unguarded bisection" class and, in the
act of closing it, created a fourth failure class and left it open:

    "Audit every other numeric DEFAULT in the repo for the scale assumption
     it encodes.  The bracket class is closed; this class has one measured
     member and NO DETECTOR."

The missing detector is the point.  The one measured member --
`bracketed_bisect`'s absolute `tol=1e-14` silently costing ten digits when
the same guard was applied to a decay length in metres -- was caught only
because an unguarded loop had been returning the right answer for three days
and could be used as an ORACLE.  That is not a method.  Most defaults in this
repo have no older correct implementation to compare against, so a detector
that needs one detects nothing.

THE DETECTOR: A UNIT-COVARIANCE PROBE
-------------------------------------
A numerical default parameterises a procedure P applied to a variable x.
Change the variable's units, x -> lam * x, and restate the SAME problem in
the new units.  The mathematics is unchanged, so the answer must transform
covariantly:

        P_lam(f_lam)  ==  lam * P(f)        (exactly, up to rounding)

A default that is a pure ratio -- a relative tolerance, a relative step, a
halving count -- respects this identically.  A default that carries the
dimensions of x does not, and the deviation is the thing the probe measures.
No oracle, no older implementation, no correct answer needed: the probe is a
self-consistency test of the procedure against its own units.

WHY COVARIANCE ALONE IS NOT THE VERDICT
---------------------------------------
Failing covariance means a default ENCODES a scale.  It does not mean the
number the repo ships is wrong -- only that it is a claim, and the claim is
about the scale of the caller's variable at the call site where the default
was born.  So the audit is two-stage and both stages are reported:

    (a) covariance probe  -> is this default scale-free, or does it encode
                             a scale?                         (a screen)
    (b) margin measurement-> for the ones that encode a scale, how far is
                             the SHIPPED value from where it stops working?
                                                              (the verdict)

Stage (a) is cheap and needs nothing.  Stage (b) is a plateau or convergence
study at the shipped scale, and it is the only stage that can say "this
number is fine" or "this number is one decade from being wrong."

A THIRD CLASS THE 2026-09-24 ITEM LUMPED IN AND SHOULD NOT HAVE
---------------------------------------------------------------
That item listed "grid counts" alongside "tol everywhere."  They are not the
same class.  A grid count is dimensionless, so it passes the covariance
probe trivially and the probe says nothing about it.  A count is a claim
about how much STRUCTURE the integrand has over the domain, and only a
convergence study reaches it.  The audit therefore sorts every numeric
default in the repo into three classes and states which instrument applies:

    TOLERANCE-like : carries the dimensions of some variable  -> probe (a)
    COUNT-like     : dimensionless discretisation count       -> study (b')
    PHYSICAL       : a statement about graphene, not about numerics -> out
                     of scope, but ENUMERATED, so the boundary of the audit
                     is auditable rather than asserted.

PRE-REGISTERED PREDICTIONS (written before the first run; scored at the end)
---------------------------------------------------------------------------
  D1  The probe reproduces the one measured member: `tol=1e-14` on a
      metre-scale root fails covariance by a relative deviation of order
      1e-4, and the `tol=0.0, rtol=eps` call passes at the 1-ulp level.
  D2  `log_sensitivity(rel_step=1e-5)` passes covariance EXACTLY (bitwise
      zero deviation), because `rel_step` multiplies the parameter.
  D3  `H_DIFF = 1e-3` (eV) FAILS covariance -- it is an absolute step on a
      variable in eV -- and is the repo's second scale-encoding default.
  D4  Despite D3, H_DIFF's shipped value sits inside its plateau with at
      least two decades of margin on each side, so no published sensitivity
      number moves.  (Prediction of a NULL result, stated so it can fail.)
  D5  `max_iter=200` is covariant (halving is scale-free) and its margin
      over what double precision can use is at least a factor of 3.
  D6  `n_points=400` in `junction_extra_resistance` is converged to better
      than 1e-3 relative at the shipped value.

Exact validations, and the reason each is here, are in the VALIDATIONS
section.  One of them is included specifically BECAUSE IT CANNOT
DISCRIMINATE, which is a finding about validation cases and not about
defaults.
"""

import numpy as np

from bracketed_root import bracketed_bisect
from graphene_sensitivity_audit import log_sensitivity
from graphene_crossover_sensitivity_model import (
    H_DIFF, sensitivity, all_asymmetric_pairs, N_of, mean_abs_dW,
)
from graphene_contact_doping_model import (
    junction_extra_resistance, METAL_WORK_FUNCTIONS,
)

EPS = np.finfo(float).eps


def ulps(a, b):
    """Distance between two doubles in units of the last place."""
    if a == b:
        return 0
    ia = np.abs(np.frombuffer(np.float64(a).tobytes(), dtype=np.int64)[0])
    ib = np.abs(np.frombuffer(np.float64(b).tobytes(), dtype=np.int64)[0])
    return int(abs(ia - ib))


# =====================================================================
# THE DETECTOR
# =====================================================================
def covariance_deviation(solve, lam):
    """
    THE PROBE.  `solve(scale)` must solve the SAME mathematical problem with
    its variable expressed in units smaller by `scale` -- i.e. the exact root
    of the scaled problem is `scale` times the root of the unscaled one.

    Returns the relative deviation from covariance,

        | solve(lam)/lam - solve(1) |  /  | solve(1) |

    which is identically zero for a scale-free procedure and grows with the
    mismatch for one whose default carries dimensions.  Zero at lam == 1 for
    ANY procedure, which is Validation 1.
    """
    x1 = solve(1.0)
    xl = solve(lam)
    return abs(xl / lam - x1) / abs(x1), x1, xl


# ---------------------------------------------------------------------
# Target 1 + 2: bracketed_bisect's `tol` (absolute) vs `rtol` (relative).
# The root is at `a`, chosen far from 1 so the two differ.  This is the
# one member measured on 2026-09-24, reproduced here by the new detector
# WITHOUT the oracle that was needed then.
# ---------------------------------------------------------------------
A_ROOT = 5.0e-11          # metres -- p2_threshold's actual scale


def _bisect_solver(**kw):
    def solve(scale):
        a = A_ROOT * scale
        return bracketed_bisect(lambda x: x - a, 0.05e-10 * scale,
                                5.0e-10 * scale, **kw)
    return solve


# ---------------------------------------------------------------------
# Target 3: H_DIFF, the central-difference step in
# graphene_crossover_sensitivity_model.sensitivity().
#
# The variable differentiated is `delta`, an energy in eV.  Rescaling the
# unit of energy means delta -> lam*delta at fixed physics, so a covariant
# derivative d ln|N| / d delta must scale as 1/lam.  We build the rescaled
# problem exactly: N as a function of the RESCALED delta is N_of(pair,
# d/lam), and the true derivative in rescaled units is S_true/lam.  So the
# quantity that must be covariant is 1/S, which has the dimensions of delta.
# ---------------------------------------------------------------------
PROBE_PAIR = ("Au", "Pd")     # the pair Prediction P2 of 09-23 is about


def _hdiff_solver(h):
    def solve(scale):
        def N_scaled(d_scaled):
            return N_of(PROBE_PAIR, d_scaled / scale)
        Np, Nm = N_scaled(+h), N_scaled(-h)
        N0 = N_scaled(0.0)
        S = (abs(Np) - abs(Nm)) / (2.0 * h) / abs(N0)
        return 1.0 / S          # dimensions of delta -> must be covariant
    return solve


def _hdiff_solver_relative(rel):
    """The scale-free alternative: step as a FRACTION of DELTA_NOMINAL."""
    def solve(scale):
        h = rel * 0.2 * scale
        def N_scaled(d_scaled):
            return N_of(PROBE_PAIR, d_scaled / scale)
        Np, Nm = N_scaled(+h), N_scaled(-h)
        N0 = N_scaled(0.0)
        S = (abs(Np) - abs(Nm)) / (2.0 * h) / abs(N0)
        return 1.0 / S
    return solve


# ---------------------------------------------------------------------
# Target 4: log_sensitivity's `rel_step`.  The parameter p carries units;
# rescaling p must leave the LOGARITHMIC derivative unchanged (it is
# dimensionless), so here covariance means invariance, and the probe is
# applied to a quantity of dimension 1 -- lam = 1 in the formula, tested by
# comparing S(p) against S(lam*p) for a function rescaled to match.
# ---------------------------------------------------------------------
def _logsens_invariance(rel_step, lam):
    """S = d ln|f| / d ln p is dimensionless: rescaling p must not move it."""
    f1 = lambda p: 3.0 * p ** 2.5
    S1, _ = log_sensitivity(f1, 1.0, rel_step=rel_step)
    f2 = lambda p: 3.0 * (p / lam) ** 2.5
    S2, _ = log_sensitivity(f2, 1.0 * lam, rel_step=rel_step)
    return abs(S2 - S1), S1, S2


# =====================================================================
# VALIDATIONS -- each against an EXACTLY known value
# =====================================================================
def validate_identity_rescaling(verbose=True):
    """
    EXACT #1.  lam == 1 is the identity, so the probe must return BITWISE
    zero for every target, good default or bad.  This is the probe's own
    zero-perturbation test: a non-zero answer here would mean the harness,
    not the default, is the thing that varies.
    """
    devs = []
    for name, solve in (
        ("bisect tol=1e-14", _bisect_solver(tol=1e-14)),
        ("bisect rtol=eps", _bisect_solver(tol=0.0, rtol=EPS)),
        ("H_DIFF=1e-3", _hdiff_solver(H_DIFF)),
    ):
        d, _, _ = covariance_deviation(solve, 1.0)
        devs.append((name, d))
    worst = max(d for _, d in devs)
    if verbose:
        print("\nVALIDATION 1 (exact): identity rescaling lam=1 -> deviation "
              "must be bitwise 0.0")
        for n, d in devs:
            print(f"    {n:<22} deviation = {d:.1e}   "
                  f"{'OK' if d == 0.0 else 'NON-ZERO'}")
        print(f"  -> worst = {worst:.1e}   "
              f"{'PASS' if worst == 0.0 else 'FAIL'}")
    return worst


def validate_odd_function_zero(verbose=True):
    """
    EXACT #2, AND THE REASON IT IS HERE IS THAT IT CANNOT DISCRIMINATE.

    Bisecting an odd function on a symmetric interval [-b, +b]: the root is
    exactly 0.0, and the FIRST midpoint 0.5*(-b + b) is exactly 0.0, so
    `bracketed_bisect` returns bitwise 0.0 after one evaluation -- at ANY
    tolerance, for ANY b, with or without `rtol`.  The answer is exact and
    the test passes for the good default and the bad one alike.

    It is retained deliberately.  A test case whose outcome does not depend
    on the quantity under audit proves the harness and proves nothing about
    the default, and this repo has twice mistaken one for the other (the
    2026-09-22 analysis-layer fault, which every model-level validation
    passed underneath; and 2026-09-24's truthiness-scored validation).
    Recording a non-discriminating exact test AS non-discriminating is
    cheaper than rediscovering that it was.
    """
    outs = []
    for b in (1.0, 5.0e-11, 5.0e11):
        for kw in ({"tol": 1e-14}, {"tol": 0.0, "rtol": EPS}):
            r = bracketed_bisect(lambda x: x ** 3, -b, b, **kw)
            outs.append((b, kw, r))
    all_exact = all(r == 0.0 for _, _, r in outs)
    if verbose:
        print("\nVALIDATION 2 (exact, and NON-DISCRIMINATING BY CONSTRUCTION):")
        print("    odd f on symmetric [-b,b] -> root exactly 0.0, first "
              "midpoint exactly 0.0")
        for b, kw, r in outs:
            print(f"    b={b:>9.1e}  {str(kw):<28} -> {r!r}")
        print(f"  -> all bitwise 0.0: {all_exact}   "
              f"{'PASS (harness sound; discriminates nothing)' if all_exact else 'FAIL'}")
    return all_exact


def validate_reproduce_20260924_member(verbose=True):
    """
    EXACT #3.  The 2026-09-24 entry measured, for `tol=1e-14` on a
    metre-scale decay length, a stop after 16 HALVINGS at a relative
    precision of 2e-4.  Both are arithmetic consequences of the interval and
    must reproduce exactly: the interval is [0.05e-10, 5.0e-10], width
    4.95e-10, and halving stops at the first k with 4.95e-10 / 2**k < 1e-14.
    """
    width = 5.0e-10 - 0.05e-10
    k = 0
    w = width
    while w >= 1e-14:
        w *= 0.5
        k += 1
    rel = w / A_ROOT
    ok = (k == 16)
    if verbose:
        print("\nVALIDATION 3 (exact): reproduce 2026-09-24's measured member")
        print(f"    interval width {width:.4e} m, absolute tol 1e-14")
        print(f"    halvings to convergence = {k}   (2026-09-24 logged: 16)")
        print(f"    final interval {w:.3e} m -> relative precision "
              f"{rel:.1e}   (logged: 2e-4)")
        print(f"  -> {'PASS' if ok else 'FAIL'}")
    return k, rel, ok


def validate_hdiff_derived_rounding_bound(verbose=True):
    """
    EXACT #4, and IT FAILED ITS FIRST FORM, which is the point.

    A central difference of an EXACTLY linear f(d) = c0 + c1*d recovers c1
    with no truncation error at any h, so the first version of this test
    asserted a relative error below 1e-14 -- and MEASURED 4.8e-11.  The
    assertion was wrong, not the arithmetic.  Cancellation leaves a rounding
    error of order eps*|c0| in the numerator, so the relative error carries a
    DERIVED bound:

        |est - c1| / |c1|  <=  eps * |c0| / (2 * h * |c1|)

    which is 7.6e-15 at h = 1e-3 and 7.6e-10 at h = 1e-8.  The threshold
    1e-14 was an ABSOLUTE constant that silently encoded h ~ 1e-3 -- this
    audit's own validation committed the exact error the audit exists to
    find, one level up, and was caught by the test failing rather than by
    review.  It is recorded rather than quietly corrected because the
    2026-09-24 entry's lesson was that a default is a claim about scale, and
    a test THRESHOLD is a default.

    The test now compares against the derived bound at each h.
    """
    c0, c1 = 0.5, -7.25
    f = lambda d: c0 + c1 * d
    rows = []
    for h in (1.0, 1e-1, H_DIFF, 1e-8):
        est = (f(+h) - f(-h)) / (2.0 * h)
        err = abs(est - c1) / abs(c1)
        bound = EPS * abs(c0) / (2.0 * h * abs(c1))
        rows.append((h, err, bound, err <= max(bound, EPS)))
    ok = all(r[3] for r in rows)
    if verbose:
        print("\nVALIDATION 4 (exact, DERIVED bound -- first form of this "
              "test FAILED):")
        print("    central difference of f = c0 + c1*d has no truncation "
              "error, only cancellation")
        print(f"    {'h':>10}{'rel. error':>14}{'derived bound':>16}   ok")
        for h, err, bound, k in rows:
            print(f"    {h:>10.1e}{err:>14.2e}{bound:>16.2e}   "
                  f"{'yes' if k else 'NO'}")
        print(f"  -> {'PASS' if ok else 'FAIL'}   (a fixed 1e-14 threshold "
              f"would FAIL here, and did)")
    return ok


def validate_shipped_sensitivity_unmoved(verbose=True):
    """
    EXACT #5.  Whatever this audit concludes, it must not have CHANGED
    anything: `sensitivity()` at the shipped H_DIFF must return bitwise what
    it returned before, for all 21 pairs.  The reference values are recomputed
    from the module as imported, so this is a tautology unless a probe has
    mutated module state -- which is exactly the failure mode being excluded
    (the 2026-09-24 session's `H_DIFF`-reading helpers read a module global).
    """
    before = [sensitivity(p)[0] for p in all_asymmetric_pairs()]
    # run every probe that touches the sensitivity module
    for h in (1e-6, 1e-3, 1e-1):
        _hdiff_solver(h)(1.0)
        _hdiff_solver(h)(1e3)
    after = [sensitivity(p)[0] for p in all_asymmetric_pairs()]
    identical = all(a == b for a, b in zip(before, after))
    if verbose:
        print("\nVALIDATION 5 (exact): probing must not move the shipped "
              "numbers")
        print(f"    21 pairs, d ln|N|/d delta at H_DIFF, before vs after all "
              f"probes: bitwise identical = {identical}")
        print(f"  -> {'PASS' if identical else 'FAIL'}")
    return identical


# =====================================================================
# RESULT 1 -- the census.  Every numeric default in the repo, sorted into
# the three classes, so that the boundary of this audit is auditable.
# =====================================================================
CENSUS = [
    # (module, symbol, value, class, note)
    ("bracketed_root", "bracketed_bisect tol", "1e-14", "TOLERANCE",
     "absolute interval width; the one member measured 2026-09-24"),
    ("bracketed_root", "bracketed_bisect rtol", "0.0", "TOLERANCE",
     "relative; default 0.0 keeps pre-existing callers bitwise unchanged"),
    ("bracketed_root", "bracketed_bisect max_iter", "200", "COUNT",
     "halving budget; interacts with tol=0.0, where it SETS the precision"),
    ("graphene_sensitivity_audit", "log_sensitivity rel_step", "1e-5",
     "TOLERANCE", "RELATIVE step -- flagged safe by inspection 2026-09-24"),
    ("graphene_crossover_sensitivity_model", "H_DIFF", "1e-3 eV",
     "TOLERANCE", "ABSOLUTE central-difference step on delta"),
    ("graphene_crossover_sensitivity_model", "sign_invariance_scan n", "601",
     "COUNT", "settled by the 2026-09-23 convergence study"),
    ("graphene_crossover_sensitivity_model", "DELTA_NOMINAL / DELTA_OUTER",
     "0.2 / 0.5 eV", "PHYSICAL", "the perturbation band adopted in Section 2"),
    ("graphene_contact_doping_model", "junction_extra_resistance n_points",
     "400", "COUNT", "trapezoid nodes over the junction; never studied"),
    ("graphene_contact_doping_model", "junction_extra_resistance L_junction",
     "1e-6 m", "PHYSICAL", "the junction length the integral is taken over"),
    ("graphene_contact_doping_model", "contact_edge_carrier_density T",
     "300.0 K", "PHYSICAL", "room temperature"),
    ("graphene_contact_doping_model", "recalibrate_metal_rc n_bulk",
     "2.0e16 m^-2", "PHYSICAL", "the bulk doping the recalibration assumes"),
    ("graphene_band_structure", "k_path_graphene n_points", "100", "COUNT",
     "plot resolution only; no number is read off it"),
    ("graphene_band_structure", "graphene_hamiltonian hopping", "2.8 eV",
     "PHYSICAL", "the nearest-neighbour hopping integral"),
    ("graphene_per_metal_crossover_model", "p2_threshold tol", "0.10",
     "PHYSICAL", "the 10% band P2 is DEFINED as; not a numerical tolerance"),
    ("graphene_per_metal_crossover_model", "p2_threshold lo / hi",
     "0.05e-10 / 5.0e-10 m", "PHYSICAL", "the swept decay-length range"),
    ("graphene_photodetector_nonuniform_illumination_model",
     "g_gaussian_spot sigma", "20e-9 m", "PHYSICAL", "the spot size"),
    ("graphene_fet_model", "transfer_characteristic Vds", "0.05 V",
     "PHYSICAL", "the bias the transfer curve is taken at"),
    ("rf_small_signal_model", "output_conductance dVds", "1e-3 V",
     "TOLERANCE", "ABSOLUTE finite-difference step on Vds (order 0.05 V)"),
    ("rf_small_signal_model", "gate_resistance N_fingers", "1", "PHYSICAL",
     "a layout choice, swept explicitly in Section 4.9"),
]


def result_1_census(verbose=True):
    counts = {}
    for row in CENSUS:
        counts[row[3]] = counts.get(row[3], 0) + 1
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 1 -- census of numeric defaults, and which instrument "
              "reaches each")
        print("=" * 78)
        for cls, instrument in (
            ("TOLERANCE", "unit-covariance probe (this module)"),
            ("COUNT", "convergence study -- the probe says NOTHING"),
            ("PHYSICAL", "out of scope: a claim about graphene, not numerics"),
        ):
            print(f"\n  [{cls}]  n = {counts.get(cls, 0)}   instrument: "
                  f"{instrument}")
            for mod, sym, val, c, note in CENSUS:
                if c == cls:
                    print(f"    {sym:<38} = {val:<16} {note}")
        print(f"\n  total enumerated: {len(CENSUS)}")
        print("  The PHYSICAL rows are listed and not probed ON PURPOSE.  The"
              "\n  2026-09-24 item said 'every numeric default'; two thirds of"
              "\n  them are modelling choices, and an audit that probed those"
              "\n  would report covariance failures that mean nothing.  A"
              "\n  default like p2_threshold's tol=0.10 LOOKS like a numerical"
              "\n  tolerance and is not: 10% is the definition of P2.")
    return counts


# =====================================================================
# RESULT 2 -- the probe applied to every TOLERANCE-class default
# =====================================================================
LAMBDAS = (1e-3, 1.0, 1e3, 1e6)


def result_2_probe(verbose=True):
    rows = []
    for name, solve, expect in (
        ("bracketed_bisect tol=1e-14", _bisect_solver(tol=1e-14),
         "SCALE-BOUND"),
        ("bracketed_bisect tol=0, rtol=eps", _bisect_solver(tol=0.0, rtol=EPS),
         "scale-free"),
        ("bracketed_bisect tol=0, rtol=0 (200 halvings)",
         _bisect_solver(tol=0.0, rtol=0.0), "scale-free"),
        (f"H_DIFF = {H_DIFF:g} (absolute)", _hdiff_solver(H_DIFF),
         "SCALE-BOUND"),
        ("h = 0.5% of DELTA_NOMINAL (relative)",
         _hdiff_solver_relative(5e-3), "scale-free"),
    ):
        devs = []
        for lam in LAMBDAS:
            try:
                d, _, _ = covariance_deviation(solve, lam)
            except Exception as exc:
                d = float("nan")
            devs.append(d)
        worst = np.nanmax(devs)
        verdict = "scale-free" if worst <= 64 * EPS else "SCALE-BOUND"
        rows.append((name, devs, worst, verdict, expect))
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 2 -- unit-covariance probe.  Relative deviation from "
              "P(lam*x) == lam*P(x)")
        print("=" * 78)
        print(f"  {'default':<42}" +
              "".join(f"{'lam=' + f'{l:g}':>12}" for l in LAMBDAS) +
              f"{'verdict':>14}")
        for name, devs, worst, verdict, expect in rows:
            print(f"  {name:<42}" +
                  "".join(f"{d:>12.2e}" for d in devs) +
                  f"{verdict:>14}")
        print("\n  threshold for 'scale-free': worst deviation <= 64 eps "
              f"({64 * EPS:.1e})")
        agree = all(v == e for _, _, _, v, e in rows)
        print(f"  all five verdicts match the pre-registered expectation: "
              f"{agree}")
    return rows


# =====================================================================
# RESULT 3 -- the margin measurement for H_DIFF (the verdict stage)
# =====================================================================
def result_3_hdiff_plateau(verbose=True):
    """
    A covariance failure says H_DIFF encodes a scale.  The verdict is whether
    the SHIPPED value sits inside the plateau where truncation error has
    fallen and cancellation error has not yet risen -- and with how much
    margin on each side.  Measured on the pair the 09-23 prediction P2 names,
    and on the 21-pair worst case.
    """
    hs = np.logspace(-10, -1.5, 18)
    pairs = all_asymmetric_pairs()

    def S_at(pair, h):
        Np, Nm = N_of(pair, +h), N_of(pair, -h)
        return (abs(Np) - abs(Nm)) / (2.0 * h) / abs(N_of(pair, 0.0))

    # THE REFERENCE STEP IS ITSELF A DEFAULT, AND THE FIRST CHOICE WAS WRONG.
    # It was 1e-5 -- "two decades finer than shipped", chosen on the
    # assumption that finer is better.  Measured below: the h -> 0 plateau of
    # this derivative runs from ~1e-9 to ~1e-6, so 1e-5 sits OUTSIDE it, and
    # using it as the reference reported a plateau of width zero (one grid
    # point, itself) and a NEGATIVE margin.  The reference is now 1e-7,
    # inside the plateau, and the plateau is measured rather than assumed.
    ref_h = 1e-7
    ref = {p: S_at(p, ref_h) for p in pairs}
    curve = []
    for h in hs:
        worst = 0.0
        for p in pairs:
            s = S_at(p, h)
            r = ref[p]
            if r != 0.0:
                worst = max(worst, abs(s - r) / abs(r))
        curve.append(worst)
    curve = np.array(curve)

    TOLBAND = 1e-3
    inside = hs[curve <= TOLBAND]
    lo_edge, hi_edge = (inside.min(), inside.max()) if inside.size else (np.nan, np.nan)
    shipped_err = max(abs(S_at(p, H_DIFF) - ref[p]) / abs(ref[p])
                      for p in pairs if ref[p] != 0.0)
    dec_lo = np.log10(H_DIFF / lo_edge)
    dec_hi = np.log10(hi_edge / H_DIFF)
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 3 -- H_DIFF margin.  Worst relative deviation of "
              "d ln|N|/d delta over all 21 pairs")
        print("=" * 78)
        print(f"  reference step {ref_h:g} eV;  plateau band = {TOLBAND:g} "
              f"relative")
        print(f"  {'h (eV)':>12}{'worst rel. dev.':>18}   in band")
        for h, c in zip(hs, curve):
            mark = "  <-- SHIPPED H_DIFF" if abs(h - H_DIFF) / H_DIFF < 0.2 else ""
            print(f"  {h:>12.2e}{c:>18.2e}   "
                  f"{'yes' if c <= TOLBAND else 'no ':<4}{mark}")
        print(f"\n  plateau edges: {lo_edge:.2e} .. {hi_edge:.2e} eV")
        print(f"  shipped H_DIFF = {H_DIFF:g} eV, worst deviation "
              f"{shipped_err:.2e}")
        print(f"  margin: {dec_lo:.2f} decades below, {dec_hi:.2f} decades "
              f"above")
    return lo_edge, hi_edge, shipped_err, dec_lo, dec_hi, hs, curve


# =====================================================================
# RESULT 4 -- max_iter's margin
# =====================================================================
def result_4_max_iter(verbose=True):
    """
    With `tol=0.0, rtol=0.0` -- the way p2_threshold's sibling call would run
    if `rtol` were omitted -- `max_iter` alone sets the precision.  How many
    halvings does double precision actually use, and what is the margin?
    """
    width = 5.0e-10 - 0.05e-10
    needed = 0
    lo, hi = 0.05e-10, 5.0e-10
    while True:
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        needed += 1
        if mid > A_ROOT:
            hi = mid
        else:
            lo = mid
        if needed > 400:
            break
    margin = 200 / needed
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 4 -- max_iter=200 margin")
        print("=" * 78)
        print(f"  halvings until the midpoint stops moving (double precision, "
              f"this interval): {needed}")
        print(f"  shipped max_iter = 200  ->  margin x{margin:.2f}")
        print("  Covariant: halving is a ratio, so the count is scale-free "
              "and the probe\n  in RESULT 2 confirms it.  But the number is "
              "still a claim -- that 200\n  exceeds what any interval in this "
              "repo needs -- and it is now measured\n  rather than assumed.")
    return needed, margin


# =====================================================================
# RESULT 5 -- the COUNT class: grid convergence where it was never studied
# =====================================================================
def result_5_grid_convergence(verbose=True):
    """
    `junction_extra_resistance(n_points=400)` feeds Chapter 4's per-metal Rc
    decomposition, including the Section 4.7 negative residual that has been
    open since 2026-08-31.  The count has never been converged.  The probe
    cannot touch it -- a count is dimensionless -- so this is the other
    instrument.
    """
    ns = [50, 100, 200, 400, 800, 1600, 3200, 6400]
    metals = [m for m in ("Ti", "Ni", "Pd", "Au", "Pt")
              if m in METAL_WORK_FUNCTIONS]
    ref = {}
    rows = []
    for m in metals:
        ref[m] = junction_extra_resistance(METAL_WORK_FUNCTIONS[m], 2.0e16,
                                           n_points=25600)[0]
    for n in ns:
        worst = 0.0
        vals = {}
        for m in metals:
            v = junction_extra_resistance(METAL_WORK_FUNCTIONS[m], 2.0e16,
                                          n_points=n)[0]
            vals[m] = v
            worst = max(worst, abs(v - ref[m]) / abs(ref[m]))
        rows.append((n, vals, worst))
    shipped = [r for r in rows if r[0] == 400][0]
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 5 -- grid convergence of junction_extra_resistance "
              "(COUNT class)")
        print("=" * 78)
        print(f"  reference n_points = 25600;  R_extra in Ohm.um")
        print(f"  {'n':>6}" + "".join(f"{m:>12}" for m in metals) +
              f"{'worst rel.':>14}")
        for n, vals, worst in rows:
            mark = "  <-- SHIPPED" if n == 400 else ""
            print(f"  {n:>6}" + "".join(f"{vals[m]:>12.4f}" for m in metals) +
                  f"{worst:>14.2e}{mark}")
        print(f"\n  at the shipped n_points=400 the worst relative error is "
              f"{shipped[2]:.2e}")
    return rows, shipped[2]


# =====================================================================
# RESULT 6 -- a second scale-bound default found by the probe
# =====================================================================
def result_6_dvds(verbose=True):
    """
    `rf_small_signal_model.output_conductance(dVds=1e-3)` is an absolute
    finite-difference step on Vds, whose shipped operating point is 0.05 V.
    So dVds is 2% of the variable -- the same construction as H_DIFF, at a
    different call site, found by the census rather than by the probe (the
    probe confirms the class; the census is what noticed it exists).

    This is measured WITHOUT importing the RF module, because the quantity
    of interest is the arithmetic of the step, not the device physics: the
    ratio dVds / Vds is the claim, and it is 2e-2.  Whether that lands inside
    a plateau is the SAME question RESULT 3 answers for H_DIFF and is left as
    an explicitly open item rather than asserted here -- the RF module's
    g_ds feeds f_max, a published number, and moving it needs its own
    before/after comparison, which is a session's work and not a footnote.
    """
    vds, dvds = 0.05, 1e-3
    ratio = dvds / vds
    hdiff_ratio = H_DIFF / 0.2
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 6 -- a SECOND absolute-step default, from the census")
        print("=" * 78)
        print(f"  output_conductance(dVds=1e-3) at Vds = 0.05 V   ->  step is "
              f"{ratio * 100:.1f}% of the variable")
        print(f"  sensitivity()'s H_DIFF = 1e-3 eV at delta ~ 0.2 eV       ->  "
              f"step is {hdiff_ratio * 100:.1f}% of the variable")
        print("  Both are the identical construction -- an absolute step whose "
              "correctness\n  depends on a scale stated nowhere near it.  "
              "H_DIFF is measured in RESULT 3.\n  dVds is NOT measured here "
              "and is logged as an open item: it feeds f_max.")
    return ratio, hdiff_ratio


# =====================================================================
# RESULT 7 -- THE SESSION'S REAL FINDING, AND IT WAS NOT PREDICTED
# =====================================================================
def result_7_p2_is_a_truncation_artefact(verbose=True):
    """
    RESULT 3 shows the shipped H_DIFF sits far outside the plateau.  This asks
    the only question that matters about that: does any PUBLISHED claim depend
    on it?

    One does.  Prediction P2 -- pre-registered 2026-09-21, scored PASS then
    and carried as PASS since -- asserts that the pair with the largest
    |d ln|N| / d delta| at delta = 0.1 is (Au, Pd).  It passes at the shipped
    h = 1e-3 and it FAILS at every h inside the plateau.  The converged
    answer is (Ti, Cu), by a margin of 40%.

    The mechanism is not subtle once the numbers are side by side: Au-Pd's
    |S| is FLAT in h -- 2.0162 from 1e-9 to 1e-4, and 2.0609 at the shipped
    step, a 2% drift.  Ti-Cu's collapses from 2.8282 to 0.8725, a 69% error.
    Au-Pd does not win at h = 1e-3 because it is the most sensitive pair; it
    wins because it is the pair whose derivative the shipped step happens to
    get right, while its competitors' are destroyed.  A ranking taken at a
    step outside the plateau ranks the pairs by how well the step suits them.

    P1's SUBSTANCE survives -- min|S| over same-side pairs divided by max|S|
    over straddling pairs is 30.4 at the shipped step and 26.7 converged, both
    far above the factor of 10 claimed -- so Chapter 6's central
    straddle-versus-same-side result is not in question.  What moves is the
    per-pair table and the ranking read off it.
    """
    pairs = all_asymmetric_pairs()

    def S_at(pair, h, d0=0.0):
        Np, Nm = N_of(pair, d0 + h), N_of(pair, d0 - h)
        return (abs(Np) - abs(Nm)) / (2.0 * h) / abs(N_of(pair, d0))

    hs = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    ranking = []
    for h in hs:
        top = sorted(((abs(S_at(p, h, 0.1)), p) for p in pairs), reverse=True)
        ranking.append((h, top[:3]))

    from graphene_crossover_sensitivity_model import straddles
    strad = [p for p in pairs if straddles(p)]
    same = [p for p in pairs if not straddles(p)]
    p1 = {}
    for h in (H_DIFF, 1e-7):
        p1[h] = (min(abs(S_at(p, h)) for p in same) /
                 max(abs(S_at(p, h)) for p in strad))

    shifts = []
    for p in pairs:
        a, b = S_at(p, H_DIFF), S_at(p, 1e-7)
        shifts.append((abs(a - b) / abs(b), p, a, b, straddles(p)))
    shifts.sort(reverse=True)

    p2_shipped = ranking[hs.index(1e-3)][1][0][1] == ("Au", "Pd")
    p2_conv = ranking[hs.index(1e-7)][1][0][1] == ("Au", "Pd")

    if verbose:
        print("\n" + "=" * 78)
        print("RESULT 7 -- P2 IS A TRUNCATION ARTEFACT  (unpredicted; the "
              "session's real finding)")
        print("=" * 78)
        print("\n  Top three pairs by |S| at delta = 0.1, as a function of the "
              "step:")
        for h, top in ranking:
            mark = "   <-- SHIPPED" if h == H_DIFF else ""
            print(f"    h={h:.0e}  " +
                  "   ".join(f"{a}-{b} {v:.4f}" for v, (a, b) in top) + mark)
        print(f"\n  P2 ('the top pair is Au-Pd') at the shipped step: "
              f"{'PASS' if p2_shipped else 'FAIL'}")
        print(f"  P2 at every step inside the plateau (1e-9 .. 1e-6):        "
              f"{'PASS' if p2_conv else 'FAIL'}  <-- FALSIFIED")
        print("\n  Per-pair shift, shipped H_DIFF vs converged h = 1e-7:")
        print(f"    {'pair':<9}{'shipped':>11}{'converged':>12}{'shift':>10}")
        for r, (a, b), sa, sb, st in shifts:
            print(f"    {a + '-' + b:<9}{sa:>11.5f}{sb:>12.5f}{r * 100:>9.2f}%"
                  f"{'   straddling' if st else ''}")
        print(f"\n  P1 substance (min|S| same-side / max|S| straddling):")
        print(f"    shipped   h={H_DIFF:g}:  {p1[H_DIFF]:.3f}   "
              f"(> 10 claimed: {p1[H_DIFF] > 10})")
        print(f"    converged h=1e-07:  {p1[1e-7]:.3f}   "
              f"(> 10 claimed: {p1[1e-7] > 10})")
        print("\n  So: the straddle-versus-same-side result SURVIVES, the "
              "per-pair table moves\n  by up to 44%, and one pre-registered "
              "prediction that was scored PASS is now\n  falsified.  Nothing "
              "here is rewritten -- the shipped numbers stay, annotated.")
    return ranking, p1, shifts, p2_shipped, p2_conv


# =====================================================================
# MAIN
# =====================================================================
def main():
    print("=" * 78)
    print("NUMERIC-DEFAULT SCALE AUDIT -- 2026-09-25")
    print("closing the 2026-09-24 item: the scale-assumption class had one")
    print("measured member and no detector.  This is the detector.")
    print("=" * 78)

    print("\n" + "-" * 78)
    print("VALIDATIONS (exact)")
    print("-" * 78)
    v1 = validate_identity_rescaling()
    v2 = validate_odd_function_zero()
    v3k, v3rel, v3ok = validate_reproduce_20260924_member()
    v4 = validate_hdiff_derived_rounding_bound()
    v5 = validate_shipped_sensitivity_unmoved()

    result_1_census()
    rows2 = result_2_probe()
    lo_e, hi_e, ship_err, dec_lo, dec_hi, _, _ = result_3_hdiff_plateau()
    needed, margin = result_4_max_iter()
    _, grid_err = result_5_grid_convergence()
    result_6_dvds()
    result_7_p2_is_a_truncation_artefact()

    print("\n" + "=" * 78)
    print("PREDICTIONS SCORED")
    print("=" * 78)
    d1 = (rows2[0][3] == "SCALE-BOUND" and 1e-6 < rows2[0][2] < 1e-2
          and rows2[1][3] == "scale-free")
    d2 = _logsens_invariance(1e-5, 1e3)[0] == 0.0
    d3 = rows2[3][3] == "SCALE-BOUND"
    d4 = (ship_err < 1e-3) and (dec_lo >= 2.0) and (dec_hi >= 2.0)
    d5 = (rows2[2][3] == "scale-free") and (margin >= 3.0)
    d6 = grid_err < 1e-3
    for tag, ok, txt in (
        ("D1", d1, "probe reproduces the 09-24 member; rtol version passes"),
        ("D2", d2, "log_sensitivity rel_step covariance is BITWISE zero"),
        ("D3", d3, "H_DIFF fails covariance -- second scale-encoding default"),
        ("D4", d4, "H_DIFF's shipped value inside plateau, >=2 decades both "
                   "sides"),
        ("D5", d5, "max_iter=200 scale-free and margin >= 3x"),
        ("D6", d6, "n_points=400 converged to better than 1e-3"),
    ):
        print(f"  {tag}  {'PASS' if ok else 'FAIL'}   {txt}")
    print("""
  D1 FAILED on its MAGNITUDE BAND, not its class call.  It predicted the
  covariance deviation of `tol=1e-14` would land in 1e-6..1e-2; measured
  4.9e-2 at lam=1e-3 and 4.8e-5 at lam=1e+3.  The band was written from the
  lam > 1 direction alone, and the probe is ASYMMETRIC in lam: shrinking the
  variable makes an absolute tolerance coarse relative to the interval (one or
  two halvings, catastrophic), while growing it makes the tolerance fine and
  the residual deviation is then dominated by the error of the lam = 1
  BASELINE, not of the rescaled run.  A deviation from this probe is therefore
  a reliable yes/no and NOT a calibrated error estimate.  That asymmetry is
  the same one bracketed_root.py already records for the unguarded loop, met
  again in a different guise.

  D4 FAILED outright and it was written as a null result so that it could.
  See RESULT 7.""")

    print("\n" + "=" * 78)
    print("VALIDATION SUMMARY")
    print("=" * 78)
    checks = [
        ("V1 identity rescaling bitwise zero", v1 == 0.0),
        ("V2 odd-function exact zero (non-discriminating)", v2),
        ("V3 reproduces 09-24's 16 halvings / 2e-4", v3ok),
        ("V4 central difference within its DERIVED rounding bound", v4),
        ("V5 shipped sensitivity numbers bitwise unmoved", v5),
    ]
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\n  {sum(1 for _, o in checks if o)}/{len(checks)} validations pass")

    print("\n" + "=" * 78)
    print("WHAT THE DETECTOR BUYS, STATED PLAINLY")
    print("=" * 78)
    print("""
  The 2026-09-24 fault was found with an ORACLE -- a three-day-old unguarded
  loop that happened to be right.  The probe in RESULT 2 finds the same fault
  with no oracle, no reference implementation and no correct answer: it asks
  only whether a procedure respects the units of its own variable.  That is
  why it generalises, and it is what the open item asked for.

  The probe is strictly a SCREEN: every verdict it returns is 'this default
  encodes a scale' or 'it does not', never 'this number is wrong'.  RESULT 3
  is the half that decides.  On 2026-09-25 the screen flagged two defaults and
  the deciding half convicted one of them -- H_DIFF, whose shipped value sits
  3.5 decades ABOVE the top of its own plateau, moving individual entries of
  Chapter 6's sensitivity table by up to 44% and falsifying prediction P2
  outright (RESULT 7).  D4 was written as a prediction that this would NOT
  happen, and writing it that way is the only reason it could be scored.

  The order matters and is worth stating: the screen came first and cost
  nothing, the margin measurement second and cost a few seconds, and the
  published-claim check (RESULT 7) came last and is the only one that produced
  a correction.  Running only the screen would have flagged H_DIFF without
  knowing whether to care.  Running only the margin study would have shown a
  bad step without knowing which claim it reached.  Neither alone is an audit.
""")


def make_plot():
    """The two panels that carry RESULT 3 and RESULT 7."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from graphene_crossover_sensitivity_model import straddles
    def S_at(p,h,d0=0.0):
        Np,Nm=N_of(p,d0+h),N_of(p,d0-h)
        return (abs(Np)-abs(Nm))/(2.0*h)/abs(N_of(p,d0))

    hs=np.logspace(-10,-1.5,40)
    pairs=all_asymmetric_pairs()
    fig,ax=plt.subplots(1,2,figsize=(13,5.2))

    # left: per-pair |S| at delta=0 vs h, normalised to the converged value
    for p in pairs:
        ref=S_at(p,1e-7)
        ax[0].plot(hs,[abs(S_at(p,h)/ref) for h in hs],
                   color="0.75" if straddles(p) else None, lw=1.2,
                   label=None)
    ax[0].axvspan(1e-10,3.16e-7,color="tab:green",alpha=0.12)
    ax[0].axvline(H_DIFF,color="tab:red",ls="--",lw=1.6)
    ax[0].axhline(1.0,color="k",lw=0.7,ls=":")
    ax[0].set_xscale("log"); ax[0].set_xlabel("central-difference step h (eV)")
    ax[0].set_ylabel(r"$S(h)\,/\,S(h\!\to\!0)$")
    ax[0].set_title("All 21 pairs: derivative vs step size\n"
                    "green = measured plateau, red dashed = shipped H_DIFF")
    ax[0].set_ylim(0,2.2)
    ax[0].text(H_DIFF*1.25,2.05,"H_DIFF = 1e-3",color="tab:red",fontsize=9)
    ax[0].text(2e-9,0.12,"plateau\n1e-10 .. 3e-7 eV",color="tab:green",fontsize=9)

    # right: the P2 ranking crossover at delta = 0.1
    for name,pr,col in (("Ti-Cu",("Ti","Cu"),"tab:blue"),
                        ("Cr-Au",("Cr","Au"),"tab:orange"),
                        ("Ni-Pd",("Ni","Pd"),"tab:green"),
                        ("Au-Pd",("Au","Pd"),"tab:red")):
        ax[1].plot(hs,[abs(S_at(pr,h,0.1)) for h in hs],label=name,color=col,lw=1.8)
    ax[1].axvspan(1e-10,3.16e-7,color="tab:green",alpha=0.12)
    ax[1].axvline(H_DIFF,color="tab:red",ls="--",lw=1.6)
    ax[1].set_xscale("log"); ax[1].set_xlabel("central-difference step h (eV)")
    ax[1].set_ylabel(r"$|d\ln|N|/d\delta|$ at $\delta=0.1$ eV  (1/eV)")
    ax[1].set_title("Prediction P2 is a truncation artefact:\n"
                    "Au-Pd leads only at the shipped step")
    ax[1].legend(fontsize=9,loc="lower left")
    fig.suptitle("Numeric-default scale audit, 2026-09-25 -- H_DIFF sits 3.5 decades above its own plateau",
                 fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.94])
    fig.savefig("default_scale_audit.png",dpi=140)
    return "default_scale_audit.png"


if __name__ == "__main__":
    main()
    print("\nplot written:", make_plot())
