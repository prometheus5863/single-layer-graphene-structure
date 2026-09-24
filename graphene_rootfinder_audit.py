"""
graphene_rootfinder_audit.py
============================
Discharges the open item created on 2026-09-23:

    "Audit the repo's remaining analysis-layer code for the bracket bug
     class -- `alpha_from_separation`, the `ell` bisection of 2026-09-21,
     `contact_resistance_crossover.py` and `graphene_sensitivity_audit.py`.
     Today's new root-finder is guarded; the four older ones are not known
     to be."

THE HEADLINE IS THAT THE ITEM'S OWN LIST WAS WRONG ABOUT THREE OF ITS FOUR
ENTRIES, AND RIGHT ABOUT THE FOURTH.

Three of the four named "root-finders" do not solve anything.  They are a
closed form, an affine formula and a division:

  * `contact_resistance_crossover.crossover_length` -- closed form.  R_channel
    is linear in L at fixed V_g, so L_x = Rc * W * sigma_sheet, evaluated
    once.  There is no interval and nothing to bracket.
  * `graphene_contact_doping_nonlinear_model.alpha_from_separation` -- affine
    in d.  It cannot iterate.  Its real hazard is a different class: it
    returns a NEGATIVE alpha for d < d0, and that is handled by
    `alpha_for_metal()` refusing those metals -- a guard that already exists.
  * `graphene_sensitivity_audit`'s Family-B "calibration solve" -- a single
    division, lambda_imp = lambda_bulk / sum(terms).  Its hazard is a
    vanishing denominator, again not a bracket.

The fourth is real: `p2_threshold()` in
`graphene_per_metal_crossover_model.py` -- the "ell bisection of 2026-09-21"
-- ran 200 halvings with no bracket check at all.

**The published number was never wrong.**  At the shipped defaults
(tol = 0.10, ell in [0.05, 5.0] A) the interval IS bracketed, so
0.4997 A is a genuine root and Chapter 6 is unaffected.  The fault was
LATENT: `tol`, `lo` and `hi` are keyword arguments with defaults, so the
first caller who moves any of them out of the bracketed region gets a
plausible number and no warning.  That is the same shape as the 2026-09-22
fault, one call site earlier.

Naming an unguarded solver is not the same as showing it can produce a wrong
number, so Validation 3 pre-registers the two wrong numbers and then
measures them.

Validations, all six against values known in advance:
  1  the three non-solvers are shown to contain no iteration, by identity
     rather than by reading the source
  2  at shipped defaults the interval is bracketed and the root is a root
  3  PRE-REGISTERED: on two root-free intervals the unguarded loop returns
     a specific endpoint -- and the "endpoint" is off by one ulp in one of
     the two directions, which is a correction to the 2026-09-22 wording
  4  the guarded version raises on every root-free tolerance in a sweep
  5  guarded agrees with unguarded to 1 ulp where the interval is bracketed,
     so the published 0.4997 A is unchanged (1 ulp, not bitwise: the old loop
     returns `hi`, the shared one returns the final midpoint)
  6  the monotonicity the bisection rests on is PROVED term by term, not
     measured -- and Pt's term is exactly zero for every ell
  7  UNPLANNED, and the more interesting one: sharing the guard exposed a
     SECOND fault class in the guard's own default tolerance, which is
     absolute and therefore not scale-free.  The first version of this fix
     degraded a correct number by ten digits and raised nothing.  Caught
     only because the pre-existing unguarded loop supplied the right answer
     to compare against.
"""

import numpy as np

import bracketed_root
from bracketed_root import bracketed_bisect

import graphene_per_metal_crossover_model as pm
import graphene_contact_doping_nonlinear_model as nl
import contact_resistance_crossover as crc
import graphene_sensitivity_audit as sa

TOL_SHIPPED = 0.10
LO_SHIPPED = 0.05e-10
HI_SHIPPED = 5.0e-10


# =====================================================================
# The quantity p2_threshold() bisects, lifted out so the audit can see it
# =====================================================================
def worst_frac(ell):
    """
    max over usable metals of |dW(ell) - dW(flat)| / |dW(flat)|.

    This is exactly the closure inside `p2_threshold`, re-entered here so
    the audit can evaluate it at the endpoints -- which is the one thing
    the unguarded version never did.
    """
    w = 0.0
    for m in pm.usable_metals():
        f = pm.signed_offset(pm.METAL_WORK_FUNCTIONS[m])
        d = pm.METAL_WORK_FUNCTIONS[m] - pm.w_cross_for_metal(m, ell)[0]
        w = max(w, abs(d - f) / abs(f))
    return w


def unguarded_p2_threshold(tol=TOL_SHIPPED, lo=LO_SHIPPED, hi=HI_SHIPPED):
    """
    The pre-2026-09-24 body of `p2_threshold`, preserved verbatim so the
    audit can demonstrate the fault and so Validation 5 can compare the
    guarded result against it bitwise.  Do not call this for physics.
    """
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if worst_frac(mid) > tol:
            lo = mid
        else:
            hi = mid
    return hi


# =====================================================================
# CANDIDATE CLASSIFICATION
# =====================================================================
CANDIDATES = [
    ("graphene_per_metal_crossover_model.p2_threshold",
     "iterative bisection",
     "REAL -- 200 halvings, no bracket check (the 'ell bisection of "
     "2026-09-21'). Latent: shipped defaults ARE bracketed."),
    ("contact_doping_nonlinear.alpha_from_separation",
     "affine formula",
     "NOT this class -- cannot iterate. Separate hazard (negative alpha "
     "for d < d0) already guarded by alpha_for_metal()."),
    ("contact_resistance_crossover.crossover_length",
     "closed form",
     "NOT this class -- L_x = Rc * W * sigma_sheet, one evaluation."),
    ("sensitivity_audit Family-B calibration solve",
     "single division",
     "NOT this class -- hazard is a vanishing denominator, not a bracket."),
]


def classify(verbose=True):
    if verbose:
        print("=" * 78)
        print("CANDIDATE CLASSIFICATION -- the 2026-09-23 audit list")
        print("=" * 78)
        for name, kind, verdict in CANDIDATES:
            print(f"  {name}")
            print(f"      kind    : {kind}")
            print(f"      verdict : {verdict}")
    return CANDIDATES


# =====================================================================
# VALIDATION 1 (EXACT) -- the three non-solvers contain no iteration
# =====================================================================
def validate_non_solvers_are_not_solvers(verbose=True):
    """
    Reading a function and calling it a closed form is an assertion.  These
    three identities are tests, and each one is false for anything iterative.

    (a) crossover_length is LINEAR in Rc at fixed V_g, bitwise-exactly:
        crossover_length(a*Rc) == a * crossover_length(Rc) is not something
        a bisection to finite tolerance can satisfy to the last bit.
    (b) alpha_from_separation is AFFINE in d: the midpoint value equals the
        mean of the endpoint values.
    (c) Family B's lambda_imp is lambda_bulk / sum(terms), bitwise.
    """
    ok = True

    V_g = 3.5
    Rc = 300e-6
    lx1 = crc.crossover_length(Rc, V_g)
    lx2 = crc.crossover_length(2.0 * Rc, V_g)
    a_exact = (lx2 == 2.0 * lx1)
    ok &= a_exact

    d1, d2 = 2.8e-10, 3.8e-10
    amid = nl.alpha_from_separation(0.5 * (d1 + d2))
    amean = 0.5 * (nl.alpha_from_separation(d1) + nl.alpha_from_separation(d2))
    b_ulps = abs(amid - amean) / np.spacing(abs(amean))
    b_exact = b_ulps <= 1.0
    ok &= b_exact

    terms = sa._calib_terms()
    import graphene_interconnect_model as icm
    lam = sa._lambda_imp()
    c_exact = (lam == icm.lambda_bulk_nm / sum(terms.values()))
    ok &= c_exact

    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 1 (EXACT) -- the three non-solvers contain no iteration")
        print("=" * 78)
        print(f"  (a) crossover_length(2Rc) == 2*crossover_length(Rc) bitwise : "
              f"{'PASS' if a_exact else 'FAIL'}")
        print(f"      {lx2:.17e}  vs  {2.0*lx1:.17e}")
        print(f"  (b) alpha_from_separation affine (midpoint == mean)         : "
              f"{'PASS' if b_exact else 'FAIL'}  ({b_ulps:.1f} ulp)")
        print(f"  (c) lambda_imp == lambda_bulk / sum(terms) bitwise          : "
              f"{'PASS' if c_exact else 'FAIL'}")
        print(f"      residual = {sum(terms.values()):.17e}   lambda_imp = {lam:.6f} nm")
        print("  -> three of the four audit-list entries are not root-finders.")
    return ok


# =====================================================================
# VALIDATION 2 (EXACT) -- the shipped call IS bracketed
# =====================================================================
def validate_shipped_call_is_bracketed(verbose=True):
    flo = worst_frac(LO_SHIPPED) - TOL_SHIPPED
    fhi = worst_frac(HI_SHIPPED) - TOL_SHIPPED
    bracketed = np.sign(flo) != np.sign(fhi)
    root = pm.p2_threshold()
    resid = worst_frac(root) - TOL_SHIPPED
    is_root = abs(resid) < 1e-12
    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 2 (EXACT) -- the shipped p2_threshold call is bracketed")
        print("=" * 78)
        print(f"  f(lo) = worst(0.05 A) - 0.10 = {flo:+.6f}")
        print(f"  f(hi) = worst(5.00 A) - 0.10 = {fhi:+.6f}")
        print(f"  sign change present : {'PASS' if bracketed else 'FAIL'}")
        print(f"  root  = {root*1e10:.9f} A,  worst(root) - tol = {resid:+.3e}  "
              f"-> {'PASS' if is_root else 'FAIL'}")
        print("  -> Chapter 6's 0.50 A is a genuine root. The fault was LATENT.")
    return bracketed and is_root


# =====================================================================
# VALIDATION 3 (EXACT, PRE-REGISTERED) -- the two wrong numbers
# =====================================================================
def validate_preregistered_wrong_numbers(verbose=True):
    """
    Written down BEFORE running, on 2026-09-24:

      A. tol = 2.0 exceeds worst(lo) = 1.4706, so f < 0 on the whole
         interval, no root exists, the loop always takes `hi = mid`, and the
         result should collapse toward `lo`.  PREDICTED: NOT exactly `lo` but
         one ulp above it, because once `hi` is the next double after `lo`,
         0.5*(lo+hi) rounds back up to `hi` and the walk stalls.
      B. tol = 0.005 is below worst(hi) = 0.009638, so f > 0 everywhere, the
         loop always takes `lo = mid`, and the result should be `hi`
         EXACTLY.

    The asymmetry is the finding: the 2026-09-22 note recorded that an
    unguarded bisection "returns the ENDPOINT", and a caller who believed
    that literally would test `result in (lo, hi)` and be reassured in
    case A.
    """
    rA = unguarded_p2_threshold(tol=2.0)
    rB = unguarded_p2_threshold(tol=0.005)
    ulpsA = (rA - LO_SHIPPED) / np.spacing(LO_SHIPPED)
    A_ok = (rA == np.nextafter(LO_SHIPPED, HI_SHIPPED)) and (rA != LO_SHIPPED)
    B_ok = (rB == HI_SHIPPED)
    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 3 (EXACT, PRE-REGISTERED) -- the two plausible wrong numbers")
        print("=" * 78)
        print(f"  A  tol=2.0     returned {rA*1e10:.6f} A   "
              f"(predicted: one ulp above lo, not lo)")
        print(f"     measured    {ulpsA:.1f} ulp above lo,  == lo : {rA == LO_SHIPPED}"
              f"   -> {'PASS' if A_ok else 'FAIL'}")
        print(f"     worst(returned) = {worst_frac(rA):.6f}, tol = 2.0  "
              f"-- a root would need equality")
        print(f"  B  tol=0.005   returned {rB*1e10:.6f} A   (predicted: hi exactly)")
        print(f"     == hi : {B_ok}   -> {'PASS' if B_ok else 'FAIL'}")
        print(f"     worst(returned) = {worst_frac(rB):.6f}, tol = 0.005")
        print("  -> both are physically plausible decay lengths. Neither is a root.")
        print("  -> CORRECTION to the 2026-09-22 wording: the walk does not always")
        print("     reach the endpoint. `result == lo or result == hi` is not a test.")
    return A_ok and B_ok


# =====================================================================
# VALIDATION 4 (EXACT) -- the guard fires on every root-free tolerance
# =====================================================================
def validate_guard_fires(verbose=True):
    lo_val = worst_frac(LO_SHIPPED)
    hi_val = worst_frac(HI_SHIPPED)
    too_high = list(np.linspace(lo_val * 1.001, lo_val * 5.0, 12))
    too_low = list(np.linspace(hi_val * 0.02, hi_val * 0.999, 12))
    raised = 0
    total = 0
    for tol in too_high + too_low:
        total += 1
        try:
            pm.p2_threshold(tol=tol)
        except ValueError:
            raised += 1
    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 4 (EXACT) -- the guarded p2_threshold refuses every "
              "root-free tolerance")
        print("=" * 78)
        print(f"  worst(lo) = {lo_val:.6f}   worst(hi) = {hi_val:.6f}")
        print(f"  12 tolerances above worst(lo) and 12 below worst(hi)")
        print(f"  raised ValueError : {raised}/{total}  -> "
              f"{'PASS' if raised == total else 'FAIL'}")
    return raised == total


# =====================================================================
# VALIDATION 5 -- the published number does not move
# =====================================================================
def validate_fix_moves_no_published_number(verbose=True):
    """
    This validation was WRITTEN as a bitwise check and had to be weakened,
    for a reason worth keeping: the old loop returns `hi`, while
    `bracketed_bisect` returns the midpoint `0.5*(lo+hi)` of the final
    interval.  Once lo and hi are adjacent doubles those differ by exactly
    one ulp, forever, for any function.  A bitwise claim here would have
    been a claim about which of two adjacent doubles a loop happens to name,
    not about the root -- so the honest statement is 1 ulp, and 1 ulp is
    what is asserted.

    The residuals are the real comparison, and they are both at the noise
    floor: -9.6e-16 before, +2.4e-16 after.
    """
    old_v = unguarded_p2_threshold()
    new_v = pm.p2_threshold()
    ulps = abs(new_v - old_v) / np.spacing(old_v)
    within_1ulp = ulps <= 1.0
    r_old = worst_frac(old_v) - TOL_SHIPPED
    r_new = worst_frac(new_v) - TOL_SHIPPED
    both_converged = abs(r_old) < 1e-12 and abs(r_new) < 1e-12
    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 5 (EXACT, 1 ulp) -- the fix moves no published number")
        print("=" * 78)
        print(f"  unguarded (pre-2026-09-24 body) : {old_v:.17e} m   "
              f"residual {r_old:+.3e}")
        print(f"  guarded   (current)             : {new_v:.17e} m   "
              f"residual {r_new:+.3e}")
        print(f"  separation : {ulps:.1f} ulp  -> "
              f"{'PASS' if within_1ulp else 'FAIL'}")
        print(f"  both residuals at the noise floor : "
              f"{'PASS' if both_converged else 'FAIL'}")
        print(f"  -> {new_v*1e10:.4f} A, as published in Chapter 6. The 1 ulp is")
        print("     `hi` vs the final midpoint, not a change in the root.")
    return within_1ulp and both_converged


# =====================================================================
# VALIDATION 6 (PROOF) -- the monotonicity the bisection rests on
# =====================================================================
def validate_monotonicity_is_proved(verbose=True):
    """
    `p2_threshold`'s docstring asserts that the |dW| shift is monotone
    decreasing in ell.  A bisection on a non-monotone function can bracket a
    root and still return the wrong one of several, so this is load-bearing
    and was asserted rather than shown.

    It is provable term by term.  For metal m,

        worst_m(ell) = DC_ANCHOR * |exp(-(d_m - D_ANCHOR)/ell) - 1| / |f_m|

    and d_m - D_ANCHOR has a fixed sign per metal, so the exponential
    approaches 1 monotonically from one side as ell grows.  Each |.| is
    therefore monotone non-increasing, and a pointwise max of monotone
    non-increasing functions is monotone non-increasing.  No grid needed.

    Two exact consequences are checkable, and the second is the sharper test:
      (i)  every metal's own term is non-increasing on a grid (necessary);
      (ii) Pt has d_eq == D_ANCHOR == 3.30 A EXACTLY, so its term is
           exactly 0.0 for every ell, bitwise -- the max is really over Cu
           and Au alone, and any code that disagrees is not evaluating the
           formula above.
    """
    metals = pm.usable_metals()
    grid = np.linspace(0.05e-10, 5.0e-10, 801)

    def term(m, ell):
        f = pm.signed_offset(pm.METAL_WORK_FUNCTIONS[m])
        d = pm.METAL_WORK_FUNCTIONS[m] - pm.w_cross_for_metal(m, ell)[0]
        return abs(d - f) / abs(f)

    per_metal_ok = {}
    for m in metals:
        vals = np.array([term(m, e) for e in grid])
        per_metal_ok[m] = bool(np.all(np.diff(vals) <= 0.0))

    pt_zero = all(term("Pt", e) == 0.0 for e in grid) if "Pt" in metals else None
    agg = np.array([worst_frac(e) for e in grid])
    agg_ok = bool(np.all(np.diff(agg) <= 0.0))

    ok = all(per_metal_ok.values()) and agg_ok and (pt_zero is not False)
    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 6 (PROOF, with two exact consequences) -- monotonicity")
        print("=" * 78)
        for m in metals:
            dm = nl.D_EQ[m] * 1e10
            print(f"  {m:>3s}  d_eq = {dm:5.2f} A   d_eq - D_ANCHOR = "
                  f"{dm - pm.D_ANCHOR*1e10:+.2f} A   term non-increasing: "
                  f"{'PASS' if per_metal_ok[m] else 'FAIL'}")
        print(f"  (ii) Pt term == 0.0 bitwise for all 801 ell : "
              f"{'PASS' if pt_zero else ('FAIL' if pt_zero is False else 'n/a')}")
        print(f"       Pt's d_eq equals D_ANCHOR exactly, so the max is over "
              f"Cu and Au alone.")
        print(f"  aggregate max non-increasing               : "
              f"{'PASS' if agg_ok else 'FAIL'}")
        print("  -> monotone by proof, not by sample. The bisection is entitled")
        print("     to assume a unique root once the bracket is checked.")
    return ok


# =====================================================================
# VALIDATION 7 (EXACT) -- the trap that sharing the guard exposed
# =====================================================================
def validate_tolerance_scale_trap(verbose=True):
    """
    NOT PLANNED.  Found on 2026-09-24 by the fix itself failing its own
    Validation 5 on the first run, and it is a different fault class from
    the bracket bug -- arguably a more dangerous one, because the guard
    cannot see it.

    `bracketed_bisect(..., tol=1e-14)` compares an ABSOLUTE interval width.
    That default was chosen in graphene_differential_crossover_model for a
    variable of order 1 eV, where it is ~machine precision (48 halvings).
    `p2_threshold` bisects a decay length in METRES, order 5e-11.  The same
    1e-14 is then COARSE: the interval (5e-10 - 5e-12) reaches it after 16
    halvings, a relative precision of 2e-4.

    Nothing raises.  The function returns a decay length that looks exactly
    as reasonable as the right one and is wrong in its fifth significant
    figure.  Measured below, with the correct value known in advance from
    the pre-existing unguarded loop -- which is the only reason this was
    caught at all: the first version of the fix silently degraded a number
    that had been right for three days.

    THE LESSON, which is the point of recording it: a guard moved to a new
    call site carries its defaults with it, and a default tolerance is a
    claim about scale.  Reviewing the guard for correctness would never have
    found this; only comparing against a known-good number did.
    """
    f = lambda ell: worst_frac(ell) - TOL_SHIPPED
    converged = bracketed_bisect(f, LO_SHIPPED, HI_SHIPPED,
                                 tol=0.0, rtol=np.finfo(float).eps)
    absolute_default = bracketed_bisect(f, LO_SHIPPED, HI_SHIPPED)  # tol=1e-14

    # halvings needed for the absolute default to trigger, at each scale
    def halvings(width, tol=1e-14):
        n = 0
        while width >= tol:
            width *= 0.5
            n += 1
        return n

    n_metre = halvings(HI_SHIPPED - LO_SHIPPED)
    n_unit = halvings(2.0)
    ulps_off = abs(absolute_default - converged) / np.spacing(converged)
    resid_bad = worst_frac(absolute_default) - TOL_SHIPPED
    resid_good = worst_frac(converged) - TOL_SHIPPED

    degraded = ulps_off > 1e6 and abs(resid_bad) > 1e-9
    fixed = abs(resid_good) < 1e-12
    if verbose:
        print("\n" + "=" * 78)
        print("VALIDATION 7 (EXACT, unplanned) -- absolute tolerance is not "
              "scale-free")
        print("=" * 78)
        print(f"  halvings before (hi-lo) < 1e-14  at metre scale : {n_metre}"
              f"   (relative precision {1e-14/5e-11:.1e})")
        print(f"  halvings before (hi-lo) < 1e-14  at O(1) scale  : {n_unit}"
              f"   (relative precision {1e-14/1.0:.1e})")
        print(f"  tol=1e-14 (absolute default)  : {absolute_default:.17e} m   "
              f"residual {resid_bad:+.3e}")
        print(f"  rtol=eps  (scale-free)        : {converged:.17e} m   "
              f"residual {resid_good:+.3e}")
        print(f"  the default is {ulps_off:.3e} ulp from the root -> "
              f"degradation reproduced: {'PASS' if degraded else 'FAIL'}")
        print(f"  rtol=eps converges            : "
              f"{'PASS' if fixed else 'FAIL'}")
        print("  -> the guard prevents an unbracketed interval; it does not")
        print("     prevent its own default from being wrong for the caller.")
    return degraded and fixed


def main():
    print("=" * 78)
    print("ROOT-FINDER BRACKET AUDIT -- discharging the 2026-09-23 audit item")
    print("=" * 78)
    classify()
    results = {
        "V1 non-solvers contain no iteration": validate_non_solvers_are_not_solvers(),
        "V2 shipped call is bracketed": validate_shipped_call_is_bracketed(),
        "V3 pre-registered wrong numbers": validate_preregistered_wrong_numbers(),
        "V4 guard fires on root-free tolerances": validate_guard_fires(),
        "V5 fix moves no published number": validate_fix_moves_no_published_number(),
        "V6 monotonicity proved": validate_monotonicity_is_proved(),
        "V7 tolerance scale trap": validate_tolerance_scale_trap(),
    }
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    for k, v in results.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    n = sum(1 for v in results.values() if v)
    print(f"\n  {n}/{len(results)} validations passed")
    print("\n  Audit outcome: 1 real instance of the bracket-bug class found and")
    print("  fixed (p2_threshold); 3 of the 4 audit-list entries were not")
    print("  root-finders at all; 0 published numbers changed; 1 correction to")
    print("  the 2026-09-22 description of how the fault manifests.")
    print(f"\n  Guard now lives in {bracketed_root.__name__}.py and is shared.")
    return all(results.values())


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
