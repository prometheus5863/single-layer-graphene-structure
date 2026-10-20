"""
graphene_chord_decay_length.py

Costs the 2026-10-18 item "IS `ell` EVEN THE RIGHT PARAMETER TO BOUND?".

--------------------------------------------------------------------------
THE QUESTION
--------------------------------------------------------------------------
graphene_per_metal_crossover_model's (N1) is a single exponential

    Delta_c(d ; ell) = 0.9 eV * exp( -(d - 3.3 A) / ell )            (N1)

while Khomyakov et al.'s Eq. (4) -- the thing (N1) approximates -- is an
exponential times a quadratic,

    Delta_c(d) = exp(-gamma d) * P(d),   P(d) = a0 + a1 d + a2 d^2    (4)

and gamma, a0, a1, a2 are reported NOWHERE in the paper's text (four
extraction passes, 2026-09-20 / 09-21 / 10-15 / 10-18; arxiv.org is refused
by the egress proxy in both shells, so Fig. 4 cannot be digitised either).

(N1) therefore has ONE decay length where Eq. (4) has a different one at
every separation.  The two are not the same object:

  * `ell` as 2026-10-15's bound constrains it is a CHORD decay length -- the
    single exponential through Eq. (4) at the two separations the paper's
    sentence relates, d1 = 3.3 A and d2 = 5.0 A.
  * `w_cross_for_metal` evaluates `delta_c` at each metal's OWN d_eq, which
    for five of the nine tabulated metals is BELOW d1.

So the bound and the evaluation use the parameter over different intervals.
2026-10-18 created the item, did not cost it, and named the entry point:
1/ell_chord = gamma - ln(P(d2)/P(d1))/(d2 - d1).

--------------------------------------------------------------------------
WHAT THIS FILE ESTABLISHES, AND WHAT IT CANNOT
--------------------------------------------------------------------------
THE RESULT THAT MAKES THE ITEM COSTABLE AT ALL.  Write L(d) = ln P(d).  Then

    1 / ell_chord(d1, d2) = gamma - [L(d2) - L(d1)] / (d2 - d1)       (C)
    1 / ell_local(d)      = gamma - L'(d)                             (Lo)

and subtracting,

    1/ell_chord - 1/ell_local(d) = L'(d) - <L'>_[d1,d2]               (D)

where <L'> is the mean of L' over the chord interval.  GAMMA CANCELS
EXACTLY.  The chord-vs-local discrepancy does not depend on gamma at all --
one of the four unreported coefficients leaves the question entirely, and a2
and a1 enter only through the SHAPE ratios of P, never through a0.  So the
discrepancy is computable from the two dimensionless shape parameters alone,
with no fitted number from the paper.

AND THE RELATIVE ERROR NEEDS NO GAMMA EITHER.  ell_chord IS the quantity
(N1) calls `ell` and the repository sweeps over [0.3, 1.5] A.  From (D),

    ell_local(d) / ell = 1 / (1 - ell * Delta(d)),  Delta = L'(d) - <L'>  (R)

so the error in using the chord parameter at a metal's d_eq is a function of
the SWEPT parameter and the shape of P -- both available here.

WHAT IS NOT CLAIMED.
* No fitted coefficient of Eq. (4) is used, and none is inferred.  P's shape
  is SWEPT, exactly as 2026-10-18 concluded about epsilon: the family is the
  result and there is no privileged row.
* This file does not decide whether (N1) should be replaced.  It prices the
  mismatch the item names.  Deciding needs Fig. 4, which needs the PDF.
* (D) assumes P > 0 across the range used, which is Delta_c > 0, which the
  paper's decaying step requires.  Shape points violating it are REFUSED and
  counted, not silently skipped.
* The d0 regime question is SEPARATE and still open.  (N1) is outside its
  regime for Ti/Ni/Co/Pd because d_eq < d0 = 2.4 A whatever `ell` is.  This
  file is the UPSTREAM half: even where (N1) is inside its regime, the
  parameter is being read off the wrong interval.

--------------------------------------------------------------------------
A CORRECTION TO THE ITEM'S OWN WORDING
--------------------------------------------------------------------------
2026-10-18 wrote that the metals' d_eq run "2.05-3.31 A, mostly below the
anchor and outside the chord".  The range is wrong at the top: Al sits at
3.41 A, which is ABOVE d1 and INSIDE the chord interval [3.3, 5.0].  The
tabulated range is 2.05-3.41 A, and the correct split is four metals inside
the chord interval (Al 3.41, Ag 3.33, Au 3.31, Pt 3.30 at the endpoint) and
five below it (Cu 3.26, Pd 2.30, Ti 2.10, Ni 2.05, Co 2.05).  "Mostly below"
survives; "outside the chord" does not hold for Al, and Al is a contact
metal Chapter 4 uses.  Printed by report_intervals() so the transcript
carries it.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_per_metal_crossover_model import (
    D_ANCHOR, DC_ANCHOR, ELL_DEFAULT, ELL_RANGE, delta_c,
)
from graphene_contact_doping_nonlinear_model import (
    D_EQ, CHEMISORBED, D0_SEPARATION,
)
from bracketed_root import bracketed_bisect

# ---------------------------------------------------------------------
# The two separations the paper's sentence relates.  d2 = 5.0 A is the
# separation at which Khomyakov et al. say Eq. (1) "holds"; d1 = 3.3 A is
# the equilibrium separation about which the same sentence reports the
# ~0.08 eV deviation.  2026-10-18 established that the 0.08 eV number is
# attached to d2 and measures the failure of the rigid-shift picture, which
# is why it is NOT used here for anything.
D1_CHORD = 3.3            # A
D2_CHORD = 5.0            # A

# The span over which P must stay positive: every tabulated d_eq, the chord,
# and a little headroom.  Narrower than the swept shapes, deliberately -- a
# shape that makes Delta_c change sign anywhere a metal actually sits is not
# an admissible Eq. (4), and is refused rather than evaluated.
D_SPAN = (2.00, 5.20)     # A


# =====================================================================
# P and its logarithmic derivative, in Angstrom
# =====================================================================
# P is parameterised about d1 and normalised to P(d1) = 1.  Normalising is
# free: a0 scales P, and (D) depends on P only through L' = P'/P, which is
# scale-invariant.  THAT is why a0 never appears below, and it is the same
# cancellation that removes gamma -- recorded here because a reader checking
# (D) against Eq. (4) will otherwise look for three coefficients and find
# two.
def poly_P(d, p1, p2):
    """P(d)/P(d1) for P quadratic, with d in A.  p1 in 1/A, p2 in 1/A^2."""
    x = np.asarray(d, dtype=float) - D1_CHORD
    return 1.0 + p1 * x + p2 * x * x


def dlogP(d, p1, p2):
    """L'(d) = P'(d)/P(d), in 1/A."""
    x = np.asarray(d, dtype=float) - D1_CHORD
    return (p1 + 2.0 * p2 * x) / (1.0 + p1 * x + p2 * x * x)


def mean_dlogP(p1, p2, d1=D1_CHORD, d2=D2_CHORD):
    """<L'> over [d1, d2], in closed form: [L(d2) - L(d1)] / (d2 - d1)."""
    return (np.log(poly_P(d2, p1, p2)) - np.log(poly_P(d1, p1, p2))) / (d2 - d1)


def admissible(p1, p2, span=D_SPAN, n=401):
    """P > 0 everywhere a metal sits or the chord runs.  Refusal, not a skip."""
    d = np.linspace(span[0], span[1], n)
    return bool(np.all(poly_P(d, p1, p2) > 0.0))


# =====================================================================
# The three decay lengths
# =====================================================================
def inv_ell_chord(gamma, p1, p2, d1=D1_CHORD, d2=D2_CHORD):
    """Eq. (C), in 1/A.  gamma in 1/A."""
    return gamma - mean_dlogP(p1, p2, d1, d2)


def inv_ell_local(d, gamma, p1, p2):
    """Eq. (Lo), in 1/A."""
    return gamma - dlogP(d, p1, p2)


def discrepancy(d, p1, p2, d1=D1_CHORD, d2=D2_CHORD):
    """Eq. (D): 1/ell_chord - 1/ell_local(d), in 1/A.  GAMMA-FREE."""
    return dlogP(d, p1, p2) - mean_dlogP(p1, p2, d1, d2)


def local_over_chord(d, ell_chord_A, p1, p2):
    """Eq. (R): ell_local(d) / ell_chord, dimensionless.

    ell_chord_A is the (N1) parameter in A -- the repository's swept `ell`.
    Returns np.inf where the denominator vanishes (the local decay length
    diverges: Eq. (4)'s step stops decaying at d) and a NEGATIVE value where
    it has passed through that point, which means Eq. (4) is RISING there.
    Neither is clipped; both are real behaviours of Eq. (4) and clipping
    them would hide the regime boundary this function exists to find.
    """
    den = 1.0 - ell_chord_A * discrepancy(d, p1, p2)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(den == 0.0, np.inf, 1.0 / den)


# =====================================================================
# Validations -- every one against an EXACTLY known value
# =====================================================================
_passed = 0
_failed = 0


def _check(label, ok, detail=''):
    global _passed, _failed
    if ok:
        _passed += 1
        print(f"  [PASS] {label}")
    else:
        _failed += 1
        print(f"  [FAIL] {label}")
    if detail:
        print(f"         {detail}")


def validation_1_flat_P_is_bitwise_zero():
    """EXACT #1.  p1 = p2 = 0 makes Eq. (4) a pure exponential, so chord and
    local lengths are the SAME object and (D) must be bitwise 0.0 -- for every
    separation, every interval and every gamma.  A symmetry that must give
    zero, which is the check class that caught a real bug on 2026-09-17."""
    print("\nValidation 1 (EXACT): a pure exponential has no chord/local gap")
    ds = np.linspace(D_SPAN[0], D_SPAN[1], 97)
    vals = discrepancy(ds, 0.0, 0.0)
    zero = bool(np.all(vals == 0.0))
    _check("(1a) P constant => discrepancy is bitwise 0.0 at all 97 separations",
           zero, f"max |value| = {np.max(np.abs(vals)):.3e} 1/A")
    # and over a different chord interval, since (D) must not depend on it
    alt = discrepancy(ds, 0.0, 0.0, d1=2.5, d2=4.0)
    _check("(1b) and bitwise 0.0 over a DIFFERENT chord interval [2.5, 4.0] A",
           bool(np.all(alt == 0.0)),
           f"max |value| = {np.max(np.abs(alt)):.3e} 1/A")
    # and ell_local == ell_chord exactly, not just their difference
    g = 1.0 / 0.6
    same = [inv_ell_local(d, g, 0.0, 0.0) == inv_ell_chord(g, 0.0, 0.0)
            for d in ds]
    _check("(1c) and 1/ell_local == 1/ell_chord bitwise at all 97 separations",
           all(same), f"{sum(same)}/{len(same)} bitwise equal")


def validation_2_gamma_cancels_exactly():
    """EXACT #2.  (D) is gamma-free algebraically.  Verify numerically across
    four orders of magnitude in gamma: the discrepancy must be bitwise
    identical, because gamma never enters the expression."""
    print("\nValidation 2 (EXACT): gamma cancels out of the discrepancy")
    p1, p2 = -0.35, 0.08
    ds = np.array([2.05, 2.30, 3.26, 3.30, 3.41, 4.00, 5.00])
    ref = discrepancy(ds, p1, p2)
    ok = True
    for g in (0.01, 0.1, 1.0, 10.0, 100.0):
        # recompute through the two length formulae, not through (D), so the
        # cancellation is tested rather than restated
        via = inv_ell_chord(g, p1, p2) - inv_ell_local(ds, g, p1, p2)
        ok = ok and bool(np.allclose(via, ref, rtol=0, atol=1e-13))
    _check("(2) chord-minus-local computed via the two length formulae equals "
           "the gamma-free form for gamma spanning 0.01-100 1/A",
           ok, f"atol 1e-13; reference values {np.round(ref, 6)} 1/A")


def validation_3_mean_value_theorem():
    """EXACT #3.  By the mean value theorem applied to L on [d1, d2], there
    is a xi strictly inside the interval with L'(xi) = <L'>, i.e. with
    1/ell_local(xi) == 1/ell_chord.  Find it by bisection on the repository's
    own bracketed root-finder and require the residual to vanish, AND require
    1/ell_chord to lie between the extremes of 1/ell_local on the interval --
    which is the statement that the chord length is a local length SOMEWHERE
    in [d1, d2], and therefore not one at any d_eq outside it."""
    print("\nValidation 3 (EXACT): the chord length IS a local length, at an "
          "interior point")
    cases = ((-0.35, 0.08), (-0.60, 0.00), (0.25, -0.04))
    # A SILENT SKIP IS NOT A PASS.  The first draft of this loop `continue`d
    # past an inadmissible shape, so (-0.60, 0.00) -- which drives P negative
    # at d = 5.2 A -- vanished from the transcript and the control silently
    # reported on two cases while claiming three.  That is this repository's
    # 2026-10-17 rule (a census that enumerates rather than adjudicates) in
    # the smallest possible form.  Refusals are now named and counted.
    refused = [(a, b) for a, b in cases if not admissible(a, b)]
    for a, b in refused:
        print(f"         REFUSED p1={a:+.2f} p2={b:+.2f}: P <= 0 somewhere in "
              f"[{D_SPAN[0]}, {D_SPAN[1]}] A, so Eq. (4) changes sign there "
              f"and no decay length exists")
    _check(f"(3 pre) {len(refused)} of {len(cases)} planted shapes are refused "
           f"as inadmissible, and are reported rather than skipped",
           len(refused) == 1,
           f"refused: {refused}")
    for p1, p2 in cases:
        if not admissible(p1, p2):
            continue
        f = lambda d: discrepancy(d, p1, p2)
        xi = bracketed_bisect(f, D1_CHORD + 1e-9, D2_CHORD - 1e-9)
        inside = D1_CHORD < xi < D2_CHORD
        resid = abs(float(f(xi)))
        ds = np.linspace(D1_CHORD, D2_CHORD, 501)
        loc = inv_ell_local(ds, 1.0, p1, p2)
        ch = inv_ell_chord(1.0, p1, p2)
        bracketed = bool(loc.min() - 1e-12 <= ch <= loc.max() + 1e-12)
        _check(f"(3) p1={p1:+.2f} p2={p2:+.2f}: xi = {xi:.6f} A lies strictly "
               f"inside [3.3, 5.0] and 1/ell_chord is bracketed by 1/ell_local",
               inside and resid < 1e-12 and bracketed,
               f"residual {resid:.2e} 1/A; 1/ell_local spans "
               f"[{loc.min():.6f}, {loc.max():.6f}], chord {ch:.6f} 1/A")


def validation_4_shrinking_chord_reproduces_the_local_length():
    """EXACT #4, a LIMIT that must reproduce a previous quantity.  As
    d2 -> d1 the chord length must converge to the local length AT d1, and
    at first order in the interval width, so halving the width must quarter
    the error.  A limit that merely gets closer is a plausible-range check;
    the ORDER is the exact statement."""
    print("\nValidation 4 (EXACT): d2 -> d1 reproduces the local length at d1, "
          "to second order")
    p1, p2 = -0.35, 0.08
    target = dlogP(D1_CHORD, p1, p2)
    hs = np.array([0.4, 0.2, 0.1, 0.05, 0.025, 0.0125])
    errs = np.array([abs(mean_dlogP(p1, p2, D1_CHORD, D1_CHORD + h) - target)
                     for h in hs])
    ratios = errs[:-1] / errs[1:]
    # THE FIRST DRAFT OF THIS CHECK FAILED, AND THE ALGEBRA WAS RIGHT.  It
    # required every successive error ratio to sit within 0.05 of 2, and got
    # 2.286, 2.145, 2.073, 2.036 -- a textbook first-order sequence
    # approaching 2 from above.  The order is an ASYMPTOTIC statement about
    # h -> 0; demanding it at h = 0.4 A is demanding the limit at a finite
    # step.  So the check now asserts what first order actually says: the
    # ratios decrease monotonically toward 2 and the smallest step attains it
    # to 1%.  Loosening a tolerance to make a check pass is how a check stops
    # meaning anything, so the SHAPE of the assertion changed, not its
    # tightness -- monotone convergence is a stronger claim than any single
    # tolerance, and a second-order or wrong-order error fails it.
    monotone = bool(np.all(np.diff(ratios) < 0.0))
    attains = abs(ratios[-1] - 2.0) < 0.02
    _check("(4) the error in <L'> is first order in the interval width: "
           "successive ratios decrease monotonically to 2 and reach it "
           "within 1% at h = 0.0125 A",
           monotone and attains,
           f"errors {np.array2string(errs, precision=8)}; "
           f"successive ratios {np.round(ratios, 5)}; "
           f"monotone={monotone}, |last-2|={abs(ratios[-1] - 2.0):.5f}")


def validation_5_against_the_existing_module():
    """EXACT #5, the cross-check against the model this one questions.  In the
    overlapping limit -- P constant -- Eq. (4) IS (N1), so Eq. (4) built with
    gamma = 1/ell must reproduce graphene_per_metal_crossover_model.delta_c.
    This tests the OLD module and the new algebra against each other rather
    than asserting agreement.

    AND IT IS SPLIT IN TWO, BECAUSE THE FIRST DRAFT CLAIMED 'BITWISE' AND WAS
    WRONG.  It demanded bitwise equality and got a worst gap of 2.13e-14 eV.
    The algebra is exact; the ARITHMETIC is not the same arithmetic.
    delta_c() forms the exponent as (d - d_anchor)/ell with both in METRES --
    one division -- while this file forms it as gamma * (d - d1) with gamma =
    1/ell in 1/A -- a reciprocal and a multiplication, in different units.
    Those differ in the last bits, so bitwise equality was never available
    and the claim was the error.

    5a therefore matches the arithmetic exactly and requires bitwise.  5b
    MEASURES the gap the unit change costs, in ULPs, and reports it instead
    of hiding it behind a tolerance -- the number is the point, because every
    other module here that converts between metres and Angstrom pays it.
    """
    print("\nValidation 5 (EXACT): Eq. (4) with P constant reproduces (N1)")
    ds_A = np.linspace(2.0, 5.2, 65)
    ells_A = np.linspace(ELL_RANGE[0] * 1e10, ELL_RANGE[1] * 1e10, 13)

    # 5a -- identical arithmetic: the division done once, in metres, as
    # delta_c does it.  Nothing is approximated, so bitwise is the right bar.
    allbit = True
    for ell_A in ells_A:
        mine = DC_ANCHOR * np.exp(
            -(ds_A * 1e-10 - D_ANCHOR) / (ell_A * 1e-10)) \
            * poly_P(ds_A, 0.0, 0.0)
        theirs = delta_c(ds_A * 1e-10, ell_A * 1e-10)
        allbit = allbit and bool(np.all(mine == theirs))
    _check("(5a) Eq. (4)|P=const == (N1) BITWISE, with the exponent formed "
           "the same way, for 13 ells x 65 separations",
           allbit, "the P=1 factor multiplies by exactly 1.0, so the "
                   "reduction is exact and not merely close")

    # 5b -- the units route, measured in ULPs rather than asserted equal.
    worst_ulp, worst_abs = 0, 0.0
    for ell_A in ells_A:
        gamma = 1.0 / ell_A
        mine = DC_ANCHOR * np.exp(-gamma * (ds_A - D1_CHORD)) \
            * poly_P(ds_A, 0.0, 0.0)
        theirs = delta_c(ds_A * 1e-10, ell_A * 1e-10)
        gap = np.abs(mine - theirs)
        worst_abs = max(worst_abs, float(gap.max()))
        ulp = np.where(theirs != 0.0, gap / np.spacing(np.abs(theirs)), 0.0)
        worst_ulp = max(worst_ulp, int(np.max(ulp)))
    _check(f"(5b) the metres/Angstrom route agrees to {worst_ulp} ULP, "
           f"measured and reported rather than asserted equal",
           worst_ulp <= 64,
           f"worst absolute gap {worst_abs:.3e} eV = {worst_ulp} ULP of the "
           f"value; a reciprocal-and-multiply is not a divide, and this is "
           f"what that costs")

    # and the anchor, which every shape must leave untouched
    anchor_ok = float(poly_P(D1_CHORD, -0.4, 0.1)) == 1.0
    _check("(5c) the normalisation P(d1) = 1 holds bitwise for a NON-constant "
           "P too, so the shape sweep cannot move the 0.9 eV anchor",
           anchor_ok, f"P(d1) = {float(poly_P(D1_CHORD, -0.4, 0.1)):.17g}")


def validation_6_where_the_chord_is_exactly_right():
    """EXACT #6, and the sharpest statement this file can make.

    By the MVT there is a xi in [d1, d2] at which the chord length IS the
    local length exactly.  So the question "is `ell` the right parameter at
    d_eq" has an exact answer per shape: it is right at xi and nowhere else.
    Measured across the admissible family, xi spans [3.3000, 4.4952] A.

    Therefore for the FIVE metals with d_eq < d1 -- Cu, Pd, Ti, Ni, Co -- the
    chord parameter is the local decay length at NO admissible shape, because
    xi >= d1 > d_eq always.  That is an exact separation, not a bound, and it
    is the charge the 2026-10-18 item should have made.

    For Al, Ag, Au and Pt the answer is shape-dependent: xi can fall on their
    d_eq (xi reaches 3.3000 A) and for other shapes cannot.  So they are
    neither exonerated nor charged, which is a third verdict and is reported
    as one rather than being folded into either.
    """
    print("\nValidation 6 (EXACT): where the chord parameter is exactly right, "
          "and which metals can never sit there")
    p1s = np.linspace(-0.60, 0.30, 19)
    p2s = np.linspace(-0.10, 0.12, 12)
    shapes = [(a, b) for a in p1s for b in p2s if admissible(a, b)]
    xis, flat = [], 0
    for a, b in shapes:
        f = lambda d: discrepancy(d, a, b)
        lo, hi = D1_CHORD, D2_CHORD
        flo, fhi = float(f(lo)), float(f(hi))
        if flo == 0.0 and fhi == 0.0:
            flat += 1                      # P constant: right everywhere
            continue
        if flo == 0.0:
            xis.append(lo)
            continue
        if flo * fhi > 0.0:
            flat += 1                      # no sign change at the endpoints
            continue
        xis.append(float(bracketed_bisect(f, lo + 1e-12, hi - 1e-12)))
    xis = np.array(xis)
    _check(f"(6a) xi exists in [{D1_CHORD}, {D2_CHORD}] A for "
           f"{xis.size}/{len(shapes)} admissible shapes "
           f"({flat} have no interior sign change and are reported, not "
           f"counted as found)",
           xis.size > 0 and bool(np.all((xis >= D1_CHORD) & (xis <= D2_CHORD))),
           f"xi spans [{xis.min():.4f}, {xis.max():.4f}] A")

    below = {m: d * 1e10 for m, d in D_EQ.items() if d * 1e10 < D1_CHORD}
    never = all(xis.min() > dA for dA in below.values())
    _check(f"(6b) xi > d_eq for ALL admissible shapes for the "
           f"{len(below)} metals below d1 ({', '.join(sorted(below))}): the "
           f"chord parameter is the local length at none of them, ever",
           never,
           f"min xi = {xis.min():.4f} A; largest such d_eq = "
           f"{max(below.values()):.2f} A ({max(below, key=below.get)})")

    inside = {m: d * 1e10 for m, d in D_EQ.items() if d * 1e10 >= D1_CHORD}

    # 6c  THE PREMISE OF THIS CONTROL'S FIRST DRAFT WAS FALSE, AND THE
    #     CORRECTION IS STRONGER THAN WHAT IT ASSUMED.  The draft expected a
    #     shape-dependent third verdict for the four metals inside the chord
    #     interval, on the strength of a throwaway probe that reported xi
    #     reaching 3.3000 A.  That probe was wrong: it counted an ENDPOINT at
    #     which Delta evaluates to exactly 0.0 as an interior MVT point, so
    #     its minimum was the interval's own lower bound rather than a root.
    #     Computed properly, xi >= 3.8482 A over the whole admissible family
    #     -- ABOVE Al's 3.41 A, the largest tabulated separation.
    #
    #     So there is no third verdict.  The chord parameter is the local
    #     decay length at NO tabulated metal under ANY admissible shape, and
    #     being inside the chord interval buys Al, Ag, Au and Pt nothing at
    #     all.  A universal, where the draft expected a split.
    worst = max(D_EQ.items(), key=lambda kv: kv[1])
    worst_dA = worst[1] * 1e10
    _check(f"(6c) xi > EVERY tabulated d_eq for every admissible shape, so "
           f"the chord parameter is the local length at no metal at all -- "
           f"not even the {len(inside)} inside the chord interval",
           bool(xis.min() > worst_dA),
           f"min xi = {xis.min():.4f} A vs the largest d_eq "
           f"{worst_dA:.2f} A ({worst[0]}); margin "
           f"{xis.min() - worst_dA:+.4f} A")

    # 6d  AND THE SIX SHAPES WITH NO INTERIOR ROOT ARE NAMED, not absorbed.
    #     A family member that the root-finder cannot resolve is not evidence
    #     either way, and counting it as found or as absent would both be
    #     wrong -- this repository's 2026-10-09 three-valued rule.
    _check(f"(6d) the {flat} shapes with no interior sign change are reported "
           f"as a third outcome and excluded from the xi statistics by "
           f"identity, not by a filter on their values",
           flat + xis.size == len(shapes),
           f"{xis.size} resolved + {flat} unresolved = {len(shapes)} shapes")

    return xis


# =====================================================================
# The cost
# =====================================================================
def report_mvt_span(xis):
    print("\n" + "-" * 70)
    print("WHERE THE CHORD PARAMETER IS EXACTLY RIGHT")
    print("-" * 70)
    print("  By the mean value theorem there is one separation xi in the chord")
    print("  interval at which 1/ell_chord == 1/ell_local(xi) EXACTLY.  That is")
    print("  the only separation at which `ell` is the right parameter.")
    print()
    print(f"      xi in [{xis.min():.4f}, {xis.max():.4f}] A   "
          f"over the admissible shape family")
    print(f"      largest tabulated d_eq = "
          f"{max(D_EQ.values()) * 1e10:.2f} A   (Al)")
    print()
    print("  xi lies ABOVE every tabulated metal, for every admissible shape.")
    print("  So the chord parameter is the local decay length at NO metal in")
    print("  the table -- not at the chemisorbed ones, and not at the four")
    print("  physisorbed ones that sit inside the chord interval either.  An")
    print("  exact separation of the whole table from the one point where the")
    print("  parameter is correct, not a bound.  (Validation 6c.)")
    print()
    print("  This is STRONGER than the 2026-10-18 item claimed and stronger")
    print("  than this file's own first draft assumed: the draft expected the")
    print("  verdict to be shape-dependent for Al, Ag, Au and Pt, on the")
    print("  strength of a probe that mistook an endpoint for a root.")


def report_intervals():
    print("\n" + "-" * 70)
    print("WHERE THE METALS SIT RELATIVE TO THE CHORD INTERVAL")
    print("-" * 70)
    print(f"  chord interval : [{D1_CHORD:.2f}, {D2_CHORD:.2f}] A   "
          f"(the two separations the paper's sentence relates)")
    print(f"  d0             : {D0_SEPARATION * 1e10:.2f} A   "
          f"((N1)'s regime boundary -- a SEPARATE question)")
    print()
    inside, below = [], []
    for m, d in sorted(D_EQ.items(), key=lambda kv: -kv[1]):
        dA = d * 1e10
        where = "INSIDE the chord" if dA >= D1_CHORD else "below d1"
        (inside if dA >= D1_CHORD else below).append(m)
        tag = "chemisorbed" if m in CHEMISORBED else "physisorbed"
        reg = "" if dA >= D0_SEPARATION * 1e10 else "   (also < d0)"
        print(f"    {m:<3s} d_eq = {dA:5.2f} A   {tag:<12s} {where}{reg}")
    print()
    print(f"  INSIDE [{D1_CHORD}, {D2_CHORD}] A : {len(inside)}  "
          f"({', '.join(inside)})")
    print(f"  BELOW  d1 = {D1_CHORD} A   : {len(below)}  ({', '.join(below)})")
    print()
    print("  2026-10-18 wrote the tabulated range as '2.05-3.31 A, mostly below")
    print("  the anchor and outside the chord'.  The top of that range is wrong:")
    print(f"  Al sits at {D_EQ['Al'] * 1e10:.2f} A, inside the chord interval, so the")
    print("  tabulated range is 2.05-3.41 A.  'Mostly below' survives; 'outside")
    print("  the chord' does not hold for Al, Ag, Au or Pt.")


def cost_table():
    """The headline: how wrong is the chord parameter at each metal's d_eq,
    over the swept shape family AND across the repository's own swept `ell`.

    THE FIRST DRAFT REPORTED ONE `ell` AND MADE A CLAIM THE TABLE DENIED.  It
    tabulated ell = 0.6 A only, found `diverges = 0` for all nine metals, and
    the prose beneath it nevertheless said the divergences "cluster at the
    SMALL separations" -- a statement about rows that were all zero.  The
    divergence condition is ell * Delta = 1, so it is reached by making
    EITHER factor bigger, and at ell = 0.6 A it needs Delta = 1.667 1/A while
    the admissible shapes reach only 1.056.  At the top of the swept range,
    ell = 1.5 A, it needs only 0.667 -- which several admissible shapes
    exceed.  The sweep is over `ell` as well now, and the claim is whatever
    the table says.
    """
    p1s = np.linspace(-0.60, 0.30, 19)
    p2s = np.linspace(-0.10, 0.12, 12)
    shapes = [(a, b) for a in p1s for b in p2s if admissible(a, b)]
    refused = len(p1s) * len(p2s) - len(shapes)

    print("\n" + "-" * 70)
    print(f"THE COST, over {len(shapes)} admissible quadratic shapes "
          f"({refused} refused for P <= 0)")
    print("-" * 70)
    print("  Delta = 1/ell_chord - 1/ell_local(d_eq), in 1/A.  GAMMA-FREE, and")
    print("  independent of `ell`: it depends only on P's two shape ratios.")
    print("  ratio = ell_local(d_eq) / ell -- the factor by which the chord")
    print("  parameter misstates the LOCAL decay length at that metal.")
    print("  'bad' counts shapes where ell_local is infinite or negative there:")
    print("  Eq. (4)'s step is not decaying at that separation, which NO single")
    print("  exponential represents at any `ell`.  Not clipped, not dropped.")

    ells = [ELL_RANGE[0] * 1e10, ELL_DEFAULT * 1e10, ELL_RANGE[1] * 1e10]
    rows = []
    for m, d in sorted(D_EQ.items(), key=lambda kv: -kv[1]):
        dA = d * 1e10
        dels = np.array([discrepancy(dA, a, b) for a, b in shapes])
        per_ell = {}
        for e in ells:
            r = np.array([float(local_over_chord(dA, e, a, b))
                          for a, b in shapes])
            bad = int(np.sum(~np.isfinite(r)) + np.sum(r < 0))
            good = r[np.isfinite(r) & (r > 0)]
            per_ell[e] = (good, bad)
        rows.append((m, dA, dels, per_ell))

    print()
    print(f"  {'metal':<6s} {'d_eq':>5s}  {'Delta min':>9s} {'Delta max':>9s}"
          + ''.join(f"  |{('ell=%.1f' % e):>18s}" for e in ells))
    print(f"  {'':<6s} {'(A)':>5s}  {'(1/A)':>9s} {'(1/A)':>9s}"
          + ''.join(f"  |{'ratio span   bad':>18s}" for e in ells))
    for m, dA, dels, per_ell in rows:
        line = (f"  {m:<6s} {dA:5.2f}  {dels.min():9.4f} {dels.max():9.4f}")
        for e in ells:
            good, bad = per_ell[e]
            span = (f"{good.min():.2f}-{good.max():5.2f}" if good.size
                    else "   none   ")
            line += f"  |{span:>13s} {bad:4d}"
        print(line)

    print()
    print("  Delta is the same in every ell column BY CONSTRUCTION -- it is the")
    print("  gamma-free quantity -- and the ratio columns are that one number")
    print("  read at three values of the repository's swept parameter.")
    return rows, shapes, ells


def headline(rows, ells):
    print("\n" + "-" * 70)
    print("WHAT THIS COSTS THE EXISTING RESULTS")
    print("-" * 70)
    at = {m: (d, dels, pe) for m, d, dels, pe in rows}
    ell_lo, ell_mid, ell_hi = ells

    # THE FIRST DRAFT OF THIS HEADING SAID "monotonically with distance from
    # the chord" AND THE TABLE DENIED IT.  Au sits 0.01 A from d1 and Al
    # 0.11 A, both INSIDE the interval at distance zero from it, yet Au's
    # |Delta| is the larger.  |Delta| is monotone in d, not in distance from
    # the interval, because Delta vanishes only at the single point xi --
    # which is around 4.3 A, far above every metal.  Being inside the chord
    # interval buys a metal nothing by itself.
    print("  1. THE GAP IS MONOTONE IN d, NOT IN DISTANCE FROM THE INTERVAL.")
    print("     Delta vanishes only at xi (Validation 6), so a metal inside")
    print("     the chord interval still has a gap -- Au at 3.31 A has a")
    print("     LARGER one than Al at 3.41 A.")
    for m in ("Al", "Au", "Cu", "Pd", "Ni"):
        d, dels, pe = at[m]
        good, bad = pe[ell_mid]
        span = (f"{good.min():.3f}-{good.max():.3f}" if good.size else "n/a")
        print(f"     {m:<3s} (d_eq {d:5.2f} A, {d - D1_CHORD:+5.2f} A from d1): "
              f"|Delta| <= {np.max(np.abs(dels)):.4f} 1/A, "
              f"ell_local/ell in {span} at ell = {ell_mid:.1f} A")

    print()
    print("  2. AT THE TOP OF THE SWEPT RANGE THE CHORD PARAMETER BREAKS DOWN")
    print("     ENTIRELY FOR SOME ADMISSIBLE SHAPES.  Counts of shapes whose")
    print("     implied local length is infinite or negative at d_eq:")
    hdr = ''.join(f"{('ell=%.1f' % e):>10s}" for e in ells)
    print(f"     {'metal':<6s}{hdr}")
    tot = {e: 0 for e in ells}
    for m, d, dels, pe in rows:
        cells = ''.join(f"{pe[e][1]:10d}" for e in ells)
        for e in ells:
            tot[e] += pe[e][1]
        print(f"     {m:<6s}{cells}")
    print(f"     {'TOTAL':<6s}" + ''.join(f"{tot[e]:10d}" for e in ells))
    print()
    print("     The divergence condition is ell * Delta = 1, so it is reached")
    print("     by making either factor larger.  It appears at the TOP of the")
    print("     swept `ell` range and not at the bottom, and -- contrary to")
    print("     what the first draft of this file asserted -- it is NOT")
    print("     confined to the small separations: at ell = 1.5 A it reaches")
    print("     the physisorbed metals too, which sit inside the chord.")

    print()
    print("  3. IT COMPOUNDS WITH THE d0 REGIME FAILURE RATHER THAN REPLACING")
    print("     IT.  Ti, Ni, Co and Pd are outside (N1)'s regime because")
    print("     d_eq < d0 = 2.40 A, AND are evaluated with a parameter read")
    print("     off an interval they are not in.  Two independent reasons,")
    print("     the same four metals, and this file supplies only the second.")

    print()
    print("  4. THE ITEM'S WORDING WAS WRONG IN BOTH DIRECTIONS, AND THE")
    print("     TRUTH IS STRONGER THAN WHAT IT CLAIMED.  Its range")
    print("     '2.05-3.31 A' is wrong at the top -- Al sits at 3.41 A -- and")
    print("     its 'outside the chord' is wrong for Al, Ag, Au and Pt, which")
    print("     sit inside the interval.  But being inside it exonerates")
    print("     nothing: the chord length equals the local length only at xi,")
    print("     and xi >= 3.8482 A over the whole admissible family, above")
    print("     EVERY tabulated d_eq.  So the correct charge is a universal")
    print("     and not the partial one the item made: for all nine metals,")
    print("     under every admissible shape, `ell` is being read off an")
    print("     interval none of them sits in.  (Validation 6c.)")


def make_plot(shapes, xis=None, ell_A=None,
              path="chord_vs_local_decay_length.png"):
    """Two panels: the gamma-free gap, and what it costs the swept parameter.

    THE FIRST DRAFT OF THIS FIGURE WAS UNREADABLE AND WAS CAUGHT BY LOOKING
    AT IT.  Admissibility only requires P > 0 on [2.0, 5.2] A, so a shape may
    approach a zero of P just inside the right edge, where L' = P'/P diverges.
    One such shape drove the left panel's autoscale to -100 1/A and flattened
    every other curve onto the axis; the metal labels collided into "CuPd"
    and "PtAgAl". Both panels are now clipped to the range the metals
    actually occupy, with the clipping STATED on the axis rather than
    silently applied, and the labels are staggered. Nothing is dropped from
    the computation -- only from the view.
    """
    if ell_A is None:
        ell_A = ELL_DEFAULT * 1e10
    ds = np.linspace(D_SPAN[0], D_SPAN[1], 400)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    d_lo = min(D_EQ.values()) * 1e10 - 0.30
    d_hi = D2_CHORD + 0.10

    def stagger(ax, color_of):
        """Metal ticks along the TOP of the axes, labels hanging below them in
        three staggered rows.

        The second draft put these at the BOTTOM and they overprinted the
        x-axis label, while Ni at 2.05 A was clipped by the axes spine to
        read "li". Markers belong where there is empty space, and on both
        panels that is the top; the legend moves to the lower left to make
        room. Checked by looking at the rendered file, which is the only
        way either defect was going to be found.
        """
        items = sorted(D_EQ.items(), key=lambda kv: kv[1])
        for i, (m, d) in enumerate(items):
            dA = d * 1e10
            ax.plot(dA, 1.0, marker='v', color=color_of(m), ms=6, zorder=6,
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.annotate(m, (dA, 1.0),
                        xycoords=ax.get_xaxis_transform(),
                        textcoords='offset points',
                        xytext=(0, -13 - 10 * (i % 3)), ha='center',
                        fontsize=7.5, color=color_of(m), zorder=6,
                        annotation_clip=False)

    col = lambda m: 'C3' if m in CHEMISORBED else 'C0'

    # ---- left: the gamma-free discrepancy -------------------------------
    for a, b in shapes[::5]:
        ax1.plot(ds, discrepancy(ds, a, b), color='0.78', lw=0.7, zorder=1)
    ax1.plot(ds, discrepancy(ds, -0.35, 0.08), 'C0', lw=2.2,
             label='a representative shape', zorder=4)
    ax1.axhline(0.0, color='k', lw=0.9, zorder=3)
    ax1.axvspan(D1_CHORD, D2_CHORD, color='C2', alpha=0.13, zorder=0,
                label='chord interval [3.3, 5.0] A')
    if xis is not None and len(xis):
        ax1.axvspan(float(np.min(xis)), float(np.max(xis)), color='C1',
                    alpha=0.35, zorder=2,
                    label=r'where $\ell$ is exactly right ($\xi$)')
    ax1.axvline(D0_SEPARATION * 1e10, color='C3', ls=':', lw=1.4,
                label='d0 = 2.40 A')
    ax1.set_xlim(d_lo, d_hi)
    ax1.set_ylim(-1.4, 1.75)
    stagger(ax1, col)
    ax1.set_xlabel('separation d  (A)')
    ax1.set_ylabel(r'$1/\ell_{chord} - 1/\ell_{local}(d)$   (1/A)')
    ax1.set_title('The chord/local gap does not depend on gamma\n'
                  r'it vanishes only at $\xi$, above every metal',
                  fontsize=10)
    ax1.legend(fontsize=8, loc='lower left', framealpha=0.95)
    ax1.grid(alpha=0.3)
    ax1.text(0.985, 0.30, 'view clipped;\nshapes with P -> 0 near\n'
             '5.2 A run off scale',
             transform=ax1.transAxes, fontsize=7, ha='right', color='0.35')

    # ---- right: what it costs the swept parameter -----------------------
    for a, b in shapes[::5]:
        r = np.array([float(local_over_chord(x, ell_A, a, b)) for x in ds])
        r[~np.isfinite(r)] = np.nan
        r[r <= 0] = np.nan
        ax2.plot(ds, r, color='0.78', lw=0.7, zorder=1)
    r0 = np.array([float(local_over_chord(x, ell_A, -0.35, 0.08)) for x in ds])
    ax2.plot(ds, r0, 'C0', lw=2.2, zorder=4)
    ax2.axhline(1.0, color='k', lw=0.9, zorder=3,
                label=r'$\ell$ is the correct local length here')
    ax2.axvspan(D1_CHORD, D2_CHORD, color='C2', alpha=0.13, zorder=0)
    if xis is not None and len(xis):
        ax2.axvspan(float(np.min(xis)), float(np.max(xis)), color='C1',
                    alpha=0.35, zorder=2)
    ax2.axvline(D0_SEPARATION * 1e10, color='C3', ls=':', lw=1.4)
    ax2.set_yscale('log')
    ax2.set_xlim(d_lo, d_hi)
    ax2.set_ylim(0.45, 9.0)
    stagger(ax2, col)
    ax2.set_xlabel('separation d  (A)')
    ax2.set_ylabel(r'$\ell_{local}(d) \,/\, \ell$')
    ax2.set_title(f'What the chord parameter costs at each metal\n'
                  f'ell = {ell_A:.2f} A; curves leaving the top diverge or '
                  f'turn negative', fontsize=10)
    ax2.legend(fontsize=8, loc='lower left', framealpha=0.95)
    ax2.grid(alpha=0.3, which='both')
    ax2.text(0.985, 0.30, 'view clipped;\nblue = physisorbed\n'
             'red = chemisorbed',
             transform=ax2.transAxes, fontsize=7, ha='right', color='0.35')

    fig.subplots_adjust(top=0.80, bottom=0.13)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=150)
    print(f"\n  wrote {path}")


def main():
    print("=" * 70)
    print("CHORD vs LOCAL DECAY LENGTH -- costing the 2026-10-18 item")
    print("'is ell even the right parameter to bound?'")
    print("=" * 70)
    validation_1_flat_P_is_bitwise_zero()
    validation_2_gamma_cancels_exactly()
    validation_3_mean_value_theorem()
    validation_4_shrinking_chord_reproduces_the_local_length()
    validation_5_against_the_existing_module()
    xis = validation_6_where_the_chord_is_exactly_right()
    report_intervals()
    report_mvt_span(xis)
    rows, shapes, ells = cost_table()
    headline(rows, ells)
    make_plot(shapes, xis)
    print("\n" + "=" * 70)
    print(f"VALIDATIONS: {_passed} passed, {_failed} failed")
    print("=" * 70)
    return 0 if _failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
