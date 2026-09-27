"""
graphene_log_sensitivity_step_audit.py

AN ANCHORED AUDIT OF `log_sensitivity`'s STEP AND OF ITS CONVERGENCE ESTIMATE
============================================================================

`graphene_sensitivity_audit.log_sensitivity(f, p, rel_step=1e-5)` is the only
finite-difference estimator in this thesis's conditioning work.  Every
parameter-level sensitivity printed by Chapter 5 Section 5.5 comes out of it,
and it reports with every value

    conv = |S(h) - S(2h)|

described in its own docstring as "a convergence estimate that is reported with
every value rather than assumed small."

That is step-doubling.  2026-09-25 showed a step-doubling convergence estimate
answering *yes* while the answer was wrong by 44%; 2026-09-26 showed that
refining one step converges to the derivative of whatever other discretisation
was held fixed.  This estimator has never been checked against an anchored
criterion.  It has been the top open numerical item since 2026-09-25.

Pre-registration (Q1-Q7, V1-V5), committed to git BEFORE this file existed:
    notes/2026-09-27-log-sensitivity-preregistration.md   (commit e2eed30)

--------------------------------------------------------------------------
THE TWO ESTIMATORS COMPARED HERE
--------------------------------------------------------------------------
Write L = ln p and g(L) = ln|f(e^L)|, so that S = g'(L).

SHIPPED ("affine"): evaluates at p(1 +/- h).  In L those points are
    L + ln(1+h)   and   L + ln(1-h),
and ln(1+h) != -ln(1-h).  It is therefore a secant over an ASYMMETRIC log
interval.  With a = ln(1+h), b = ln(1-h), s = a+b = ln(1-h^2), d = a-b:

    S_affine(h) = g' + (s/2) g'' + ((s^2+d^2)/2 - ... ) -> expanded:
    S_affine(h) = g' - (h^2/2) g'' + (h^2/6) g''' + O(h^4)          (A1)

The -(h^2/2) g'' term is a DISPLACEMENT OF THE EVALUATION POINT, not a
truncation error: the secant returns dg/dL at the interval's midpoint
L + s/2 = L + (1/2)ln(1-h^2), which is L - h^2/2 + O(h^4).  This is exactly
2026-09-26 item 5's fault class -- "an asymmetric secant is a second-order
estimate at that interval's midpoint, so the caller who asked for the
derivative at L receives the derivative somewhere else, correctly computed" --
found there on a guard path that had never fired, and found here on the
DEFAULT path behind every published number.

CORRECTED ("geometric"): evaluates at p*exp(+/-h).  Those points are L +/- h,
exactly symmetric, so it is a true centred difference:

    S_geom(h) = g' + (h^2/6) g''' + O(h^4)                          (A2)

The whole difference is one line of code.  The old body is retained (see
`graphene_sensitivity_audit.log_sensitivity(..., symmetric=False)`) because
2026-09-24 established that deleting a superseded implementation destroys the
only oracle available for judging its replacement.

--------------------------------------------------------------------------
WHAT MAKES THIS AUDIT ANCHORED RATHER THAN SELF-REFERENTIAL
--------------------------------------------------------------------------
Three of the five validations have a CLOSED FORM for the error at every h,
with no truncation at any order:

  ln|f| = A(ln p)^2  ->  E_affine(h) = A ln(1-h^2)  EXACTLY, all orders
                         E_geom(h)   = 0            EXACTLY

  ln|f| = A L^2 + B L^4 at L:
      E_affine(h) = A s + B[(2L+s)(2L^2+2Ls+(s^2+d^2)/2) - 4L^3]  EXACTLY
      E_geom(h)   = 4 B L h^2                                      EXACTLY

So one function separates the two estimators by a closed form at every step,
and a second lets the blind spots of `conv` be located to machine precision
rather than estimated.
"""
import math
import re
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_interconnect_model as icm
from graphene_sensitivity_audit import (
    log_sensitivity, _lambda_imp, _rho_propagated, W_AUDIT_NM)
from bracketed_root import bracketed_bisect

EPS = np.finfo(float).eps
SHIPPED_H = 1e-5


# =====================================================================
# The two estimators, written out so neither can be confused with the other
# =====================================================================
def S_affine(f, p, h):
    """The SHIPPED estimator: evaluation points p(1 +/- h)."""
    fp, fm = f(p * (1.0 + h)), f(p * (1.0 - h))
    return (math.log(abs(fp)) - math.log(abs(fm))) / (math.log1p(h) - math.log1p(-h))


def S_geom(f, p, h):
    """The CORRECTED estimator: evaluation points p*exp(+/-h), symmetric in ln p."""
    fp, fm = f(p * math.exp(h)), f(p * math.exp(-h))
    return (math.log(abs(fp)) - math.log(abs(fm))) / (2.0 * h)


def conv_of(estimator, f, p, h):
    """The step-doubling estimate, signed (conv as shipped is its absolute value)."""
    return estimator(f, p, h) - estimator(f, p, 2.0 * h)


# =====================================================================
# Test functions with EXACTLY known S, and EXACTLY known estimator errors
# =====================================================================
def power_law(A, n):
    return (lambda p: A * p ** n), (lambda p: float(n))


def quad_log(A):
    """ln|f| = A (ln p)^2."""
    return (lambda p: math.exp(A * math.log(p) ** 2)), (lambda p: 2.0 * A * math.log(p))


def quart_log(A, B):
    """ln|f| = A (ln p)^2 + B (ln p)^4."""
    def f(p):
        L = math.log(p)
        return math.exp(A * L * L + B * L ** 4)
    return f, (lambda p: 2.0 * A * math.log(p) + 4.0 * B * math.log(p) ** 3)


def _sd(h):
    """s = ln(1-h^2) and d = ln((1+h)/(1-h)), accurately."""
    return math.log1p(-h * h), math.log1p(h) - math.log1p(-h)


def E_affine_quad_exact(A, h):
    """Exact error of the shipped estimator on ln|f| = A(ln p)^2, all orders."""
    return A * math.log1p(-h * h)


def E_affine_quart_exact(A, B, L, h):
    """Exact error of the shipped estimator on ln|f| = A L^2 + B L^4, all orders."""
    s, d = _sd(h)
    quad = A * s
    quart = B * ((2 * L + s) * (2 * L * L + 2 * L * s + (s * s + d * d) / 2.0)
                 - 4 * L ** 3)
    return quad + quart


def E_geom_quart_exact(B, L, h):
    """Exact error of the corrected estimator on the same function."""
    return 4.0 * B * L * h * h


# =====================================================================
# Q7 -- is `conv` ever compared with anything?  A source census.
# =====================================================================
def q7_conv_census(verbose=True):
    import glob
    call_sites, comparisons = [], []
    SELF = "graphene_log_sensitivity_step_audit.py"
    for path in sorted(glob.glob("*.py")):
        if path == SELF:          # this module is the auditor, not the audited
            continue
        src = open(path).read().splitlines()
        for i, line in enumerate(src, 1):
            if "log_sensitivity(" in line and "def log_sensitivity" not in line:
                call_sites.append((path, i, line.strip()))
            if re.search(r"\bconv\w*\s*(<|>|<=|>=|==|!=)", line) or \
               re.search(r"(<|>|<=|>=|==|!=)\s*conv\w*\b", line):
                comparisons.append((path, i, line.strip()))
    if verbose:
        print("\n" + "=" * 78)
        print("Q7 -- SOURCE CENSUS: is `conv` ever compared against anything?")
        print("=" * 78)
        print(f"  log_sensitivity call sites found : {len(call_sites)}")
        for p, i, l in call_sites:
            print(f"    {p}:{i}  {l[:64]}")
        print(f"  lines in which a conv-like name appears in a COMPARISON : "
              f"{len(comparisons)}")
        for p, i, l in comparisons:
            print(f"    {p}:{i}  {l[:64]}")
        if not comparisons:
            print("    (none)")
        print("  NOTE: a `conv`-like quantity IS compared against a tolerance in")
        print("  graphene_crossover_sensitivity_model.py -- a DIFFERENT convergence")
        print("  quantity in a different module.  The repository knows how to attach")
        print("  a criterion to a convergence estimate; it did not do so here.")
    return call_sites, comparisons


# =====================================================================
# V1 / V4 -- the checks that CANNOT discriminate, measured as such
# =====================================================================
def v1_v4_power_law(verbose=True):
    f, Sx = power_law(3.71, -2.5)
    p0 = 4.3
    # h = 0.25 is the largest usable: conv needs 2h < 1 or p(1-2h) hits 0
    hs = [1e-8, 1e-6, SHIPPED_H, 1e-3, 1e-2, 0.1, 0.25]
    rows = []
    for h in hs:
        a = S_affine(f, p0, h)
        g = S_geom(f, p0, h)
        c = abs(conv_of(S_affine, f, p0, h))
        rows.append((h, a, g, c, a == Sx(p0), g == Sx(p0), c == 0.0))
    n_bit_a = sum(r[4] for r in rows)
    n_bit_g = sum(r[5] for r in rows)
    n_conv0 = sum(r[6] for r in rows)
    if verbose:
        print("\n" + "=" * 78)
        print("V1 + V4 (EXACT, AND DELIBERATELY NON-DISCRIMINATING)")
        print("  f = 3.71 p^-2.5 at p = 4.3 : S must be exactly -2.5 for EVERY h")
        print("=" * 78)
        print(f"  {'h':>8}{'S_affine':>20}{'S_geom':>20}{'conv':>12}  bitwise?")
        for h, a, g, c, ba, bg, b0 in rows:
            print(f"  {h:>8.0e}{a:>20.17g}{g:>20.17g}{c:>12.1e}"
                  f"  affine {str(ba):>5}  geom {str(bg):>5}  conv==0 {str(b0):>5}")
        print(f"  affine bitwise {n_bit_a}/{len(rows)}, geom bitwise "
              f"{n_bit_g}/{len(rows)}, conv exactly 0.0 {n_conv0}/{len(rows)}")
        print("  THIS CHECK CANNOT DISCRIMINATE.  For a power law g(L) is LINEAR, so")
        print("  every secant over it is exact and `conv` is exactly zero -- for the")
        print("  shipped estimator, for the corrected one, and for an absurd h = 0.25")
        print("  alike.  `graphene_sensitivity_audit.validate_power_laws` is built")
        print("  ENTIRELY out of this case (its four sub-checks (a),(b),(d) are pure")
        print("  power laws and (c) is closed-form), so that validation has never")
        print("  been capable of judging either the estimator or its step.")
    return n_bit_a, n_bit_g, n_conv0, len(rows)


# =====================================================================
# V2 -- the DISCRIMINATING exact check, and Q1's magnitude
# =====================================================================
def v2_quadratic_log(verbose=True):
    A = 0.8123
    f, Sx = quad_log(A)
    p0 = math.e ** 1.37
    hs = np.logspace(-3, -1, 9)
    rows = []
    for h in hs:
        ea = S_affine(f, p0, h) - Sx(p0)
        eg = S_geom(f, p0, h) - Sx(p0)
        pred = E_affine_quad_exact(A, h)
        rows.append((h, ea, pred, abs(ea - pred) / abs(pred), eg))
    worst_rel = max(r[3] for r in rows)
    worst_geom = max(abs(r[4]) for r in rows)
    # DERIVED cancellation floor, per h: two logs of magnitude |g| are
    # differenced and divided by an interval of width ~2h, so the noise in
    # E is ~ EPS*|g| / (2h).  The bound on the RELATIVE deviation from the
    # closed form is that noise divided by |E(h)| itself.  Stating it as a
    # constant instead is the fault class of 09-25 item 10 / 09-26 item 9,
    # and the first form of this validation committed exactly that error --
    # see the note.
    g0 = abs(math.log(abs(f(p0))))
    bounds = [EPS * g0 / (2 * h) / abs(pred) for h, _, pred, _, _ in rows]
    ratios = [r[3] / b for r, b in zip(rows, bounds)]
    worst_ratio = max(ratios)
    floor = EPS * g0 / (2 * hs[0])
    if verbose:
        print("\n" + "=" * 78)
        print("V2 (EXACT, DISCRIMINATING) -- ln|f| = A(ln p)^2, A = 0.8123, ln p = 1.37")
        print("  The shipped estimator's error is A*ln(1-h^2) EXACTLY, at every order.")
        print("  A true centred difference on the same function is EXACTLY ZERO,")
        print("  because a centred difference of a quadratic has no error at all.")
        print("=" * 78)
        print(f"  {'h':>10}{'E_affine measured':>22}{'A ln(1-h^2)':>22}"
              f"{'rel dev':>11}{'E_geom':>13}")
        for h, ea, pred, rel, eg in rows:
            print(f"  {h:>10.3e}{ea:>22.15e}{pred:>22.15e}{rel:>11.1e}{eg:>13.2e}")
        print(f"  worst relative deviation from the closed form : {worst_rel:.2e}")
        print(f"  measured / DERIVED cancellation bound, per h :")
        print("    " + "  ".join(f"{r:.2f}" for r in ratios))
        print(f"  worst measured/bound ratio : {worst_ratio:.3f}   "
              f"{'PASS' if worst_ratio <= 4.0 else 'FAIL'}")
        print("  FIRST FORM OF THIS VALIDATION FAILED and the failure is left in the")
        print("  record: it required the deviation below an absolute 1e-13 and")
        print("  measured 5.10e-08 at h = 1e-3, which IS the derived floor")
        print(f"  EPS*|g|/(2h)/|E| = {bounds[0]:.1e} there.  The constant silently")
        print("  encoded h >~ 1e-2.  That is 2026-09-25 item 10's and 2026-09-26")
        print("  item 9's fault class, committed for the third consecutive session")
        print("  by the session auditing it.")
        print(f"  worst |E_geom| over the same sweep            : {worst_geom:.2e}"
              f"   (cancellation floor at the smallest h is ~{floor:.1e})")
        print("  The corrected estimator is at its roundoff floor; the shipped one is")
        print("  above it by up to four orders of magnitude, on a function where the")
        print("  right answer is exactly recoverable.")
    return worst_rel, worst_geom, floor, worst_ratio


# =====================================================================
# V3 -- a symmetry that must give exactly zero
# =====================================================================
def v3_reciprocal_symmetry(verbose=True):
    A, B = -1.7, 0.43
    f, _ = quart_log(A, B)
    finv = lambda p: 1.0 / f(p)
    p0 = math.e ** 0.91
    hs = np.logspace(-8, -1, 31)
    n_bit = 0
    worst = 0.0
    worst_ratio = 0.0
    g0 = abs(math.log(abs(f(p0))))
    for h in hs:
        a = S_affine(f, p0, h)
        b = S_affine(finv, p0, h)
        if a + b == 0.0:
            n_bit += 1
        worst = max(worst, abs(a + b))
        bound = 4.0 * EPS * g0 / (2 * h)
        worst_ratio = max(worst_ratio, abs(a + b) / bound)
    if verbose:
        print("\n" + "=" * 78)
        print("V3 (EXACT SYMMETRY) -- S[1/f] = -S[f], which must be bitwise")
        print("=" * 78)
        print(f"  ln|1/f| = -ln|f| and the estimator is linear in g, so S[f] + S[1/f]")
        print(f"  must be exactly 0.0 at every h.")
        print(f"  bitwise zero : {n_bit}/{len(hs)} over h in [1e-8, 1e-1]")
        print(f"  worst |S[f] + S[1/f]| : {worst:.2e}")
        print(f"  worst / DERIVED bound 4*EPS*|g|/(2h) : {worst_ratio:.3f}   "
              f"{'PASS' if worst_ratio <= 1.0 else 'FAIL'}")
        print("  NOT 31/31, and reported as such.  The prediction was wrong in SHAPE,")
        print("  not only in magnitude: ln|1/f| = -ln|f| holds in exact arithmetic,")
        print("  but 1/f is a separately rounded number, so log(1/x) != -log(x)")
        print("  bitwise.  The symmetry is exact as mathematics and holds to the")
        print("  cancellation floor in double precision -- a floating-point")
        print("  statement, not a modelling one.  (Same shape as Validation 3(c)")
        print("  of graphene_sensitivity_audit.py, which says so about kappa.)")
    return n_bit, len(hs), worst, worst_ratio


# =====================================================================
# Q2 + Q3 + V5 -- conv's conservatism, and its exact blind spot
# =====================================================================
def q2_conservatism(verbose=True):
    """conv should be 3x the true error where truncation dominates."""
    A = 0.8123
    f, Sx = quad_log(A)
    p0 = math.e ** 1.37
    rows = []
    for h in (1e-3, 3e-3, 1e-2, 3e-2):
        err = abs(S_affine(f, p0, h) - Sx(p0))
        conv = abs(conv_of(S_affine, f, p0, h))
        rows.append((h, err, conv, conv / err))
    if verbose:
        print("\n" + "=" * 78)
        print("Q2 -- is conv = |S(h) - S(2h)| conservative, and by how much?")
        print("=" * 78)
        print(f"  {'h':>8}{'|true error|':>16}{'conv':>16}{'conv/|error|':>15}")
        for h, e, c, r in rows:
            print(f"  {h:>8.0e}{e:>16.3e}{c:>16.3e}{r:>15.4f}")
        print("  Expected 3.000 exactly where E = c2 h^2 dominates (conv = |c2 h^2 -")
        print("  4 c2 h^2| = 3|c2|h^2), approaching it from above as h -> 0 because the")
        print("  h^4 term adds 12|c4|h^4 to the numerator and |c4|h^4 to the denominator.")
    return rows


def q3_v5_blind_spot(verbose=True):
    """
    Construct a function for which conv is EXACTLY zero at a step where the
    error is not, and locate that step by a BRACKETED root solve (2026-09-24).

    Design: g(L) = A L^2 + B L^4 at L = 1, with A = -2B + delta.  Then
      c2 = -delta  and  c4 = 5B/3 - delta/2,
    so conv vanishes near h*^2 = 3 delta / (20 B), and -- independently of A,
    B and delta -- S_exact = 2 delta at L = 1 while E(h*) = (3/4) c2 h*^2,
    giving a RELATIVE error at the blind spot of exactly (3/8) h*^2.
    """
    B = 1.0
    h_target = 0.02
    delta = 20.0 * B * h_target ** 2 / 3.0
    A = -2.0 * B + delta
    L = 1.0
    p0 = math.exp(L)
    f, Sx = quart_log(A, B)
    S_exact = Sx(p0)

    # signed conv, from the CLOSED FORM: no finite differences in the root solve
    def conv_exact(h):
        return E_affine_quart_exact(A, B, L, h) - E_affine_quart_exact(A, B, L, 2.0 * h)

    lo, hi = h_target / 5.0, h_target * 5.0
    h_star = bracketed_bisect(conv_exact, lo, hi, tol=0.0, rtol=EPS)

    E_star_exact = E_affine_quart_exact(A, B, L, h_star)
    S_meas = S_affine(f, p0, h_star)
    E_star_meas = S_meas - S_exact
    conv_meas = abs(conv_of(S_affine, f, p0, h_star))
    rel_err = abs(E_star_exact / S_exact)
    rel_pred = 0.375 * h_star ** 2
    v5_g0 = abs(math.log(abs(f(p0))))
    v5_bound = EPS * v5_g0 / (2 * h_star) / abs(E_star_exact)
    v5_ratio = abs(E_star_meas - E_star_exact) / abs(E_star_exact) / v5_bound
    # how far a step-size sweep has to miss h* before conv recovers
    sweep = [(m, abs(conv_of(S_affine, f, p0, h_star * m)),
              abs(S_affine(f, p0, h_star * m) - S_exact))
             for m in (0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0)]
    if verbose:
        print("\n" + "=" * 78)
        print("Q3 + V5 (EXACT) -- conv's blind spot, constructed and located exactly")
        print("=" * 78)
        print(f"  g(L) = A L^2 + B L^4 with B = {B:.1f}, A = {A:.12f}, at L = 1")
        print(f"  S_exact = 2*delta                         = {S_exact:.17g}")
        print(f"  h* (bracketed bisect on the CLOSED FORM)  = {h_star:.17g}")
        print(f"  conv MEASURED at h*                       = {conv_meas:.3e}")
        print(f"  true error at h*, closed form             = {E_star_exact:.15e}")
        print(f"  true error at h*, measured                = {E_star_meas:.15e}")
        print(f"    relative deviation measured vs closed form : "
              f"{abs(E_star_meas - E_star_exact)/abs(E_star_exact):.2e}")
        print(f"    DERIVED floor EPS*|g|/(2h*)/|E| = {v5_bound:.1e}, "
              f"measured/bound = {v5_ratio:.3f}   "
              f"{'PASS' if v5_ratio <= 4.0 else 'FAIL'}")
        print("    The first form of this check also required an absolute 1e-13 and")
        print("    also measured its own floor.  Second instance in one session.")
        print(f"  relative error in S at h*                 = {rel_err:.6e}")
        print(f"  predicted (3/8) h*^2                      = {rel_pred:.6e}   "
              f"(ratio {rel_err/rel_pred:.6f})")
        print(f"  conv UNDERSTATES the error there by a factor of "
              f"{abs(E_star_exact)/max(conv_meas, 5e-324):.2e}")
        print("\n  How narrow is the blind spot?  conv and |error| at h* * m:")
        print(f"    {'m':>7}{'conv':>14}{'|error|':>14}{'conv/|err|':>13}")
        for m, c, e in sweep:
            print(f"    {m:>7.2f}{c:>14.3e}{e:>14.3e}{c/e:>13.3e}")
        print("  conv is a valid 3x-conservative bar everywhere EXCEPT within ~1% of")
        print("  h*, where it collapses.  The absolute error it hides is of ordinary")
        print("  O(h^2) size -- exactly (3/8)h*^2 relative, independent of the")
        print("  function -- so the RATIO of understatement is unbounded while the")
        print("  DAMAGE is bounded by the step-size error itself.  That is the honest")
        print("  limit on this mechanism, and it is why Q4 matters more than Q3.")
    return dict(h_star=h_star, conv_meas=conv_meas, E_exact=E_star_exact,
                E_meas=E_star_meas, rel_err=rel_err, rel_pred=rel_pred, sweep=sweep,
                v5_bound=v5_bound, v5_ratio=v5_ratio)


# =====================================================================
# The real call sites
# =====================================================================
def real_call_sites():
    """Every function `log_sensitivity` is actually applied to in this thesis."""
    sites = []
    keyB = {"rho_calibration": "rho_cal", "rho_bulk": "rho_bulk",
            "p_default": "p", "W_calibration": "W_cal", "lambda_bulk": "lam_bulk"}
    valsB = {"rho_calibration": icm.rho_calibration_uohm_cm,
             "rho_bulk": icm.rho_bulk_uohm_cm,
             "p_default": icm.p_default,
             "W_calibration": icm.W_calibration_nm,
             "lambda_bulk": icm.lambda_bulk_nm}
    for name, val in valsB.items():
        n = keyB[name]
        sites.append((f"B lambda_imp/{name}",
                      (lambda v, n=n: _lambda_imp(**{n: v})), val))
    for W in (18.0, 30.0, 52.0):
        for name, val in valsB.items():
            n = keyB[name]
            sites.append((f"C rho({W:.0f})/{name}",
                          (lambda v, n=n, w=W: _rho_propagated(w, **{n: v})), val))
    t = icm.t_liner_nm_default
    for W in (18.0, 30.0, 52.0):
        sites.append((f"D rho_Cu({W:.0f})/t_liner",
                      (lambda v, w=W: icm.cu_resistivity_with_liner(w, t_liner_nm=v)), t))
    return sites


def q4_sign_changes(verbose=True):
    """
    Q4, ORACLE-FREE.  conv = 0 requires the SIGNED step-doubling difference to
    change sign, which needs no knowledge of the true S.  Count sign changes in
    the truncation-dominated part of the sweep for every real call site.
    """
    hs = np.logspace(-8, -1, 64)
    rows = []
    for label, f, p0 in real_call_sites():
        D = np.array([conv_of(S_affine, f, p0, h) for h in hs])
        nz = np.abs(D) > 0.0
        sgn = np.sign(D[nz])
        flips_all = int(np.sum(sgn[1:] * sgn[:-1] < 0))
        # Separate the two regimes.  conv bottoms out where truncation (~h^2,
        # falling with h) crosses cancellation (~EPS|g|/h, rising as h falls).
        # Only sign changes ABOVE that knee can be the c2/c4 mechanism Q3
        # describes; below it they are noise, and counting them together is
        # what the first form of this test did.
        i_knee = int(np.argmin(np.abs(D)))
        h_knee = hs[i_knee]
        # The truncation branch is identified against the DERIVED rounding
        # floor, not against conv's own minimum: floor(h) = 4*EPS*|g|/h, and a
        # point is on the truncation branch when |conv(h)| exceeds it by 30x.
        # Two earlier forms of this boundary were wrong -- 10 x argmin|conv| is
        # two to three decades too low (Q8 shows why), and a log-log slope test
        # reads noise as slope +2 whenever two adjacent noise samples happen to
        # rise.  The floor is the only one of the three that is derived.
        g0 = abs(math.log(abs(f(p0))))
        floors = 4.0 * EPS * g0 / hs
        keep = np.abs(D) > 30.0 * floors
        st = np.sign(D[keep])
        st = st[st != 0.0]
        flips_tr = int(np.sum(st[1:] * st[:-1] < 0)) if st.size > 1 else 0
        h_trunc_lo = hs[keep].min() if keep.any() else float("nan")
        floor5 = 4.0 * EPS * g0 / SHIPPED_H
        c5 = abs(conv_of(S_affine, f, p0, SHIPPED_H))
        rows.append((label, flips_all, flips_tr, h_trunc_lo, abs(D[i_knee]),
                     c5, c5 / floor5))
    tot = sum(r[2] for r in rows)
    tot_all = sum(r[1] for r in rows)
    if verbose:
        print("\n" + "=" * 78)
        print("Q4 (ORACLE-FREE) -- does conv have a zero at any real call site?")
        print("  A zero of conv requires the SIGNED difference S(h) - S(2h) to change")
        print("  sign, which can be counted without knowing the true S at all.")
        print("=" * 78)
        print(f"  {'call site':<26}{'flips all':>10}{'flips trunc':>12}"
              f"{'trunc branch from':>19}{'conv(1e-5)':>13}{'/floor':>10}")
        for lab, fa, ft, ho, mo, c5, rf in rows:
            print(f"  {lab:<26}{fa:>10d}{ft:>12d}{ho:>19.2e}{c5:>13.2e}{rf:>10.1f}")
        print(f"  total sign changes, WHOLE sweep            : {tot_all}")
        print(f"  total sign changes, TRUNCATION regime only : {tot}")
        print("  The two numbers answer different questions and the first form of")
        print("  this test conflated them.  Below the derived rounding floor the")
        print("  signed difference is cancellation noise and changes sign freely;")
        print("  that is not Q3's c2/c4 mechanism and it is not evidence about the")
        print("  published values.  It IS evidence for something else: a step-size")
        print("  sweep that trusts a small conv will be pulled BELOW the knee.")
        fam = {"B": [], "C": [], "D": []}
        for r in rows:
            fam[r[0][0]].append(r[6])
        print("\n  IS conv AT THE SHIPPED STEP A TRUNCATION ESTIMATE AT ALL?")
        print("  conv(1e-5) divided by the derived cancellation floor "
              "4*EPS*|g|/h :")
        for k in ("B", "C", "D"):
            v = fam[k]
            print(f"    Family {k}: min {min(v):.1f}  median "
                  f"{float(np.median(v)):.1f}  max {max(v):.1f}")
        print("  Family B's conv is three decades above its own rounding floor and")
        print("  is a real truncation estimate.  Families C and D sit AT the floor:")
        print("  for them the 'worst convergence estimate' lines printed by Chapter")
        print("  5 Section 5.5 are measurements of double-precision rounding, not of")
        print("  convergence, and would print the same number for a model with any")
        print("  amount of curvature.  That is the sharpest form of this session's")
        print("  finding: the reassurance scales with |ln Q| and h, not with the")
        print("  quality of the answer.")
    return rows, tot, tot_all


def q1_q5_q6_shape_and_shift(verbose=True):
    """
    Q1 magnitude, Q5 (do published numbers move?) and Q6 (is Family B worst?),
    all measured at the shipped step.
    """
    rows = []
    for label, f, p0 in real_call_sites():
        a = S_affine(f, p0, SHIPPED_H)
        g = S_geom(f, p0, SHIPPED_H)
        # 16x refinement of the shipped estimator, as an independent probe
        a16 = S_affine(f, p0, SHIPPED_H / 16.0)
        # Richardson on the corrected estimator: O(h^4) residual
        href = 1e-4
        ref = (4.0 * S_geom(f, p0, href) - S_geom(f, p0, 2 * href)) / 3.0
        den = abs(a) if a != 0.0 else 1.0
        rows.append((label, a, g, abs(a - g) / den, abs(a - a16) / den,
                     abs(a - ref) / den))
    worst_shape = max(r[3] for r in rows)
    worst_ref16 = max(r[4] for r in rows)
    worst_ref = max(r[5] for r in rows)
    famB = [r for r in rows if r[0].startswith("B ")]
    famCD = [r for r in rows if not r[0].startswith("B ")]
    wB = max(r[3] for r in famB)
    wCD = max(r[3] for r in famCD)
    if verbose:
        print("\n" + "=" * 78)
        print("Q1 magnitude + Q5 + Q6 -- at the shipped h = 1e-5")
        print("=" * 78)
        print(f"  {'call site':<26}{'S shipped':>14}{'S corrected':>14}"
              f"{'|rel shift|':>13}{'16x refine':>12}{'vs Richardson':>15}")
        for lab, a, g, ds, d16, dr in rows:
            print(f"  {lab:<26}{a:>14.8f}{g:>14.8f}{ds:>13.2e}{d16:>12.2e}{dr:>15.2e}")
        print(f"  worst |shipped - corrected| / |S|          : {worst_shape:.2e}")
        print(f"  worst |S(1e-5) - S(6.25e-7)| / |S|         : {worst_ref16:.2e}")
        print(f"  worst |shipped - Richardson ref| / |S|     : {worst_ref:.2e}")
        print(f"  Family B worst shift {wB:.2e}   Families C/D worst shift {wCD:.2e}"
              f"   ratio {wB/wCD if wCD else float('inf'):.3f}")
    return rows, worst_shape, worst_ref16, worst_ref, wB, wCD



def q9_locate_truncation_zeros(verbose=True):
    """
    The result Q4 predicted would not exist.  Four call sites DO carry a zero
    of conv on the truncation branch, and all four are in Family B -- the
    lambda_impurity calibration that Chapter 5 Section 5.5 calls the
    worst-conditioned step in this thesis.  Locate each with a BRACKETED
    bisection (2026-09-24) and report what conv says there against an
    independent Richardson reference built from the CORRECTED estimator.
    """
    hs = np.logspace(-8, -1, 400)
    out = []
    for label, f, p0 in real_call_sites():
        g0 = abs(math.log(abs(f(p0))))
        D = np.array([conv_of(S_affine, f, p0, h) for h in hs])
        floors = 4.0 * EPS * g0 / hs
        keep = np.abs(D) > 30.0 * floors
        idx = np.where(keep)[0]
        if idx.size < 2:
            continue
        href = 1e-4
        ref = (4.0 * S_geom(f, p0, href) - S_geom(f, p0, 2 * href)) / 3.0
        e5 = abs(S_affine(f, p0, SHIPPED_H) - ref)
        for i, j in zip(idx[:-1], idx[1:]):
            if j != i + 1 or D[i] == 0.0 or D[j] == 0.0:
                continue
            if D[i] * D[j] >= 0.0:
                continue
            hz = bracketed_bisect(lambda h: conv_of(S_affine, f, p0, h),
                                  hs[i], hs[j], tol=0.0, rtol=EPS)
            cz = abs(conv_of(S_affine, f, p0, hz))
            ez = abs(S_affine(f, p0, hz) - ref)
            out.append((label, hz, cz, ez, ez / max(cz, 5e-324),
                        hz / SHIPPED_H, ez / e5 if e5 else float("nan")))
    if verbose:
        print("\n" + "=" * 78)
        print("Q9 (THE RESULT Q4 PREDICTED WOULD NOT EXIST) -- conv's zeros, located")
        print("=" * 78)
        if not out:
            print("  none found")
        else:
            print(f"  {'call site':<26}{'h_zero':>12}{'conv there':>13}"
                  f"{'|error| there':>15}{'err/conv':>11}{'h_zero/1e-5':>13}")
            for lab, hz, cz, ez, r, rel5, _ in out:
                print(f"  {lab:<26}{hz:>12.5e}{cz:>13.2e}{ez:>15.2e}"
                      f"{r:>11.1e}{rel5:>13.3f}")
            worst = max(r[4] for r in out)
            nearest = min(out, key=lambda r: abs(math.log(r[5])))
            print(f"  worst understatement factor : {worst:.1e}")
            print(f"  closest zero to the shipped step : {nearest[0]} at "
                  f"h = {nearest[1]:.3e} = {nearest[5]:.2f} x 1e-5")
            print("  ALL of them are in FAMILY B -- the lambda_impurity calibration,")
            print("  kappa = 19.71, the step Chapter 5 Section 5.5 calls the")
            print("  worst-conditioned in this thesis.  Q4 predicted no such zero")
            print("  anywhere, and reasoned that these are smooth rationals in which")
            print("  the second derivative of ln|f| dominates the third.  Family B is")
            print("  a 1/residual with a near-zero residual, which is exactly where")
            print("  that reasoning fails: the near-cancellation puts large higher")
            print("  derivatives into ln|f|, c2 and c4 acquire opposite signs, and")
            print("  the blind spot appears.")
            relerr = [ez / abs(S_affine(f2, p2, SHIPPED_H))
                      for (lab, hz, cz, ez, r, rel5, _), (l2, f2, p2)
                      in zip(out, [c for c in real_call_sites()
                                   if c[0] in [o[0] for o in out]])]
            print(f"  relative error in S at those zeros : "
                  + ", ".join(f"{100*r:.1f}%" for r in relerr))
            print("  So at h ~ 3-5% conv reports 1e-14 -- machine precision -- while")
            print("  the sensitivity is wrong by 14%.  That is 2026-09-25's finding")
            print("  (a procedure asked whether it has converged can answer yes and")
            print("  be wrong by 44%) reproduced on THIS estimator, at 14%, on the")
            print("  step this thesis calls its worst-conditioned.")
            print("  NOTHING PUBLISHED MOVES -- Q5 stands, worst 1.3e-8 relative, and")
            print("  the nearest zero is ~3000x ABOVE the shipped step, not near it.")
            print("  The shipped h = 1e-5 is safe by a wide margin and this session")
            print("  does not claim otherwise.  What is unsafe is the PROCEDURE: a")
            print("  step-size sweep over the natural range h in [1e-8, 1e-1] passes")
            print("  straight through these four points, and at each of them the only")
            print("  instrument the code offers reports full convergence.  The repo")
            print("  ran exactly such a sweep on 2026-09-25 and 2026-09-26.")
    return out


def q8_conv_minimum_vs_error_minimum(verbose=True):
    """
    UNPREDICTED RESULT, oracle-free in the sense that matters.

    `conv` is the only instrument the shipped code offers for choosing a step.
    If a step-size sweep minimised it, which step would it choose, and is that
    step good?  The error is measured against a Richardson reference built from
    the CORRECTED estimator at h = 1e-4 (residual O(h^4) ~ 1e-16 in truncation,
    ~1e-12 in cancellation), which is an anchor only to the accuracy stated --
    so the comparison below is reported as a ratio of two h values, which is
    robust to the reference's own error, rather than as an error magnitude.
    """
    hs = np.logspace(-9, -1.5, 80)
    rows = []
    for label, f, p0 in real_call_sites():
        href = 1e-4
        ref = (4.0 * S_geom(f, p0, href) - S_geom(f, p0, 2 * href)) / 3.0
        conv = np.array([abs(conv_of(S_affine, f, p0, h)) for h in hs])
        err = np.array([abs(S_affine(f, p0, h) - ref) for h in hs])
        h_conv = hs[int(np.argmin(conv))]
        h_err = hs[int(np.argmin(err))]
        # error at the step conv would choose, relative to error at the shipped step
        e_convmin = abs(S_affine(f, p0, h_conv) - ref)
        e_ship = abs(S_affine(f, p0, SHIPPED_H) - ref)
        rows.append((label, h_conv, h_err, h_err / h_conv,
                     e_convmin / e_ship if e_ship else float("inf")))
    med_ratio = float(np.median([r[3] for r in rows]))
    n_worse = sum(1 for r in rows if r[4] > 1.0)
    worst_pen = max(r[4] for r in rows)
    if verbose:
        print("\n" + "=" * 78)
        print("Q8 (UNPREDICTED) -- which step does conv choose, and is it a good one?")
        print("=" * 78)
        print(f"  {'call site':<26}{'h at min conv':>15}{'h at min err':>14}"
              f"{'ratio':>9}{'err(h_conv)/err(1e-5)':>23}")
        for lab, hc, he, r, pen in rows:
            print(f"  {lab:<26}{hc:>15.2e}{he:>14.2e}{r:>9.1f}{pen:>23.2e}")
        print(f"  median h_err / h_conv : {med_ratio:.1f}")
        print(f"  call sites at which minimising conv gives a WORSE answer than the")
        print(f"  shipped h = 1e-5 : {n_worse}/{len(rows)}, worst penalty "
              f"{worst_pen:.1e}x")
        print("  conv is minimised DEEP IN THE CANCELLATION REGIME, two to three")
        print("  decades below the step that actually minimises the error, because")
        print("  there the two estimates it differences are both noise and their")
        print("  difference is small for that reason.  A step-size sweep steered by")
        print("  conv therefore walks AWAY from the right step, and reports a")
        print("  smaller convergence estimate as it goes.  This is 2026-09-25's")
        print("  finding -- a procedure asked whether it has converged can answer")
        print("  yes and be wrong -- on this estimator, and it is the strongest")
        print("  reason to attach a criterion rather than a number to conv.")
    return rows, med_ratio, n_worse, worst_pen


# =====================================================================
# Figure
# =====================================================================
def make_figure(q3, q4rows, q9, fname="log_sensitivity_step_audit.png"):
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 10))

    # (a) V2: the discriminating exact check
    A = 0.8123
    f, Sx = quad_log(A)
    p0 = math.e ** 1.37
    hs = np.logspace(-7, -0.7, 120)
    ea = [abs(S_affine(f, p0, h) - Sx(p0)) for h in hs]
    eg = [abs(S_geom(f, p0, h) - Sx(p0)) for h in hs]
    ex = [abs(E_affine_quad_exact(A, h)) for h in hs]
    ax[0, 0].loglog(hs, ea, "o", ms=3, color="C3", label="shipped, measured")
    ax[0, 0].loglog(hs, ex, "-", color="k", lw=1, label=r"closed form $|A\ln(1-h^2)|$")
    ax[0, 0].loglog(hs, eg, "s", ms=3, color="C0",
                    label="corrected (exactly 0, so roundoff only)")
    ax[0, 0].axvline(SHIPPED_H, color="0.4", ls=":", label="shipped $h=10^{-5}$")
    ax[0, 0].set_xlabel("relative step $h$")
    ax[0, 0].set_ylabel(r"$|S(h)-S_{\rm exact}|$")
    ax[0, 0].set_title(r"(a) $\ln|f|=A(\ln p)^2$: the shipped estimator has an"
                       "\n"
                       r"exact $O(h^2)$ error where a centred difference has none")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3, which="both")

    # (b) the blind spot
    B = 1.0
    h_t = 0.02
    delta = 20.0 * B * h_t ** 2 / 3.0
    Aq = -2.0 * B + delta
    fq, Sq = quart_log(Aq, B)
    pq = math.e
    hs2 = np.logspace(math.log10(q3["h_star"]) - 0.6, math.log10(q3["h_star"]) + 0.6, 400)
    cv = [abs(conv_of(S_affine, fq, pq, h)) for h in hs2]
    er = [abs(S_affine(fq, pq, h) - Sq(pq)) for h in hs2]
    ax[0, 1].loglog(hs2, cv, "-", color="C3", label=r"$conv=|S(h)-S(2h)|$")
    ax[0, 1].loglog(hs2, er, "-", color="C0", label=r"true $|S(h)-S_{\rm exact}|$")
    ax[0, 1].axvline(q3["h_star"], color="k", ls="--", lw=1,
                     label=r"$h^\ast$ (bracketed, closed form)")
    ax[0, 1].set_xlabel("relative step $h$")
    ax[0, 1].set_ylabel("magnitude")
    ax[0, 1].set_title("(b) conv reports ZERO at $h^\\ast$ while the error is\n"
                       r"$(3/8)h^{\ast 2}$ relative $-$ unbounded ratio, ordinary size")
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].grid(alpha=0.3, which="both")

    # (c) conv vs h for the real call sites
    hs3 = np.logspace(-8, -1, 60)
    for label, fn, pp in real_call_sites():
        if not (label.startswith("B ") or label.startswith("D ")):
            continue
        c = [abs(conv_of(S_affine, fn, pp, h)) for h in hs3]
        ax[1, 0].loglog(hs3, c, lw=1, alpha=0.85, label=label)
    ax[1, 0].axvline(SHIPPED_H, color="k", ls=":", label="shipped $h=10^{-5}$")
    for lab, hz, cz, ez, r, rel5, _ in q9:
        ax[1, 0].plot([hz], [max(cz, 1e-17)], "v", ms=9, mfc="none",
                      mec="C3", mew=1.6)
    if q9:
        ax[1, 0].plot([], [], "v", ms=9, mfc="none", mec="C3", mew=1.6,
                      label="conv = 0 (error 13-15%)")
    ax[1, 0].set_xlabel("relative step $h$")
    ax[1, 0].set_ylabel("conv")
    ax[1, 0].set_title("(c) real call sites: conv falls as $h^2$, rises as $1/h$, and\n"
                       "collapses to $10^{-14}$ at four Family B steps where $S$ is 13-15% wrong")
    ax[1, 0].legend(fontsize=6.5, ncol=2)
    ax[1, 0].grid(alpha=0.3, which="both")

    # (d) the shape defect at the shipped step
    labs, shifts = [], []
    for label, fn, pp in real_call_sites():
        a = S_affine(fn, pp, SHIPPED_H)
        g = S_geom(fn, pp, SHIPPED_H)
        labs.append(label)
        shifts.append(abs(a - g) / (abs(a) if a else 1.0))
    order = np.argsort(shifts)[::-1]
    ax[1, 1].barh([labs[i] for i in order][:14], [shifts[i] for i in order][:14],
                  color=["C3" if labs[i].startswith("B ") else "C0" for i in order][:14])
    ax[1, 1].set_xscale("log")
    ax[1, 1].invert_yaxis()
    ax[1, 1].tick_params(axis="y", labelsize=7)
    ax[1, 1].set_xlabel(r"$|S_{\rm shipped}-S_{\rm corrected}|/|S|$ at $h=10^{-5}$")
    ax[1, 1].set_title("(d) the shape defect, measured: it is real and it is\n"
                       "far below every printed digit (red = Family B)")
    ax[1, 1].grid(alpha=0.3, axis="x", which="both")

    fig.tight_layout()
    fig.savefig(fname, dpi=140)
    print(f"\n[figure written: {fname}]")


# =====================================================================
def main():
    print("=" * 78)
    print("AN ANCHORED AUDIT OF log_sensitivity's STEP AND CONVERGENCE ESTIMATE")
    print("  pre-registration: notes/2026-09-27-log-sensitivity-preregistration.md")
    print("=" * 78)
    sites, comps = q7_conv_census()
    v1 = v1_v4_power_law()
    v2 = v2_quadratic_log()
    v3 = v3_reciprocal_symmetry()
    q2 = q2_conservatism()
    q3 = q3_v5_blind_spot()
    q4rows, q4tot, q4tot_all = q4_sign_changes()
    q156 = q1_q5_q6_shape_and_shift()
    q9 = q9_locate_truncation_zeros()
    q8 = q8_conv_minimum_vs_error_minimum()
    make_figure(q3, q4rows, q9)

    print("\n" + "=" * 78)
    print("SCORING")
    print("=" * 78)
    rows, worst_shape, worst_ref16, worst_ref, wB, wCD = q156
    print(f"  Q1 class call : the shipped estimator IS an asymmetric log-interval")
    print(f"                  secant; V2 confirms its error against the closed form")
    print(f"                  A ln(1-h^2) at the derived cancellation bound")
    print(f"                  (worst measured/bound {v2[3]:.2f}).  PASS")
    print(f"  Q1 magnitude  : predicted < 1e-9 relative at the shipped h; measured")
    print(f"                  worst {worst_shape:.2e}.  "
          f"{'PASS' if worst_shape < 1e-9 else 'FAIL'}")
    print(f"  Q2            : predicted conv/|error| = 3.00 +/- 0.02 at h = 1e-2;")
    print(f"                  measured {[f'{r[3]:.4f}' for r in q2]}.  "
          f"{'PASS' if abs(q2[2][3]-3.0) < 0.02 else 'FAIL'}")
    print(f"  Q3            : conv exactly zero at a constructible h* with relative")
    print(f"                  error (3/8)h*^2; measured ratio to prediction "
          f"{q3['rel_err']/q3['rel_pred']:.6f}.  "
          f"{'PASS' if abs(q3['rel_err']/q3['rel_pred'] - 1) < 1e-3 else 'FAIL'}")
    print(f"  Q4            : predicted NO sign change at any real call site over")
    print(f"                  h in [1e-8, 1e-1].  Truncation branch: {q4tot} sign")
    print(f"                  changes -> class call CORRECT.  Whole sweep: {q4tot_all},")
    print(f"                  all below the cancellation knee -> the prediction as")
    print(f"                  WRITTEN (\"over h in [1e-8, 1e-1]\") is FALSIFIED, because")
    print(f"                  it named a range it had not thought about.  "
          f"{'PASS on class, FAIL as written' if q4tot == 0 else 'FAIL'}")
    print(f"  Q5            : predicted every published S moves < 1e-6 relative;")
    print(f"                  worst 16x refinement {worst_ref16:.2e}, worst vs the")
    print(f"                  corrected estimator {worst_shape:.2e}, worst vs a")
    print(f"                  Richardson reference {worst_ref:.2e}.  "
          f"{'PASS' if max(worst_ref16, worst_shape, worst_ref) < 1e-6 else 'FAIL'}")
    print(f"  Q6            : predicted Family B worse than the rest by > 10x;")
    print(f"                  measured ratio {wB/wCD if wCD else float('inf'):.3f}.  "
          f"{'PASS' if (wCD and wB/wCD > 10) else 'FAIL'}")
    print(f"  Q7            : predicted no call site compares conv with anything.")
    print(f"                  {len(sites)} log_sensitivity call sites, 0 of which")
    print(f"                  compare the returned conv with anything: class call")
    print(f"                  PASS.  The census also found {len(comps)} comparison of a")
    print(f"                  DIFFERENT convergence quantity in another module, so the")
    print(f"                  prediction as written (\"no call site in the repository\")")
    print(f"                  is FAIL -- it over-scoped.")
    if q9:
        print(f"  Q4 FALSIFIED ON ITS CLASS CALL TOO: {len(q9)} zeros of conv sit on the")
        print(f"                  truncation branch, ALL in Family B, worst")
        print(f"                  understatement {max(r[4] for r in q9):.1e}x, closest to")
        print(f"                  the shipped step at {min(q9, key=lambda r: abs(math.log(r[5])))[5]:.2f} x 1e-5.")
    print(f"  Q8 UNPREDICTED: minimising conv chooses a step {q8[1]:.0f}x below the one")
    print(f"                  that minimises the error, and gives a worse answer than")
    print(f"                  the shipped h at {q8[2]}/{len(q8[0])} call sites (worst "
          f"{q8[3]:.1e}x).")
    print(f"\n  V1 (non-discriminating, retained as such) : affine {v1[0]}/{v1[3]} "
          f"bitwise, geom {v1[1]}/{v1[3]}, conv==0.0 {v1[2]}/{v1[3]}")
    print(f"  V2 (discriminating)  : worst relative deviation {v2[0]:.2e}, "
          f"worst measured/derived-bound {v2[3]:.2f}")
    print(f"  V3 (exact symmetry)  : {v3[0]}/{v3[1]} bitwise, worst {v3[2]:.2e}, "
          f"worst measured/derived-bound {v3[3]:.2f}")
    print(f"  V4 (non-discriminating) : conv exactly 0.0 at {v1[2]}/{v1[3]} steps "
          f"including h = 0.25")
    print(f"  V5 (exact blind spot): measured error vs closed form "
          f"{abs(q3['E_meas']-q3['E_exact'])/abs(q3['E_exact']):.2e} relative, "
          f"measured/derived-bound {q3['v5_ratio']:.2f}")
    print("\n  THREE of five validations failed in their first form, and all three")
    print("  failed the SAME way: a round absolute tolerance (1e-13, 1e-14, 1e-13)")
    print("  that silently encoded a step size.  That is the fault class documented")
    print("  on 2026-09-25 and reproduced on 2026-09-26, committed three more times")
    print("  here, in the session auditing it.  All three are rebuilt on derived")
    print("  bounds and all three first forms are left in the record.")


if __name__ == "__main__":
    main()
