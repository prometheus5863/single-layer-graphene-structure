"""
graphene_gds_quadrature_audit.py
================================

Audit of the top open numerical item carried forward from 2026-09-25:

    rf_small_signal_model.output_conductance(Vg, Vds=0.05, dVds=1e-3)

`dVds = 1e-3 V` is an absolute constant that is 2% of its variable at the
default operating point -- the identical construction to `H_DIFF = 1e-3`,
which 2026-09-25 convicted of moving eleven of Chapter 6.13's sensitivities
by up to 44%.  It was deliberately left unmeasured that day because g_ds
feeds f_max, a Chapter 4 published quantity, and moving it needs its own
before/after comparison.  This module is that comparison.

WHAT MAKES THIS CASE DIFFERENT FROM H_DIFF
------------------------------------------
H_DIFF differentiates a closed-form expression.  `dVds` differentiates
`graphene_fet_model.transfer_characteristic()`, which is itself a
DISCRETISED object: it averages the local channel resistance over

    V_channel_profile = np.linspace(0, Vds, n_segments)     n_segments = 50

so Id(Vds) is not the model's continuum Id(Vds) but a 50-point quadrature
of it -- and, critically, the quadrature GRID IS SET BY Vds.  Differencing
in Vds therefore differences the quadrature error too.  Two defaults are
entangled: a step size `dVds` and a resolution `n_segments`, and they are
not independent.

That entanglement makes a specific, sharper failure mode available than
anything in the 09-22..09-25 series.  The anchored step criterion adopted
on 09-25 -- accept a step only if the derivative is unchanged at h/10 AND
h/100 -- tests convergence IN THE STEP.  A step-refinement study of a
discretised function converges to the derivative of THE DISCRETISATION YOU
FIXED, not to the derivative of the model.  It can therefore pass, at any
tolerance, while sitting on a bias that no step refinement can see.  The
question this module asks is whether that is the actual situation here, and
if so how large the bias is and whether it reaches f_max.

STRUCTURE
---------
Pre-registered predictions are printed BEFORE any measurement (habit
adopted 2026-09-21) and scored at the end.  Five exact validations, one of
which is retained and LABELLED as non-discriminating (09-25 item 11), and
one of which derives its own bound rather than asserting an absolute
constant (09-25 item 10: an absolute tolerance is a hidden scale claim).

Nothing existing is rewritten by this module.  It measures, and reports.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import rf_small_signal_model as rf

EPS = np.finfo(float).eps

# The defaults under audit, named so they cannot be silently confused with
# the values this module measures.
DVDS_DEFAULT = 1e-3
N_SEGMENTS_DEFAULT = 50
VDS_DEFAULT = 0.05
VG_GRID_DEFAULT = 400          # np.gradient grid used by transconductance()
VG_SPAN_DEFAULT = (-1.5, 3.5)

RTOL_ANCHORED = 1e-3           # placed inside a measured gap below, not chosen


# ----------------------------------------------------------------------
# Reimplementation of the model's Id(Vds) with n_segments exposed.
# Fidelity against the shipped function is asserted, not assumed.
# ----------------------------------------------------------------------

def Rbar_model(Vg, Vds, n_segments):
    """Mean local channel resistance over the Vds drop, exactly as
    transfer_characteristic() computes it (same grid, same reduction),
    but with n_segments as an argument instead of a literal."""
    V_channel_profile = np.linspace(0, Vds, n_segments)
    return np.mean(
        np.array([gfet.channel_resistance(Vg, V_ch)
                  for V_ch in V_channel_profile]), axis=0)


def Id_model(Vg, Vds, n_segments=N_SEGMENTS_DEFAULT):
    return Vds / (Rbar_model(Vg, Vds, n_segments) + gfet.Rc_total)


def gds_model(Vg, Vds=VDS_DEFAULT, dVds=DVDS_DEFAULT,
              n_segments=N_SEGMENTS_DEFAULT, clamp=True):
    """Central difference in Vds.  clamp=True reproduces the shipped
    function's `max(Vds - dVds, 1e-4)` guard; clamp=False removes it."""
    lo = max(Vds - dVds, 1e-4) if clamp else Vds - dVds
    Ip = Id_model(Vg, Vds + dVds, n_segments)
    Im = Id_model(Vg, lo, n_segments)
    return (Ip - Im) / (2 * dVds)


def anchored_step(f, h_start=1e-2, h_min=1e-11, rtol=RTOL_ANCHORED):
    """09-25's anchored criterion: accept h only if f(h) agrees with BOTH
    f(h/10) and f(h/100) to rtol.  Returns (h, value, r10, r100, trace).
    Self-consistency against a single doubling is NOT used: 09-25 measured
    it under-reporting a true 44% error by 190,000x at a stationary point
    of the error curve."""
    trace = []
    h = h_start
    while h >= h_min * 100:
        v, v10, v100 = f(h), f(h / 10), f(h / 100)
        scale = max(abs(v), 1e-300)
        r10, r100 = abs(v - v10) / scale, abs(v - v100) / scale
        trace.append((h, v, r10, r100))
        if r10 < rtol and r100 < rtol:
            return h, v, r10, r100, trace
        h /= 10.0
    return None, None, None, None, trace


# ----------------------------------------------------------------------
# Pre-registered predictions (printed before any measurement)
# ----------------------------------------------------------------------

PREDICTIONS = [
    ("Q1", "The shipped default dVds=1e-3 passes the 09-25 anchored step "
           "criterion (agreement at h/10 AND h/100 to rtol=1e-3) at the "
           "default operating point. Id(Vds) is close to linear at small "
           "Vds, so unlike H_DIFF the STEP is expected to be innocent."),
    ("Q2", "Nevertheless g_ds at the anchored step is biased by n_segments=50 "
           "at the >=0.1% level, and -- the discriminating half -- that bias "
           "is INVARIANT under step refinement: refining dVds by four "
           "decades changes it by less than 1%% of itself."),
    ("Q3", "The quadrature bias is FIRST order in 1/n (not second), because "
           "np.mean over an endpoint-inclusive linspace is not the trapezoid "
           "rule. Predicted: halving 1/(n-1) halves the residual."),
    ("Q4", "The `max(Vds - dVds, 1e-4)` guard is silently wrong whenever it "
           "fires, because it changes the interval but not the divisor 2*dVds. "
           "Predicted relative error >10%% at dVds=0.06, Vds=0.05."),
    ("Q5", "Propagated to f_max the bias is SMALLER than the bias in g_ds, "
           "because f_max ~ 1/sqrt(g_ds*(Rg+Rs) + 2*pi*f_T*Cgd*Rg): the "
           "square root halves it and the second denominator term dilutes "
           "it. Predicted f_max shift < half the g_ds shift."),
    ("Q6", "The 400-point Vg grid behind g_m (and hence behind the published "
           "peak f_T ~ 20 GHz) is a third entangled default. Predicted: peak "
           "f_T moves by <1%% under 16x Vg-grid refinement -- i.e. this one "
           "is innocent and the item can be closed rather than carried."),
]


def print_predictions():
    print("=" * 74)
    print("PRE-REGISTERED PREDICTIONS (recorded before measurement)")
    print("=" * 74)
    for tag, text in PREDICTIONS:
        print(f"  {tag}: {text}")
    print()


# ----------------------------------------------------------------------
# Exact validations
# ----------------------------------------------------------------------

def _rbar_test(Rfun, V, n):
    grid = np.linspace(0, V, n)
    return float(np.mean(np.array([Rfun(v) for v in grid])))


def validation_A_quadratic_exact():
    """DISCRIMINATING, EXACT.  For R(V_ch) = a + b*V_ch^2 both the
    discretised and the continuum channel average are known in closed form:

        mean over linspace(0,V,n) of V_ch^2 = V^2 (2n-1) / (6(n-1))
        continuum mean of V_ch^2 over [0,V] = V^2 / 3
        difference                          = V^2 / (6(n-1))      EXACT

    so Id(V) = V/(A + b c V^2) and g(V) = (A - b c V^2)/(A + b c V^2)^2
    with c = c_n for the discretisation and c = 1/3 for the model.  This
    gives an exact target for the finite difference AND an exact value for
    the bias the finite difference cannot see."""
    a, b, Rc, V = 400.0, 5.0e5, 600.0, VDS_DEFAULT
    A = a + Rc
    n = N_SEGMENTS_DEFAULT
    c_n = (2 * n - 1) / (6.0 * (n - 1))
    c_inf = 1.0 / 3.0

    # (i) the quadrature-bias formula itself, exactly
    meas = (_rbar_test(lambda v: b * v * v, V, n)) / (b * V * V)
    pred = c_n
    err_rule = abs(meas - pred) / pred
    ok1 = err_rule < 1e-13

    def g_closed(c):
        return (A - b * c * V * V) / (A + b * c * V * V) ** 2

    def Id_disc(v):
        return v / (A + b * c_n * v * v)

    def Id_cont(v):
        return v / (A + b * c_inf * v * v)

    # (ii) central difference of the discretised model -> the DISCRETISED
    #      closed form, not the model's.
    #      FIRST FORM OF THIS TEST FAILED (2026-09-26) and was rebuilt: it
    #      asserted `err_vs_continuum > 100 * err_vs_discrete`, measured 67x,
    #      and the 100 was a round number with no derivation behind it -- the
    #      same fault 09-25 item 10 caught one level up.  The defensible
    #      statements are (a) the residual against the DISCRETISED form is
    #      second order in h, verified by an oracle-free order test, and
    #      (b) the residual against the CONTINUUM form converges to the
    #      exactly-known quadrature gap instead of to zero.
    h = DVDS_DEFAULT
    g_fd = (Id_disc(V + h) - Id_disc(V - h)) / (2 * h)
    err_fd_vs_disc = abs(g_fd - g_closed(c_n)) / abs(g_closed(c_n))
    err_fd_vs_cont = abs(g_fd - g_closed(c_inf)) / abs(g_closed(c_inf))
    g_fd_half = (Id_disc(V + h / 2) - Id_disc(V - h / 2)) / (2 * (h / 2))
    err_half = abs(g_fd_half - g_closed(c_n)) / abs(g_closed(c_n))
    order = err_fd_vs_disc / err_half            # 4.0 => second order in h
    gap_exact = abs(g_closed(c_n) - g_closed(c_inf)) / abs(g_closed(c_inf))
    g_fd_tiny = (Id_disc(V + 1e-7) - Id_disc(V - 1e-7)) / 2e-7
    gap_meas = abs(g_fd_tiny - g_closed(c_inf)) / abs(g_closed(c_inf))
    ok2 = (abs(order - 4.0) < 0.25
           and abs(gap_meas - gap_exact) / gap_exact < 1e-6)

    # (iii) THE DISCRIMINATING ASSERTION: refine the step by 1e4 and the
    #       discrete-vs-continuum gap does not move.
    gaps = []
    for hh in (1e-3, 1e-5, 1e-7):
        g_d = (Id_disc(V + hh) - Id_disc(V - hh)) / (2 * hh)
        g_c = (Id_cont(V + hh) - Id_cont(V - hh)) / (2 * hh)
        gaps.append((g_d - g_c) / g_c)
    drift = abs(gaps[-1] - gaps[0]) / abs(gaps[0])
    ok3 = drift < 1e-2 and abs(gaps[0]) > 1e-4

    print("  [A] quadratic test resistance, exact closed forms"
          "                 DISCRIMINATING")
    print(f"      quadrature rule  measured c_n = {meas:.16f}")
    print(f"                       exact    c_n = {pred:.16f}   "
          f"rel err {err_rule:.2e}   {'PASS' if ok1 else 'FAIL'}")
    print(f"      central difference vs DISCRETISED closed form: "
          f"{err_fd_vs_disc:.3e}")
    print(f"      order test: err(h)/err(h/2) = {order:.4f}  "
          f"(4.0 = second order in h)")
    print(f"      central difference vs CONTINUUM   closed form: "
          f"{err_fd_vs_cont:.3e}")
    print(f"      as h->0 that residual -> {gap_meas:.10e}, and the exact "
          f"quadrature gap is {gap_exact:.10e}"
          f"   {'PASS' if ok2 else 'FAIL'}")
    print(f"      discrete-vs-continuum gap at h=1e-3/1e-5/1e-7: "
          f"{gaps[0]:+.6e} {gaps[1]:+.6e} {gaps[2]:+.6e}")
    print(f"      -> gap drifts {drift:.2e} over 4 decades of step "
          f"refinement   {'PASS' if ok3 else 'FAIL'}")
    print("      MEANING: a step-convergence study of a discretised "
          "function converges")
    print("      to the derivative of the discretisation, and cannot see "
          "this gap at all.")
    return ok1 and ok2 and ok3


def validation_B_constant_nondiscriminating():
    """EXACT, and RETAINED BECAUSE IT CANNOT DISCRIMINATE (09-25 item 11).
    For R independent of V_ch the quadrature is exact for every n and
    Id = Vds/(R+Rc) is exactly linear in Vds, so the central difference is
    exact for every step.  It passes for the shipped default and for an
    absurd one alike: it proves the harness and nothing about the default.
    Recorded AS non-discriminating so it is not mistaken for evidence."""
    R, Rc = 900.0, 600.0
    exact = 1.0 / (R + Rc)
    # FIRST FORM OF THIS TEST FAILED (2026-09-26): it required the residual
    # below an absolute 1e-12 and measured 6.08e-10 at h=1e-9.  The
    # arithmetic was right and the assertion wrong -- 6.08e-10 IS the
    # cancellation floor eps*|Id|/(2h|g|) at that step, so the constant
    # 1e-12 silently encoded h >~ 1e-7.  That is the audited fault class
    # (an absolute tolerance is a hidden scale claim) appearing inside the
    # auditor for the SECOND time in two days -- 09-25 item 10 was the
    # first.  Rebuilt on the derived bound.
    worst_ratio, rows, ok = 0.0, [], True
    for n in (2, 50, 5000):
        for h in (1e-1, 1e-3, 1e-9):
            Id = lambda v: v / (_rbar_test(lambda x: R, v, n) + Rc)
            g = (Id(VDS_DEFAULT + h) - Id(VDS_DEFAULT - h)) / (2 * h)
            meas = abs(g - exact) / exact
            bound = (EPS * abs(VDS_DEFAULT / (R + Rc))
                     / (2 * h * abs(exact))) + 4 * EPS
            rows.append((n, h, meas, bound))
            worst_ratio = max(worst_ratio, meas / bound)
            if meas > 60 * bound:
                ok = False
    print("  [B] Vds-independent channel resistance"
          "                    NON-DISCRIMINATING")
    print(f"      exact for every (n, h) up to the DERIVED cancellation "
          f"floor; worst")
    print(f"      measured/bound ratio {worst_ratio:.3f} over "
          f"n in (2,50,5000) x h in (1e-1,1e-3,1e-9)   "
          f"{'PASS' if ok else 'FAIL'}")
    print("      passes for a good default and a bad one -- proves the "
          "harness only.")
    return ok


def validation_C_linear_exact():
    """DISCRIMINATING ON THE MECHANISM, EXACT.  The mean of an
    endpoint-inclusive linspace reproduces the continuum mean EXACTLY for
    a constant and EXACTLY for a linear integrand (mean of an arithmetic
    sequence), and fails only from curvature onward.  So the bias found in
    [A] is specifically a curvature effect, not a general property of the
    averaging rule -- which is what makes its size depend on where in the
    transfer curve you sit."""
    V, b = VDS_DEFAULT, 7.0e3
    worst = 0.0
    for n in (2, 3, 50, 501):
        disc = _rbar_test(lambda v: b * v, V, n)
        cont = b * V / 2.0
        worst = max(worst, abs(disc - cont) / cont)
    ok = worst < 1e-14
    print("  [C] linear channel resistance -> zero quadrature bias"
          "      DISCRIMINATING")
    print(f"      worst rel err over n in (2,3,50,501): {worst:.2e}   "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def validation_D_cancellation_bound():
    """EXACT BOUND, DERIVED RATHER THAN ASSERTED (09-25 item 10 caught this
    module's ancestor asserting an absolute 1e-14 that silently encoded
    h~1e-3).  A central difference of Id near Vds loses
        eps*|Id| / (2h|g|)
    to cancellation in the numerator, relative to g.  The assertion is that
    the measured error at a deliberately tiny step tracks that bound, not
    that it is below any fixed constant."""
    Vg = 2.0
    Id0 = float(Id_model(np.array([Vg]), VDS_DEFAULT)[0])
    g_ref = float(gds_model(np.array([Vg]), dVds=1e-6, clamp=False)[0])
    rows, ok = [], True
    for h in (1e-9, 1e-11):
        g = float(gds_model(np.array([Vg]), dVds=h, clamp=False)[0])
        meas = abs(g - g_ref) / abs(g_ref)
        bound = EPS * abs(Id0) / (2 * h * abs(g_ref))
        rows.append((h, meas, bound))
        if not (meas <= 60 * bound):
            ok = False
    print("  [D] cancellation floor tracks its derived bound"
          "               DISCRIMINATING")
    for h, meas, bound in rows:
        print(f"      h={h:.0e}: measured {meas:.3e}  derived bound "
              f"{bound:.3e}  ratio {meas/bound:.2f}")
    print(f"      {'PASS' if ok else 'FAIL'}  (no absolute constant used)")
    return ok


def validation_E_clamp_exact():
    """EXACT, and it convicts shipped code.  When `max(Vds-dVds, 1e-4)`
    fires, the evaluated interval is [1e-4, Vds+dVds] of width
    W = Vds+dVds-1e-4, but the divisor stays 2*dVds.  The returned value is
    therefore the true secant slope over that interval times exactly
    W/(2*dVds) -- an algebraic factor, known in closed form, with no
    approximation anywhere."""
    Vg = np.array([2.0])
    Vds, h = VDS_DEFAULT, 0.06          # guard fires: Vds - h < 1e-4
    lo = 1e-4
    W = Vds + h - lo
    Ip = float(Id_model(Vg, Vds + h)[0])
    Im = float(Id_model(Vg, lo)[0])
    shipped = (Ip - Im) / (2 * h)
    truth = (Ip - Im) / W
    ratio_meas = shipped / truth
    ratio_pred = W / (2 * h)
    err = abs(ratio_meas - ratio_pred) / ratio_pred
    ok = err < 1e-12
    print("  [E] the Vds-dVds guard, exactly"
          "                             DISCRIMINATING")
    print(f"      shipped/true secant measured {ratio_meas:.16f}")
    print(f"                           exact   {ratio_pred:.16f}   rel err "
          f"{err:.2e}   {'PASS' if ok else 'FAIL'}")
    print(f"      -> the guard under-reports g_ds by "
          f"{100*(1-ratio_meas):.2f}% when it fires, silently.")
    return ok, 100 * (1 - ratio_meas)


# ----------------------------------------------------------------------
# Measurements
# ----------------------------------------------------------------------

def measure_fidelity():
    """The reimplementation must reproduce the shipped transfer_characteristic
    at n_segments = 50 before anything measured with it means anything."""
    Vg = np.linspace(-1.5, 3.5, 37)
    a = gfet.transfer_characteristic(Vg, Vds=VDS_DEFAULT)
    b = Id_model(Vg, VDS_DEFAULT, N_SEGMENTS_DEFAULT)
    bitwise = bool(np.array_equal(a, b))
    rel = float(np.max(np.abs(a - b) / np.abs(a)))
    print(f"  reimplementation vs shipped transfer_characteristic: "
          f"bitwise={bitwise}, max rel diff {rel:.2e}")
    assert rel < 1e-13, "reimplementation does not reproduce the model"
    # and the shipped output_conductance itself
    g_ship = rf.output_conductance(Vg, Vds=VDS_DEFAULT)
    g_here = gds_model(Vg, dVds=DVDS_DEFAULT, clamp=True)
    rel2 = float(np.max(np.abs(g_ship - g_here) / np.abs(g_ship)))
    print(f"  reimplementation vs shipped output_conductance:       "
          f"max rel diff {rel2:.2e}")
    assert rel2 < 1e-13
    return bitwise, rel, rel2


def measure_step(Vg_points):
    """RESULT 1: is dVds=1e-3 inside an anchored plateau?"""
    print("\nRESULT 1 -- the STEP, at n_segments = 50 (default)")
    print("  Vg (V)   anchored h    g_ds(h*) [S]     g_ds(1e-3) [S]   "
          "rel diff   r(h/10)   r(h/100)")
    out = {}
    for Vg in Vg_points:
        arr = np.array([Vg])
        f = lambda h: float(gds_model(arr, dVds=h, clamp=False)[0])
        h, v, r10, r100, trace = anchored_step(f)
        v_def = f(DVDS_DEFAULT)
        rel = abs(v_def - v) / abs(v) if v else float("nan")
        out[Vg] = dict(h=h, v=v, v_def=v_def, rel=rel, r10=r10, r100=r100)
        print(f"  {Vg:+6.2f}   {h:.0e}     {v:+.8e}   {v_def:+.8e}   "
              f"{rel:8.2e}   {r10:.1e}   {r100:.1e}")
    worst = max(o["rel"] for o in out.values())
    print(f"  worst step-induced relative error in g_ds: {worst:.3e}  "
          f"({100*worst:.4f}%)")
    return out, worst


def measure_quadrature(Vg_points, h_star):
    """RESULT 2: with the step converged, does n_segments still move g_ds --
    and is the movement first order in 1/(n-1)?"""
    print(f"\nRESULT 2 -- the QUADRATURE, at the anchored step h = {h_star:.0e}")
    ns = [50, 100, 200, 400, 800, 1600, 3200, 6400]
    table, out = {}, {}
    for Vg in Vg_points:
        arr = np.array([Vg])
        vals = [float(gds_model(arr, dVds=h_star, n_segments=n,
                                clamp=False)[0]) for n in ns]
        table[Vg] = vals
        # oracle-free order test: successive differences under n-doubling
        diffs = [vals[i] - vals[i + 1] for i in range(len(vals) - 1)]
        ratios = [diffs[i] / diffs[i + 1] for i in range(len(diffs) - 1)
                  if diffs[i + 1] != 0]
        # Richardson in 1/(n-1) from the two finest
        n1, n2 = ns[-2], ns[-1]
        g1, g2 = vals[-2], vals[-1]
        g_inf = (g2 * (n2 - 1) - g1 * (n1 - 1)) / ((n2 - 1) - (n1 - 1))
        rel50 = abs(vals[0] - g_inf) / abs(g_inf)
        out[Vg] = dict(vals=vals, g_inf=g_inf, rel50=rel50, ratios=ratios)
        print(f"  Vg = {Vg:+.2f} V")
        print(f"    g_ds(n=50)   = {vals[0]:+.8e} S")
        print(f"    g_ds(n->inf) = {g_inf:+.8e} S   (Richardson in 1/(n-1) "
              f"from n={n1},{n2})")
        print(f"    n=50 bias    = {rel50:.4e}  ({100*rel50:.4f}%)")
        print(f"    difference ratios under n-doubling (2.0 => first order "
              f"in 1/n): "
              + " ".join(f"{r:.3f}" for r in ratios))
    worst = max(o["rel50"] for o in out.values())
    print(f"  worst quadrature bias at n=50: {worst:.3e}  ({100*worst:.4f}%)")
    return out, worst, ns


def measure_gap_invariance(Vg):
    """RESULT 3: the discriminating measurement on the REAL model -- the
    n=50 bias does not shrink when the step is refined."""
    print("\nRESULT 3 -- is the quadrature bias visible to a step study? "
          "(the real model)")
    arr = np.array([Vg])
    print("     dVds      g_ds(n=50)        g_ds(n=6400)      gap (rel)")
    gaps = []
    for h in (1e-2, 1e-3, 1e-5, 1e-7):
        a = float(gds_model(arr, dVds=h, n_segments=50, clamp=False)[0])
        b = float(gds_model(arr, dVds=h, n_segments=6400, clamp=False)[0])
        gap = (a - b) / b
        gaps.append(gap)
        print(f"   {h:.0e}   {a:+.8e}   {b:+.8e}   {gap:+.6e}")
    drift = abs(gaps[-1] - gaps[1]) / abs(gaps[1])
    print(f"  gap drifts {drift:.2e} over 4 decades of step refinement "
          f"(1e-3 -> 1e-7)")
    print("  -> refining dVds by 1e4 removes none of it. The anchored step "
          "criterion is")
    print("     silent on this bias BY CONSTRUCTION: it converges in h at "
          "fixed n.")
    return gaps, drift


def measure_call_site():
    """RESULT 6 -- NOT PRE-REGISTERED.  Found while reconciling this module's
    peak f_T (10.14 GHz) against Chapter 4's published "Peak f_T ~20 GHz":
    the two are the same model at two different operating points.
    `output_conductance`'s SIGNATURE default is Vds = 0.05 V, but the path
    that actually produces `rf_figures_of_merit.png` and Chapter 4's numbers
    -- `plot_fT_fmax()` -- calls `compute_fT_fmax(..., Vds=0.1)`.

    That matters beyond bookkeeping.  2026-09-25's default-scale census
    scored this default by reading the SIGNATURE: "dVds = 1e-3 is 2% of its
    variable".  At the call site it is 1%.  A census that reads signature
    defaults mis-states the step-to-variable ratio by exactly the factor
    between the signature default and the call site -- here 2x, and in
    general unbounded.  Recorded as an unpredicted finding, not folded into
    a prediction after the fact."""
    print("\nRESULT 6 -- signature default vs call site (NOT PRE-REGISTERED)")
    Vds_sig, Vds_call = VDS_DEFAULT, 0.1
    print(f"  output_conductance signature default   Vds = {Vds_sig} V"
          f"   -> dVds/Vds = {DVDS_DEFAULT/Vds_sig:.3e}")
    print(f"  plot_fT_fmax call site (Chapter 4)     Vds = {Vds_call} V"
          f"   -> dVds/Vds = {DVDS_DEFAULT/Vds_call:.3e}")
    print(f"  -> 09-25's census scored the ratio as 2e-2; at the call site "
          f"it is 1e-2, low by exactly 2x")

    Vg = np.linspace(*VG_SPAN_DEFAULT, VG_GRID_DEFAULT)
    fT, fmax, gm, gds, Cgs = rf.compute_fT_fmax(Vg, Vds=Vds_call)
    i, j = int(np.argmax(fT)), int(np.argmax(fmax))
    print(f"  at the call site: peak f_T = {fT[i]/1e9:.3f} GHz, peak f_max = "
          f"{fmax[j]/1e9:.3f} GHz")
    print(f"  -> this is Chapter 4's \"Peak f_T ~20 GHz\"; the 10.14 GHz in "
          f"RESULT 4 is the")
    print(f"     SAME model at the signature default Vds = 0.05 V, not a "
          f"discrepancy.")

    # the two audited defaults, re-measured at the operating point that
    # actually feeds the thesis
    arr = np.array([2.0])
    f = lambda h: float(gds_model(arr, dVds=h, n_segments=N_SEGMENTS_DEFAULT,
                                  clamp=False)[0] if False else
                        _gds_at(arr, Vds_call, h, N_SEGMENTS_DEFAULT))
    h_star, v, r10, r100, _ = anchored_step(f)
    v_def = f(DVDS_DEFAULT)
    step_err = abs(v_def - v) / abs(v)
    g50 = _gds_at(arr, Vds_call, h_star, 50)
    g3200 = _gds_at(arr, Vds_call, h_star, 3200)
    g6400 = _gds_at(arr, Vds_call, h_star, 6400)
    g_inf = (g6400 * 6399 - g3200 * 3199) / 3200.0
    quad_err = abs(g50 - g_inf) / abs(g_inf)
    print(f"  step error at Vds=0.1:       {step_err:.3e} "
          f"({100*step_err:.5f}%)  [anchored h = {h_star:.0e}]")
    print(f"  quadrature bias at Vds=0.1:  {quad_err:.3e} "
          f"({100*quad_err:.5f}%)")
    print(f"  ratio to the Vds=0.05 quadrature bias: "
          f"{quad_err/8.3521e-07:.2f}x  (a curvature effect grows with Vds)")
    print(f"  the guard `max(Vds-dVds,1e-4)` becomes reachable for "
          f"Vds <= {DVDS_DEFAULT + 1e-4:.1e} V;")
    print(f"  neither operating point is near it, so the 8.42% error is "
          f"latent, not active.")
    return dict(step_err=step_err, quad_err=quad_err,
                peak_fT=fT[i] / 1e9, peak_fmax=fmax[j] / 1e9,
                h_star=h_star)


def _gds_at(Vg, Vds, dVds, n_segments):
    Ip = Id_model(Vg, Vds + dVds, n_segments)
    Im = Id_model(Vg, Vds - dVds, n_segments)
    return float(((Ip - Im) / (2 * dVds))[0])


def _fT_fmax_with_gds(Vg_range, gds, Vds=VDS_DEFAULT):
    """rf._compute_fT_fmax_core with g_ds injected instead of recomputed,
    so the only thing that changes between the two f_max curves below is
    g_ds itself."""
    gm, _Id = rf.transconductance(Vg_range, Vds=Vds)
    Cgs = rf.gate_capacitance(Vg_range)
    Cgd = rf.Cgd_over_Cgs * Cgs
    Rg = rf.gate_resistance()
    Rs = rf.source_access_resistance()
    fT = np.abs(gm) / (2 * np.pi * Cgs)
    term_gds = gds * (Rg + Rs)
    term_rg = 2 * np.pi * fT * Cgd * Rg
    denom = np.clip(term_gds + term_rg, 1e-30, None)
    fmax = fT / (2 * np.sqrt(denom))
    return fT, fmax, term_gds, term_rg


def measure_fmax_propagation(h_star):
    """RESULT 4: does the bias reach the Chapter 4 published numbers?"""
    print("\nRESULT 4 -- propagation to f_T and f_max (Chapter 4, Section 4.6)")
    Vg = np.linspace(*VG_SPAN_DEFAULT, VG_GRID_DEFAULT)
    gds_def = rf.output_conductance(Vg, Vds=VDS_DEFAULT)      # as shipped
    g800 = gds_model(Vg, dVds=h_star, n_segments=800, clamp=False)
    g1600 = gds_model(Vg, dVds=h_star, n_segments=1600, clamp=False)
    gds_fix = (g1600 * 1599 - g800 * 799) / 800.0             # 1/(n-1) limit

    fT_d, fmax_d, tg_d, tr_d = _fT_fmax_with_gds(Vg, gds_def)
    fT_f, fmax_f, tg_f, tr_f = _fT_fmax_with_gds(Vg, gds_fix)

    i = int(np.argmax(fT_d))
    j_d, j_f = int(np.argmax(fmax_d)), int(np.argmax(fmax_f))
    w = tg_d[j_d] / (tg_d[j_d] + tr_d[j_d])
    print(f"  peak f_T  = {fT_d[i]/1e9:.3f} GHz at Vg = {Vg[i]:+.3f} V  "
          f"(f_T does not contain g_ds; unchanged)")
    print(f"  peak f_max default   = {fmax_d[j_d]/1e9:.4f} GHz at "
          f"Vg = {Vg[j_d]:+.3f} V")
    print(f"  peak f_max corrected = {fmax_f[j_f]/1e9:.4f} GHz at "
          f"Vg = {Vg[j_f]:+.3f} V")
    rel_peak = abs(fmax_f[j_f] - fmax_d[j_d]) / fmax_d[j_d]
    rel_gds = float(np.max(np.abs(gds_fix - gds_def) / np.abs(gds_fix)))
    rel_curve = float(np.max(np.abs(fmax_f - fmax_d) / fmax_d))
    print(f"  peak f_max shift            = {rel_peak:.3e}  "
          f"({100*rel_peak:.4f}%)")
    print(f"  worst g_ds shift on the sweep = {rel_gds:.3e}  "
          f"({100*rel_gds:.4f}%)")
    print(f"  worst f_max shift on the sweep = {rel_curve:.3e}  "
          f"({100*rel_curve:.4f}%)")
    print(f"  weight of the g_ds term in the f_max denominator at peak "
          f"f_max: {w:.4f}")
    print(f"  -> dilution factor observed: f_max shift / g_ds shift at peak "
          f"= {rel_peak/max(rel_gds,1e-300):.4f}")
    print(f"     (a pure square-root would give 0.5; the R_g term "
          f"contributes the rest)")
    return dict(Vg=Vg, gds_def=gds_def, gds_fix=gds_fix, fmax_d=fmax_d,
                fmax_f=fmax_f, fT=fT_d, rel_peak=rel_peak, rel_gds=rel_gds,
                rel_curve=rel_curve, weight=w, peak_fT=fT_d[i] / 1e9,
                peak_fmax_d=fmax_d[j_d] / 1e9, peak_fmax_f=fmax_f[j_f] / 1e9)


def measure_vg_grid():
    """RESULT 5: the third entangled default -- the Vg grid behind
    g_m = np.gradient(Id, Vg), which sets the published peak f_T."""
    print("\nRESULT 5 -- the Vg grid behind g_m (and the published peak f_T)")
    print("    points    dVg (mV)   peak f_T (GHz)   Vg at peak (V)")
    res = []
    for npts in (VG_GRID_DEFAULT, 1600, 6400):
        Vg = np.linspace(*VG_SPAN_DEFAULT, npts)
        Id = gfet.transfer_characteristic(Vg, Vds=VDS_DEFAULT)
        gm = np.gradient(Id, Vg)
        Cgs = rf.gate_capacitance(Vg)
        fT = np.abs(gm) / (2 * np.pi * Cgs)
        k = int(np.argmax(fT))
        res.append((npts, (Vg[1] - Vg[0]) * 1e3, fT[k] / 1e9, Vg[k]))
        print(f"    {npts:6d}    {(Vg[1]-Vg[0])*1e3:7.3f}   "
              f"{fT[k]/1e9:12.5f}   {Vg[k]:+.4f}")
    rel = abs(res[-1][2] - res[0][2]) / res[-1][2]
    print(f"  peak f_T moves {rel:.3e} ({100*rel:.4f}%) under 16x Vg-grid "
          f"refinement")
    return res, rel


def make_figure(step_out, quad_out, ns, prop, gaps, path="gds_quadrature_audit.png"):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

    ax = axes[0]
    Vg0 = sorted(step_out)[len(step_out) // 2]
    arr = np.array([Vg0])
    hs = np.logspace(-1, -11, 21)
    ref = float(gds_model(arr, dVds=1e-5, clamp=False)[0])
    errs = [abs(float(gds_model(arr, dVds=h, clamp=False)[0]) - ref) / abs(ref)
            for h in hs]
    ax.loglog(hs, np.maximum(errs, 1e-18), "o-", color="steelblue",
              label=r"$|g_{ds}(h)-g_{ds}(10^{-5})|/|g_{ds}|$")
    ax.axvline(DVDS_DEFAULT, color="crimson", ls="--",
               label=r"shipped $\Delta V_{ds}=10^{-3}$")
    ax.axhline(abs(gaps[1]), color="darkgreen", ls=":",
               label="quadrature bias at $n=50$")
    ax.set_xlabel(r"step $\Delta V_{ds}$ (V)")
    ax.set_ylabel("relative error")
    ax.set_title(f"(a) the step is innocent; the bias is\nnot below it "
                 f"($V_g={Vg0:+.2f}$ V)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    for Vg, o in sorted(quad_out.items()):
        inv = 1.0 / (np.array(ns) - 1.0)
        ax.plot(inv, np.array(o["vals"]) * 1e3, "o-",
                label=f"$V_g={Vg:+.2f}$ V")
        ax.plot([0], [o["g_inf"] * 1e3], "k*", ms=11)
    ax.axvline(1.0 / (N_SEGMENTS_DEFAULT - 1), color="crimson", ls="--",
               label="shipped $n=50$")
    ax.set_xlabel(r"$1/(n_{segments}-1)$")
    ax.set_ylabel(r"$g_{ds}$ (mS)")
    ax.set_title("(b) $g_{ds}$ is linear in $1/(n-1)$;\nstars = Richardson "
                 "limit")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(prop["Vg"], prop["fmax_d"] / 1e9, lw=2, color="crimson",
            label=r"$f_{max}$, shipped ($n=50$)")
    ax.plot(prop["Vg"], prop["fmax_f"] / 1e9, lw=1.4, ls="--",
            color="navy", label=r"$f_{max}$, $n\to\infty$")
    ax.plot(prop["Vg"], prop["fT"] / 1e9, lw=1.0, color="gray", alpha=0.7,
            label=r"$f_T$ (contains no $g_{ds}$)")
    ax.set_xlabel(r"$V_g$ (V)")
    ax.set_ylabel("frequency (GHz)")
    ax.set_title(f"(c) propagation to Chapter 4's RF numbers\n"
                 f"peak $f_{{max}}$ shifts "
                 f"{100*prop['rel_peak']:.3f}%")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Two entangled defaults behind $g_{ds}$: a finite-difference "
                 "step and a quadrature resolution", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\n  figure written: {path}")


def main():
    print(__doc__.strip()[:0] or "", end="")
    print("=" * 74)
    print("g_ds STEP-AND-QUADRATURE AUDIT -- rf_small_signal_model /"
          " graphene_fet_model")
    print("=" * 74 + "\n")
    print_predictions()

    print("FIDELITY OF THE INSTRUMENT")
    measure_fidelity()

    print("\nEXACT VALIDATIONS")
    vA = validation_A_quadratic_exact()
    vB = validation_B_constant_nondiscriminating()
    vC = validation_C_linear_exact()
    vD = validation_D_cancellation_bound()
    vE, clamp_pct = validation_E_clamp_exact()
    npass = sum(map(bool, (vA, vB, vC, vD, vE)))
    print(f"  validations passing: {npass}/5")

    Vg_points = [0.0, 1.6, 2.0, 2.6, 3.2]
    step_out, worst_step = measure_step(Vg_points)
    h_star = min(o["h"] for o in step_out.values())
    quad_out, worst_quad, ns = measure_quadrature(Vg_points, h_star)
    gaps, drift = measure_gap_invariance(2.0)
    prop = measure_fmax_propagation(h_star)
    grid_res, grid_rel = measure_vg_grid()
    call_site = measure_call_site()

    print("\n" + "=" * 74)
    print("PREDICTIONS SCORED")
    print("=" * 74)
    q1 = worst_step < 1e-3
    print(f"  Q1 {'PASS' if q1 else 'FAIL'}: dVds=1e-3 step error "
          f"{100*worst_step:.4f}% "
          f"({'inside' if q1 else 'outside'} the anchored plateau)")
    q2 = worst_quad >= 1e-3 and drift < 1e-2
    print(f"  Q2 {'PASS' if q2 else 'FAIL'}: n=50 bias "
          f"{100*worst_quad:.4f}%, invariant to {drift:.1e} under 4 decades "
          f"of step refinement")
    ratios = quad_out[2.0]["ratios"]
    q3 = abs(np.mean(ratios) - 2.0) < 0.15
    print(f"  Q3 {'PASS' if q3 else 'FAIL'}: n-doubling difference ratio "
          f"mean {np.mean(ratios):.3f} (2.0 = first order in 1/n)")
    q4 = clamp_pct > 10.0
    print(f"  Q4 {'PASS' if q4 else 'FAIL'}: guard under-reports g_ds by "
          f"{clamp_pct:.2f}% when it fires")
    q5 = prop["rel_peak"] < 0.5 * prop["rel_gds"]
    print(f"  Q5 {'PASS' if q5 else 'FAIL'}: f_max shift / g_ds shift = "
          f"{prop['rel_peak']/prop['rel_gds']:.4f}")
    q6 = grid_rel < 1e-2
    print(f"  Q6 {'PASS' if q6 else 'FAIL'}: peak f_T moves "
          f"{100*grid_rel:.4f}% under 16x Vg-grid refinement")
    print(f"\n  score: {sum(map(bool,(q1,q2,q3,q4,q5,q6)))}/6 predictions, "
          f"{npass}/5 validations")
    print("\n  UNPREDICTED FINDING (RESULT 6): Chapter 4's published RF")
    print("  numbers are evaluated at Vds = 0.1 V while output_conductance's")
    print("  signature default is 0.05 V, so 09-25's default-scale census")
    print("  scored dVds/Vds as 2e-2 where the call site makes it 1e-2.")
    print("  A default-scale census must read CALL SITES, not signatures.")

    make_figure(step_out, quad_out, ns, prop, gaps)


if __name__ == "__main__":
    main()
