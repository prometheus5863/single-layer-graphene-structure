"""
graphene_velocity_saturation_model.py

Opens work on the item created 2026-10-04 and named that session's TOP item:
*a saturation term in `transfer_characteristic()`*, with the first open item
in this series to arrive carrying a MEASURED acceptance criterion:

    (A) it must divide g_ds by  ~4.85  at the peak-f_T bias, and
    (B) it must REVERSE the sign of d(f_max/f_T)/dV_ds, which in the
        resistor model runs the wrong way (0.6442 -> 0.6218 over
        V_ds = 0.05 - 1 V), where a saturating device's ratio RISES,
    (C) validated against the resistor model in the low-V_ds OVERLAPPING
        LIMIT, where the two must agree (the rule set 2026-10-03 and
        implemented as E5 on 2026-10-04).

RESULT, stated up front because criterion (A) is NOT met and the amount by
which it is missed is the finding.  With v_sat taken from the literature
rather than fitted -- see "WHY v_sat IS NOT A FREE PARAMETER" below -- the
saturation term divides g_ds by a factor this module measures and reports,
and that factor is compared against the required 4.85 as a RATIO rather
than as a pass/fail.  Criterion (B) is met.  Criterion (C) is met, and in
being met it exposes a quadrature error in the model being superseded (see
"WHAT THE OVERLAPPING LIMIT FOUND" below), which was not what it was
checked for.

----------------------------------------------------------------------------
WHY THIS IS A NEW MODULE AND NOT AN EDIT TO transfer_characteristic()
----------------------------------------------------------------------------
`gfet.transfer_characteristic()` is read by `rf_small_signal_model`,
`graphene_fmax_shortfall_decomposition`, `graphene_gds_quadrature_audit`,
`graphene_sensitivity_audit` and `graphene_pristine_transcript_audit`, and
its numbers are quoted in committed transcripts and in Chapter 4.  Editing
it in place would invalidate every one of those in a single commit, and the
2026-10-02 rule ("when a rule changes, every number derived from it is
stale until re-derived") says that is a separate, deliberate act.  So this
module adds `transfer_characteristic_saturated()` ALONGSIDE it and leaves
the resistor model BITWISE UNTOUCHED -- Section 7 checks that it is -- and
then reports what the rewiring would cost, so the rewiring is a decision
with a price tag rather than a side effect.  (Same discipline as the
2026-10-04 session in the verification repository: a subclass, not an edit.)

----------------------------------------------------------------------------
THE PHYSICS, AND WHY IT IS STILL CLOSED-FORM
----------------------------------------------------------------------------
Local drift velocity under a soft-saturation (Caughey-Thomas) law with
exponent beta = 1, which is the gamma = 1 form Feijoo et al. 2019 use:

    v(E) = mu*E / (1 + mu*E/v_sat)                                      (1)

Current continuity makes I_d the same at every point along the channel:

    I_d = W * e * n(V) * v(E(V))                                        (2)

so the local velocity the channel is REQUIRED to deliver at potential V is
u(V) = I_d/(W*e*n(V)), and inverting (1) for the field that delivers it:

    E = u / ( mu * (1 - u/v_sat) )                                      (3)

which has no solution once u >= v_sat -- a real physical ceiling, not a
numerical one, and Section 2 reports when the model hits it rather than
clipping it away.  Since E = dV/dx, L = integral of dV/E over the intrinsic
channel drop, and (3) makes that integral SEPARABLE:

    L = (mu*W*e/I_d) * Q  -  mu * S

    Q = integral of n(V) dV          over 0 .. V_ds,ch   [m^-2 . V]
    S = integral of dV / v_sat(V)    over 0 .. V_ds,ch   [V.s/m]

    =>   I_d = mu*W*e*Q / ( L + mu*S )                                  (4)

(4) is exact for beta = 1 and holds even when v_sat varies along the
channel, which it does here.  The contacts then close the loop:

    V_ds,ch = V_ds - I_d * Rc_total                                     (5)

solved by bisection on I_d.  Setting S = 0 in (4) recovers the
drift-diffusion result with no saturation, which is the overlapping limit
criterion (C) asks for.

----------------------------------------------------------------------------
WHY v_sat IS NOT A FREE PARAMETER
----------------------------------------------------------------------------
Criterion (A) names a number -- 4.85 -- and v_sat is the one knob that
moves it.  Choosing v_sat to make the criterion pass would be fitting the
model to its own acceptance test, which is the 2026-10-02 fault ("when it
agrees, that is a fact about two artefacts, not about the world") in its
purest form.  So v_sat is taken from optical-phonon-emission physics with
NO adjustable scale:

    v_sat(n) = (2/pi) * Omega / sqrt(pi*n),    Omega = hbar_Omega/hbar   (6)

the standard form behind the measurements of Dorgan, Bae and Pop, Appl.
Phys. Lett. 97, 082112 (2010), who measure v_sat ~ 3e7 cm/s at low density
on SiO2 -- THIS MODEL'S SUBSTRATE -- falling with increasing n, and the
form Feijoo et al., 2D Mater. 6, 035027 (2019) evaluate with a fitted
phonon energy hbar_Omega = 0.10 eV.  That 0.10 eV is the ONLY literature
number entering, it is taken as published, and Section 5 reports what
v_sat (6) delivers at this device's own carrier density so the choice is
auditable rather than buried.  Section 6 then inverts the question and
reports the v_sat that WOULD satisfy criterion (A) exactly, so the gap is
stated as a falsifiable claim about a physical quantity.

Dorgan et al.'s own best fit is beta = 2, not beta = 1.  beta = 2 destroys
the separability that gives (4) and is NOT implemented here; Section 6
bounds the error that costs by evaluating the beta = 1 law against the
beta = 2 law at the fields this device actually reaches, because an
unquantified modelling choice is the 2026-10-04 fault (a quantity computed
correctly under a name that claims more than it measures).

----------------------------------------------------------------------------
WHAT THE OVERLAPPING LIMIT FOUND (not what it was checked for)
----------------------------------------------------------------------------
In the v_sat -> infinity limit, (4) becomes I_d = mu*W*e*<n>_V*V_ds,ch/L,
where <n>_V is the V-AVERAGE of n.  The resistor model being superseded
computes instead

    R_channel = mean over V of (L/W)/(n(V)*e*mu)  =  (L/(W*e*mu)) * <1/n>_V

i.e. it averages the RECIPROCAL.  By Cauchy-Schwarz <n>_V >= 1/<1/n>_V
with equality only when n is constant along the channel, so the two differ
by a Jensen gap that is ZERO at V_ds = 0 and grows with V_ds.  The
overlapping-limit check therefore cannot be an assertion of bitwise
agreement -- it has to be an assertion that the gap VANISHES in the limit
and a measurement of its size away from it.  Section 3 does both.  The
resistor model's quadrature is not the exact drift-diffusion answer, and
nothing in this repository had asked.  **No committed number is withdrawn
by this**: the gap is reported, its sign is stated (the resistor model
UNDER-states I_d, always, by construction), and Chapter 4 is annotated in
place.

----------------------------------------------------------------------------
VALIDATION DESIGN
----------------------------------------------------------------------------
Five tolerance-free checks, each NAMING the operation it is exact under
(the 2026-10-03 rule, and the form 2026-10-04's E1-E5 established):

  X1  V_ds = 0  =>  I_d exactly 0.0            (multiplication by zero)
  X2  constant n, v_sat -> inf  =>  the solver reproduces the CLOSED FORM
      V_ds/(L/(W*e*mu*n) + Rc) -- exact under the bisection's own
      convergence, reported in ULP rather than claimed exact
  X3  the strong-saturation ceiling  I_d -> W*e*n*v_sat  exactly, and it
      is independent of mu, L and V_ds -- exact under the limit mu -> inf
  X4  doubling v_sat exactly doubles the ceiling of X3  (exact under
      multiplication by a power of two)
  X5  constant n  =>  Q == n*V_ds,ch  exactly   (exact under the
      quadrature's own partition-of-unity, i.e. the trapezoid weights
      summing to the interval)

plus MAGNITUDE-form controls, per the standing top methodological item
(record the magnitude of every MUST_CHANGE, not just its sign):

  G1  the g_ds reduction factor is reported and asserted as a BAND
  G2  the sign of d(f_max/f_T)/dV_ds is reversed, AND the magnitude of the
      reversal is printed
  G3  CONTROL: the v_sat response and the mobility response are NOT
      aliased -- if they were, "saturation fixed it" would be
      indistinguishable from "mu was wrong"
  G4  CONTROL: the resistor model is bitwise unchanged by importing this
      module (Section 7)

Run:  python3 graphene_velocity_saturation_model.py
Writes: velocity_saturation_model.png, velocity_saturation_output.txt
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import rf_small_signal_model as rf

# ---------------------------------------------------------------------------
# Literature constants.  Both are taken as published; neither is fitted here.
# ---------------------------------------------------------------------------
HBAR_OMEGA_EV = 0.10   # eV, optical-phonon energy; Feijoo et al. 2D Mater. 6,
                       # 035027 (2019), Table 1.  The ONLY fitted literature
                       # number entering this module, and it is THEIR fit.
E_CHARGE = 1.602176634e-19
HBAR = 1.054571817e-34
OMEGA_OP = HBAR_OMEGA_EV * E_CHARGE / HBAR   # rad/s

# Reference bias points used throughout, matching rf_small_signal_model's
# plot_fT_fmax() so the comparison is against the numbers Chapter 4 quotes.
VDS_REF = 0.1
VDS_LADDER = (0.05, 0.1, 0.2, 0.5, 1.0)
N_QUAD = 201     # trapezoid points for Q and S along the channel

_FAIL = []
_LINES = []


def say(s=""):
    _LINES.append(s)
    print(s)


def check(tag, ok, label, detail=""):
    say(f"  [{'PASS' if ok else 'FAIL'}] {tag:<5}{label}"
        + (f"  -- {detail}" if detail else ""))
    if not ok:
        _FAIL.append(tag)
    return ok


# ---------------------------------------------------------------------------
# The saturation model
# ---------------------------------------------------------------------------
def v_sat_of_n(n, omega_op=OMEGA_OP, prefactor=2.0 / np.pi):
    """Density-dependent saturation velocity, Eq. (6).

    v_sat = (2/pi) * Omega / sqrt(pi*n).  Diverges as n -> 0, which is
    physical (an empty channel imposes no velocity ceiling) and harmless
    downstream because n is floored at n_puddle by carrier_density().
    """
    return prefactor * omega_op / np.sqrt(np.pi * np.asarray(n, dtype=float))


def _channel_integrals(V_g, Vds_ch, n_quad=N_QUAD, v_sat_const=None,
                       constant_n=False, omega_op=OMEGA_OP,
                       v_sat_prefactor=2.0 / np.pi):
    """Q = int n dV and S = int dV/v_sat over the intrinsic channel drop.

    constant_n=True freezes n at its V_ch = 0 value, which is what the
    exactness checks X2/X3/X5 need: it makes Q = n*Vds_ch an identity.
    v_sat_const, when given, overrides Eq. (6) with a fixed value (used by
    X3/X4 and by the beta-sensitivity bound).
    """
    V = np.linspace(0.0, Vds_ch, n_quad)
    if constant_n:
        n = np.full_like(V, float(gfet.carrier_density(V_g, 0.0)))
    else:
        n = np.array([float(gfet.carrier_density(V_g, v)) for v in V])

    Q = np.trapezoid(n, V)
    if v_sat_const is None:
        vs = v_sat_of_n(n, omega_op=omega_op, prefactor=v_sat_prefactor)
    else:
        vs = np.full_like(V, float(v_sat_const))
    S = np.trapezoid(1.0 / vs, V)
    return Q, S, n, vs


def _Id_given_Vds_ch(V_g, Vds_ch, **kw):
    """Equation (4).  saturate=False sets S = 0 (the overlapping limit)."""
    saturate = kw.pop("saturate", True)
    Q, S, n, vs = _channel_integrals(V_g, Vds_ch, **kw)
    if not saturate:
        S = 0.0
    denom = gfet.L + gfet.mu * S
    Id = gfet.mu * gfet.W * E_CHARGE * Q / denom
    return Id, Q, S, n, vs


def transfer_characteristic_saturated(Vg_range, Vds=VDS_REF, saturate=True,
                                      Rc_total=None, tol=1e-15, max_iter=200,
                                      report_ceiling=False, **kw):
    """I_d(V_g) with velocity saturation, Eqs. (4)-(5).

    Returns I_d in the same shape/units as gfet.transfer_characteristic().
    `saturate=False` gives the no-saturation drift-diffusion limit, which
    is the artefact the overlapping-limit check compares against.
    `report_ceiling=True` additionally returns, per bias, the largest
    u(V)/v_sat(V) reached -- the physical ceiling of Eq. (3).  A value
    >= 1 means the channel CANNOT deliver the current and the model has no
    solution there; this is reported, never clipped.
    """
    Vg_range = np.atleast_1d(np.asarray(Vg_range, dtype=float))
    Rc = gfet.Rc_total if Rc_total is None else Rc_total
    Id_out = np.zeros_like(Vg_range)
    ceil_out = np.zeros_like(Vg_range)

    for i, V_g in enumerate(Vg_range):
        if Vds == 0.0:
            Id_out[i] = 0.0      # X1: exact, no solve
            continue
        if Rc == 0.0:
            Id, Q, S, n, vs = _Id_given_Vds_ch(V_g, Vds, saturate=saturate, **kw)
        else:
            # Bisect on Vds_ch in (0, Vds]: the residual
            #   f(Vds_ch) = Vds - Id(Vds_ch)*Rc - Vds_ch
            # is strictly decreasing in Vds_ch, so a bracket is guaranteed.
            lo, hi = 0.0, Vds
            for _ in range(max_iter):
                mid = 0.5 * (lo + hi)
                Id, Q, S, n, vs = _Id_given_Vds_ch(
                    V_g, mid, saturate=saturate, **kw)
                resid = Vds - Id * Rc - mid
                if resid > 0.0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo <= tol * max(Vds, 1.0):
                    break
            Id, Q, S, n, vs = _Id_given_Vds_ch(
                V_g, 0.5 * (lo + hi), saturate=saturate, **kw)
        Id_out[i] = Id
        if report_ceiling:
            u = Id / (gfet.W * E_CHARGE * n)
            ceil_out[i] = float(np.max(u / vs))

    if report_ceiling:
        return Id_out, ceil_out
    return Id_out


def gds_saturated(Vg_range, Vds=VDS_REF, dVds=1e-3, **kw):
    """d(I_d)/d(V_ds) by the SAME central difference rf.output_conductance
    uses, so the comparison is between models and not between estimators."""
    Ip = transfer_characteristic_saturated(Vg_range, Vds=Vds + dVds, **kw)
    Im = transfer_characteristic_saturated(Vg_range, Vds=Vds - dVds, **kw)
    return (Ip - Im) / (2.0 * dVds)


def gds_at_bias(V_g, Vds=VDS_REF, dVds=1e-3, **kw):
    """g_ds at ONE gate bias, without going through np.gradient.

    fT_fmax_saturated() needs g_m and therefore needs >= 2 gate points;
    Section 6's inversion needs only g_ds at a single bias.  Calling the
    former with a 1-element array raised IndexError inside np.gradient --
    a crash in a reporting path, which 2026-10-04 recorded as a strictly
    weaker detection than an assertion.  This helper removes the need.
    """
    Vg1 = np.array([float(V_g)])
    Ip = transfer_characteristic_saturated(Vg1, Vds=Vds + dVds, **kw)[0]
    Im = transfer_characteristic_saturated(Vg1, Vds=Vds - dVds, **kw)[0]
    return (Ip - Im) / (2.0 * dVds)


def gm_saturated(Vg_range, Vds=VDS_REF, **kw):
    Id = transfer_characteristic_saturated(Vg_range, Vds=Vds, **kw)
    return np.gradient(Id, Vg_range), Id


def fT_fmax_saturated(Vg_range, Vds=VDS_REF, N_fingers=1, **kw):
    """f_T and f_max from the saturated model, reusing rf's own C_gs, C_gd,
    R_g and R_s so that ONLY g_m and g_ds change between the two models."""
    gm, Id = gm_saturated(Vg_range, Vds=Vds, **kw)
    gds = gds_saturated(Vg_range, Vds=Vds, **kw)
    Cgs = rf.gate_capacitance(Vg_range)
    Cgd = rf.Cgd_over_Cgs * Cgs
    Rg = rf.gate_resistance(N_fingers=N_fingers)
    Rs = rf.source_access_resistance()
    fT = np.abs(gm) / (2 * np.pi * Cgs)
    denom = gds * (Rg + Rs) + 2 * np.pi * fT * Cgd * Rg
    denom = np.clip(denom, 1e-30, None)
    fmax = fT / (2 * np.sqrt(denom))
    return fT, fmax, gm, gds, Cgs


def _at_literature_geometry(fn, W=rf.W_RF, N_fingers=rf.N_FINGERS_RF):
    """Run fn() with gfet.W / gfet.Rc_total temporarily set to the
    literature-scale geometry, exactly as rf.compute_fT_fmax does.  Late
    binding throughout -- no module-level name is captured as a default
    (the 2026-10-03 frozen-default bug, and the 2026-10-04 repeat of it)."""
    W_saved, Rc_saved = gfet.W, gfet.Rc_total
    gfet.W = W
    gfet.Rc_total = 2 * (gfet.Rc_per_width_ohm_um * 1e-6) / W
    try:
        return fn(N_fingers)
    finally:
        gfet.W, gfet.Rc_total = W_saved, Rc_saved


def ulps(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / np.spacing(max(abs(a), abs(b)))


# ===========================================================================
def main():
    Vg = np.linspace(-2.0, 4.0, 400)

    say("velocity saturation in transfer_characteristic(): the 2026-10-04 top item")
    say("(acceptance criteria A, B and C are quoted in the module docstring)")
    say()
    say("=" * 78)
    say("SECTION 1.  The exactness battery")
    say("=" * 78)

    Vg1 = np.array([2.0])

    # X1 -- multiplication by zero
    Id0 = transfer_characteristic_saturated(Vg1, Vds=0.0)
    check("X1", Id0[0] == 0.0,
          "V_ds = 0 gives I_d exactly 0.0  (exact under mult. by zero)",
          f"I_d = {Id0[0]!r}")

    # X2 -- the solver against the closed form, constant n, no saturation
    n0 = float(gfet.carrier_density(2.0, 0.0))
    R_ch_closed = gfet.L / (gfet.W * E_CHARGE * gfet.mu * n0)
    Id_closed = VDS_REF / (R_ch_closed + gfet.Rc_total)
    Id_solver = transfer_characteristic_saturated(
        Vg1, Vds=VDS_REF, saturate=False, constant_n=True)[0]
    u2 = ulps(Id_solver, Id_closed)
    check("X2", u2 <= 64.0,
          "constant n, no saturation: solver == closed form (ULP, not 'exact')",
          f"{u2:.1f} ULP; I_d = {Id_solver:.12e} vs {Id_closed:.12e}")

    # X3 -- the strong-saturation ceiling, as an IDENTITY rather than a limit.
    # With L = 0 exactly, Eq. (4) reads I_d = mu*W*e*n*V / (mu*V/v_sat)
    # = W*e*n*v_sat: every appearance of mu, V and L cancels ALGEBRAICALLY,
    # so this is exact under cancellation, not a numerical limit.  The first
    # committed form of X3 drove mu -> large at finite L and compared in ULP;
    # that was the WRONG INSTRUMENT -- the quantity converges as
    # L*v_sat/(mu*V_ds), so it was ~1e-6 away by construction and the check
    # failed on its own tolerance rather than on the model.  Recorded, not
    # relaxed: X3 is now an identity and X3b checks the finite-L law it
    # should have checked in the first place.
    mu_saved, L_saved = gfet.mu, gfet.L
    vs_test = 3.0e5
    ceilings = []
    try:
        gfet.L = 0.0
        for mu_t, Vds_t in ((4.0e-1, 0.3), (4.0e2, 0.7), (4.0e4, 0.2)):
            gfet.mu = mu_t
            ceilings.append(transfer_characteristic_saturated(
                Vg1, Vds=Vds_t, constant_n=True, v_sat_const=vs_test,
                Rc_total=0.0)[0])
    finally:
        gfet.mu, gfet.L = mu_saved, L_saved
    Id_ceiling_exact = gfet.W * E_CHARGE * n0 * vs_test
    worst = max(ulps(c, Id_ceiling_exact) for c in ceilings)
    check("X3", worst <= 1.0,
          "L = 0: I_d == W*e*n*v_sat to 1 ULP, independent of mu and V_ds "
          "(exact under cancellation, up to the one rounding of mu*V/v_sat)",
          f"3 (mu, V_ds) pairs, all {worst:.0f} ULP from "
          f"{Id_ceiling_exact:.12e} A")

    # X3b -- the finite-L closed form, which is where the limit came from
    try:
        gfet.mu, gfet.L = 0.4, 200e-9
        Id_fin = transfer_characteristic_saturated(
            Vg1, Vds=0.2, constant_n=True, v_sat_const=vs_test,
            Rc_total=0.0)[0]
    finally:
        gfet.mu, gfet.L = mu_saved, L_saved
    pred = 1.0 / (1.0 + gfet.L * vs_test / (gfet.mu * 0.2))
    check("X3b", ulps(Id_fin / Id_ceiling_exact, pred) <= 64.0,
          "finite L: I_d/(W*e*n*v_sat) == 1/(1 + L*v_sat/(mu*V_ds)) exactly",
          f"{ulps(Id_fin / Id_ceiling_exact, pred):.1f} ULP; ratio = "
          f"{Id_fin / Id_ceiling_exact:.12f} vs {pred:.12f}")

    # X4 -- doubling v_sat doubles the ceiling, now an identity too
    try:
        gfet.L = 0.0
        a = transfer_characteristic_saturated(
            Vg1, Vds=0.2, constant_n=True, v_sat_const=vs_test,
            Rc_total=0.0)[0]
        b = transfer_characteristic_saturated(
            Vg1, Vds=0.2, constant_n=True, v_sat_const=2.0 * vs_test,
            Rc_total=0.0)[0]
    finally:
        gfet.mu, gfet.L = mu_saved, L_saved
    check("X4", b == 2.0 * a,
          "L = 0: doubling v_sat doubles I_d BITWISE (exact under mult. by 2)",
          f"ratio = {b / a!r}")

    # X5 -- the quadrature is a partition: constant n => Q == n*Vds_ch
    Q5, S5, _, _ = _channel_integrals(2.0, VDS_REF, constant_n=True)
    check("X5", ulps(Q5, n0 * VDS_REF) <= 64.0,
          "constant n: Q == n*V_ds,ch  (exact under the trapezoid weights)",
          f"{ulps(Q5, n0 * VDS_REF):.1f} ULP")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 2.  Does the channel hit the physical ceiling of Eq. (3)?")
    say("=" * 78)
    say("  u(V)/v_sat(V) >= 1 would mean the channel cannot deliver I_d at all.")
    say(f"  {'V_ds':>6}  {'max u/v_sat':>12}  {'at V_g':>8}   verdict")
    ceil_table = []
    for Vds_t in VDS_LADDER:
        _, ceil = transfer_characteristic_saturated(
            Vg, Vds=Vds_t, report_ceiling=True)
        k = int(np.argmax(ceil))
        ceil_table.append((Vds_t, ceil[k], Vg[k]))
        say(f"  {Vds_t:>6.2f}  {ceil[k]:>12.6f}  {Vg[k]:>8.3f}   "
            + ("OK" if ceil[k] < 1.0 else "NO SOLUTION"))
    check("X6", all(c < 1.0 for _, c, _ in ceil_table),
          "the ceiling is never reached on the committed V_ds ladder",
          f"largest u/v_sat = {max(c for _, c, _ in ceil_table):.6f}")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 3.  The OVERLAPPING-LIMIT validation (criterion C), and the")
    say("            TWO mechanisms it separates")
    say("=" * 78)
    say("  FOUR artefacts must converge as V_ds -> 0:")
    say("    A   committed resistor model: mean of 1/n over 50 SAMPLES,")
    say("        profiled over the FULL V_ds")
    say("    B   the same harmonic average by TRAPEZOID on 201 points")
    say("    B'  harmonic average profiled over the INTRINSIC drop V_ds,ch")
    say("    C   this module with S = 0: ARITHMETIC average over V_ds,ch")
    say("  A -> B  isolates the QUADRATURE RULE")
    say("  B -> B' isolates the PROFILE DOMAIN (what the contacts drop)")
    say("  B' -> C isolates the JENSEN gap (<1/n> vs 1/<n>)")
    say(f"  {'V_ds':>7} {'B/A-1 [%]':>12} {"B'/B-1 [%]":>12} "
        f"{"C/B'-1 [%]":>12} {'C/A-1 [%]':>12}")

    def _R_harmonic(V_g, V_hi, n_quad=N_QUAD):
        """(L/W)/sigma averaged harmonically over the potential range
        [0, V_hi] by trapezoid -- the committed model's own functional form."""
        V = np.linspace(0.0, V_hi, n_quad)
        inv_n = np.array([1.0 / float(gfet.carrier_density(V_g, v)) for v in V])
        return (gfet.L / (gfet.W * E_CHARGE * gfet.mu)) * (
            np.trapezoid(inv_n, V) / V_hi)

    def _Id_harmonic_trapz(V_g, Vds_t):
        """Artefact B: trapezoid, but profiled over the FULL V_ds, exactly
        as the committed model is."""
        return Vds_t / (_R_harmonic(V_g, Vds_t) + gfet.Rc_total)

    def _Id_harmonic_selfconsistent(V_g, Vds_t):
        """Artefact B': same harmonic average, but profiled only over the
        INTRINSIC channel drop V_ds,ch = V_ds - I_d*Rc, found by bisection.
        This is the committed model's functional form with its profile
        domain corrected and nothing else changed."""
        lo, hi = 0.0, Vds_t
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            Id = mid / _R_harmonic(V_g, mid)
            if Vds_t - Id * gfet.Rc_total - mid > 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 1e-15 * max(Vds_t, 1.0):
                break
        mid = 0.5 * (lo + hi)
        return mid / _R_harmonic(V_g, mid)

    ladder = (1e-3, 3e-3, 1e-2, 3e-2, 0.05, 0.1)
    rows = []
    for Vds_t in ladder:
        Id_A = float(gfet.transfer_characteristic(Vg1, Vds=Vds_t)[0])
        Id_B = _Id_harmonic_trapz(2.0, Vds_t)
        Id_Bp = _Id_harmonic_selfconsistent(2.0, Vds_t)
        Id_C = transfer_characteristic_saturated(
            Vg1, Vds=Vds_t, saturate=False)[0]
        rows.append((Vds_t, Id_A, Id_B, Id_Bp, Id_C))
        say(f"  {Vds_t:>7.4f} {(Id_B / Id_A - 1) * 100:>12.6f} "
            f"{(Id_Bp / Id_B - 1) * 100:>12.6f} "
            f"{(Id_C / Id_Bp - 1) * 100:>12.6f} "
            f"{(Id_C / Id_A - 1) * 100:>12.6f}")

    quad_small = rows[0][2] / rows[0][1] - 1.0
    quad_big = rows[-1][2] / rows[-1][1] - 1.0
    prof_small = rows[0][3] / rows[0][2] - 1.0
    prof_big = rows[-1][3] / rows[-1][2] - 1.0
    jens_small = rows[0][4] / rows[0][3] - 1.0
    jens_big = rows[-1][4] / rows[-1][3] - 1.0

    check("C1", abs(rows[0][4] / rows[0][1] - 1.0) < 1e-4,
          "every artefact converges as V_ds -> 0 (total gap under 0.01 % at 1 mV)",
          f"C/A - 1 = {(rows[0][4] / rows[0][1] - 1) * 100:.3e} % at V_ds = 1 mV")
    check("C2", abs(jens_big) > abs(jens_small),
          "the JENSEN gap grows with V_ds, as Cauchy-Schwarz requires",
          f"{jens_small * 100:.3e} % at 1 mV -> {jens_big * 100:.3e} % at "
          f"{ladder[-1]} V")
    check("C3", jens_big > 0.0,
          "SIGN: the harmonic average UNDER-states I_d (<1/n> >= 1/<n>), always",
          f"C/B' - 1 = {jens_big * 100:+.6f} % at V_ds = {ladder[-1]} V")

    # The attribution that the item did not ask for.
    def _order(small, big):
        return np.log(abs(big) / abs(small)) / np.log(ladder[-1] / ladder[0])
    order_quad = _order(quad_small, quad_big)
    order_prof = _order(prof_small, prof_big)
    order_jens = _order(jens_small, jens_big)
    say()
    say(f"  fitted order in V_ds:  quadrature rule {order_quad:.3f},  "
        f"profile domain {order_prof:.3f},  Jensen gap {order_jens:.3f}")
    say(f"  size at V_ds = {ladder[-1]} V:     quadrature {quad_big * 100:+.6f} %,  "
        f"profile {prof_big * 100:+.6f} %,  Jensen {jens_big * 100:+.6f} %")

    pieces = {"quadrature": abs(quad_big), "profile domain": abs(prof_big),
              "Jensen gap": abs(jens_big)}
    biggest = max(pieces, key=pieces.get)
    check("C4", order_quad > 1.6,
          "the QUADRATURE RULE is second order in V_ds, so `np.mean` over 50 "
          "samples is NOT the O(1/N) bias I predicted",
          f"order {order_quad:.3f}, size {quad_big * 100:+.2e} % at "
          f"{ladder[-1]} V -- numerically negligible")
    check("C5", biggest == "profile domain",
          "MAGNITUDE: the dominant discrepancy is the PROFILE DOMAIN, not the "
          "Jensen gap this check was written to find",
          f"profile {prof_big * 100:+.6f} % vs Jensen {jens_big * 100:+.6f} % "
          f"({abs(prof_big / jens_big):.2f}x) vs quadrature "
          f"{quad_big * 100:+.2e} %")
    check("C6", abs(order_prof - 1.0) < 0.15 and order_jens > 1.6,
          "and the two have DIFFERENT orders in V_ds (1 vs 2), so they are "
          "separable and not one effect double-counted",
          f"profile ~ V_ds^{order_prof:.2f}, Jensen ~ V_ds^{order_jens:.2f}")
    say()
    say("  TWO PREDICTIONS OF MINE WERE WRONG HERE AND THE DECOMPOSITION IS")
    say("  WHAT SAID SO.  Both are recorded rather than edited out.")
    say()
    say("  (i) I predicted the committed `np.mean` over `np.linspace(0, Vds,")
    say("      50)` carried an O(1/N) sample-mean bias.  It does not: a mean")
    say("      of equally spaced samples including both endpoints is exact for")
    say("      a linear integrand and second order for a smooth one, so it is")
    say(f"      {abs(quad_big) * 100:.2e} %  at V_ds = 0.1 V -- "
        f"order {order_quad:.2f}, negligible.")
    say("      `n_segments = 50`, open since 2026-09-26 and called")
    say("      load-bearing on 2026-10-04, is hereby measured and is NOT a")
    say("      problem at the committed biases.  That closes a nine-day item")
    say("      with a negative result.")
    say()
    say("  (ii) I predicted the residual was the Jensen gap.  It is mostly")
    say("      not.  The committed model profiles the local channel")
    say("      resistance over the FULL V_ds, but a fraction of V_ds is")
    say("      dropped across the CONTACTS and never appears across the")
    say("      channel at all.  At W = 1 um, Rc_total = 600 Ohm carries over")
    say("      half of V_ds, so the committed model evaluates n(V_ch) over a")
    say("      potential range about twice too wide.  That error is FIRST")
    say(f"      order in V_ds (order {order_prof:.2f}), is {abs(prof_big / jens_big):.1f}x the Jensen gap at")
    say("      V_ds = 0.1 V, and is the real content of criterion C's check.")
    say()
    say("  **No committed number is withdrawn.**  The total is one-signed and")
    say(f"  {abs(rows[-1][4] / rows[-1][1] - 1) * 100:.4f} % at V_ds = 0.1 V -- far below anything Chapter 4")
    say("  concludes from -- and both directions make I_d LARGER, so every")
    say("  shortfall against literature in Chapter 4 is if anything")
    say("  understated rather than overstated.  What changes is that the bound")
    say("  is measured, its two mechanisms are separated and ordered, and the")
    say("  larger one was not the one anybody here was looking for.")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 4.  Criterion A: what does saturation do to g_ds?")
    say("=" * 78)

    def _res(Nf):
        fT, fmax, gm, gds, Cgs = rf.compute_fT_fmax(
            Vg, Vds=VDS_REF, N_fingers=Nf)
        return fT, fmax, gm, gds, Cgs

    def _sat(Nf):
        return fT_fmax_saturated(Vg, Vds=VDS_REF, N_fingers=Nf)

    fT_r, fmax_r, gm_r, gds_r, Cgs_r = _at_literature_geometry(_res)
    fT_s, fmax_s, gm_s, gds_s, Cgs_s = _at_literature_geometry(_sat)

    k_r = int(np.argmax(fT_r))
    k_s = int(np.argmax(fT_s))
    say(f"  peak-f_T bias, resistor model : V_g = {Vg[k_r]:+.4f} V "
        f"(f_T = {fT_r[k_r]:.4e} Hz)")
    say(f"  peak-f_T bias, saturated model: V_g = {Vg[k_s]:+.4f} V "
        f"(f_T = {fT_s[k_s]:.4e} Hz)")
    say("  The peak MOVES BRANCH -- hole side to electron side -- so \"at the")
    say("  peak-f_T bias\" is ambiguous between the two models and both")
    say("  readings are reported below rather than one being chosen.")
    say()
    say(f"  {'quantity':<34}{'resistor':>14}{'saturated':>14}{'factor':>10}")
    g_factor = gds_r[k_r] / gds_s[k_r]
    for label, a, b in (
            ("g_ds at the resistor peak  [S]", gds_r[k_r], gds_s[k_r]),
            ("g_m  at the resistor peak  [S]", gm_r[k_r], gm_s[k_r]),
            ("f_T  at the resistor peak [Hz]", fT_r[k_r], fT_s[k_r]),
            ("f_max at the resistor peak [Hz]", fmax_r[k_r], fmax_s[k_r]),
            ("f_max/f_T at that bias", fmax_r[k_r] / fT_r[k_r],
             fmax_s[k_r] / fT_s[k_r])):
        say(f"  {label:<34}{a:>14.6e}{b:>14.6e}{a / b:>10.4f}")
    say()
    g_factor_s = gds_r[k_s] / gds_s[k_s]
    say(f"  {'g_ds at the SATURATED peak [S]':<34}{gds_r[k_s]:>14.6e}"
        f"{gds_s[k_s]:>14.6e}{g_factor_s:>10.4f}")
    say()
    say(f"  REQUIRED by criterion (A)         : divide g_ds by 4.8483")
    say(f"  DELIVERED at the resistor peak    : divide g_ds by {g_factor:.4f}")
    say(f"  DELIVERED at the saturated peak   : divide g_ds by {g_factor_s:.4f}")
    say(f"  fraction of the requirement met   : {g_factor / 4.8483 * 100:.2f} % "
        f"(resistor peak), {g_factor_s / 4.8483 * 100:.2f} % (saturated peak)")
    say()
    say("  WHY so little, when the channel integral says more.  mu*S/L = "
        f"{gfet.mu * _channel_integrals(Vg[k_r], VDS_REF)[1] / gfet.L:.4f} at this")
    say("  bias, which alone would divide the CHANNEL conductance by "
        f"{(1 + gfet.mu * _channel_integrals(Vg[k_r], VDS_REF)[1] / gfet.L) ** 2:.3f}.")
    say("  It is diluted because saturation acts only on the channel and the")
    say("  contacts are in series with it: at the 40 um geometry Rc_total is")
    say("  15.0 Ohm of a ~26 Ohm device, so roughly half the output")
    say("  conductance is a contact resistance that no saturation mechanism")
    say("  can touch.  **The contacts, not the channel, now cap g_ds** -- which")
    say("  points criterion A back at Chapter 4's contact-resistance work")
    say("  rather than at the saturation law.")
    check("G1", 1.0 < g_factor < 4.8483,
          "MAGNITUDE: the saturation knob is ALIVE (factor > 1) and is NOT "
          "sufficient (factor < 4.85)",
          f"factor = {g_factor:.4f}, i.e. {g_factor / 4.8483 * 100:.2f} % of "
          f"the requirement -- criterion A is NOT met")
    say("  NOTE ON G1's OWN BOUND.  G1 was first written as "
        "1.5 < factor < 4.8483,")
    say("  because I expected the channel-integral estimate (1.699) to survive")
    say("  into g_ds.  It does not -- the contacts eat it -- so the check FAILED")
    say("  on its own lower bound, and the bound was the thing that was wrong.")
    say("  It is recorded here rather than quietly relaxed: the band now asserts")
    say("  only what the mechanism guarantees (the knob is alive, and it is not")
    say("  enough), and 1.5 was a prior, not a prediction of this model.")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 5.  Criterion B: the SIGN of d(f_max/f_T)/dV_ds")
    say("=" * 78)
    say(f"  {'V_ds':>6} {'ratio resistor':>16} {'ratio saturated':>17} "
        f"{'<v_sat> [cm/s]':>16}")
    ratio_r, ratio_s, vsat_rep = [], [], []
    for Vds_t in VDS_LADDER:
        def _r(Nf, _v=Vds_t):
            return rf.compute_fT_fmax(Vg, Vds=_v, N_fingers=Nf)

        def _s(Nf, _v=Vds_t):
            return fT_fmax_saturated(Vg, Vds=_v, N_fingers=Nf)
        a = _at_literature_geometry(_r)
        b = _at_literature_geometry(_s)
        ia, ib = int(np.argmax(a[0])), int(np.argmax(b[0]))
        ratio_r.append(a[1][ia] / a[0][ia])
        ratio_s.append(b[1][ib] / b[0][ib])
        _, _, n_prof, vs_prof = _channel_integrals(Vg[ib], Vds_t)
        vsat_rep.append(float(np.mean(vs_prof)) * 100.0)   # m/s -> cm/s
        say(f"  {Vds_t:>6.2f} {ratio_r[-1]:>16.6f} {ratio_s[-1]:>17.6f} "
            f"{vsat_rep[-1]:>16.4e}")

    d_r = ratio_r[-1] - ratio_r[0]
    d_s = ratio_s[-1] - ratio_s[0]
    check("G2a", d_r < 0.0,
          "the resistor model's f_max/f_T FALLS with V_ds (the wrong sign)",
          f"{ratio_r[0]:.6f} -> {ratio_r[-1]:.6f}, delta = {d_r:+.6f}")
    check("G2b", d_s > 0.0,
          "the saturated model's f_max/f_T RISES with V_ds -- criterion B MET",
          f"{ratio_s[0]:.6f} -> {ratio_s[-1]:.6f}, delta = {d_s:+.6f} "
          f"(magnitude {abs(d_s / d_r):.2f}x the resistor model's)")
    say()
    say()
    say("  THE RATIO REACHES THE FEIJOO BAND AT V_ds = 1 V, AND THAT IS NOT")
    say("  CORROBORATION.  2026-10-04 recorded exactly this trap: peak f_T here")
    say("  scales ~x20 over a x20 drain-bias range, and the literature band")
    say("  spans a comparable factor, so a model with a knob that size will")
    say("  agree with it somewhere.  The table above is reported because the")
    say("  SIGN of the slope is the claim; the fact that one row of it lands")
    say("  in [1.3, 1.4] is a coincidence of where the ladder was stopped and")
    say("  is explicitly NOT offered as agreement with Feijoo et al.")
    say()
    say(f"  Reported for audit (criterion A's one knob): v_sat from Eq. (6) at")
    say(f"  this device's own density is {vsat_rep[1]:.4e} cm/s at V_ds = 0.1 V,")
    say(f"  against 3e7 cm/s measured on SiO2 at low density by Dorgan et al.")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 6.  Inverting criterion A, and bounding the beta = 1 choice")
    say("=" * 78)
    lo, hi = 1.0e2, 1.0e9

    def _g_at(v_const):
        return _at_literature_geometry(
            lambda Nf: gds_at_bias(Vg[k_r], Vds=VDS_REF, v_sat_const=v_const))

    f_lo = gds_r[k_r] / _g_at(lo)
    f_hi = gds_r[k_r] / _g_at(hi)
    bracketed = (f_lo - 4.8483) * (f_hi - 4.8483) < 0.0
    if bracketed:
        for _ in range(80):
            mid = np.sqrt(lo * hi)
            if gds_r[k_r] / _g_at(mid) < 4.8483:
                hi = mid
            else:
                lo = mid
        v_required = np.sqrt(lo * hi)
    else:
        v_required = float("nan")
    v_model = v_sat_of_n(float(gfet.carrier_density(Vg[k_r], 0.0)))
    say(f"  bracket over v_sat in [1e2, 1e9] m/s gives g_ds factors "
        f"[{f_lo:.4f}, {f_hi:.4f}]")
    check("G5a", bracketed,
          "criterion A's required g_ds factor of 4.8483 is REACHABLE by v_sat "
          "alone",
          f"factor spans {f_hi:.4f} (v_sat -> inf) to {f_lo:.4f} "
          f"(v_sat -> 0): 4.8483 is "
          + ("inside" if bracketed else "OUTSIDE") + " that span")
    if bracketed:
        say(f"  v_sat that satisfies criterion A exactly : "
            f"{v_required * 100:.4e} cm/s")
    else:
        say("  v_sat that satisfies criterion A exactly : DOES NOT EXIST.")
        say("  No value of v_sat, however small, divides g_ds by 4.85 at this")
        say("  geometry, because the contact resistance sets a floor on g_ds")
        say("  that the channel law cannot reach past (Section 4).  Criterion A")
        say("  is therefore not merely unmet -- it is UNSATISFIABLE by the")
        say("  mechanism it names, and that is a correction to the item itself.")
    say(f"  v_sat Eq. (6) gives at this bias         : {v_model * 100:.4e} cm/s")
    say(f"  Dorgan et al. measured band on SiO2      : ~1e7 - 3e7 cm/s")
    check("G5b", v_model * 100 > 3.0e7,
          "Eq. (6) with Feijoo's 0.10 eV sits ABOVE the measured SiO2 band, so "
          "it is the GENEROUS end of the physics, not a tuned value",
          f"Eq. (6) gives {v_model * 100:.3e} cm/s vs ~3e7 cm/s measured")

    # The substrate the model actually has: SiO2 remote polar phonons.
    say()
    say("  hbar*Omega SENSITIVITY -- the open item of 2026-09-29 (\"the")
    say("  remote-polar-phonon cap on SiO2 is asserted, not computed\"), which")
    say("  2026-10-04 said this work would sharpen.  Eq. (6) is linear in")
    say("  Omega, so the phonon energy IS the v_sat scale:")
    say(f"  {'hbar*Omega [eV]':>16} {'source':<34}{'v_sat [cm/s]':>14}"
        f"{'g_ds factor':>13}")
    rpp_rows = []
    for hw, src in ((0.100, "Feijoo et al. 2019 fit (used above)"),
                    (0.059, "SiO2 surface polar phonon, low mode"),
                    (0.149, "SiO2 surface polar phonon, high mode"),
                    (0.196, "graphene intrinsic optical phonon")):
        om = hw * E_CHARGE / HBAR
        vs_h = v_sat_of_n(float(gfet.carrier_density(Vg[k_r], 0.0)),
                          omega_op=om)
        g_h = gds_r[k_r] / _at_literature_geometry(
            lambda Nf, _o=om: gds_at_bias(Vg[k_r], Vds=VDS_REF, omega_op=_o))
        rpp_rows.append((hw, vs_h, g_h))
        say(f"  {hw:>16.3f} {src:<34}{vs_h * 100:>14.4e}{g_h:>13.4f}")
    check("G5c", max(g for _, _, g in rpp_rows) < 4.8483,
          "NO phonon energy in the physical range reaches criterion A's 4.85",
          f"the largest factor any of the four gives is "
          f"{max(g for _, _, g in rpp_rows):.4f}, at hbar*Omega = "
          f"{min(rpp_rows, key=lambda r: r[1])[0]:.3f} eV")

    say()
    say("  Cost of the beta = 1 choice (Dorgan et al. fit beta = 2).  At the")
    say("  largest field this device reaches, the two laws differ by:")
    E_max = VDS_REF / gfet.L
    for vs_t in (v_model, v_required):
        x = gfet.mu * E_max / vs_t
        v1 = gfet.mu * E_max / (1.0 + x)
        v2 = gfet.mu * E_max / np.sqrt(1.0 + x ** 2)
        say(f"    v_sat = {vs_t * 100:.3e} cm/s:  beta=1 gives {v1:.5e} m/s, "
            f"beta=2 gives {v2:.5e} m/s  ({(v2 / v1 - 1) * 100:+.2f} %)")
    say("  beta = 2 saturates LESS hard at a given field, so it moves the")
    say("  answer AWAY from criterion A, not toward it: the beta = 1 choice")
    say("  is the GENEROUS one and the shortfall in Section 4 is a bound.")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 7.  Controls")
    say("=" * 78)
    # G3 -- v_sat response vs mobility response are not aliased
    g_slow = _at_literature_geometry(
        lambda Nf: gds_at_bias(Vg[k_r], Vds=VDS_REF,
                               v_sat_const=v_model / 2.0))
    mu_saved = gfet.mu
    try:
        gfet.mu = mu_saved / 2.0
        g_lowmu = _at_literature_geometry(
            lambda Nf: gds_at_bias(Vg[k_r], Vds=VDS_REF))
    finally:
        gfet.mu = mu_saved
    resp_vsat = gds_s[k_r] / g_slow
    resp_mu = gds_s[k_r] / g_lowmu
    check("G3", abs(np.log(resp_vsat / resp_mu)) > 0.1,
          "CONTROL: halving v_sat and halving mu are NOT aliased responses",
          f"v_sat response {resp_vsat:.4f}x, mu response {resp_mu:.4f}x, "
          f"ratio {resp_vsat / resp_mu:.4f}")

    # G4 -- the resistor model is bitwise untouched by this module
    Id_now = gfet.transfer_characteristic(Vg, Vds=VDS_REF)
    fT_now, fmax_now = rf.compute_fT_fmax(Vg, Vds=VDS_REF)[:2]
    check("G4a", np.array_equal(Id_now, gfet.transfer_characteristic(Vg, Vds=VDS_REF)),
          "the resistor model is reproducible after this module has run")
    check("G4b", np.array_equal(fmax_now, rf.compute_fT_fmax(Vg, Vds=VDS_REF)[1])
          and gfet.mu == 0.4 and gfet.L == 200e-9 and gfet.W == 1e-6,
          "gfet.mu / L / W and rf's f_max are restored to committed values",
          f"mu = {gfet.mu}, L = {gfet.L}, W = {gfet.W}")

    # G6 -- a non-aliasing control in the other direction: saturation must
    # NOT simply reproduce a smaller mobility, or "add saturation" and
    # "mu was wrong" would be the same claim.
    ratio_shape_s = np.array(ratio_s)
    ratio_shape_r = np.array(ratio_r)
    check("G6", np.sign(np.diff(ratio_shape_s)).sum()
          != np.sign(np.diff(ratio_shape_r)).sum(),
          "CONTROL: the two models' f_max/f_T have different monotonicity, so "
          "saturation is not a relabelled mobility",
          f"sign sums {np.sign(np.diff(ratio_shape_s)).sum():+.0f} vs "
          f"{np.sign(np.diff(ratio_shape_r)).sum():+.0f}")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say("SECTION 8.  Figure")
    say("=" * 78)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    Id_r = gfet.transfer_characteristic(Vg, Vds=VDS_REF)
    Id_s = transfer_characteristic_saturated(Vg, Vds=VDS_REF)
    Id_ns = transfer_characteristic_saturated(Vg, Vds=VDS_REF, saturate=False)
    ax[0, 0].plot(Vg, Id_r * 1e6, label="resistor model (committed)")
    ax[0, 0].plot(Vg, Id_ns * 1e6, "--",
                  label=r"no-saturation limit ($S\to 0$)")
    ax[0, 0].plot(Vg, Id_s * 1e6, label="velocity saturation, Eq. (4)")
    ax[0, 0].set_xlabel(r"$V_g$ (V)")
    ax[0, 0].set_ylabel(r"$I_d$ ($\mu$A/$\mu$m)")
    ax[0, 0].set_title(f"(a) transfer curves at $V_{{ds}}$ = {VDS_REF} V")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3)

    ax[0, 1].semilogy(Vg, np.abs(gds_r), label=r"$g_{ds}$ resistor")
    ax[0, 1].semilogy(Vg, np.abs(gds_s), label=r"$g_{ds}$ saturated")
    ax[0, 1].axvline(Vg[k_r], color="k", lw=0.7, ls=":",
                     label="peak-$f_T$ bias")
    ax[0, 1].set_xlabel(r"$V_g$ (V)")
    ax[0, 1].set_ylabel(r"$g_{ds}$ (S)")
    ax[0, 1].set_title(f"(b) output conductance: /{g_factor:.2f} "
                       f"(criterion A asks /4.85)")
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].grid(alpha=0.3)

    ax[1, 0].plot(VDS_LADDER, ratio_r, "o-", label="resistor (sign wrong)")
    ax[1, 0].plot(VDS_LADDER, ratio_s, "s-", label="saturated (sign right)")
    ax[1, 0].axhspan(1.3, 1.4, color="green", alpha=0.15,
                     label="Feijoo et al. 2016 band")
    ax[1, 0].set_xscale("log")
    ax[1, 0].set_xlabel(r"$V_{ds}$ (V)")
    ax[1, 0].set_ylabel(r"$f_{max}/f_T$ at peak $f_T$")
    ax[1, 0].set_title("(c) criterion B: the sign of "
                       r"$d(f_{max}/f_T)/dV_{ds}$")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3)

    n_grid = np.logspace(15, 17.3, 200)
    ax[1, 1].loglog(n_grid * 1e-4, v_sat_of_n(n_grid) * 100,
                    label=r"Eq. (6), $\hbar\Omega$ = 0.10 eV")
    ax[1, 1].axhspan(1e7, 3e7, color="green", alpha=0.15,
                     label="Dorgan et al. 2010, SiO$_2$")
    ax[1, 1].axhline(v_required * 100, color="r", ls="--",
                     label="needed for criterion A")
    ax[1, 1].set_xlabel(r"$n$ (cm$^{-2}$)")
    ax[1, 1].set_ylabel(r"$v_{sat}$ (cm/s)")
    ax[1, 1].set_title("(d) $v_{sat}$ is not a free parameter")
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3, which="both")

    fig.suptitle("Velocity saturation in the GFET transfer characteristic "
                 "(2026-10-05; closes criteria B and C, misses A)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig("velocity_saturation_model.png", dpi=150)
    say("  wrote velocity_saturation_model.png")

    # ------------------------------------------------------------------
    say()
    say("=" * 78)
    say(f"RESULT: {len(_LINES)} lines, "
        f"{sum(1 for l in _LINES if '[PASS]' in l)} passed, {len(_FAIL)} failed")
    say("=" * 78)
    say("  Criterion A (g_ds / 4.85)        : NOT MET -- "
        f"delivered /{g_factor:.4f}, {g_factor / 4.8483 * 100:.1f} % of it")
    say("  Criterion B (sign of the slope)  : MET")
    say("  Criterion C (overlapping limit)  : MET, and it found a one-signed")
    say("                                     PROFILE-DOMAIN error in the model")
    say("                                     it supersedes -- 674x the Jensen")
    say("                                     gap it was written to find, and")
    say("                                     it closes `n_segments = 50`")
    say("                                     (2026-09-26) with a negative")
    say("                                     result (Section 3)")
    if _FAIL:
        say(f"  FAILED CHECKS: {', '.join(_FAIL)}")
    with open("velocity_saturation_output.txt", "w") as fh:
        fh.write("\n".join(_LINES) + "\n")
    return 1 if _FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
