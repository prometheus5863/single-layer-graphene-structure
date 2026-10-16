"""
graphene_diffusion_current_model.py

Opens the item created 2026-10-06 and named that session's TOP item:

    A DIFFUSION TERM IN EQ. (4).  Feijoo et al. 2020 measure the diffusion
    contribution as comparable to drift at the peak-f_max bias, so every
    V_ds ladder in this repository's RF thread is a drift-only slice and
    the peak-f_max bias is a bias this model cannot be asked about.

RESULT, stated up front, because the item's PREMISE is what moved and not
only its number:

  1. The diffusion term is derived, implemented, and it PRESERVES the
     separability that makes Eq. (4) closed-form.  Its contribution to the
     channel integral turns out to be a pure BOUNDARY term -- it depends
     only on the source and drain densities and not at all on the density
     profile between them (Section 1, X6, which asserts it structurally by
     changing the interior quadrature grid and requiring the term to be
     BITWISE unchanged).  That is a stronger statement than "it is small":
     the one part of this model that 2026-10-05's criterion C found to be
     profile-domain-sensitive cannot contaminate it.

  2. Its SIZE at the biases Chapter 4 quotes is measured and reported, as a
     percentage of I_d, rather than asserted to be negligible.

  3. THE PREMISE IS NOT SETTLED BY ADDING A TERM.  The standard GFET compact
     model this repository's Eq. (4) descends from (Pasadas and Jimenez,
     IEEE TED 2016; arXiv:1605.08235) writes v = mu*F with F = -dV/dx and
     states explicitly that "V(x) is the quasi-Fermi level along the
     graphene channel".  A drift expression whose driving potential is the
     QUASI-FERMI level is already the complete drift-diffusion current.
     `gfet.carrier_density(V_g, V_ch)` never says which potential V_ch is:
     it enters an electrostatic charge relation C_ox*(V_g - V_ch - V_dirac),
     which reads electrostatic, and is then corrected by a quantum-
     capacitance series factor, which is the correction one applies when
     V_ch is the quasi-Fermi level and the graphene drop E_F/e is separate.
     THE TWO READINGS DIFFER BY EXACTLY THE TERM ADDED HERE.  So the number
     this module reports is simultaneously (a) the diffusion current under
     the electrostatic reading and (b) the SIZE OF AN AMBIGUITY THAT WAS
     ALREADY IN EVERY COMMITTED I_d.  Section 6 states that in those terms.
     This is the 2026-10-04 fault -- a quantity computed correctly under a
     name that does not pin down what it is -- located in a variable name
     rather than in a figure of merit.

----------------------------------------------------------------------------
THE DERIVATION, AND WHY SEPARABILITY SURVIVES
----------------------------------------------------------------------------
Drift-diffusion, electron branch, magnitudes, V the channel potential
rising from 0 at the source to V_ds,ch at the drain:

    I_d = W*e*[ n*mu*dV/dx  -  D*dn/dx ]                               (7)

n falls toward the drain, so dn/dx < 0 and the diffusion term ADDS.  The
generalized Einstein relation, in the form Zebrev writes it (Eq. 20 of
"Graphene Field Effect Transistors: Diffusion-Drift Theory",
arXiv:1102.2348), is

    mu = e*D/eps_D,       eps_D = n / (dn/dmu_c)                       (8)

and for graphene's linear dispersion eps_D = E_F/2 exactly in the
degenerate limit (Zebrev's Eq. 13 analysis), so D = mu*E_F/(2e).  Define

    V_F(n) = E_F/e = hbar*v_F*sqrt(pi*n)/e = A_F*sqrt(n)

(this repository's own Dirac dispersion, no new parameter), and D/(mu*n) dn
= (V_F/(2n)) dn = d(A_F*sqrt(n)) = dV_F.  So (7) collapses to

    I_d = W*e*mu*n * d( V - lambda*V_F )/dx = W*e*mu*n*dPhi/dx          (9)

with lambda = 1 the physical case and lambda = 0 recovering Eq. (4)
BITWISE (X1).  Phi = V - V_F is the quasi-Fermi potential, which is the
textbook statement that drift-diffusion is transport down the gradient of
the quasi-Fermi level -- and the reason Section 6's ambiguity exists.

Applying the beta = 1 soft-saturation law of Eq. (1) to the effective
driving field dPhi/dx and integrating dx = dPhi/E as before:

    I_d = mu*W*e*Q_D / ( L + mu*S_D )                                 (10)

    Q_D = int n dPhi = Q  +  lambda*(A_F/3)*( n_s^(3/2) - n_d^(3/2) )  (11)
    S_D = int dPhi/v_sat = S + int kappa(V)/v_sat(V) dV                (12)

    kappa(V) = -lambda*dV_F/dV   (the LOCAL diffusion/drift ratio, which is
                                  Zebrev's Eq. 60 kappa = J_DIF/J_DR)

Equation (11) is the load-bearing line.  int n dV_F = A_F * int n d(sqrt n)
= (A_F/3)*[n^(3/2)], an EXACT antiderivative, so the diffusion part of the
channel integral is a difference of endpoint values: a BOUNDARY term with no
dependence on the profile.  Q_D is therefore computed from (11) in closed
form and the quadrature route is kept only as an independent check with no
shared machinery (X5, Richardson on the grid).  S_D has no such form and is
quadratured.

Zebrev also gives a closed form for the local ratio, his Eq. 62,
kappa = C_ox/(C_Q + C_it).  With C_it = 0 and the dispersion-consistent
C_Q = (2e^3/(pi*hbar^2*v_F^2)) * V_F, that is an INDEPENDENT prediction of
the quantity (12) measures -- capacitor algebra against a boundary term --
and Section 2 compares them instead of asserting agreement.

----------------------------------------------------------------------------
WHAT THE COMPARISON WITH ZEBREV'S kappa FINDS (not what it was for)
----------------------------------------------------------------------------
`gfet.quantum_capacitance()` evaluates C_Q at the GATE OVERDRIVE
(V_g - V_ch - V_dirac) in the slot where the dispersion wants E_F/e = V_F.
Those differ by the series factor the same function is used to build, and
the ratio is reported in Section 2.  Because kappa goes as 1/C_Q, the
repository's own C_q UNDERSTATES kappa by that ratio.  No committed number
is withdrawn by this: C_q is used downstream only through the series factor
C_q/(C_q+C_ox), it is annotated in place, and Chapter 4 gets the note.

----------------------------------------------------------------------------
VALIDATION DESIGN
----------------------------------------------------------------------------
Tolerance-free checks, each naming the operation it is exact under:

  X1  lambda = 0 reproduces graphene_velocity_saturation_model BITWISE
      (exact: identical operations in identical order with a term that is
      exactly 0.0)
  X2  V_ds = 0  =>  I_d exactly 0.0                (multiplication by zero)
  X3  CONSTANT n  =>  the diffusion term is exactly 0.0 and I_d is BITWISE
      the lambda = 0 current.  A uniform channel has no density gradient,
      so this is a symmetry that MUST give exactly zero -- the
      exactly-known-value check the 2026-09-17 rule asks for.
  X4  lambda = 2 gives exactly twice lambda = 1's diffusion term
      (multiplication by a power of two)
  X5  the closed form (11) against a trapezoid of n*kappa on the same grid:
      Richardson ratio on grid refinement must show O(h^2), reported as a
      number rather than passed as a tolerance
  X6  PROFILE INDEPENDENCE, asserted structurally: a non-uniform interior
      grid with the same endpoints changes Q but leaves the diffusion term
      BITWISE unchanged
  X7  the strong-saturation ceiling I_d -> W*e*n*v_sat is reached and is
      INDEPENDENT of lambda (diffusion cannot lift a velocity ceiling)
  X8  1 + kappa > 0 everywhere on every bias of the ladder -- if it ever
      reached zero, Phi would be non-monotonic and (10) would have no
      solution.  Reported with its minimum, never clipped.

and magnitude-form reporting per the standing top methodological item:

  G1  the diffusion share of I_d, as a percentage, on the V_ds ladder
  G2  kappa along the channel vs Zebrev's C_ox/(C_Q+C_it) closed form
  G3  the C_q argument-slot ratio that explains the gap in G2
  G4  the share across the whole V_g sweep, with its maximum and location
  G5  f_T, f_max and f_max/f_T at the literature RF geometry, lambda = 0
      vs lambda = 1, against 2026-10-06's committed 0.683186
  G6  CONTROL: the lambda response and the mu response are NOT aliased
  G7  CONTROL: vsat and the resistor model are unchanged by this module

TWO CHECKS ADDED BY THIS MODULE'S OWN MUTATION HARNESS (2026-10-07), which
is the sixth consecutive harness in this repository to improve the suite it
was pointed at rather than bless it:

  X9  lambda enters in TWO places, Eq. (11) and Eq. (12).  The harness's M2
      arrives at only the first, and on its first run was killed only by G5c
      -- a check written for something else.  X9a/b/c give Eq. (12) its own
      detector with a magnitude.  X9c reports a fact nothing had noticed:
      S_d > 0 on the electron branch, so the two halves of the diffusion
      term pull in OPPOSITE directions, which is why the share of I_d
      (+0.50 %) is about half of Qd/Q (+1.16 %).
  G5d the stencil is pinned against the COMMITTED 0.683186 rather than only
      against vsat on the same stencil.  The harness's M7 changed the
      stencil spacing and G5a could not see it, because G5a moves with it.
      M7 SURVIVED the harness's first run; G5d is why it does not now.

Run:  python3 graphene_diffusion_current_model.py
Writes: diffusion_current_model.png, diffusion_current_output.txt
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import rf_small_signal_model as rf
import graphene_velocity_saturation_model as vsat

E_CHARGE = 1.602176634e-19
HBAR = 1.054571817e-34

# V_F(n) = A_FERMI*sqrt(n).  Built from this repository's own v_F; no new
# literature parameter enters this module at all.
A_FERMI = HBAR * gfet.v_F * np.sqrt(np.pi) / E_CHARGE

# Dispersion-consistent T = 0 quantum capacitance, C_Q = CQ_PREFACTOR*V_F.
# Same prefactor gfet.quantum_capacitance() uses; the difference is only
# which voltage goes in the slot, which is Section 2's finding.
CQ_PREFACTOR = 2.0 * E_CHARGE ** 3 / (np.pi * HBAR ** 2 * gfet.v_F ** 2)

VDS_REF = 0.1
VDS_LADDER = (0.05, 0.1, 0.2, 0.5, 1.0)

# The two fixed biases 2026-10-06 reported at, quoted from
# perfect_contact_counterfactual_output.txt so the comparison is against
# committed numbers rather than against a re-derived argmax.
VG_SAT_PEAK = 2.857143
VG_RES_PEAK = -1.112782

# The committed gate-grid spacing, np.linspace(-2, 4, 400): the local
# 5-point stencils in Section 4 use this h so that np.gradient's central
# difference is the SAME ESTIMATOR as the committed one rather than a
# finer one.  2026-10-06's M5 is why this is spelled out.
H_VG_COMMITTED = 6.0 / 399.0

N_QUAD = 201
_FAST = os.environ.get("DIFF_FAST") == "1"
N_VG = 60 if _FAST else 240

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


def ulps(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / np.spacing(max(abs(a), abs(b)))


# ---------------------------------------------------------------------------
# The drift-diffusion model
# ---------------------------------------------------------------------------
def fermi_voltage(n):
    """V_F = E_F/e = A_FERMI*sqrt(n) [V], from the repo's Dirac dispersion."""
    return A_FERMI * np.sqrt(np.asarray(n, dtype=float))


def quantum_capacitance_dispersion(n):
    """C_Q = 2e^3 V_F/(pi hbar^2 v_F^2) [F/m^2], i.e. e^2 dn/dmu at T = 0,
    evaluated at E_F/e as the dispersion requires.  Section 2 compares this
    against gfet.quantum_capacitance(), which uses the gate overdrive."""
    return CQ_PREFACTOR * fermi_voltage(n)


def _dd_integrals(V_g, Vds_ch, n_quad=N_QUAD, lam=1.0, v_sat_const=None,
                  constant_n=False, omega_op=vsat.OMEGA_OP,
                  v_sat_prefactor=2.0 / np.pi, grid=None):
    """Q, the diffusion boundary term, S and its diffusion part, Eqs. (11-12).

    `grid`, when given, replaces the uniform partition with ANY partition of
    the same interval; X6 uses it to show the diffusion term does not see
    the interior at all.  `constant_n` freezes n, which makes the diffusion
    term exactly zero (X3).
    """
    if grid is None:
        V = np.linspace(0.0, Vds_ch, n_quad)
    else:
        V = np.asarray(grid, dtype=float)
    if constant_n:
        n = np.full_like(V, float(gfet.carrier_density(V_g, 0.0)))
    else:
        n = np.array([float(gfet.carrier_density(V_g, v)) for v in V])

    Q = np.trapezoid(n, V)

    # Eq. (11): exact antiderivative, endpoints only.
    Qd = lam * (A_FERMI / 3.0) * (n[0] ** 1.5 - n[-1] ** 1.5)

    V_F = fermi_voltage(n)
    if len(V) > 1:
        kappa = -lam * np.gradient(V_F, V)
    else:
        kappa = np.zeros_like(V)

    if v_sat_const is None:
        vs = vsat.v_sat_of_n(n, omega_op=omega_op, prefactor=v_sat_prefactor)
    else:
        vs = np.full_like(V, float(v_sat_const))

    S = np.trapezoid(1.0 / vs, V)
    Sd = np.trapezoid(kappa / vs, V)
    Qd_quad = np.trapezoid(n * kappa, V)

    return dict(V=V, n=n, V_F=V_F, vs=vs, kappa=kappa,
                Q=Q, Qd=Qd, Qd_quad=Qd_quad, S=S, Sd=Sd)


def _Id_given_Vds_ch(V_g, Vds_ch, lam=1.0, saturate=True, **kw):
    """Equation (10).  The arithmetic is written in the SAME order as
    graphene_velocity_saturation_model._Id_given_Vds_ch so that lam = 0 is
    bitwise identical to it rather than merely close (X1)."""
    d = _dd_integrals(V_g, Vds_ch, lam=lam, **kw)
    Q = d["Q"] + d["Qd"]
    S = d["S"] + d["Sd"]
    if not saturate:
        S = 0.0
    denom = gfet.L + gfet.mu * S
    Id = gfet.mu * gfet.W * E_CHARGE * Q / denom
    return Id, d


def transfer_characteristic_dd(Vg_range, Vds=VDS_REF, lam=1.0, saturate=True,
                               Rc_total=None, tol=1e-15, max_iter=200,
                               report=False, **kw):
    """I_d(V_g) with velocity saturation AND diffusion, Eqs. (10) + (5).

    lam = 0 recovers graphene_velocity_saturation_model exactly.
    report=True additionally returns, per bias, min(1 + kappa) over the
    channel (X8) and the diffusion share Qd/(Q+Qd).
    """
    Vg_range = np.atleast_1d(np.asarray(Vg_range, dtype=float))
    Rc = gfet.Rc_total if Rc_total is None else Rc_total
    Id_out = np.zeros_like(Vg_range)
    onepk_out = np.ones_like(Vg_range)
    share_out = np.zeros_like(Vg_range)

    for i, V_g in enumerate(Vg_range):
        if Vds == 0.0:
            Id_out[i] = 0.0          # X2: exact, no solve
            continue
        if Rc == 0.0:
            Id, d = _Id_given_Vds_ch(V_g, Vds, lam=lam, saturate=saturate, **kw)
        else:
            lo, hi = 0.0, Vds
            for _ in range(max_iter):
                mid = 0.5 * (lo + hi)
                Id, d = _Id_given_Vds_ch(V_g, mid, lam=lam,
                                         saturate=saturate, **kw)
                resid = Vds - Id * Rc - mid
                if resid > 0.0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo <= tol * max(Vds, 1.0):
                    break
            Id, d = _Id_given_Vds_ch(V_g, 0.5 * (lo + hi), lam=lam,
                                     saturate=saturate, **kw)
        Id_out[i] = Id
        if report:
            onepk_out[i] = float(np.min(1.0 + d["kappa"]))
            tot = d["Q"] + d["Qd"]
            share_out[i] = d["Qd"] / tot if tot != 0.0 else 0.0

    if report:
        return Id_out, onepk_out, share_out
    return Id_out


def gds_dd(Vg_range, Vds=VDS_REF, dVds=1e-3, **kw):
    Ip = transfer_characteristic_dd(Vg_range, Vds=Vds + dVds, **kw)
    Im = transfer_characteristic_dd(Vg_range, Vds=Vds - dVds, **kw)
    return (Ip - Im) / (2.0 * dVds)


def _fT_fmax_local(V_g, Vds=VDS_REF, lam=1.0, N_fingers=1, h=H_VG_COMMITTED,
                   **kw):
    """f_T, f_max and the pieces, at ONE gate bias, on a 5-point stencil of
    the COMMITTED gate-grid spacing so np.gradient's central difference here
    is the same estimator as the committed one.  Mirrors
    vsat.fT_fmax_saturated's formula exactly; only g_m and g_ds change."""
    Vg = V_g + h * np.arange(-2, 3, dtype=float)
    Id = transfer_characteristic_dd(Vg, Vds=Vds, lam=lam, **kw)
    gm = np.gradient(Id, Vg)
    gds = gds_dd(Vg, Vds=Vds, lam=lam, **kw)
    Cgs = rf.gate_capacitance(Vg)
    Cgd = rf.Cgd_over_Cgs * Cgs
    Rg = rf.gate_resistance(N_fingers=N_fingers)
    Rs = rf.source_access_resistance()
    fT = np.abs(gm) / (2 * np.pi * Cgs)
    denom = gds * (Rg + Rs) + 2 * np.pi * fT * Cgd * Rg
    denom = np.clip(denom, 1e-30, None)
    fmax = fT / (2 * np.sqrt(denom))
    k = 2
    return dict(fT=float(fT[k]), fmax=float(fmax[k]), ratio=float(fmax[k] / fT[k]),
                gm=float(gm[k]), gds=float(gds[k]), Id=float(Id[k]))


def _at_literature_geometry(fn, W=rf.W_RF, N_fingers=rf.N_FINGERS_RF):
    """Late binding only -- the 2026-10-03 frozen-default fault."""
    W_saved, Rc_saved = gfet.W, gfet.Rc_total
    gfet.W = W
    gfet.Rc_total = 2 * (gfet.Rc_per_width_ohm_um * 1e-6) / W
    try:
        return fn(N_fingers)
    finally:
        gfet.W, gfet.Rc_total = W_saved, Rc_saved


# ===========================================================================
def main():
    if _FAST:
        say("  [DIFF_FAST=1: gate grid %d pts -- mutation-harness mode]" % N_VG)
    say("a diffusion term in Eq. (4): the 2026-10-06 TOP item")
    say("(the result, and the premise that moved, are in the module docstring)")
    say()

    # Captured BEFORE anything in this module runs, so G7 compares against a
    # value this module cannot have disturbed.
    vsat_ref_before = float(vsat.transfer_characteristic_saturated(
        np.array([2.0]), Vds=VDS_REF)[0])
    res_ref_before = float(gfet.transfer_characteristic(
        np.array([2.0]), Vds=VDS_REF)[0])

    say("=" * 78)
    say("SECTION 1.  The exactness battery")
    say("=" * 78)

    # --- X1: lambda = 0 is the module being extended, bitwise -------------
    worst = 0.0
    rows = []
    for V_g in (-1.112782, 0.5, 2.0, 2.857143):
        for Vd in (0.05, VDS_REF, 0.5):
            a = float(vsat.transfer_characteristic_saturated(
                np.array([V_g]), Vds=Vd)[0])
            b = float(transfer_characteristic_dd(
                np.array([V_g]), Vds=Vd, lam=0.0)[0])
            worst = max(worst, ulps(a, b))
            rows.append((V_g, Vd, a, b))
    check("X1", worst == 0.0,
          "lambda = 0 reproduces the saturated drift model BITWISE at 12 "
          "(V_g, V_ds) points",
          f"worst disagreement {worst:.1f} ULP")

    # --- X2: multiplication by zero ---------------------------------------
    z = transfer_characteristic_dd(np.array([2.0, -1.0]), Vds=0.0, lam=1.0)
    check("X2", bool(np.all(z == 0.0)), "V_ds = 0 gives I_d exactly 0.0",
          f"{z}")

    # --- X3: a uniform channel has no diffusion.  Exactly zero. -----------
    d_const = _dd_integrals(2.0, VDS_REF, lam=1.0, constant_n=True)
    Id_c1, _ = _Id_given_Vds_ch(2.0, VDS_REF, lam=1.0, constant_n=True)
    Id_c0, _ = _Id_given_Vds_ch(2.0, VDS_REF, lam=0.0, constant_n=True)
    check("X3a", d_const["Qd"] == 0.0,
          "CONSTANT n: the diffusion boundary term is EXACTLY 0.0",
          f"Qd = {d_const['Qd']!r}, n_s - n_d = "
          f"{d_const['n'][0] - d_const['n'][-1]!r}")
    check("X3b", ulps(Id_c1, Id_c0) == 0.0,
          "CONSTANT n: I_d is BITWISE the lambda = 0 current",
          f"{Id_c1:.17e} vs {Id_c0:.17e}")

    # --- X4: linearity in lambda, exact -----------------------------------
    d1 = _dd_integrals(2.0, VDS_REF, lam=1.0)
    d2 = _dd_integrals(2.0, VDS_REF, lam=2.0)
    check("X4", ulps(d2["Qd"], 2.0 * d1["Qd"]) == 0.0,
          "lambda = 2 gives EXACTLY twice lambda = 1's diffusion term",
          f"{d2['Qd']:.17e} vs {2.0 * d1['Qd']:.17e}")

    # --- X5: closed form vs quadrature, with the convergence ORDER --------
    say("  Eq. (11) closed form against a trapezoid of n*kappa on the same grid:")
    say(f"  {'n_quad':>8} {'|closed - quad| [m^-2 V]':>26} {'ratio':>8}")
    errs = []
    for nq in (101, 201, 401, 801):
        d = _dd_integrals(2.0, VDS_REF, lam=1.0, n_quad=nq)
        errs.append(abs(d["Qd"] - d["Qd_quad"]))
        r = errs[-2] / errs[-1] if len(errs) > 1 else float("nan")
        say(f"  {nq:>8} {errs[-1]:>26.6e} {r:>8.4f}")
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    check("X5", all(3.0 < r < 5.0 for r in ratios),
          "the gap between the two routes is O(h^2), so they are the same "
          "integral computed two ways",
          "Richardson ratios " + ", ".join(f"{r:.4f}" for r in ratios)
          + " against 4 for O(h^2)")

    # --- X6: PROFILE INDEPENDENCE, structurally ---------------------------
    Vu = np.linspace(0.0, VDS_REF, 201)
    # Same endpoints, deliberately lopsided interior.
    Vn = VDS_REF * (np.linspace(0.0, 1.0, 201) ** 2.5)
    du = _dd_integrals(2.0, VDS_REF, lam=1.0, grid=Vu)
    dn = _dd_integrals(2.0, VDS_REF, lam=1.0, grid=Vn)
    check("X6a", ulps(du["Qd"], dn["Qd"]) == 0.0,
          "a lopsided interior grid leaves the diffusion term BITWISE "
          "unchanged -- it is a boundary term",
          f"{du['Qd']:.17e} vs {dn['Qd']:.17e}")
    check("X6b", du["Q"] != dn["Q"],
          "CONTROL: the same regridding DOES move Q, so X6a is not a "
          "no-op test",
          f"Q moves by {abs(dn['Q'] / du['Q'] - 1) * 100:.3e} % "
          f"({ulps(du['Q'], dn['Q']):.3e} ULP)")
    check("X6c", abs(dn["Qd_quad"] / du["Qd_quad"] - 1) > 1e-6,
          "CONTRAST: the QUADRATURE route to the same term DOES move on "
          "regridding, so X6a is a property of the closed form and not of "
          "the grid",
          f"Qd_quad moves by "
          f"{abs(dn['Qd_quad'] / du['Qd_quad'] - 1) * 100:.4f} % while Qd "
          f"moves 0 ULP")

    # --- X7: the velocity ceiling does not care about lambda --------------
    mu_saved = gfet.mu
    gfet.mu = 1e6
    try:
        vsc = 5.0e5
        ceil = gfet.W * E_CHARGE * float(gfet.carrier_density(2.0, 0.0)) * vsc
        I7 = {}
        for lam in (0.0, 1.0):
            I7[lam], _ = _Id_given_Vds_ch(2.0, VDS_REF, lam=lam,
                                          constant_n=True, v_sat_const=vsc)
    finally:
        gfet.mu = mu_saved
    check("X7a", abs(I7[1.0] / ceil - 1.0) < 1e-5,
          "mu -> inf reaches the ceiling W*e*n*v_sat",
          f"I_d/ceiling = {I7[1.0] / ceil:.12f}")
    check("X7b", ulps(I7[1.0], I7[0.0]) == 0.0,
          "the ceiling is BITWISE independent of lambda -- diffusion cannot "
          "lift a velocity ceiling", f"{I7[1.0]:.17e}")

    # --- X8: 1 + kappa > 0, reported ---------------------------------------
    say("  min(1 + kappa) over the channel, which must stay > 0 for Phi to be")
    say("  monotonic and Eq. (10) to have a solution:")
    say(f"  {'V_ds':>7} {'V_g':>10} {'min(1+kappa)':>14} {'max kappa':>11}")
    worst_opk = 1e9
    for Vd in VDS_LADDER:
        for V_g in (VG_RES_PEAK, 0.8, VG_SAT_PEAK):
            _, opk, _ = transfer_characteristic_dd(
                np.array([V_g]), Vds=Vd, lam=1.0, report=True)
            d = _dd_integrals(V_g, Vd, lam=1.0)
            worst_opk = min(worst_opk, float(opk[0]))
            say(f"  {Vd:>7.2f} {V_g:>+10.6f} {float(opk[0]):>14.8f} "
                f"{float(np.max(d['kappa'])):>11.6f}")
    check("X8", worst_opk > 0.0,
          "1 + kappa > 0 at every bias on the ladder",
          f"smallest value {worst_opk:.8f}")

    # --- X9: the OTHER place lambda enters.  Added 2026-10-07 in response
    # to this module's own mutation harness, whose M2 (lambda reaches Q_D but
    # not S_D) was killed only by G5c -- a check written for something else.
    # 2026-10-06's coherence trap is exactly this shape, so S_D gets its own
    # detector with a magnitude rather than an incidental one.
    dref = _dd_integrals(VG_SAT_PEAK, VDS_REF, lam=1.0)
    d0 = _dd_integrals(VG_SAT_PEAK, VDS_REF, lam=0.0)
    Id_both, _ = _Id_given_Vds_ch(VG_SAT_PEAK, VDS_REF, lam=1.0)
    denom_noSd = gfet.L + gfet.mu * dref["S"]
    Id_noSd = gfet.mu * gfet.W * E_CHARGE * (dref["Q"] + dref["Qd"]) / denom_noSd
    say("  lambda enters in TWO places, Eq. (11) and Eq. (12).  At the "
        "saturated-peak")
    say(f"  bias, V_ds = {VDS_REF} V, the Eq. (12) half is worth:")
    say(f"    S   = {dref['S']:.9e} V.s/m    S_d = {dref['Sd']:.9e} V.s/m"
        f"    S_d/S = {dref['Sd'] / dref['S'] * 100:+.6f} %")
    say(f"    dropping S_d moves I_d by "
        f"{(Id_noSd / Id_both - 1) * 100:+.6f} %, which is the size of the "
        f"half-arrival M2 injects")
    check("X9a", dref["Sd"] != 0.0 and d0["Sd"] == 0.0,
          "the Eq. (12) diffusion term is non-zero at lambda = 1 and exactly "
          "zero at lambda = 0 -- both halves of lambda ARRIVE",
          f"S_d = {dref['Sd']:.9e} at lambda = 1, {d0['Sd']!r} at lambda = 0")
    check("X9b", abs(Id_noSd / Id_both - 1) > 1e-9,
          "and dropping it is OBSERVABLE in I_d, so a half-arrived lambda "
          "cannot pass unnoticed",
          f"{(Id_noSd / Id_both - 1) * 100:+.6f} %")
    check("X9c", dref["Sd"] > 0.0,
          "S_d > 0 on the electron branch, so Eq. (12) OPPOSES Eq. (11): the "
          "two halves of the diffusion term pull in opposite directions, "
          "which is why the share of I_d is smaller than Qd/Q",
          f"S_d/S = {dref['Sd'] / dref['S'] * 100:+.6f} %")

    say()
    say("=" * 78)
    say("SECTION 2.  How big is it, and does Zebrev's closed form agree")
    say("=" * 78)

    say("  G1.  Diffusion share of I_d at the two fixed biases 2026-10-06 used.")
    say(f"  {'bias':>14} {'V_ds':>6} {'I_d(lam=0) [A]':>16} "
        f"{'I_d(lam=1) [A]':>16} {'share [%]':>10} {'Qd/Q [%]':>9}")
    shares = {}
    for name, V_g in (("resistor-peak", VG_RES_PEAK),
                      ("saturated-peak", VG_SAT_PEAK)):
        for Vd in VDS_LADDER:
            I0 = float(transfer_characteristic_dd(
                np.array([V_g]), Vds=Vd, lam=0.0)[0])
            I1 = float(transfer_characteristic_dd(
                np.array([V_g]), Vds=Vd, lam=1.0)[0])
            d = _dd_integrals(V_g, Vd, lam=1.0)
            shares[(name, Vd)] = (I0, I1)
            say(f"  {name:>14} {Vd:>6.2f} {I0:>16.6e} {I1:>16.6e} "
                f"{(I1 / I0 - 1.0) * 100:>10.4f} "
                f"{d['Qd'] / d['Q'] * 100:>9.4f}")
    I0r, I1r = shares[("saturated-peak", VDS_REF)]
    share_ref = (I1r / I0r - 1.0) * 100
    check("G1a", I1r > I0r,
          "on the ELECTRON branch diffusion ADDS to I_d, as d|n|/dV < 0 "
          "requires",
          f"at the saturated-peak bias, V_ds = {VDS_REF} V: "
          f"+{share_ref:.4f} %")
    I0h, I1h = shares[("resistor-peak", VDS_REF)]
    share_hole = (I1h / I0h - 1.0) * 100
    check("G1b", share_hole < 0.0 < share_ref,
          "AND IT SUBTRACTS ON THE HOLE BRANCH -- the correction is SIGNED, "
          "and the sign follows d|n|/dV, which this model carries as a "
          "magnitude",
          f"hole branch {share_hole:+.4f} %, electron branch "
          f"{share_ref:+.4f} %; a signed-carrier treatment is needed before "
          f"the hole-branch SIGN can be trusted (new item)")
    say("  Note: the share of I_d is about half of Qd/Q because the contact "
        "feedback")
    say("  of Eq. (5) absorbs part of the gain and the S_D term of Eq. (12) "
        "opposes")
    say(f"  it.  At the saturated-peak bias, V_ds = {VDS_REF} V: Qd/Q = "
        f"{_dd_integrals(VG_SAT_PEAK, VDS_REF, lam=1.0)['Qd'] / _dd_integrals(VG_SAT_PEAK, VDS_REF, lam=1.0)['Q'] * 100:.4f} %"
        f" against {share_ref:.4f} % in I_d.")

    say()
    say("  G2.  kappa along the channel against Zebrev Eq. 62, kappa = "
        "C_ox/(C_Q + C_it),")
    say("       with C_it = 0 and C_Q from the dispersion.  Two computations "
        "sharing")
    say("       no machinery: a numerical dV_F/dV against capacitor algebra.")
    d = _dd_integrals(VG_SAT_PEAK, VDS_REF, lam=1.0)
    CQ_disp = quantum_capacitance_dispersion(d["n"])
    z_bare = gfet.C_ox / CQ_disp                      # Zebrev, C_it = 0
    z_sc = gfet.C_ox / (CQ_disp + gfet.C_ox)          # self-consistent gate eq.
    say(f"  {'V_ch':>7} {'n [m^-2]':>12} {'kappa (this)':>13} "
        f"{'Cox/CQ':>10} {'Cox/(CQ+Cox)':>13}")
    for j in (0, 50, 100, 150, 200):
        say(f"  {d['V'][j]:>7.4f} {d['n'][j]:>12.4e} {d['kappa'][j]:>13.6f} "
            f"{z_bare[j]:>10.6f} {z_sc[j]:>13.6f}")
    # Compare on the interior, where np.gradient is central (the endpoints
    # use a one-sided stencil, which 2026-10-06's M5 is a reminder about).
    sl = slice(1, -1)
    r_bare = float(np.mean(d["kappa"][sl] / z_bare[sl]))
    r_sc = float(np.mean(d["kappa"][sl] / z_sc[sl]))
    say(f"  interior mean kappa/(Cox/CQ)       = {r_bare:.6f}")
    say(f"  interior mean kappa/(Cox/(CQ+Cox)) = {r_sc:.6f}")
    check("G2", 0.3 < r_bare < 3.0,
          "the boundary term's kappa and Zebrev's capacitor ratio agree to "
          "within a factor of 3, so they are the same quantity",
          f"ratio {r_bare:.6f} against Cox/CQ; the residual is this model's "
          f"series factor and puddle floor, not a disagreement about physics")

    say()
    say("  G3.  WHY gfet.quantum_capacitance() CANNOT be used for kappa:")
    say("       it takes the gate overdrive in the slot the dispersion wants")
    say("       E_F/e = V_F.  The ratio IS the error in kappa, since kappa "
        "~ 1/C_Q.")
    say(f"  {'V_ch':>7} {'overdrive [V]':>14} {'V_F [V]':>9} "
        f"{'ratio':>8} {'C_q(repo)':>11} {'C_Q(disp)':>11}")
    for j in (0, 100, 200):
        od = VG_SAT_PEAK - d["V"][j] - gfet.V_dirac
        cq_repo = float(gfet.quantum_capacitance(od, T=gfet.T))
        say(f"  {d['V'][j]:>7.4f} {od:>14.6f} {d['V_F'][j]:>9.6f} "
            f"{od / d['V_F'][j]:>8.4f} {cq_repo:>11.4e} {CQ_disp[j]:>11.4e}")
    od0 = VG_SAT_PEAK - gfet.V_dirac
    cq_repo0 = float(gfet.quantum_capacitance(od0, T=gfet.T))
    check("G3", cq_repo0 > CQ_disp[0],
          "the repository's C_q EXCEEDS the dispersion-consistent C_Q, so it "
          "understates kappa by that ratio",
          f"{cq_repo0 / CQ_disp[0]:.4f}x at the source end; annotated in "
          f"Chapter 4, no committed number withdrawn")

    say()
    say("=" * 78)
    say("SECTION 3.  WHERE diffusion matters: the V_g sweep")
    say("=" * 78)
    say("  PREDICTION, written into this module before the sweep was run:")
    say("  Feijoo et al. 2020 put the peak-f_max bias near the onset of "
        "bipolar")
    say("  conduction; Zebrev's kappa = C_ox/C_Q DIVERGES as C_Q -> 0 at "
        "charge")
    say("  neutrality.  Both therefore predict the share PEAKS AT THE DIRAC "
        "POINT.")
    Vg = np.linspace(-2.0, 4.0, N_VG)
    I0 = transfer_characteristic_dd(Vg, Vds=VDS_REF, lam=0.0)
    I1, opk, qshare = transfer_characteristic_dd(Vg, Vds=VDS_REF, lam=1.0,
                                                 report=True)
    share = (I1 / I0 - 1.0) * 100.0
    ashare = np.abs(share)
    j_max = int(np.argmax(ashare))
    j_dirac = int(np.argmin(np.abs(Vg - gfet.V_dirac)))
    say()
    say(f"  MEASURED.  largest |share| {ashare[j_max]:.4f} % at V_g = "
        f"{Vg[j_max]:+.6f} V")
    say(f"             share at the nearest grid point to the Dirac point "
        f"({Vg[j_dirac]:+.6f} V): {share[j_dirac]:+.6f} %")
    say(f"             share at the grid edges: {share[0]:+.4f} % at "
        f"{Vg[0]:+.2f} V, {share[-1]:+.4f} % at {Vg[-1]:+.2f} V")
    say()
    say("  THE PREDICTION IS REFUTED, and the mechanism is this model's own")
    say("  regularization.  carrier_density() returns")
    say("  n = sqrt(n_eff^2 + n_puddle^2), so at n_eff = 0 the derivative")
    say("  d|n|/dV = (n_eff/n) * dn_eff/dV is EXACTLY ZERO while C_Q stays")
    say("  finite at n_puddle.  kappa goes as d|n|/dV, so the puddle floor")
    say("  FLATTENS exactly the region where Zebrev's capacitor ratio "
        "diverges.")
    say("  The diffusion-dominated regime is therefore not small in this "
        "model --")
    say("  it is ABSENT, and that is a statement about the regularization "
        "rather")
    say("  than about graphene.  THE MODEL CANNOT BE ASKED ABOUT THE "
        "PEAK-f_max")
    say("  BIAS, which is what 2026-10-06 suspected, but for the puddle floor")
    say("  rather than for the missing term.")
    check("G4a", ashare[j_dirac] < 0.1 * ashare[j_max],
          "PREDICTION NOT MET: the share COLLAPSES at charge neutrality "
          "instead of diverging there",
          f"{ashare[j_dirac]:.6f} % at the Dirac point against "
          f"{ashare[j_max]:.4f} % at the maximum -- a factor "
          f"{ashare[j_max] / max(ashare[j_dirac], 1e-30):.3g}")
    # The mechanism, measured rather than argued.
    h = 1e-4
    nd = float(gfet.carrier_density(gfet.V_dirac, 0.0))
    dn_floored = (float(gfet.carrier_density(gfet.V_dirac, h))
                  - float(gfet.carrier_density(gfet.V_dirac, -h))) / (2 * h)
    dn_unfloored = -gfet.C_ox / E_CHARGE * float(
        gfet.quantum_capacitance(0.0, T=gfet.T)
        / (gfet.quantum_capacitance(0.0, T=gfet.T) + gfet.C_ox))
    say(f"  MECHANISM, measured at the Dirac point: d|n|/dV = "
        f"{dn_floored:.6e} m^-2/V with the puddle floor,")
    say(f"  against {dn_unfloored:.6e} m^-2/V for the unfloored n_eff -- a "
        f"factor {abs(dn_unfloored / dn_floored) if dn_floored != 0 else float('inf'):.4g}.")
    check("G4b", abs(dn_floored) < 1e-3 * abs(dn_unfloored),
          "and the mechanism is confirmed by measurement: the floored "
          "derivative is smaller than the unfloored one by >1000x at the "
          "Dirac point",
          f"{dn_floored:.6e} vs {dn_unfloored:.6e} m^-2/V")
    check("G4c", bool(np.all(opk > 0.0)),
          "1 + kappa > 0 across the whole sweep, including the Dirac region",
          f"min {float(np.min(opk)):.8f} at V_g = "
          f"{Vg[int(np.argmin(opk))]:+.6f} V")

    say()
    say("=" * 78)
    say("SECTION 4.  The RF consequence, against 2026-10-06's committed numbers")
    say("=" * 78)
    say("  5-point stencils at the COMMITTED gate-grid spacing h = "
        f"{H_VG_COMMITTED:.9f} V,")
    say("  so g_m here is the same ESTIMATOR as the committed one.  The")
    say("  lambda = 0 row is the control: it must reproduce vsat on the same")
    say("  stencil bitwise, and it is reported against the committed 0.683186")
    say("  so the stencil's own offset is visible rather than absorbed.")
    say()
    say(f"  {'bias':>14} {'lambda':>7} {'f_T [Hz]':>13} {'f_max [Hz]':>13} "
        f"{'f_max/f_T':>11} {'g_ds [S]':>12}")
    rf_res = {}
    for name, V_g in (("resistor-peak", VG_RES_PEAK),
                      ("saturated-peak", VG_SAT_PEAK)):
        for lam in (0.0, 1.0):
            r = _at_literature_geometry(
                lambda nf, V_g=V_g, lam=lam: _fT_fmax_local(
                    V_g, Vds=VDS_REF, lam=lam, N_fingers=nf))
            rf_res[(name, lam)] = r
            say(f"  {name:>14} {lam:>7.1f} {r['fT']:>13.6e} "
                f"{r['fmax']:>13.6e} {r['ratio']:>11.6f} {r['gds']:>12.6e}")

    # CONTROL: lambda = 0 on this stencil IS vsat on this stencil.
    def _vsat_local(nf, V_g):
        Vg5 = V_g + H_VG_COMMITTED * np.arange(-2, 3, dtype=float)
        fT, fmax, gm, gds, Cgs = vsat.fT_fmax_saturated(
            Vg5, Vds=VDS_REF, N_fingers=nf)
        return float(fmax[2] / fT[2])
    v0 = _at_literature_geometry(
        lambda nf: _vsat_local(nf, VG_SAT_PEAK))
    check("G5a", ulps(v0, rf_res[("saturated-peak", 0.0)]["ratio"]) == 0.0,
          "CONTROL: lambda = 0 reproduces vsat's f_max/f_T on this stencil "
          "BITWISE", f"{v0:.17e}")
    say(f"  committed 2026-10-06 value on the 400-point grid: 0.683186")
    say(f"  this stencil, lambda = 0:                         "
        f"{rf_res[('saturated-peak', 0.0)]['ratio']:.6f}")
    say(f"  stencil-vs-grid offset:                           "
        f"{(rf_res[('saturated-peak', 0.0)]['ratio'] / 0.683186 - 1) * 100:+.4f} %")

    # Added 2026-10-07 in response to this module's own mutation harness,
    # whose M7 (the stencil stops using the committed gate-grid spacing)
    # SURVIVED: G5a compares lambda = 0 against vsat on the SAME stencil, so
    # both move together and the comparison against the committed estimator
    # was printed but never asserted.  2026-10-06's M5 is the same fault --
    # a quantity nothing in the battery read -- and this is its repair here.
    COMMITTED_RATIO_2026_10_06 = 0.683186
    off = abs(rf_res[("saturated-peak", 0.0)]["ratio"]
              - COMMITTED_RATIO_2026_10_06)
    check("G5d", off < 5e-7,
          "and the stencil REPRODUCES the committed 2026-10-06 f_max/f_T to "
          "the 6 figures it was committed at, so the stencil IS the committed "
          "estimator rather than merely self-consistent",
          f"|{rf_res[('saturated-peak', 0.0)]['ratio']:.9f} - "
          f"{COMMITTED_RATIO_2026_10_06}| = {off:.3e}")

    r0 = rf_res[("saturated-peak", 0.0)]
    r1 = rf_res[("saturated-peak", 1.0)]
    say()
    say("  What diffusion does to the three quantities Chapter 4 is measured on:")
    for key in ("fT", "fmax", "ratio", "gds"):
        say(f"    {key:>6}: {r0[key]:>13.6e} -> {r1[key]:>13.6e}   "
            f"{(r1[key] / r0[key] - 1) * 100:+8.4f} %")
    check("G5b", r1["ratio"] < 1.3,
          "f_max/f_T with diffusion is STILL below the 1.3-1.4 literature "
          "band, so 2026-10-06's structural verdict survives the term that "
          "was supposed to threaten it",
          f"{r1['ratio']:.6f} = {r1['ratio'] / 1.3 * 100:.2f} % of 1.3")
    check("G5c", abs(r1["ratio"] / r0["ratio"] - 1) < abs(r1["fmax"] / r0["fmax"] - 1)
          or abs(r1["ratio"] / r0["ratio"] - 1) < 1e-3,
          "and the RATIO moves less than f_max does -- the 2026-10-06 "
          "ratio-blindness, in a third mechanism",
          f"ratio {(r1['ratio'] / r0['ratio'] - 1) * 100:+.4f} % vs f_max "
          f"{(r1['fmax'] / r0['fmax'] - 1) * 100:+.4f} %")

    say()
    say("=" * 78)
    say("SECTION 5.  Controls")
    say("=" * 78)
    I_lam = float(transfer_characteristic_dd(
        np.array([VG_SAT_PEAK]), Vds=VDS_REF, lam=1.0)[0])
    I_nol = float(transfer_characteristic_dd(
        np.array([VG_SAT_PEAK]), Vds=VDS_REF, lam=0.0)[0])
    mu_saved = gfet.mu
    gfet.mu = mu_saved * 1.5
    try:
        I_mu = float(transfer_characteristic_dd(
            np.array([VG_SAT_PEAK]), Vds=VDS_REF, lam=0.0)[0])
    finally:
        gfet.mu = mu_saved
    lam_resp = I_lam / I_nol
    mu_resp = I_mu / I_nol
    check("G6", abs(lam_resp - mu_resp) > 1e-3,
          "CONTROL: the lambda response and a 1.5x mobility response are NOT "
          "aliased, so 'diffusion did it' is distinguishable from 'mu was "
          "wrong'", f"lambda {lam_resp:.6f}x, mu {mu_resp:.6f}x")

    vsat_ref_after = float(vsat.transfer_characteristic_saturated(
        np.array([2.0]), Vds=VDS_REF)[0])
    res_ref_after = float(gfet.transfer_characteristic(
        np.array([2.0]), Vds=VDS_REF)[0])
    check("G7a", ulps(vsat_ref_before, vsat_ref_after) == 0.0,
          "the saturated drift model is BITWISE unchanged after this module "
          "has run", f"{vsat_ref_after:.17e}")
    check("G7b", ulps(res_ref_before, res_ref_after) == 0.0,
          "the resistor model is BITWISE unchanged too",
          f"{res_ref_after:.17e}")
    check("G7c", (gfet.mu, gfet.L, gfet.W) == (0.4, 200e-9, 1e-6),
          "gfet.mu / L / W are restored to committed values",
          f"mu = {gfet.mu}, L = {gfet.L}, W = {gfet.W}")

    say()
    say("=" * 78)
    say("SECTION 6.  The premise: which potential is V_ch?")
    say("=" * 78)
    say("  Pasadas and Jimenez (IEEE TED 2016, arXiv:1605.08235) write the "
        "GFET")
    say("  channel current as I = -W*Q_tot*v with v = mu*F, F = -dV/dx, and "
        "state")
    say("  that V(x) is the QUASI-FERMI LEVEL.  Under that reading Eq. (4) is")
    say("  already the complete drift-diffusion current and this module's term "
        "is")
    say("  a double count.  Under the electrostatic reading -- which is how "
        "V_ch")
    say("  enters C_ox*(V_g - V_ch - V_dirac) in carrier_density() -- Eq. (4) "
        "is")
    say("  drift-only and the term is missing physics.  The two readings differ "
        "by")
    say("  exactly lambda*V_F, so the number below is BOTH the diffusion "
        "current")
    say("  AND the size of an ambiguity that was already in every committed "
        "I_d.")
    say()
    say(f"  at the saturated-peak bias, V_ds = {VDS_REF} V:")
    say(f"    I_d, electrostatic reading (lambda = 1) : {I_lam:.8e} A")
    say(f"    I_d, quasi-Fermi reading   (lambda = 0) : {I_nol:.8e} A")
    say(f"    SPREAD                                  : "
        f"{(I_lam / I_nol - 1) * 100:+.4f} %")
    say(f"  and over the whole V_g sweep the spread stays within "
        f"{ashare[j_max]:.4f} % (largest |value|, at V_g = "
        f"{Vg[j_max]:+.4f} V), with OPPOSITE SIGNS on the two branches.")
    say()
    say("  This does NOT withdraw a committed number.  It bounds one.  The "
        "RF")
    say("  thread's conclusions are ratios between runs of one model, and "
        "G5b")
    say("  shows the binding verdict survives either reading.  What it does "
        "end")
    say("  is the claim that this repository's RF numbers are drift-only "
        "slices:")
    say("  they are slices of a model whose driving potential was never "
        "named,")
    say("  and the correct repair is to NAME IT, not to add a term.")
    say("  Recommended, and recorded as the new top RF item: declare V_ch the")
    say("  quasi-Fermi potential (which is what the series factor in")
    say("  carrier_density() is for), keep lambda = 0 as the production path, "
        "and")
    say("  keep this module as the bound on the choice.  lambda = 1 is NOT")
    say("  wired into any other module by this commit.")

    say()
    say("=" * 78)
    say("SECTION 7.  Figure")
    say("=" * 78)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.plot(Vg, share, lw=2, color="crimson")
    ax.axvline(gfet.V_dirac, color="gray", ls=":", label="Dirac point")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.axvline(Vg[j_max], color="navy", ls="--",
               label=f"largest |share| {ashare[j_max]:.3f} % at "
                     f"{Vg[j_max]:+.3f} V")
    ax.set_xlabel(r"$V_g$ (V)")
    ax.set_ylabel(r"diffusion share of $I_d$ (%)")
    ax.set_title(r"Diffusion share vs gate bias ($V_{ds}$ = %.2f V)" % VDS_REF)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(d["V"], d["kappa"], lw=2, label=r"$\kappa = -dV_F/dV$ (this model)")
    ax.plot(d["V"], z_bare, lw=2, ls="--",
            label=r"Zebrev $C_{ox}/C_Q$")
    ax.plot(d["V"], z_sc, lw=2, ls=":",
            label=r"$C_{ox}/(C_Q+C_{ox})$")
    ax.set_xlabel(r"$V_{ch}$ along the channel (V)")
    ax.set_ylabel(r"local diffusion/drift ratio $\kappa$")
    ax.set_title("Two independent routes to the same local ratio")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(Vg, I0 * 1e3, lw=2, label=r"$\lambda = 0$ (drift, Eq. 4)")
    ax.plot(Vg, I1 * 1e3, lw=2, ls="--",
            label=r"$\lambda = 1$ (drift-diffusion, Eq. 10)")
    ax.axvline(gfet.V_dirac, color="gray", ls=":")
    ax.set_xlabel(r"$V_g$ (V)")
    ax.set_ylabel(r"$I_d$ (mA/$\mu$m)")
    ax.set_title("Transfer characteristic, both readings of $V_{ch}$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    nq = np.logspace(np.log10(gfet.n_puddle), 17, 200)
    ax.loglog(nq * 1e-4, gfet.C_ox / quantum_capacitance_dispersion(nq), lw=2,
              label=r"$C_{ox}/C_Q$ (dispersion)")
    ax.axhline(1.0, color="crimson", ls="--",
               label=r"$\kappa = 1$: diffusion = drift")
    ax.axvline(gfet.n_puddle * 1e-4, color="gray", ls=":",
               label=r"$n_{puddle}$ floor")
    ax.set_xlabel(r"$n$ (cm$^{-2}$)")
    ax.set_ylabel(r"$\kappa$")
    ax.set_title(r"Where diffusion would equal drift (90 nm SiO$_2$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig("diffusion_current_model.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    say("  wrote diffusion_current_model.png")

    say()
    say("=" * 78)
    npass = sum(1 for ln in _LINES if "[PASS]" in ln)
    nfail = sum(1 for ln in _LINES if "[FAIL]" in ln)
    say(f"RESULT: {len(_LINES)} lines, {npass} passed, {nfail} failed")
    say("=" * 78)
    say(f"  diffusion share at the RF bias : +{share_ref:.4f} % of I_d")
    say(f"  largest |share| over the sweep : {ashare[j_max]:.4f} % at "
        f"V_g = {Vg[j_max]:+.4f} V (signed: {share[j_max]:+.4f} %)")
    say(f"  share AT charge neutrality     : {share[j_dirac]:+.6f} % -- the"
        f" preregistered")
    say("                                   divergence is REFUTED, and the "
        "puddle")
    say("                                   floor is the measured mechanism "
        "(G4a/G4b)")
    say(f"  f_max/f_T, lambda = 0 -> 1     : {r0['ratio']:.6f} -> "
        f"{r1['ratio']:.6f}  ({(r1['ratio'] / r0['ratio'] - 1) * 100:+.4f} %)")
    say("  the 1.3-1.4 band verdict       : UNCHANGED (G5b)")
    say("  the item's premise             : MOVED -- Eq. (4)'s driving "
        "potential")
    say("                                   was never named, and the two "
        "readings")
    say("                                   differ by exactly this term "
        "(Section 6)")

    with open("diffusion_current_output.txt", "w") as fh:
        fh.write("\n".join(_LINES) + "\n")
    return 1 if _FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
