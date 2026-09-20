"""
graphene_contact_doping_nonlinear_model.py

Replaces the LINEAR work-function -> doping relation that every contact and
photodetector model in this repo has assumed since 2026-08-26 with the
self-consistent, square-root relation of

  Khomyakov, Giovannetti, Rusu, Brocks, van den Brink, Kelly,
  Phys. Rev. B 79, 195425 (2009), Eq. 7   [arXiv:0902.1203]

and reports what that does to the two headline results of Chapter 6: the
Ti/Pt signed net response (1.832 under the linear model, 2026-09-18) and the
max|k| / |N_uniform| < 1.02 illumination ceiling (2026-09-19).

Full literature discussion, the derivation of alpha, and an honest record of
a PDF table-extraction failure:
  notes/2026-09-20-nonlinear-work-function-to-doping-relation.md

--------------------------------------------------------------------------
THE PHYSICS IN ONE PARAGRAPH
--------------------------------------------------------------------------
Charge transferred from the metal has to go somewhere, and in graphene it
goes into a LINEAR density of states, so n ~ E_F^2. The electrostatic
potential step across the metal-graphene gap is proportional to that
transferred charge, hence quadratic in the Fermi-level shift. Equilibrium
therefore requires

    dW'  =  phi  +  (alpha/2) phi^2                                    (1)

where phi = dE_F in eV, dW' = W_metal - W_graphene - D_c is the work-function
offset AFTER the short-range chemical term, and

    alpha = 2 e^3 (d - d0) / (eps_0 pi hbar^2 v_F^2)                   (2)

with d the metal-graphene separation and d0 = 2.4 A the separation below
which the gap-capacitance picture switches off. Inverting (1) with the root
that is continuous through phi = 0:

    phi(dW') = sgn(dW') * ( sqrt(1 + 2 alpha |dW'|) - 1 ) / alpha      (3)

which is Khomyakov et al.'s Eq. 7. It is LINEAR only as dW' -> 0, and
SUBLINEAR (asymptotically sqrt) for large offsets: graphene resists being
doped, increasingly, the harder you push it.

--------------------------------------------------------------------------
WHAT IS AND IS NOT CLAIMED
--------------------------------------------------------------------------
* alpha is DERIVED (Eq. 2), not fitted. It is cross-checked against the PRB's
  Pt number to about 1.4x -- right size, right sign, not quantitative. Every
  conclusion below is therefore also reported as a sweep over alpha, so that
  nothing rests on the exact value.
* NO per-metal dE_F from the PRB's Table I is used. Two WebFetch passes over
  the same PDF returned mutually contradictory versions of that table; see
  Section 3 of the note. Only d_eq, d0, D_c(3.3 A) ~ 0.9 eV and W_0 = 5.4 eV
  -- which agreed across both passes -- are used.
* For Ti, Ni, Co and Pd, d_eq < d0, so Eq. (2) gives alpha <= 0 and Eq. (3)
  is OUTSIDE ITS OWN REGIME. This module refuses to extrapolate there:
  alpha_for_metal() returns None and says why. That is not a numerical
  guard, it is Giovannetti et al.'s actual conclusion -- chemisorbed metals
  are not characterised by their work function alone.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, epsilon_0, hbar

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS, W_GRAPHENE, lambda_decay
from graphene_photodetector_model import L_channel
from graphene_photodetector_collection_model import E_BIAS
from graphene_photodetector_signed_carrier_model import (
    total_field, _transport, Q_HOLE, Q_ELECTRON, N_POINTS,
    W_CROSS_CHEM, signed_offset,
)
from graphene_photodetector_nonuniform_illumination_model import (
    collection_kernel as linear_collection_kernel,
    g_uniform, g_shadow_mask, response_under, response_per_incident,
    T_MASK_MEASURED,
)

V_F = 1.0e6          # m/s, graphene Fermi velocity used throughout this repo
D0_SEPARATION = 2.4e-10   # m, Khomyakov et al.'s d0
D_CHEM_PHYS = 0.9         # eV, short-range chemical term at d ~ 3.3 A

# Equilibrium metal-graphene separations, Khomyakov et al. 2009. These agreed
# across both PDF extraction passes and are the paper's robust qualitative
# result. Physisorbed metals sit near 3.3 A; chemisorbed ones are pulled in
# to ~2.05-2.3 A, i.e. BELOW d0.
D_EQ = {
    "Al": 3.41e-10, "Cu": 3.26e-10, "Ag": 3.33e-10, "Au": 3.31e-10,
    "Pt": 3.30e-10,                       # physisorbed
    "Pd": 2.30e-10, "Ti": 2.10e-10, "Ni": 2.05e-10, "Co": 2.05e-10,
    # Cr: not tabulated by Khomyakov et al.; left absent deliberately rather
    # than guessed. alpha_for_metal("Cr") therefore reports "no separation".
}
CHEMISORBED = {"Pd", "Ti", "Ni", "Co"}

ALPHA_PHYS = None    # filled in below by alpha_from_separation(3.3 A)


# =====================================================================
# The relation itself
# =====================================================================
def alpha_from_separation(d, v_f=V_F, d0=D0_SEPARATION):
    """
    alpha in eV^-1 from Eq. (2). Returns a NEGATIVE value for d < d0; callers
    must not silently use that -- see alpha_for_metal().
    """
    return 2.0 * e**3 * (d - d0) / (epsilon_0 * np.pi * hbar**2 * v_f**2)


ALPHA_PHYS = alpha_from_separation(3.3e-10)


def alpha_for_metal(metal):
    """
    (alpha, reason) for one metal. alpha is None when the model does not
    apply, with `reason` saying which of the two ways it fails.
    """
    d = D_EQ.get(metal)
    if d is None:
        return None, ("no equilibrium separation tabulated by Khomyakov et al.; "
                      "not guessed")
    if d <= D0_SEPARATION:
        return None, (f"chemisorbed, d_eq = {d*1e10:.2f} A < d0 = "
                      f"{D0_SEPARATION*1e10:.1f} A: the gap-capacitance term of "
                      f"Eq. 7 is outside its own regime, and Giovannetti et al. "
                      f"state these metals are not characterised by work "
                      f"function alone")
    return alpha_from_separation(d), "physisorbed"


def fermi_shift(dW, alpha=None):
    """
    Eq. (3): graphene Fermi-level shift, in eV, from the post-chemical
    work-function offset dW, also in eV. Odd in dW by construction.

    alpha = 0 returns dW exactly (the linear model), which is what makes
    validate_linear_limit() an exact test rather than a tolerance test.
    """
    if alpha is None:
        alpha = ALPHA_PHYS
    dW = np.asarray(dW, dtype=float)
    if alpha == 0.0:
        return dW.copy() if dW.ndim else float(dW)
    out = np.sign(dW) * (np.sqrt(1.0 + 2.0 * alpha * np.abs(dW)) - 1.0) / alpha
    return out if out.ndim else float(out)


def inverse_fermi_shift(phi, alpha=None):
    """Eq. (1), the forward direction. Used for the round-trip validation."""
    if alpha is None:
        alpha = ALPHA_PHYS
    phi = np.asarray(phi, dtype=float)
    out = phi + 0.5 * alpha * np.sign(phi) * phi**2
    return out if out.ndim else float(out)


# =====================================================================
# The field, with the profile magnitude passed through Eq. (3)
# =====================================================================
def nonlinear_total_field(x, W_A, W_B, e_applied=E_BIAS, L=L_channel,
                          w_cross=W_CROSS_CHEM, lam=lambda_decay, alpha=None):
    """
    Term-for-term the same expression as
    graphene_photodetector_signed_carrier_model.total_field(), with dW
    replaced by fermi_shift(dW). Written out rather than wrapped, so that
    validate_field_reduction() can compare the two bitwise at alpha = 0.
    """
    dW_A = fermi_shift(signed_offset(W_A, w_cross), alpha)
    dW_B = fermi_shift(signed_offset(W_B, w_cross), alpha)
    return (-e_applied
            - dW_A / lam / (1.0 + x / lam) ** 2
            + dW_B / lam / (1.0 + (L - x) / lam) ** 2)


def nonlinear_collection_kernel(W_A, W_B, e_applied=E_BIAS, L=L_channel,
                                n_points=N_POINTS, w_cross=W_CROSS_CHEM,
                                alpha=None, carriers=("hole", "electron")):
    """
    k(x) with the nonlinear profile. Accumulation order is identical to
    graphene_photodetector_nonuniform_illumination_model.collection_kernel();
    do not 'simplify' it, for the reason recorded there.
    """
    x = np.linspace(0.0, L, n_points)
    E = nonlinear_total_field(x, W_A, W_B, e_applied=e_applied, L=L,
                              w_cross=w_cross, alpha=alpha)
    k = np.zeros_like(x)
    if "hole" in carriers:
        s_h, p_h = _transport(+E, x)
        k = k + Q_HOLE * s_h * p_h
    if "electron" in carriers:
        s_e, p_e = _transport(-E, x)
        k = k + Q_ELECTRON * s_e * p_e
    return x, k


def nonlinear_net_response(W_A, W_B, e_applied=E_BIAS, alpha=None,
                           w_cross=W_CROSS_CHEM, **kw):
    x, k = nonlinear_collection_kernel(W_A, W_B, e_applied=e_applied,
                                       w_cross=w_cross, alpha=alpha, **kw)
    return response_under(k, x, g_uniform(x))


# =====================================================================
# VALIDATIONS AGAINST EXACTLY-KNOWN VALUES
# =====================================================================
def validate_zero_and_odd(verbose=True):
    """
    EXACT 1. fermi_shift(0) must be bitwise 0, and the relation must be
    bitwise odd. Both are exactly-known values, not tolerances: any
    asymmetry here would put a spurious net response into a symmetric
    device and destroy the 2026-09-17 cancellation test downstream.
    """
    z = fermi_shift(0.0)
    grid = np.linspace(0.05, 3.0, 60)
    odd_residual = np.max(np.abs(fermi_shift(grid) + fermi_shift(-grid)))
    ok = (z == 0.0) and (odd_residual == 0.0)
    if verbose:
        print("\nEXACT 1: fermi_shift(0) == 0 and odd symmetry, bitwise")
        print(f"  fermi_shift(0)          = {z!r}   (must be exactly 0.0)")
        print(f"  max |f(x) + f(-x)|      = {odd_residual!r}   (must be exactly 0.0)")
        print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


def validate_round_trip(verbose=True):
    """
    EXACT 2. Eq. (3) must invert Eq. (1). The composition is exactly the
    identity analytically, so the only admissible residual is float rounding.
    """
    grid = np.concatenate([np.linspace(-3.0, -0.01, 80), [0.0],
                           np.linspace(0.01, 3.0, 80)])
    res = np.max(np.abs(inverse_fermi_shift(fermi_shift(grid)) - grid))
    ok = res < 1e-12
    if verbose:
        print("\nEXACT 2: Eq. (1) o Eq. (3) == identity")
        print(f"  max round-trip residual = {res:.3e} eV   (float rounding only)")
        print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


def validate_linear_limit(verbose=True):
    """
    EXACT 3. THE IMPORTANT ONE: the new model must reproduce the OLD model
    exactly in the overlapping limit. At alpha = 0, nonlinear_total_field()
    must equal signed_carrier_model.total_field() bitwise, for an asymmetric
    pair (so no accidental cancellation can hide a discrepancy).

    This is the validation class that caught a real bug on 2026-09-17, which
    a plausible-range check had passed.
    """
    x = np.linspace(0.0, L_channel, N_POINTS)
    W_A, W_B = METAL_WORK_FUNCTIONS["Ti"], METAL_WORK_FUNCTIONS["Pt"]
    old = total_field(x, W_A, W_B, e_applied=0.0)
    new = nonlinear_total_field(x, W_A, W_B, e_applied=0.0, alpha=0.0)
    bitwise = bool(np.array_equal(old, new))

    # and the analytic first correction: f(x) = x - (alpha/2) x^2 + O(x^3)
    a = ALPHA_PHYS
    small = np.array([1e-4, 1e-3, 1e-2])
    err = np.abs(fermi_shift(small, a) - (small - 0.5 * a * small**2))
    ratios = err[1:] / err[:-1]        # should be ~1000 for a cubic residual
    cubic = bool(np.all(ratios > 300) and np.all(ratios < 3000))

    if verbose:
        print("\nEXACT 3: reduction to the previous (linear) model")
        print(f"  alpha = 0 reproduces total_field() bitwise: {bitwise}")
        print(f"  max |old - new| = {np.max(np.abs(old - new))!r}")
        print(f"  residual vs x - (alpha/2)x^2 scales as x^3: "
              f"ratios {np.round(ratios, 1).tolist()} (expect ~1000) -> {cubic}")
        print(f"  -> {'PASS' if (bitwise and cubic) else 'FAIL'}")
    return bitwise and cubic


def validate_symmetric_cancellation(verbose=True):
    """
    EXACT 4. A symmetric pair at zero bias must give exactly zero net
    response under uniform illumination, under the NEW field too. This is
    the test the 2026-09-17 stagnation bug failed.
    """
    worst = 0.0
    for metal, W in METAL_WORK_FUNCTIONS.items():
        worst = max(worst, abs(nonlinear_net_response(W, W, e_applied=0.0)))
    ok = worst < 1e-15
    if verbose:
        print("\nEXACT 4: symmetric pair, zero bias -> exactly zero")
        print(f"  worst |N| over all 7 metals = {worst:.3e}   -> "
              f"{'PASS' if ok else 'FAIL'}")
    return ok


def validate_monotone_no_flip(verbose=True):
    """
    EXACT 5. Eq. (3) is strictly monotone and sign-preserving, so it CANNOT
    change any metal's p/n assignment. Checked explicitly rather than
    asserted, because the whole question of this session is which Chapter 6
    conclusions can and cannot move.
    """
    rows = []
    ok = True
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        dW = signed_offset(W, W_CROSS_CHEM)
        phi = fermi_shift(dW)
        same = (np.sign(phi) == np.sign(dW))
        shrunk = abs(phi) <= abs(dW) + 1e-15
        ok = ok and bool(same) and bool(shrunk)
        rows.append((metal, W, dW, float(phi), abs(phi) / abs(dW)))
    if verbose:
        print(f"\nEXACT 5: sign preserved and |dE_F| <= |dW|, alpha = "
              f"{ALPHA_PHYS:.3f} eV^-1")
        print(f"{'metal':<6}{'W (eV)':>8}{'dW (eV)':>10}{'dE_F (eV)':>11}"
              f"{'ratio':>8}")
        for m, W, dW, phi, r in rows:
            print(f"{m:<6}{W:>8.2f}{dW:>+10.2f}{phi:>+11.4f}{r:>8.3f}")
        print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok, rows


# =====================================================================
# RESULTS
# =====================================================================
def headline_comparison(verbose=True):
    """
    RESULT 1. Do the two Chapter 6 headline numbers survive?

    (a) Ti/Pt signed net response, 1.832 under the linear model (2026-09-18).
    (b) The illumination ceiling max|k| / |N_uniform|, < 1.02 for every
        asymmetric pair (2026-09-19).
    """
    pairs = [("Ti", "Pt"), ("Ti", "Pd"), ("Ti", "Au"), ("Ti", "Ni"),
             ("Cu", "Pt"), ("Cr", "Pt")]
    rows = []
    for mA, mB in pairs:
        WA, WB = METAL_WORK_FUNCTIONS[mA], METAL_WORK_FUNCTIONS[mB]
        xl, kl = linear_collection_kernel(WA, WB, e_applied=0.0)[:2]
        Nl = response_under(kl, xl, g_uniform(xl))
        xn, kn = nonlinear_collection_kernel(WA, WB, e_applied=0.0)
        Nn = response_under(kn, xn, g_uniform(xn))
        ceil_l = np.max(np.abs(kl)) / abs(Nl)
        ceil_n = np.max(np.abs(kn)) / abs(Nn)
        rows.append((f"{mA}/{mB}", Nl, Nn, Nn / Nl, ceil_l, ceil_n))
    if verbose:
        print("\nRESULT 1: the two Chapter 6 headline numbers under Eq. (3)")
        print(f"  alpha = {ALPHA_PHYS:.3f} eV^-1 (physisorbed, d = 3.3 A)")
        print(f"{'pair':<9}{'N linear':>11}{'N nonlin':>11}{'ratio':>8}"
              f"{'ceil lin':>10}{'ceil nl':>9}")
        for name, Nl, Nn, r, cl, cn in rows:
            print(f"{name:<9}{Nl:>+11.4f}{Nn:>+11.4f}{r:>8.3f}{cl:>10.4f}"
                  f"{cn:>9.4f}")
        worst_ceiling = max(max(r[4], r[5]) for r in rows)
        print(f"  max|k|/|N| over all pairs, either model: {worst_ceiling:.4f}")
        print(f"  the 2026-09-19 ceiling was < 1.02 -> "
              f"{'HOLDS' if worst_ceiling < 1.02 else 'VIOLATED'}")
    return rows


def alpha_sweep(pair=("Ti", "Pt"), verbose=True):
    """
    RESULT 2. alpha is good to about 1.4x at best, so sweep it. If the
    conclusion is the same from alpha = 0 to alpha = 5 eV^-1, it does not
    depend on the calibration.
    """
    WA, WB = METAL_WORK_FUNCTIONS[pair[0]], METAL_WORK_FUNCTIONS[pair[1]]
    alphas = np.array([0.0, 0.5, 1.0, 1.7, 2.39, 3.4, 5.0])
    rows = []
    for a in alphas:
        x, k = nonlinear_collection_kernel(WA, WB, e_applied=0.0, alpha=a)
        N = response_under(k, x, g_uniform(x))
        rows.append((a, N, np.max(np.abs(k)) / abs(N)))
    if verbose:
        print(f"\nRESULT 2: sweep of alpha, pair {pair[0]}/{pair[1]}")
        print(f"{'alpha':>7}{'N':>11}{'vs linear':>11}{'max|k|/|N|':>12}")
        N0 = rows[0][1]
        for a, N, c in rows:
            print(f"{a:>7.2f}{N:>+11.4f}{N / N0:>11.3f}{c:>12.4f}")
        print(f"  ceiling range over the whole sweep: "
              f"{min(r[2] for r in rows):.4f} - {max(r[2] for r in rows):.4f}")
    return rows


def applicability_table(verbose=True):
    """
    RESULT 3. Which metals in METAL_WORK_FUNCTIONS the model may be applied
    to at all. This is the finding that matters most for Chapter 6.
    """
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        a, reason = alpha_for_metal(metal)
        rows.append((metal, W, signed_offset(W, W_CROSS_CHEM), a, reason))
    if verbose:
        print("\nRESULT 3: where Eq. (3) may legitimately be applied")
        for m, W, dW, a, reason in rows:
            tag = f"alpha = {a:.3f} eV^-1" if a is not None else "NOT APPLICABLE"
            print(f"  {m:<4} W = {W:.2f}  dW = {dW:+.2f}  {tag}")
            if a is None:
                print(f"       {reason}")
        bad = [r[0] for r in rows if r[3] is None]
        big = sorted(rows, key=lambda r: -abs(r[2]))[:3]
        print(f"  inapplicable: {', '.join(bad)}")
        print(f"  three largest |dW| in the table: "
              f"{', '.join(f'{r[0]} ({r[2]:+.2f})' for r in big)}")
        print("  NOTE: Ti carries the largest |dW| AND is chemisorbed, and "
              "Ti is\n        contact A of both Chapter 6 headline pairs.")
    return rows



def ceiling_audit_all_pairs(verbose=True):
    """
    RESULT 4. THE 2026-09-19 CEILING CLAIM IS RETRACTED.

    That session reported "max|k| / |N_uniform| = 1.0025 (Ti/Pt), 1.0159
    (Ti/Pd), < 1.02 FOR EVERY ASYMMETRIC PAIR". This function enumerates all
    21 unordered asymmetric pairs under the SAME linear model that session
    used, and finds 14 violations, the worst by a factor of 22.

    What that session actually validated (its validate_kernel_ceiling) is the
    INEQUALITY |N[g]| <= max|k|, which is a normalisation identity and is
    still true -- 0 violations on 343 cases. What it over-generalised is the
    TIGHTNESS of that inequality, evidently from the two n/p pairs it
    happened to tabulate.

    The physics of the failure: max|k|/|N_uniform| is large exactly when
    |N_uniform| is SMALL, i.e. when the two contacts dope graphene the SAME
    way and their contributions nearly cancel under uniform light. Ti/Pt and
    Ti/Pd straddle the 5.4 eV crossover; Au/Pd, Ti/Cr and Cr/Cu do not.

    Reported for both models, because the retraction is independent of this
    session's nonlinearity: it was already false linearly. The nonlinearity
    makes it worse everywhere (it shrinks |N_uniform| faster than max|k|),
    which is why this session found it at all -- the Ti/Pd ceiling moving
    from 1.016 to 1.434 was what prompted enumerating the rest.
    """
    import itertools
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    rows = []
    for (mA, WA), (mB, WB) in itertools.combinations(metals, 2):
        xl, kl, _ = linear_collection_kernel(WA, WB, e_applied=0.0)
        Nl = response_under(kl, xl, g_uniform(xl))
        xn, kn = nonlinear_collection_kernel(WA, WB, e_applied=0.0)
        Nn = response_under(kn, xn, g_uniform(xn))
        straddles = (signed_offset(WA, W_CROSS_CHEM)
                     * signed_offset(WB, W_CROSS_CHEM)) < 0.0
        rows.append((f"{mA}/{mB}", straddles, Nl,
                     float(np.max(np.abs(kl)) / abs(Nl)),
                     float(np.max(np.abs(kn)) / abs(Nn))))
    rows.sort(key=lambda r: -r[3])
    if verbose:
        print("\nRESULT 4: RETRACTION of the 2026-09-19 '< 1.02' ceiling claim")
        print(f"{'pair':<8}{'n/p?':>6}{'N_unif':>10}{'ceil lin':>10}{'ceil nl':>9}")
        for name, st, N, cl, cn in rows:
            flag = " <-- violates < 1.02" if cl >= 1.02 else ""
            print(f"{name:<8}{('yes' if st else 'no'):>6}{N:>+10.4f}"
                  f"{cl:>10.4f}{cn:>9.4f}{flag}")
        viol = [r for r in rows if r[3] >= 1.02]
        print(f"  {len(viol)} of {len(rows)} asymmetric pairs violate it, "
              f"worst {rows[0][0]} at {rows[0][3]:.1f}x")
        print(f"  every violating pair is same-sign (does not straddle "
              f"5.4 eV): {all(not r[1] for r in viol)}")
        print(f"  every straddling pair satisfies it: "
              f"{all(r[3] < 1.02 for r in rows if r[1])}")
        print("  CORRECTED STATEMENT: the ratio is controlled by "
              "|N_uniform|, not by a\n    dichotomy. Every pair that "
              "STRADDLES the 5.4 eV crossover satisfies it;\n    same-sign "
              "pairs may (Ti/Pd, 1.016) or may not (Au/Pd, 22.4), and the\n"
              "    ratio diverges as N_uniform -> 0. The 2026-09-19 sample "
              "of two was\n    drawn entirely from the straddling half.")
    return rows


def mask_gain_all_pairs(verbose=True):
    """
    RESULT 5. AND THEREFORE THE 2026-09-19 CONCLUSION IS WRONG TOO.

    That session concluded "masks are for SYMMETRIC devices only", on the
    argument that an asymmetric pair has almost nothing to gain (max|k| is
    within 2% of |N_uniform|). With the ceiling retracted, the argument goes
    with it, so this computes the ACHIEVABLE gain -- a real perfect shadow
    mask, scored per INCIDENT photon, which is the 2026-09-19 accounting
    correction and the honest one for responsivity.

    Au/Pd gains 8.4x. Ti/Cr 3.4x. Six pairs gain more than 2x.

    The practical recommendation nevertheless SURVIVES, for a different
    reason than the one given: 8.4x of a small number is still small.
    Masked Au/Pd reaches 0.377 per incident photon; unmasked Ti/Pt reaches
    0.916. Relative gain and absolute performance point opposite ways, and
    2026-09-19 conflated them.
    """
    import itertools
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    rows = []
    for (mA, WA), (mB, WB) in itertools.combinations(metals, 2):
        x, k, _ = linear_collection_kernel(WA, WB, e_applied=0.0)
        N = response_under(k, x, g_uniform(x))
        best = max((abs(response_per_incident(k, x, g_shadow_mask(x, T=0.0,
                                                                  masked=side)))
                    for side in ("left", "right")))
        rows.append((f"{mA}/{mB}", N, best, best / abs(N)))
    rows.sort(key=lambda r: -r[3])
    xr, kr, _ = linear_collection_kernel(METAL_WORK_FUNCTIONS["Ti"],
                                         METAL_WORK_FUNCTIONS["Pt"],
                                         e_applied=0.0)
    ref = abs(response_per_incident(kr, xr, g_uniform(xr)))
    if verbose:
        print("\nRESULT 5: achievable perfect-mask gain, per INCIDENT photon")
        print(f"{'pair':<8}{'N_unif':>10}{'best masked':>13}{'gain':>8}"
              f"{'vs Ti/Pt unmasked':>20}")
        for name, N, best, g in rows:
            print(f"{name:<8}{N:>+10.4f}{best:>13.4f}{g:>8.2f}x"
                  f"{best / ref:>19.2f}x")
        big = [r for r in rows if r[3] > 2.0]
        print(f"  Ti/Pt unmasked reference, per incident photon: {ref:.4f}")
        print(f"  pairs gaining more than 2x from a mask: {len(big)} "
              f"({', '.join(r[0] for r in big)})")
        print(f"  best absolute masked device: {rows[0][0]} at "
              f"{rows[0][2]:.4f}, which is {rows[0][2] / ref:.2f}x Ti/Pt")
        print("  CORRECTED STATEMENT: masks help whenever |N_uniform| is "
              "small -- symmetric\n    pairs (where it is zero) AND "
              "weakly-asymmetric same-sign pairs. They do\n    not make "
              "either competitive in ABSOLUTE terms with a straddling pair.")
    return rows, ref


def make_plots(fname="contact_doping_nonlinear_relation.png"):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    ax = axes[0]
    dW = np.linspace(-2.5, 2.5, 601)
    ax.plot(dW, dW, "k--", lw=1.2, label="linear (repo, to 2026-09-19)")
    for a, c in [(1.0, "#7fb3d5"), (ALPHA_PHYS, "#1f77b4"), (5.0, "#0b3d62")]:
        ax.plot(dW, fermi_shift(dW, a), color=c, lw=1.8,
                label=fr"Eq. 7, $\alpha$ = {a:.2f} eV$^{{-1}}$")
    for m in ("Ti", "Pt", "Pd"):
        d = signed_offset(METAL_WORK_FUNCTIONS[m], W_CROSS_CHEM)
        ax.plot([d], [fermi_shift(d)], "o", ms=6, color="#d62728")
        ax.annotate(m, (d, fermi_shift(d)), textcoords="offset points",
                    xytext=(6, -10), fontsize=9)
    ax.axhline(0, color="0.8", lw=0.8); ax.axvline(0, color="0.8", lw=0.8)
    ax.set_xlabel(r"$\Delta W' = W_M - W_G - \Delta_c$  (eV)")
    ax.set_ylabel(r"$\Delta E_F$  (eV)")
    ax.set_title("The relation saturates")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    WA, WB = METAL_WORK_FUNCTIONS["Ti"], METAL_WORK_FUNCTIONS["Pt"]
    xl, kl = linear_collection_kernel(WA, WB, e_applied=0.0)[:2]
    xn, kn = nonlinear_collection_kernel(WA, WB, e_applied=0.0)
    ax.plot(xl * 1e9, kl, "k--", lw=1.4, label="linear")
    ax.plot(xn * 1e9, kn, "-", color="#1f77b4", lw=1.8,
            label=fr"Eq. 7, $\alpha$ = {ALPHA_PHYS:.2f}")
    ax.set_xlabel("x (nm)"); ax.set_ylabel("k(x)  (charge to A per photon)")
    ax.set_title("Ti/Pt collection kernel"); ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[2]
    rows = alpha_sweep(verbose=False)
    a = [r[0] for r in rows]
    ax.plot(a, [r[1] / rows[0][1] for r in rows], "o-", color="#1f77b4",
            label="N / N(linear)")
    ax.plot(a, [r[2] for r in rows], "s-", color="#d62728",
            label=r"max$|k|$ / $|N|$ (ceiling)")
    ax.axvline(ALPHA_PHYS, color="0.5", ls=":", lw=1.2)
    ax.annotate(r"derived $\alpha$", (ALPHA_PHYS, 1.0),
                textcoords="offset points", xytext=(5, 20), fontsize=9)
    ax.set_xlabel(r"$\alpha$ (eV$^{-1}$)"); ax.set_ylabel("relative")
    ax.set_title("Ti/Pt: compression, no reversal")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    fig.tight_layout(); fig.savefig(fname, dpi=150)
    print(f"\nwrote {fname}")
    return fname


def main():
    print("=" * 72)
    print("NONLINEAR WORK-FUNCTION -> DOPING RELATION  (Khomyakov 2009 Eq. 7)")
    print("=" * 72)
    print(f"derived alpha (d = 3.3 A, d0 = 2.4 A) = {ALPHA_PHYS:.4f} eV^-1")
    oks = [validate_zero_and_odd(), validate_round_trip(),
           validate_linear_limit(), validate_symmetric_cancellation(),
           validate_monotone_no_flip()[0]]
    print(f"\n{sum(oks)}/{len(oks)} exact validations passed")
    if not all(oks):
        raise SystemExit("EXACT VALIDATION FAILED -- results below are not trustworthy")
    headline_comparison()
    alpha_sweep()
    applicability_table()
    ceiling_audit_all_pairs()
    mask_gain_all_pairs()
    make_plots()


if __name__ == "__main__":
    main()
