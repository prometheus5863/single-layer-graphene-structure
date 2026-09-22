"""
graphene_sensitivity_audit.py

A CONDITION-NUMBER AUDIT of the headline computed quantities of Chapter 4
(contact resistance) and Chapter 5 (interconnects), asking the question
created on 2026-09-21 and sharpened there:

    Are these quantities NEAR-CANCELLATIONS?

Background.  On 2026-09-21 the per-metal-crossover session measured a
65-fold difference in how two classes of quantity respond to the SAME
17.12% input perturbation: a same-sign (near-cancelling) contact pair moved
77.94% (amplification 4.55) while a straddling (reinforcing) pair moved
1.21% (amplification 0.07).  The difference is arithmetic, not physics, so
it applies to every quantity this thesis computes.  Chapters 4 and 5 had
never been asked which kind they are.

Pre-registered predictions P1-P6, committed to git BEFORE this file existed:
  notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md   (commit 968a9c0)

--------------------------------------------------------------------------
THE TWO DIAGNOSTICS
--------------------------------------------------------------------------
(1) TERM-LEVEL condition number.  Write Q as a sum of additive terms,
    Q = sum_i T_i.  Then the signed logarithmic sensitivity of Q to term
    T_i is exactly

        S_i = T_i / Q                                              (C1)

    and sum_i S_i == 1 EXACTLY, by construction.  kappa_term = max_i |S_i|.
    This is closed-form: no finite differences, no tolerance.

(2) PARAMETER-LEVEL sensitivity.  For an underlying physical parameter p,

        S_p = d ln|Q| / d ln p                                     (C2)

    by central difference in log space, with a Richardson-style two-step
    convergence check reported alongside every value.

CLASSIFICATION used throughout:
    kappa <= 1     reinforcing
    1 < kappa < 3  mildly ill-conditioned
    kappa >= 3     near-cancelling   (Chapter 6's same-sign pairs: 4.55)

--------------------------------------------------------------------------
WHAT IS AND IS NOT CLAIMED
--------------------------------------------------------------------------
* kappa is a property of a FORMULA, not of the physics.  A large kappa does
  not make a number wrong; it says the number inherits its inputs' errors
  multiplied by kappa.
* NO input value is changed anywhere in this file and nothing is
  recalibrated.  Superseded numbers are annotated in the thesis, never
  overwritten.
* Where no input uncertainty is available from the literature, the audit
  reports kappa and a HYPOTHETICAL 5% bar, labelled as hypothetical.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_interconnect_model as icm
from graphene_contact_doping_model import (
    METAL_WORK_FUNCTIONS, METAL_LITERATURE_RC, junction_extra_resistance,
)

NEAR_CANCELLING = 3.0
HYPOTHETICAL_INPUT_ERR = 0.05     # 5%, used only where literature gives none


# =====================================================================
# Diagnostic 1: term-level (closed form, exact)
# =====================================================================
def term_sensitivities(terms):
    """
    terms : dict name -> signed additive contribution T_i.
    Returns (S dict, Q, kappa).  sum(S.values()) == 1 exactly up to the
    rounding of the division itself -- this is P1's exact identity.
    """
    Q = float(sum(terms.values()))
    if Q == 0.0:
        # Exact total cancellation: every sensitivity is infinite and the
        # sum rule is undefined (0/0), which is the correct answer rather
        # than an error.  Validation 3(a) exercises this branch deliberately.
        S = {k: (np.inf if v != 0 else np.nan) for k, v in terms.items()}
        return S, Q, np.inf
    S = {k: float(v) / Q for k, v in terms.items()}
    return S, Q, max(abs(s) for s in S.values())


def classify(kappa):
    if kappa <= 1.0:
        return "reinforcing"
    if kappa < NEAR_CANCELLING:
        return "mildly ill-cond."
    return "NEAR-CANCELLING"


# =====================================================================
# Diagnostic 2: parameter-level (central difference in log space)
# =====================================================================
def log_sensitivity(f, p, rel_step=1e-5, check=True):
    """
    S = d ln|f| / d ln p by central difference.  Returns (S, conv) where
    conv is |S(h) - S(2h)|, a convergence estimate that is reported with
    every value rather than assumed small.
    """
    def _S(h):
        fp, fm = f(p * (1 + h)), f(p * (1 - h))
        return (np.log(abs(fp)) - np.log(abs(fm))) / (np.log(1 + h) - np.log(1 - h))
    S = _S(rel_step)
    conv = abs(S - _S(2 * rel_step)) if check else float("nan")
    return S, conv


# =====================================================================
# FAMILY A -- Chapter 4: R_transmission = Rc_measured - R_extra
# =====================================================================
N_BULK_DEFAULT = 2.0e16


def family_A_rows(n_bulk=N_BULK_DEFAULT):
    rows = []
    for metal, lit in METAL_LITERATURE_RC.items():
        if metal not in METAL_WORK_FUNCTIONS:
            continue
        Rc = lit["Rc_measured_ohm_um"]
        R_extra, _, _, _ = junction_extra_resistance(METAL_WORK_FUNCTIONS[metal], n_bulk)
        S, Q, k = term_sensitivities({"Rc_measured": +Rc, "R_extra": -R_extra})
        rows.append((metal, Rc, R_extra, Q, S, k))
    rows.sort(key=lambda r: -r[5])
    return rows


def family_A(verbose=True):
    rows = family_A_rows()
    if verbose:
        print("\n" + "=" * 78)
        print("FAMILY A -- Chapter 4, Section 4.7:  R_transmission = Rc_measured - R_extra")
        print("=" * 78)
        print(f"  {'metal':<6}{'Rc_meas':>10}{'R_extra':>10}{'residual':>11}"
              f"{'S(Rc)':>10}{'S(R_ex)':>10}{'sum S':>9}{'kappa':>9}  class")
        for m, Rc, Rx, Q, S, k in rows:
            print(f"  {m:<6}{Rc:>10.1f}{Rx:>10.1f}{Q:>11.1f}"
                  f"{S['Rc_measured']:>10.3f}{S['R_extra']:>10.3f}"
                  f"{sum(S.values()):>9.6f}{k:>9.2f}  {classify(k)}")
        print(f"\n  Hypothetical {HYPOTHETICAL_INPUT_ERR*100:.0f}% error on a single input"
              f" -> resulting error on the residual:")
        for m, Rc, Rx, Q, S, k in rows:
            print(f"    {m:<6}{k*HYPOTHETICAL_INPUT_ERR*100:>8.1f} %"
                  + ("   <-- exceeds 100%" if k * HYPOTHETICAL_INPUT_ERR > 1.0 else ""))
    return rows


# =====================================================================
# FAMILY B -- Chapter 5: the lambda_impurity calibration solve
#
#   residual = target_ratio - 1 - edge_term,    lambda_imp = lambda_bulk/residual
#
# The "1" is NOT an arbitrary constant: it is Matthiessen's bulk term,
# rho_bulk/rho_bulk.  Treating it as a term (rather than dropping it) is
# what makes P1's sum rule hold exactly here.
# =====================================================================
def _calib_terms(rho_cal=icm.rho_calibration_uohm_cm, rho_bulk=icm.rho_bulk_uohm_cm,
                 p=icm.p_default, W_cal=icm.W_calibration_nm,
                 lam_bulk=icm.lambda_bulk_nm):
    target = rho_cal / rho_bulk
    edge = ((1 - p) / (1 + p)) * (lam_bulk / W_cal)
    return {"target_ratio": +target, "bulk_term": -1.0, "edge_term": -edge}


def _lambda_imp(**kw):
    t = _calib_terms(**kw)
    return icm.lambda_bulk_nm / sum(t.values()) if "lam_bulk" not in kw \
        else kw["lam_bulk"] / sum(t.values())


def family_B(verbose=True):
    terms = _calib_terms()
    S, resid, k = term_sensitivities(terms)
    params = {
        "rho_calibration": icm.rho_calibration_uohm_cm,
        "rho_bulk": icm.rho_bulk_uohm_cm,
        "p_default": icm.p_default,
        "W_calibration": icm.W_calibration_nm,
        "lambda_bulk": icm.lambda_bulk_nm,
    }
    key = {"rho_calibration": "rho_cal", "rho_bulk": "rho_bulk",
           "p_default": "p", "W_calibration": "W_cal", "lambda_bulk": "lam_bulk"}
    par = {}
    for name, val in params.items():
        par[name] = log_sensitivity(lambda v, n=key[name]: _lambda_imp(**{n: v}), val)
    if verbose:
        print("\n" + "=" * 78)
        print("FAMILY B -- Chapter 5:  the lambda_impurity calibration solve")
        print("=" * 78)
        print(f"  residual = target_ratio - bulk_term - edge_term = {resid:.6f}")
        print(f"  lambda_impurity = lambda_bulk / residual        = "
              f"{icm.lambda_impurity_nm:.2f} nm")
        print(f"\n  TERM level:")
        for n in ("target_ratio", "bulk_term", "edge_term"):
            print(f"    {n:<16}{terms[n]:>10.4f}   S = {S[n]:>9.4f}")
        print(f"    {'sum S':<16}{'':>10}   S = {sum(S.values()):>9.6f}   "
              f"kappa = {k:.2f}  {classify(k)}")
        print(f"\n  PARAMETER level (d ln lambda_imp / d ln p):")
        for n, (s, c) in sorted(par.items(), key=lambda kv: -abs(kv[1][0])):
            print(f"    {n:<16}{s:>10.4f}   (conv {c:.1e})")
    return terms, S, resid, k, par


# =====================================================================
# FAMILY C -- Chapter 5: rho_GNR(W), audited two ways
#
#  (i)  lambda_impurity HELD FIXED (the formula as written): a Matthiessen
#       sum of three positive terms.
#  (ii) lambda_impurity RE-SOLVED from the calibration (the chapter's actual
#       dependency chain): the near-cancellation of Family B propagates in.
# =====================================================================
W_AUDIT_NM = (18.0, 22.0, 30.0, 40.0, 52.0)


def _rho_terms(W_nm, p=icm.p_default, lam_imp=None,
               rho_bulk=icm.rho_bulk_uohm_cm, lam_bulk=icm.lambda_bulk_nm):
    lam_imp = icm.lambda_impurity_nm if lam_imp is None else lam_imp
    edge = ((1 - p) / (1 + p)) * (lam_bulk / W_nm)
    imp = lam_bulk / lam_imp
    return {"bulk": rho_bulk * 1.0, "edge": rho_bulk * edge, "impurity": rho_bulk * imp}


def _rho_propagated(W_nm, **kw):
    """rho at width W with lambda_impurity RE-SOLVED for the given params."""
    lam_imp = _lambda_imp(**kw)
    sub = {k: v for k, v in kw.items() if k in ("p", "rho_bulk", "lam_bulk")}
    if "p" in sub:
        sub["p"] = sub.pop("p")
    return sum(_rho_terms(W_nm, lam_imp=lam_imp, **sub).values())


def family_C(verbose=True):
    fixed, prop = [], []
    key = {"rho_calibration": "rho_cal", "rho_bulk": "rho_bulk",
           "p_default": "p", "W_calibration": "W_cal", "lambda_bulk": "lam_bulk"}
    params = {
        "rho_calibration": icm.rho_calibration_uohm_cm,
        "rho_bulk": icm.rho_bulk_uohm_cm,
        "p_default": icm.p_default,
        "W_calibration": icm.W_calibration_nm,
        "lambda_bulk": icm.lambda_bulk_nm,
    }
    for W in W_AUDIT_NM:
        terms = _rho_terms(W)
        S, Q, k = term_sensitivities(terms)
        fixed.append((W, terms, S, Q, k))
        row = {}
        for name, val in params.items():
            row[name] = log_sensitivity(
                lambda v, n=key[name], w=W: _rho_propagated(w, **{n: v}), val)
        prop.append((W, row, max(abs(v[0]) for v in row.values())))
    if verbose:
        print("\n" + "=" * 78)
        print("FAMILY C -- Chapter 5:  rho_GNR(W)")
        print("=" * 78)
        print("  (i) TERM level, lambda_impurity held fixed (formula as written):")
        print(f"      {'W/nm':>6}{'rho':>9}{'S bulk':>9}{'S edge':>9}{'S imp':>9}"
              f"{'sum S':>10}{'kappa':>8}  class")
        for W, t, S, Q, k in fixed:
            print(f"      {W:>6.0f}{Q:>9.3f}{S['bulk']:>9.4f}{S['edge']:>9.4f}"
                  f"{S['impurity']:>9.4f}{sum(S.values()):>10.6f}{k:>8.3f}  {classify(k)}")
        print("\n  (ii) PARAMETER level, lambda_impurity RE-SOLVED "
              "(the chapter's real dependency chain):")
        names = list(params)
        print("      " + f"{'W/nm':>6}" + "".join(f"{n[:13]:>15}" for n in names)
              + f"{'kappa':>8}  class")
        for W, row, k in prop:
            print("      " + f"{W:>6.0f}"
                  + "".join(f"{row[n][0]:>15.4f}" for n in names)
                  + f"{k:>8.3f}  {classify(k)}")
        print(f"      worst convergence estimate over the table: "
              f"{max(c for _, r, _ in prop for _, c in r.values()):.1e}")
    return fixed, prop


# =====================================================================
# FAMILY D -- Chapter 5: the Cu liner model, W_eff = W - 2t
# =====================================================================
def family_D(verbose=True):
    t = icm.t_liner_nm_default
    rows = []
    for W in W_AUDIT_NM:
        S, Weff, k = term_sensitivities({"W_drawn": +W, "liner_2t": -2.0 * t})
        Sp, conv = log_sensitivity(
            lambda v, w=W: icm.cu_resistivity_with_liner(w, t_liner_nm=v), t)
        rows.append((W, Weff, S, k, Sp, conv))
    if verbose:
        print("\n" + "=" * 78)
        print(f"FAMILY D -- Chapter 5:  Cu liner model, W_eff = W - 2t  (t = {t:.1f} nm)")
        print("=" * 78)
        print(f"  {'W/nm':>6}{'W_eff':>9}{'S(W)':>9}{'S(2t)':>9}{'sum S':>10}"
              f"{'kappa_W_eff':>13}{'S_t(rho_eff)':>14}  class (W_eff)")
        for W, Weff, S, k, Sp, conv in rows:
            print(f"  {W:>6.0f}{Weff:>9.2f}{S['W_drawn']:>9.4f}{S['liner_2t']:>9.4f}"
                  f"{sum(S.values()):>10.6f}{k:>13.3f}{Sp:>14.4f}  {classify(k)}")
        print(f"  worst convergence estimate: {max(r[5] for r in rows):.1e}")
    return rows


# =====================================================================
# FAMILY E -- Chapter 6 cross-check (P5): does the general machinery
# reproduce 2026-09-21's published amplifications?
# =====================================================================
def family_E(verbose=True):
    from graphene_per_metal_crossover_model import (
        usable_metals, dw_for_metal, net_response_dw, ELL_RANGE,
    )
    ell_short, ell_flat = ELL_RANGE[0], np.inf
    ms = usable_metals()
    rows = []
    for mA in ms:
        for mB in ms:
            if mA == mB:
                continue
            dA0, dB0 = dw_for_metal(mA, ell_flat), dw_for_metal(mB, ell_flat)
            dA1, dB1 = dw_for_metal(mA, ell_short), dw_for_metal(mB, ell_short)
            N0 = net_response_dw(dA0, dB0, e_applied=0.0)
            N1 = net_response_dw(dA1, dB1, e_applied=0.0)
            pct_in = max(abs(dA1 - dA0) / abs(dA0), abs(dB1 - dB0) / abs(dB0)) * 100
            pct_out = abs(N1 - N0) / abs(N0) * 100
            same_sign = (dA0 * dB0) > 0
            rows.append((f"{mA}/{mB}", "same-sign" if same_sign else "straddling",
                         N0, pct_in, pct_out, pct_out / pct_in))
    rows.sort(key=lambda r: -r[5])
    if verbose:
        print("\n" + "=" * 78)
        print("FAMILY E -- Chapter 6 cross-check: same machinery, 2026-09-21's pairs")
        print("=" * 78)
        print(f"  {'pair':<9}{'class':<12}{'N(flat)':>10}{'in %':>9}{'out %':>9}"
              f"{'amplif.':>10}  published 2026-09-21")
        pub = {"Cu/Au": 4.55, "Au/Cu": 4.55, "Cu/Pt": 0.07, "Pt/Cu": 0.07,
               "Au/Pt": 0.11, "Pt/Au": 0.11}
        for p_, c, N0, i, o, a in rows:
            ref = pub.get(p_)
            d = f"{ref:.2f}   ({abs(a-ref)/ref*100:>5.1f}% off)" if ref else ""
            print(f"  {p_:<9}{c:<12}{N0:>10.4f}{i:>9.2f}{o:>9.2f}{a:>10.2f}  {d}")
    return rows


# =====================================================================
# VALIDATIONS against exactly known values
# =====================================================================
def validate_sum_rule(verbose=True):
    """
    EXACT check #1 (P1).  sum_i S_i == 1 for every term-level decomposition
    in every family.  Reported as the worst |sum - 1| over all of them.
    """
    worst, n, where = 0.0, 0, None
    for _, _, _, _, S, _ in [(m, a, b, c, S, k) for m, a, b, c, S, k in family_A_rows()]:
        n += 1
        d = abs(sum(S.values()) - 1.0)
        if d > worst:
            worst, where = d, "A"
    _, S, _, _, _ = family_B(verbose=False)
    n += 1
    if abs(sum(S.values()) - 1.0) > worst:
        worst, where = abs(sum(S.values()) - 1.0), "B"
    fixed, _ = family_C(verbose=False)
    for _, _, S, _, _ in fixed:
        n += 1
        d = abs(sum(S.values()) - 1.0)
        if d > worst:
            worst, where = d, "C"
    for _, _, S, _, _, _ in family_D(verbose=False):
        n += 1
        d = abs(sum(S.values()) - 1.0)
        if d > worst:
            worst, where = d, "D"
    if verbose:
        print("\nValidation 1 (EXACT, P1): Euler sum rule  sum_i S_i == 1")
        print(f"  decompositions checked : {n}")
        print(f"  worst |sum(S) - 1|     : {worst:.3e}   (family {where})")
        print(f"  {'PASS' if worst < 1e-12 else 'FAIL'}")
    return worst, n


def validate_power_laws(verbose=True):
    """
    EXACT check #2.  The finite-difference machinery must recover EXACTLY
    KNOWN exponents for quantities that are pure power laws in an input:

      (a) lambda_impurity  proportional to lambda_bulk^(+1) when the edge
          term is held at its own lambda_bulk -- so we test the clean case
          lambda_imp(residual fixed) -> exponent exactly +1.
      (b) rho_GNR proportional to rho_bulk^(+1) at fixed lambda_impurity
          -> exponent exactly +1.
      (c) the bulk-only limit: with edge and impurity terms zeroed,
          S(rho_bulk) == 1 and every other S == 0 EXACTLY.
      (d) a pure inverse: Q = c/x has S == -1 exactly.
    """
    checks = []
    S, c = log_sensitivity(lambda v: v / 0.1522, 55.0)
    checks.append(("(a) lambda_bulk/resid, expect +1", S, 1.0, c))
    S, c = log_sensitivity(
        lambda v: sum(_rho_terms(30.0, rho_bulk=v).values()), icm.rho_bulk_uohm_cm)
    checks.append(("(b) rho ~ rho_bulk^1 (fixed lam_imp), expect +1", S, 1.0, c))
    Sd, Qd, kd = term_sensitivities({"bulk": 1.2, "edge": 0.0, "impurity": 0.0})
    checks.append(("(c) bulk-only S(bulk), expect +1 EXACT", Sd["bulk"], 1.0, 0.0))
    checks.append(("(c) bulk-only S(edge), expect  0 EXACT", Sd["edge"], 0.0, 0.0))
    S, c = log_sensitivity(lambda v: 7.0 / v, 3.0)
    checks.append(("(d) Q = c/x, expect -1", S, -1.0, c))
    worst = max(abs(s - e) for _, s, e, _ in checks)
    if verbose:
        print("\nValidation 2 (EXACT): known power-law exponents recovered")
        for name, s, e, c in checks:
            print(f"  {name:<46}{s:>12.9f}  |err| {abs(s-e):.2e}  (conv {c:.1e})")
        print(f"  worst |error| : {worst:.3e}   {'PASS' if worst < 1e-8 else 'FAIL'}")
    return worst


def validate_degenerate_cancellation(verbose=True):
    """
    EXACT check #3.  The diagnostic must blow up where it is supposed to and
    stay finite where it is not:
      (a) Q = A - A is an exact total cancellation: 1/kappa == 0 exactly.
      (b) Q = A + A is maximally reinforcing: kappa == 0.5 exactly.
      (c) kappa is invariant under rescaling ALL terms by any c != 0
          (it is a ratio of homogeneous-degree-1 quantities) -- bitwise.
    """
    A = 517.3
    _, _, k_sum = term_sensitivities({"a": A, "b": A})
    _, Qd, k_diff = term_sensitivities({"a": A, "b": -A})
    inv = 0.0 if np.isinf(k_diff) else 1.0 / k_diff
    scales = [1e-9, 0.5, 1.0, 3.7, 1e9]
    base = term_sensitivities({"a": 584.0, "b": -533.1})[2]
    ks = [term_sensitivities({"a": 584.0 * s, "b": -533.1 * s})[2] for s in scales]
    bitwise = sum(k == base for k in ks)
    worst_rel = max(abs(k - base) / base for k in ks)
    if verbose:
        print("\nValidation 3 (EXACT): degenerate limits and scale invariance")
        print(f"  (a) Q = A - A  : Q = {Qd:.1e}, 1/kappa = {inv:.3e}  "
              f"{'PASS' if inv == 0.0 else 'FAIL'}")
        print(f"  (b) Q = A + A  : kappa = {k_sum:.17g}  "
              f"{'PASS' if k_sum == 0.5 else 'FAIL'}")
        print(f"  (c) scale invariance of kappa : {bitwise}/{len(scales)} bitwise, "
              f"worst relative deviation {worst_rel:.2e}")
        print(f"      (NOT 5/5, and reported as such: kappa is a ratio, so rescaling "
              f"all terms by 1e+-9\n       changes the rounding of the division. "
              f"The invariance is exact in exact\n       arithmetic and holds to "
              f"{worst_rel:.0e} in double precision -- a floating-point\n       "
              f"statement, not a modelling one.)")
    return inv, k_sum, bitwise, worst_rel


def validate_calibration_point_pinning(verbose=True):
    """
    EXACT check #4, NOT designed in -- it fell out of Family C and is
    recorded because it is the sharpest exact statement the audit produced.

    At W == W_calibration the model is pinned to the calibration datum by
    construction, so:
      (a) rho(W_calibration) == rho_calibration BITWISE;
      (b) S(rho_calibration) == 1 and
      (c) S(rho_bulk) == S(p) == S(lambda_bulk) == 0 EXACTLY (0.0, not
          "small") -- three independent exact zeros, the strongest form of
          check this repo uses.
    (c) is a real structural fact: at the calibration width the impurity
    term absorbs whatever the other terms do, so the prediction there
    carries no information from rho_bulk, p or lambda_bulk at all.
    """
    Wc = icm.W_calibration_nm
    bitwise = _rho_propagated(Wc) == icm.rho_calibration_uohm_cm
    S_cal, c_cal = log_sensitivity(
        lambda v: _rho_propagated(Wc, rho_cal=v), icm.rho_calibration_uohm_cm)
    zeros = {}
    for name, kw, val in (("rho_bulk", "rho_bulk", icm.rho_bulk_uohm_cm),
                          ("p_default", "p", icm.p_default),
                          ("lambda_bulk", "lam_bulk", icm.lambda_bulk_nm)):
        zeros[name] = log_sensitivity(
            lambda v, k=kw: _rho_propagated(Wc, **{k: v}), val)[0]
    n_exact = sum(v == 0.0 for v in zeros.values())
    if verbose:
        print("\nValidation 4 (EXACT): the calibration point pins rho")
        print(f"  (a) rho({Wc:.0f} nm) == rho_calibration bitwise : {bitwise}")
        print(f"  (b) S(rho_calibration) = {S_cal:.12f}  |err| {abs(S_cal-1):.1e}")
        print(f"  (c) exact zeros : {n_exact}/3  "
              + ", ".join(f"S({k}) = {v:.1f}" for k, v in zeros.items()))
        print(f"  {'PASS' if (bitwise and n_exact == 3 and abs(S_cal-1) < 1e-9) else 'FAIL'}")
    return bitwise, S_cal, zeros


def sign_flip_margins(rows, verbose=True):
    """
    RESULT.  kappa says how an input error is amplified; the operationally
    useful inverse is: HOW BIG an error in R_extra would be needed to flip
    the SIGN of the residual?  For Q = Rc - R_extra, R_extra must move by
    the fraction

        f = (R_extra - Rc) / R_extra = -Q / R_extra                  (C3)

    A small |f| means the sign of the published residual is fragile; a
    large |f| means it is robust NO MATTER how large kappa is.  These two
    statements are independent, and Section 4.7's negative residuals are
    exactly where they come apart.
    """
    out = []
    for m, Rc, Rx, Q, S, k in rows:
        out.append((m, Q, k, -Q / Rx))
    if verbose:
        print("\n" + "=" * 78)
        print("RESULT -- sign-flip margin: how wrong would R_extra have to be?")
        print("=" * 78)
        print(f"  {'metal':<6}{'residual':>11}{'kappa':>9}"
              f"{'R_extra must move by':>24}  verdict on the SIGN")
        for m, Q, k, f in sorted(out, key=lambda r: abs(r[3])):
            v = ("fragile" if abs(f) < 0.10 else
                 "robust" if abs(f) > 0.50 else "intermediate")
            print(f"  {m:<6}{Q:>11.1f}{k:>9.2f}{f*100:>23.1f} %  {v}")
    return out


# =====================================================================
def plot_audit(famA, famC_prop, famD, fname="sensitivity_audit.png"):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    names = [r[0] for r in famA]
    ks = [r[5] for r in famA]
    ax[0].bar(names, ks, color=["#c0392b" if k >= NEAR_CANCELLING else "#2980b9" for k in ks])
    ax[0].axhline(NEAR_CANCELLING, ls="--", c="k", lw=1)
    ax[0].axhline(4.55, ls=":", c="#c0392b", lw=1.2)
    ax[0].text(0.02, 4.7, "Ch.6 same-sign pair, 4.55", fontsize=7, color="#c0392b",
               transform=ax[0].get_yaxis_transform())
    ax[0].text(0.02, 3.15, "near-cancelling threshold", fontsize=7,
               transform=ax[0].get_yaxis_transform())
    ax[0].set_ylabel(r"$\kappa$")
    ax[0].set_title("(a) Ch.4 $R_{transmission}=R_c-R_{extra}$", fontsize=10)

    Ws = [r[0] for r in famC_prop]
    ax[1].plot(Ws, [r[2] for r in famC_prop], "o-", label=r"$\lambda_{imp}$ re-solved")
    ax[1].plot(Ws, [1.0] * len(Ws), "s--", c="grey",
               label=r"term level, $\lambda_{imp}$ fixed ($\kappa\leq1$)")
    ax[1].axhline(NEAR_CANCELLING, ls="--", c="k", lw=1)
    ax[1].axhline(1.0, ls=":", c="k", lw=0.8)
    ax[1].set_xlabel("linewidth W (nm)")
    ax[1].set_ylabel(r"$\kappa(\rho_{GNR})$")
    ax[1].set_title(r"(b) Ch.5 $\rho_{GNR}$: the calibration propagates in", fontsize=10)
    ax[1].legend(fontsize=7)

    ax[2].plot([r[0] for r in famD], [r[3] for r in famD], "o-", c="#8e44ad",
               label=r"$\kappa(W_{eff}=W-2t)$")
    ax[2].plot([r[0] for r in famD], [abs(r[4]) for r in famD], "^-", c="#16a085",
               label=r"$|S_t(\rho_{eff})|$")
    ax[2].axhline(NEAR_CANCELLING, ls="--", c="k", lw=1)
    ax[2].axhline(1.0, ls=":", c="k", lw=0.8)
    ax[2].set_xlabel("drawn linewidth W (nm)")
    ax[2].set_ylabel("sensitivity")
    ax[2].set_title("(c) Ch.5 Cu liner: a subtraction in the geometry", fontsize=10)
    ax[2].legend(fontsize=7)

    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle("Condition-number audit of Chapters 4 and 5 "
                 "(2026-09-22): which computed quantities are near-cancellations?",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fname, dpi=150)
    print(f"\nFigure written: {fname}")


def main():
    print("=" * 78)
    print("CONDITION-NUMBER AUDIT OF CHAPTERS 4 AND 5 -- 2026-09-22")
    print("Predictions P1-P6 pre-registered in "
          "notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md")
    print("=" * 78)

    w1, n1 = validate_sum_rule()
    w2 = validate_power_laws()
    validate_degenerate_cancellation()
    validate_calibration_point_pinning()

    A = family_A()
    famB = family_B()
    fixedC, propC = family_C()
    D = family_D()
    E = family_E()
    sign_flip_margins(A)

    print("\n" + "=" * 78)
    print("VERDICT ON THE PRE-REGISTERED PREDICTIONS")
    print("=" * 78)
    kA = max(r[5] for r in A)
    kB = famB[3]
    kCfix = max(r[4] for r in fixedC)
    kCprop = max(r[2] for r in propC)
    kD = max(r[3] for r in D)

    print(f"  P1 sum rule == 1 exactly           : "
          f"{'PASS' if w1 < 1e-12 else 'FAIL'}  (worst {w1:.1e} over {n1} decompositions)")
    p2a = kCprop > 1.0
    p2b = kCprop < NEAR_CANCELLING
    print(f"  P2 rho_GNR kappa > 1 when propagated: "
          f"{'PASS' if p2a else 'FAIL'}  (kappa = {kCprop:.3f})")
    print(f"     ...and stays below 3            : "
          f"{'PASS' if p2b else 'FAIL'}")
    order = sorted([("Ch5 calib", kB), ("Ch4 R_trans", kA),
                    ("Ch5 liner", kD), ("Ch5 rho", kCfix)], key=lambda t: -t[1])
    pred = ["Ch5 calib", "Ch4 R_trans", "Ch5 liner", "Ch5 rho"]
    got = [o[0] for o in order]
    print(f"  P3 rank order                      : "
          f"{'PASS' if got == pred else 'FAIL'}")
    print(f"     predicted : {' > '.join(pred)}")
    print(f"     measured  : {' > '.join(f'{n} ({k:.2f})' for n, k in order)}")
    nA3 = sum(1 for r in A if r[5] >= NEAR_CANCELLING)
    worstA = max(A, key=lambda r: r[5])[0]
    print(f"  P4 >=3 metals near-cancelling, Pd worst : "
          f"{'PASS' if (nA3 >= 3 and worstA == 'Pd') else 'FAIL'}  "
          f"({nA3}/4 metals, worst = {worstA})")
    pub = {"Cu/Au": 4.55, "Cu/Pt": 0.07, "Au/Pt": 0.11}
    offs = [abs(a - pub[p_]) / pub[p_] for p_, _, _, _, _, a in E if p_ in pub]
    print(f"  P5 Ch.6 amplifications within 10%  : "
          f"{'PASS' if max(offs) <= 0.10 else 'FAIL'}  (worst {max(offs)*100:.1f}% off)")
    bars = [(f"Ch4 R_trans({r[0]})", r[5] * HYPOTHETICAL_INPUT_ERR) for r in A]
    bars += [("Ch5 lambda_imp", kB * HYPOTHETICAL_INPUT_ERR),
             ("Ch5 rho_GNR", kCprop * HYPOTHETICAL_INPUT_ERR),
             ("Ch5 liner W_eff", kD * HYPOTHETICAL_INPUT_ERR)]
    over = [(n, b) for n, b in bars if b > 1.0]
    print(f"  P6 a published number with >100% bar at 5% input error : "
          f"{'PASS' if over else 'FAIL'}")
    for n, b in sorted(over, key=lambda t: -t[1]):
        print(f"       {n:<24}{b*100:>8.0f} %")

    plot_audit(A, propC, D)


if __name__ == "__main__":
    main()
