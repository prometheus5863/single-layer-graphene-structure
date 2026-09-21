"""
graphene_per_metal_crossover_model.py

Replaces the SINGLE p/n crossover work function

    W_CROSS_CHEM = 5.4 eV        (shared by all seven metals since 2026-09-18)

with a PER-METAL crossover derived from the short-range chemical interface
term of

  Khomyakov, Giovannetti, Rusu, Brocks, van den Brink, Kelly,
  Phys. Rev. B 79, 195425 (2009)   [arXiv:0902.1203]

whose Eq. 4 makes the interface potential step a function of the
metal-graphene separation d.  Charge neutrality then puts the crossover at

    w_cross(d) = W_G + Delta_c(d)                                    (*)

so 5.4 eV is (*) evaluated at the PHYSISORBED separation d ~ 3.3 A, not a
constant of nature.

Literature, the extraction failure, the anchored one-parameter form of
Delta_c, and FOUR PRE-REGISTERED PREDICTIONS (committed to git before this
file existed):
  notes/2026-09-21-per-metal-crossover-from-the-chemical-interface-term.md

--------------------------------------------------------------------------
WHAT IS AND IS NOT CLAIMED
--------------------------------------------------------------------------
* NO fitted coefficient of Khomyakov et al.'s Eq. 4 is used.  Two ar5iv
  fetches with differently-worded prompts both returned Eq. 4 in symbolic
  form with the numbers absent.  That is the SECOND extraction failure on
  this paper (2026-09-20's was Table I, contradictory rather than missing).
* Instead Delta_c is a ONE-PARAMETER family pinned to the single value both
  2026-09-20 passes and today's fetch agree on, Delta_c(3.3 A) ~ 0.9 eV:

      Delta_c(d ; ell) = 0.9 * exp( -(d - 3.3 A) / ell )             (N1)

  ell is the short-range decay length, swept over 0.3-1.5 A.  Nothing below
  rests on a single value of it.
* (N1) drops Eq. 4's polynomial prefactor.  That is a real approximation,
  and RESULT 2 below is precisely the place it breaks.
* ell -> infinity reproduces the flat 5.4 eV model BITWISE.  That is an
  exact reduction test (Validation 2), not a plausible-range check -- the
  distinction that caught real bugs on 2026-09-17 and 2026-09-19.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_contact_doping_model import (
    METAL_WORK_FUNCTIONS, W_GRAPHENE, lambda_decay,
)
from graphene_photodetector_model import L_channel
from graphene_photodetector_collection_model import E_BIAS
from graphene_photodetector_signed_carrier_model import (
    total_field, net_response, _transport,
    Q_HOLE, Q_ELECTRON, N_POINTS, W_CROSS_CHEM, signed_offset,
)
from graphene_contact_doping_nonlinear_model import D_EQ, CHEMISORBED, D0_SEPARATION

# ---------------------------------------------------------------------
# The one anchor, and the one unknown
# ---------------------------------------------------------------------
D_ANCHOR = 3.3e-10        # m   -- physisorbed equilibrium separation
DC_ANCHOR = 0.9           # eV  -- Delta_c there; agreed across all extractions
ELL_DEFAULT = 0.6e-10     # m   -- midpoint of the swept range
ELL_RANGE = (0.3e-10, 1.5e-10)

# Highest elemental work function in common tables (Pt/Se ~ 5.9 eV).  Used
# only as a reductio bound in RESULT 2 -- nothing depends on its exact value.
W_MAX_ELEMENTAL = 5.9
# Largest Fermi-level shift Khomyakov et al. report for any metal.
DEF_MAX_OBSERVED = 0.5


def delta_c(d, ell=ELL_DEFAULT):
    """
    Eq. (N1).  Exactly DC_ANCHOR at d == D_ANCHOR, for any ell.
    ell = np.inf returns DC_ANCHOR for every d, bitwise.
    """
    return DC_ANCHOR * np.exp(-(np.asarray(d, dtype=float) - D_ANCHOR) / ell)


def w_cross_for_metal(metal, ell=ELL_DEFAULT):
    """
    (w_cross_eV, status) for one metal.

    status is one of:
      "ok"             -- physisorbed, d_eq >= d0, (N1) is inside its regime
      "extrapolated"   -- chemisorbed, d_eq < d0; value returned for RESULT 2's
                          reductio ONLY, and refused by usable_metals()
      "no-separation"  -- Cr: Khomyakov et al. tabulate no d_eq and it was not
                          guessed (same stance as 2026-09-20)
    """
    d = D_EQ.get(metal)
    if d is None:
        return None, "no-separation"
    status = "ok" if d >= D0_SEPARATION else "extrapolated"
    return W_GRAPHENE + float(delta_c(d, ell)), status


def usable_metals():
    """The metals of METAL_WORK_FUNCTIONS that (N2) may be applied to."""
    out = []
    for m in METAL_WORK_FUNCTIONS:
        _, st = w_cross_for_metal(m)
        if st == "ok":
            out.append(m)
    return sorted(out, key=lambda m: METAL_WORK_FUNCTIONS[m])


def dw_for_metal(metal, ell=ELL_DEFAULT):
    """Signed offset dW = W_metal - w_cross(metal), or None if out of regime."""
    wc, st = w_cross_for_metal(metal, ell)
    if st != "ok":
        return None
    return METAL_WORK_FUNCTIONS[metal] - wc


# ---------------------------------------------------------------------
# Field and response at the dW level
#
# The existing net_response() takes ONE w_cross for both contacts and forms
# dW internally, so it cannot express a per-metal crossover.  These two
# functions are the same physics re-entered one level lower, and Validation
# 1 below pins them to the originals BITWISE.
# ---------------------------------------------------------------------
def total_field_dw(x, dW_A, dW_B, e_applied=E_BIAS, L=L_channel, lam=lambda_decay):
    """Mirror of signed_carrier_model.total_field(), taking dW directly."""
    return (-e_applied
            - dW_A / lam / (1.0 + x / lam) ** 2
            + dW_B / lam / (1.0 + (L - x) / lam) ** 2)


def net_response_dw(dW_A, dW_B, e_applied=E_BIAS, L=L_channel,
                    n_points=N_POINTS, carriers=("hole", "electron")):
    """Mirror of signed_carrier_model.net_response(), taking dW directly."""
    x = np.linspace(0.0, L, n_points)
    E = total_field_dw(x, dW_A, dW_B, e_applied=e_applied, L=L)
    integrand = np.zeros_like(x)
    if "hole" in carriers:
        s_h, p_h = _transport(+E, x)
        integrand = integrand + Q_HOLE * s_h * p_h
    if "electron" in carriers:
        s_e, p_e = _transport(-E, x)
        integrand = integrand + Q_ELECTRON * s_e * p_e
    return float(np.trapezoid(integrand, x) / L)


def net_response_per_metal(mA, mB, ell=ELL_DEFAULT, e_applied=0.0):
    """Zero-bias net response for an ordered pair under the per-metal crossover."""
    dA, dB = dw_for_metal(mA, ell), dw_for_metal(mB, ell)
    if dA is None or dB is None:
        return None
    return net_response_dw(dA, dB, e_applied=e_applied)


# =====================================================================
# VALIDATIONS -- all four against EXACTLY known values
# =====================================================================
def validate_reduction_to_signed_model(verbose=True):
    """
    EXACT check #1.  net_response_dw(W_A - 5.4, W_B - 5.4) must equal
    net_response(W_A, W_B, w_cross=5.4) BITWISE, for every ordered pair of
    all seven metals.  This is what licenses using the dW-level entry point
    for everything that follows: it is not "a similar model", it is the same
    arithmetic with the subtraction moved outward.
    """
    worst, worst_pair, n_bitwise = 0.0, None, 0
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    total = 0
    for mA, WA in metals:
        for mB, WB in metals:
            total += 1
            a = net_response(WA, WB, e_applied=0.0, w_cross=W_CROSS_CHEM)
            b = net_response_dw(signed_offset(WA), signed_offset(WB), e_applied=0.0)
            if a == b:
                n_bitwise += 1
            d = abs(a - b)
            if d > worst:
                worst, worst_pair = d, (mA, mB)
    if verbose:
        print("\nValidation 1 (EXACT): dW-level entry point == signed model")
        print(f"  ordered pairs checked : {total}")
        print(f"  bitwise identical     : {n_bitwise}/{total}")
        print(f"  worst |difference|    : {worst:.3e}   ({worst_pair})")
    return n_bitwise, total, worst


def validate_flat_limit(verbose=True):
    """
    EXACT check #2.  ell -> infinity must collapse the per-metal machinery
    onto the flat 5.4 eV convention BITWISE: Delta_c(d, inf) == 0.9 for every
    d, so w_cross == W_G + 0.9 for every metal, and every pair response
    matches the existing model exactly.

    NOTE the one place this is NOT bitwise, reported rather than hidden:
    W_GRAPHENE + DC_ANCHOR (4.5 + 0.9) need not be the same double as the
    literal 5.4.  The check prints the ULP gap and then compares responses
    against W_GRAPHENE + DC_ANCHOR, which is the quantity the model actually
    uses.
    """
    dc_exact = all(delta_c(d, np.inf) == DC_ANCHOR for d in D_EQ.values())
    wc_flat = W_GRAPHENE + DC_ANCHOR
    ulp_gap = abs(wc_flat - W_CROSS_CHEM)
    worst, n_bitwise, total = 0.0, 0, 0
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    for mA, WA in metals:
        for mB, WB in metals:
            total += 1
            a = net_response(WA, WB, e_applied=0.0, w_cross=wc_flat)
            b = net_response_dw(WA - wc_flat, WB - wc_flat, e_applied=0.0)
            if a == b:
                n_bitwise += 1
            worst = max(worst, abs(a - b))
    wc_all_flat = all(
        w_cross_for_metal(m, np.inf)[0] == wc_flat
        for m in D_EQ if m in METAL_WORK_FUNCTIONS
    )
    if verbose:
        print("\nValidation 2 (EXACT): ell -> infinity reproduces the flat model")
        print(f"  Delta_c(d, inf) == 0.9 bitwise for every tabulated d : {dc_exact}")
        print(f"  w_cross(metal, inf) == W_G + 0.9 bitwise, all metals : {wc_all_flat}")
        print(f"  (4.5 + 0.9) vs literal 5.4, |gap|                    : {ulp_gap:.3e} eV")
        print(f"  pair responses bitwise identical                     : {n_bitwise}/{total}")
        print(f"  worst |difference|                                   : {worst:.3e}")
    return dc_exact and wc_all_flat, n_bitwise, total, ulp_gap


def validate_anchor(verbose=True):
    """
    EXACT check #3.  Delta_c(3.3 A) == 0.9 eV bitwise for EVERY ell in the
    swept range -- the anchor is the one datum that is trusted, and it must
    survive the parametrisation untouched.
    """
    ells = np.linspace(ELL_RANGE[0], ELL_RANGE[1], 61)
    oks = [float(delta_c(D_ANCHOR, e)) == DC_ANCHOR for e in ells]
    if verbose:
        print("\nValidation 3 (EXACT): the anchor is untouched by ell")
        print(f"  Delta_c(3.3 A, ell) == 0.9 bitwise for {sum(oks)}/{len(oks)} "
              f"values of ell in [{ELL_RANGE[0]*1e10:.1f}, {ELL_RANGE[1]*1e10:.1f}] A")
    return all(oks)


def validate_symmetric_and_conjugation(verbose=True):
    """
    EXACT check #4, two exact symmetries that a per-metal crossover could
    plausibly have broken and must not:

    (a) SYMMETRIC PAIR.  Identical contacts share a w_cross by construction,
        so dW_A == dW_B and E(x) stays antisymmetric about L/2: N == 0 exactly,
        for every usable metal and every ell.
    (b) CHARGE CONJUGATION.  At zero bias, N(-dW_A, -dW_B) == -N(dW_A, dW_B)
        exactly.  Stated at the dW level this is sharper than the 2026-09-18
        version, which had to reflect W through a crossover to negate dW and
        so could only test it at one w_cross.
    """
    ells = [ELL_RANGE[0], ELL_DEFAULT, ELL_RANGE[1]]
    worst_sym, n_sym = 0.0, 0
    for ell in ells:
        for m in usable_metals():
            d = dw_for_metal(m, ell)
            n_sym += 1
            worst_sym = max(worst_sym, abs(net_response_dw(d, d, e_applied=0.0)))
    worst_cc, n_cc = 0.0, 0
    for ell in ells:
        ms = usable_metals()
        for mA in ms:
            for mB in ms:
                dA, dB = dw_for_metal(mA, ell), dw_for_metal(mB, ell)
                N = net_response_dw(dA, dB, e_applied=0.0)
                Nf = net_response_dw(-dA, -dB, e_applied=0.0)
                worst_cc = max(worst_cc, abs(N + Nf))
                n_cc += 1
    if verbose:
        print("\nValidation 4 (EXACT): symmetry survives the per-metal crossover")
        print(f"  (a) symmetric pair N == 0   : worst |N| = {worst_sym:.3e} "
              f"over {n_sym} (metal, ell) cases")
        print(f"  (b) charge conjugation      : worst |N + N_flip| = {worst_cc:.3e} "
              f"over {n_cc} (pair, ell) cases")
    return worst_sym, worst_cc


# =====================================================================
# RESULTS
# =====================================================================
def result_1_per_metal_offsets(verbose=True):
    """
    RESULT 1.  Per-metal crossover and dW for the metals the model admits,
    swept over ell.  Tests P1 (no sign changes) and P2 (|dW| shifts <= 10%).
    ALL metals are listed, with the reason when one cannot be evaluated --
    2026-09-20's "enumerate the space" rule.
    """
    ells = np.linspace(*ELL_RANGE, 25)
    rows, worst_frac, sign_changes = [], 0.0, []
    for m in sorted(METAL_WORK_FUNCTIONS, key=lambda k: METAL_WORK_FUNCTIONS[k]):
        W = METAL_WORK_FUNCTIONS[m]
        wc, st = w_cross_for_metal(m)
        flat_dw = signed_offset(W)
        if st != "ok":
            rows.append((m, W, None, None, None, None, st))
            continue
        dws = np.array([W - w_cross_for_metal(m, e)[0] for e in ells])
        fr = np.max(np.abs(dws - flat_dw) / abs(flat_dw))
        worst_frac = max(worst_frac, fr)
        if np.any(np.sign(dws) != np.sign(flat_dw)):
            sign_changes.append(m)
        rows.append((m, W, flat_dw, float(dws.min()), float(dws.max()), float(fr), st))
    if verbose:
        print("\nRESULT 1: per-metal crossover, all seven metals")
        print(f"  ell swept over [{ELL_RANGE[0]*1e10:.1f}, {ELL_RANGE[1]*1e10:.1f}] A")
        print(f"  {'metal':<6}{'W':>7}{'d_eq':>8}{'dW(flat)':>11}"
              f"{'dW min':>10}{'dW max':>10}{'max|Δ|/|dW|':>13}  status")
        for m, W, fdw, lo, hi, fr, st in rows:
            d = D_EQ.get(m)
            ds = f"{d*1e10:.2f}" if d else "   --"
            if st != "ok":
                print(f"  {m:<6}{W:>7.2f}{ds:>8}{signed_offset(W):>11.4f}"
                      f"{'--':>10}{'--':>10}{'--':>13}  {st}")
            else:
                print(f"  {m:<6}{W:>7.2f}{ds:>8}{fdw:>11.4f}"
                      f"{lo:>10.4f}{hi:>10.4f}{fr*100:>12.2f}%  {st}")
        print(f"  P1 (no doping-sign change): "
              f"{'PASS' if not sign_changes else 'FAIL ' + str(sign_changes)}")
        print(f"  P2 (|dW| shifts <= 10%)   : worst {worst_frac*100:.2f}%  "
              f"{'PASS' if worst_frac <= 0.10 else 'FAIL'}")
        if worst_frac > 0.10:
            print(f"      P2 holds for ell >= {p2_threshold()*1e10:.2f} A and fails "
                  f"below it; the failure is driven by the SHORT end of the "
                  f"swept range, not by the midpoint.")
    return rows, worst_frac, sign_changes


def p2_threshold(tol=0.10, lo=0.05e-10, hi=5.0e-10):
    """
    The shortest decay length ell for which EVERY usable metal's |dW| is
    still within `tol` of its flat-crossover value.  |dW| shift is monotone
    decreasing in ell (Delta_c -> the anchor as ell grows), so a bisection
    is exact to machine precision rather than a scan.
    """
    def worst(ell):
        w = 0.0
        for m in usable_metals():
            f = signed_offset(METAL_WORK_FUNCTIONS[m])
            d = METAL_WORK_FUNCTIONS[m] - w_cross_for_metal(m, ell)[0]
            w = max(w, abs(d - f) / abs(f))
        return w
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if worst(mid) > tol:
            lo = mid
        else:
            hi = mid
    return hi


def result_2_chemisorbed_reductio(verbose=True):
    """
    RESULT 2.  Extrapolating (N1) inward to the chemisorbed separations.
    Tests P4.  This is where dropping Eq. 4's polynomial prefactor breaks,
    and the breakage is informative rather than merely a limitation.
    """
    ells = [0.3e-10, 0.6e-10, 1.0e-10, 1.5e-10]
    chem = [m for m in METAL_WORK_FUNCTIONS
            if w_cross_for_metal(m)[1] == "extrapolated"]
    chem.sort(key=lambda m: D_EQ[m])
    rows, all_above = [], True
    for m in chem:
        vals = []
        for e in ells:
            wc, _ = w_cross_for_metal(m, e)
            vals.append(wc)
            if wc <= W_MAX_ELEMENTAL:
                all_above = False
        rows.append((m, D_EQ[m], vals))
    if verbose:
        print("\nRESULT 2: (N1) extrapolated inward to the chemisorbed metals")
        print("  w_cross (eV), and the implied |dW| it would force:")
        print(f"  {'metal':<6}{'d_eq/A':>8}" +
              "".join(f"{'l=' + f'{e*1e10:.1f}A':>16}" for e in ells))
        for m, d, vals in rows:
            cells = "".join(f"{v:>8.2f}{'(' + f'{abs(METAL_WORK_FUNCTIONS[m]-v):.2f}' + ')':>8}"
                            for v in vals)
            print(f"  {m:<6}{d*1e10:>8.2f}" + cells)
        print(f"  highest elemental work function anywhere : {W_MAX_ELEMENTAL:.2f} eV")
        print(f"  largest |dE_F| Khomyakov et al. report   : {DEF_MAX_OBSERVED:.2f} eV")
        print(f"  P4 (every extrapolated w_cross above every elemental W): "
              f"{'PASS' if all_above else 'FAIL'}")
    return rows, all_above


def result_3_pairs(verbose=True):
    """
    RESULT 3.  Every ordered pair of the usable metals, flat vs per-metal,
    swept over ell.  Tests P3 (same-sign amplification).
    """
    ms = usable_metals()
    ells = np.linspace(*ELL_RANGE, 25)
    rows = []
    for mA in ms:
        for mB in ms:
            if mA == mB:
                continue
            flat = net_response_dw(signed_offset(METAL_WORK_FUNCTIONS[mA]),
                                   signed_offset(METAL_WORK_FUNCTIONS[mB]),
                                   e_applied=0.0)
            vals = np.array([net_response_per_metal(mA, mB, e) for e in ells])
            pair_frac = float(np.max(np.abs(vals - flat) / abs(flat))) if flat else float("nan")
            dw_frac = 0.0
            for m in (mA, mB):
                fdw = signed_offset(METAL_WORK_FUNCTIONS[m])
                dws = np.array([METAL_WORK_FUNCTIONS[m] - w_cross_for_metal(m, e)[0]
                                for e in ells])
                dw_frac = max(dw_frac, float(np.max(np.abs(dws - fdw) / abs(fdw))))
            same_sign = (signed_offset(METAL_WORK_FUNCTIONS[mA]) *
                         signed_offset(METAL_WORK_FUNCTIONS[mB])) > 0
            rows.append((mA, mB, same_sign, flat, float(vals.min()),
                         float(vals.max()), pair_frac, dw_frac,
                         pair_frac / dw_frac if dw_frac else float("nan")))
    if verbose:
        print("\nRESULT 3: every ordered pair of the usable metals")
        print(f"  {'pair':<10}{'kind':<11}{'N(flat)':>10}{'N min':>10}{'N max':>10}"
              f"{'pair Δ%':>10}{'dW Δ%':>9}{'ampl.':>8}")
        for mA, mB, ss, flat, lo, hi, pf, df, amp in rows:
            print(f"  {mA + '/' + mB:<10}{'same-sign' if ss else 'straddling':<11}"
                  f"{flat:>10.4f}{lo:>10.4f}{hi:>10.4f}"
                  f"{pf*100:>9.2f}%{df*100:>8.2f}%{amp:>8.2f}")
        ss_amp = [r[8] for r in rows if r[2]]
        st_amp = [r[8] for r in rows if not r[2]]
        print(f"  same-sign amplification  : "
              f"{min(ss_amp):.2f}-{max(ss_amp):.2f}" if ss_amp else
              "  same-sign amplification  : (no same-sign pair admitted)")
        print(f"  straddling amplification : "
              f"{min(st_amp):.2f}-{max(st_amp):.2f}" if st_amp else
              "  straddling amplification : (none)")
        if ss_amp and st_amp:
            ok = min(ss_amp) > 1.0 and max(st_amp) <= min(ss_amp)
            print(f"  P3 (same-sign amplified, straddling not): "
                  f"{'PASS' if ok else 'FAIL'}")
        else:
            print("  P3: NOT TESTABLE on the metals this model admits -- see note")
    return rows


def make_plot(path="per_metal_crossover.png"):
    ells = np.linspace(*ELL_RANGE, 200)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))

    ax = axes[0]
    dd = np.linspace(1.9e-10, 3.8e-10, 400)
    for e, ls in [(0.3e-10, ":"), (0.6e-10, "-"), (1.5e-10, "--")]:
        ax.plot(dd * 1e10, W_GRAPHENE + delta_c(dd, e), ls,
                label=f"$\\ell$ = {e*1e10:.1f} $\\AA$")
    ax.axhline(W_CROSS_CHEM, color="k", lw=1, alpha=0.6,
               label="flat 5.4 eV (repo, 2026-09-18)")
    ax.axhline(W_MAX_ELEMENTAL, color="crimson", lw=1, ls="-.",
               label=f"highest elemental $W$ ({W_MAX_ELEMENTAL} eV)")
    ax.axvspan(1.9, D0_SEPARATION * 1e10, color="crimson", alpha=0.08)
    ax.text(2.13, 9.2, "chemisorbed:\n(N1) out of regime\n(RESULT 2)",
            color="crimson", fontsize=8, ha="center")
    for m in METAL_WORK_FUNCTIONS:
        d = D_EQ.get(m)
        if d is None:
            continue
        wc, st = w_cross_for_metal(m)
        ax.plot(d * 1e10, wc, "o" if st == "ok" else "x",
                color="tab:blue" if st == "ok" else "crimson", ms=7, zorder=5)
        ax.annotate(m, (d * 1e10, wc), textcoords="offset points",
                    xytext=(6, 4), fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("metal-graphene separation $d$ ($\\AA$)")
    ax.set_ylabel("$w_{cross} = W_G + \\Delta_c(d)$  (eV)")
    ax.set_title("Per-metal crossover, and where the anchored\nexponential stops being physical")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ms = usable_metals()
    for mA in ms:
        for mB in ms:
            if mA == mB:
                continue
            flat = net_response_dw(signed_offset(METAL_WORK_FUNCTIONS[mA]),
                                   signed_offset(METAL_WORK_FUNCTIONS[mB]),
                                   e_applied=0.0)
            vals = np.array([net_response_per_metal(mA, mB, e) for e in ells])
            ss = (signed_offset(METAL_WORK_FUNCTIONS[mA]) *
                  signed_offset(METAL_WORK_FUNCTIONS[mB])) > 0
            ax.plot(ells * 1e10, 100 * (vals - flat) / abs(flat),
                    "-" if ss else "--",
                    label=f"{mA}/{mB} ({'same-sign' if ss else 'straddling'})")
    ax.axhline(0, color="k", lw=1, alpha=0.5)
    ax.set_xlabel("$\\Delta_c$ decay length $\\ell$ ($\\AA$)")
    ax.set_ylabel("change in zero-bias net response vs flat 5.4 eV (%)")
    ax.set_title("Same-sign pairs amplify the crossover shift;\nstraddling pairs do not")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


if __name__ == "__main__":
    print("=" * 74)
    print("PER-METAL p/n CROSSOVER FROM THE CHEMICAL INTERFACE TERM")
    print("=" * 74)
    validate_reduction_to_signed_model()
    validate_flat_limit()
    validate_anchor()
    validate_symmetric_and_conjugation()
    result_1_per_metal_offsets()
    result_2_chemisorbed_reductio()
    result_3_pairs()
    print(f"\nplot -> {make_plot()}")
