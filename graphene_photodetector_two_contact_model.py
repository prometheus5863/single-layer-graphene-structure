"""
graphene_photodetector_two_contact_model.py

Self-consistent TWO-CONTACT photocarrier collection model, closing the
scope limitation that `graphene_photodetector_collection_model.py`
(2026-09-07) flagged explicitly as its own first open item and that
AUTOMATION_LOG.md has carried at the top of "not yet covered" since:

    "Self-consistent two-contact collection model (Section 6.6's single
     reinforcing-contact scope limitation -- the natural next step for
     this session's new model)"

WHAT THE PREVIOUS MODEL DID, AND WHY IT OVERSTATED THE DEVICE
-------------------------------------------------------------
`graphene_photodetector_collection_model.py` superposes ONE contact's
doping field on the uniform bias field and integrates the drift transit
time to that single contact at x=0. Every carrier is assumed to be
collected there. That is a legitimate single-junction calculation, and
its metal-dependent enhancement factors (1.00x for Cr up to 1.47x for
Pt) are correct *for one junction in isolation*.

But a real two-terminal metal-graphene-metal photodetector has a second
contact at x=L whose doping field points the other way, and the
literature is unambiguous that this matters at first order rather than
as a correction:

  Weiss & Duan, NPG Asia Materials 5, e74 (2013) -- "symmetric
  metal-graphene-metal devices generate an equal positive and negative
  flow with a net zero photocurrent"; identical contacts create
  symmetric Schottky barriers whose opposing fields cancel. "Using
  metals with asymmetric band structures breaks this equilibrium."
  https://www.nature.com/articles/am201364

  Shimomura, Imai, Nakagawa, Kawai, Hashimoto, Ideguchi & Maki,
  Carbon Trends 5, 100100 (2021) -- in symmetric
  two-electrode graphene photodetectors "the polarities of the
  photovoltages at each graphene/electrode interface ... are canceled
  out under macroscopic light irradiation"; the paper's whole device
  strategy (a shadow mask over one interface, or unequal contact areas
  via a comb-shaped electrode) exists to break that cancellation.
  https://www.sciencedirect.com/science/article/pii/S2667056921000778

So the single-contact model does not merely omit a small term: under
uniform illumination and zero bias it predicts a large collection
efficiency where the real symmetric device produces *zero* net
photocurrent. This module supersedes it for any two-terminal
prediction. See notes/2026-09-17-two-contact-self-consistent-collection.md.

THE MODEL
---------
Contact A (metal A) at x=0, contact B (metal B) at x=L. Sign convention:
a POSITIVE total field F(x) drifts a photocarrier toward contact A.

    F(x) = E_bias + g_A(x) - g_B(x)

    g_A(x) = |W_A - W_gr| / lam / (1 + x/lam)^2          (sweeps toward A)
    g_B(x) = |W_B - W_gr| / lam / (1 + (L-x)/lam)^2      (sweeps toward B)

g_A is exactly the `doping_field_magnitude()` of the single-contact
module (same saturating profile, same lambda_decay, same
METAL_WORK_FUNCTIONS table); g_B is its mirror image about x=L/2. The
minus sign on g_B is the entire physical content that was missing
before.

A carrier generated at x0 then drifts along F. Because F can change
sign inside the channel, the model does NOT assume a destination:

  * F(x0) > 0 -> heads toward A. Collected at A (signed +1) only if F
    stays positive on all of [0, x0]; the transit time is
    t = integral_0^{x0} dx / (mu F(x)).
  * F(x0) < 0 -> heads toward B, signed -1, t = integral over [x0, L].
  * If F changes sign along the path, the carrier drifts into a
    STAGNATION POINT (F=0) instead of a contact, and is counted as not
    collected. This is a real feature of the two-contact geometry, not a
    numerical guard: with two opposing doping fields and a small bias
    there is generically an interior null.

Survival to the contact is the same exp(-t/tau) as before, with the same
tau = 1 ps photocarrier lifetime (Mueller, Xia & Avouris, Nature
Photonics 4, 297 (2010)). The reported figure of merit is the SIGNED,
channel-averaged net response

    N(A,B) = (1/L) integral_0^L s(x) exp(-t(x)/tau) dx,   s in {+1,-1,0}

which is directly comparable to the single-contact module's
`mean_collection_efficiency()` -- and reduces to it exactly in the
overlapping limit (see `validate_against_single_contact()` below, run as
part of __main__).

STATED SIMPLIFICATIONS (not hidden)
-----------------------------------
1. Both doping fields are taken as sweeping carriers toward their own
   contact, via |W_metal - W_gr|, which is the same magnitude convention
   the single-contact module already uses. This is what makes identical
   contacts cancel exactly. A fully signed treatment would distinguish
   n-type (W_metal < W_gr, e.g. Ti, Cu) from p-type (W_metal > W_gr,
   e.g. Pt, Pd) contacts and track electrons and holes separately; a
   p-n pair (e.g. Ti/Pt) would then ADD rather than partially cancel for
   one carrier species. That is a genuine extension, not a refinement,
   and is listed as the next open item rather than approximated here.
2. Uniform illumination (every x generates equally). CLOSED 2026-09-19 by
   graphene_photodetector_nonuniform_illumination_model.py. The Shimomura et al.
   shadow-mask device is precisely a non-uniform-generation experiment,
   so a g(x) weight is the natural way to model it; not done here.
3. Drift only, no diffusion; single effective mobility; no photogain,
   photo-thermoelectric or bolometric contribution (the bolometric
   mechanism remains unmodeled anywhere in this repo).
"""

import numpy as np
import matplotlib.pyplot as plt

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS, W_GRAPHENE, lambda_decay
from graphene_fet_model import mu
from graphene_photodetector_model import L_channel, V_bias
from graphene_photodetector_collection_model import (
    TAU_CARRIER, E_BIAS, mean_collection_efficiency,
)

N_POINTS = 4001  # odd, so x = L/2 is on the grid (matters for symmetry tests)


def _g(x, work_function, lam=lambda_decay):
    """Doping-field magnitude at distance x from a contact of given WF."""
    dV = abs(work_function - W_GRAPHENE)  # V (1 eV/e = 1 V)
    return dV / lam / (1.0 + x / lam) ** 2


def signed_field(x, W_A, W_B, E_bias=E_BIAS, L=L_channel):
    """
    F(x): total field, positive = drifts carriers toward contact A (x=0).
    Contact A sits at x=0, contact B at x=L.
    """
    return E_bias + _g(x, W_A) - _g(L - x, W_B)


def net_response(W_A, W_B, E_bias=E_BIAS, L=L_channel, n_points=N_POINTS,
                 return_profiles=False):
    """
    Signed, channel-averaged net collection N(A,B) as defined in the
    module docstring. Positive means net carrier flow to contact A.

    Returns N, or (N, dict of profiles) if return_profiles.
    """
    x = np.linspace(0.0, L, n_points)
    F = signed_field(x, W_A, W_B, E_bias=E_bias, L=L)

    speed = mu * np.abs(F)
    # Guard only against a true zero at a grid point (infinite transit
    # time); the stagnation logic below is what physically handles nulls.
    with np.errstate(divide="ignore"):
        inv_v = np.where(speed > 0.0, 1.0 / np.maximum(speed, 1e-300), np.inf)

    dx = np.diff(x)
    seg = 0.5 * (inv_v[1:] + inv_v[:-1]) * dx
    # Two SEPARATE cumulative integrals, forward and reverse:
    #   cum_fwd[i] = integral from 0    to x[i] of dx/|v|
    #   cum_rev[i] = integral from x[i] to L    of dx/|v|
    # These must not be derived from one another (cum_rev[i] =
    # cum_fwd[-1] - cum_fwd[i] is WRONG here): 1/|v| diverges at a
    # stagnation point, so a cumulative sum that has crossed a null is
    # infinite from there on, and subtracting two such values gives
    # inf-inf for every point beyond the null. That bug silently made
    # every carrier on the far side of a null uncollectable, which broke
    # the symmetric-cancellation check below (it returned N = +0.46 for
    # Pt/Pt at zero bias instead of 0). Each direction is integrated only
    # over the interval the carrier actually traverses, which by
    # construction contains no null for collected carriers.
    cum_fwd = np.concatenate([[0.0], np.cumsum(seg)])
    cum_rev = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])

    # Running extrema used to detect a sign change between a generation
    # point and the boundary it is heading toward.
    min_from_left = np.minimum.accumulate(F)                  # min of F on [x0_0, x_i]
    max_from_right = np.maximum.accumulate(F[::-1])[::-1]     # max of F on [x_i, L]

    sign = np.zeros_like(x)
    t = np.full_like(x, np.inf)

    toward_A = F > 0.0
    toward_B = F < 0.0

    # Toward A: needs F > 0 on all of [0, x_i]
    ok_A = toward_A & (min_from_left > 0.0)
    sign[ok_A] = +1.0
    t[ok_A] = cum_fwd[ok_A]

    # Toward B: needs F < 0 on all of [x_i, L]
    ok_B = toward_B & (max_from_right < 0.0)
    sign[ok_B] = -1.0
    t[ok_B] = cum_rev[ok_B]

    p = np.where(np.isfinite(t), np.exp(-t / TAU_CARRIER), 0.0)
    integrand = sign * p
    N = float(np.trapezoid(integrand, x) / L)

    if return_profiles:
        stalled = ~(ok_A | ok_B)
        return N, {
            "x": x, "F": F, "sign": sign, "p": p,
            "collected_A_frac": float(np.trapezoid(np.where(ok_A, p, 0.0), x) / L),
            "collected_B_frac": float(np.trapezoid(np.where(ok_B, p, 0.0), x) / L),
            "stalled_frac": float(np.count_nonzero(stalled)) / len(x),
        }
    return N


# ---------------------------------------------------------------------
# Validation against the model this one supersedes
# ---------------------------------------------------------------------
def validate_against_single_contact(verbose=True):
    """
    In the overlapping limit -- contact B's doping field switched off by
    giving it graphene's own work function, so g_B == 0 -- this model must
    reproduce `graphene_photodetector_collection_model.
    mean_collection_efficiency()` exactly, because F = E_bias + g_A > 0
    everywhere, every carrier reaches contact A, and the signed average
    collapses to the old unsigned one.

    This is a real numerical check of agreement, not an assertion that
    the two models agree. Returns a list of per-metal rows and the worst
    relative deviation.
    """
    rows = []
    worst = 0.0
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        new = net_response(W, W_GRAPHENE)              # g_B = 0 exactly
        old = mean_collection_efficiency(W, include_doping=True)
        rel = abs(new - old) / old
        worst = max(worst, rel)
        rows.append((metal, W, old, new, rel))

    if verbose:
        print("Validation: two-contact model with contact B's doping field off")
        print("vs. the 2026-09-07 single-contact model (must agree)")
        print(f"{'metal':<6}{'W (eV)':>8}{'single':>12}{'two-contact':>14}{'rel. dev.':>12}")
        for metal, W, old, new, rel in rows:
            print(f"{metal:<6}{W:>8.2f}{old:>12.6f}{new:>14.6f}{rel:>12.2e}")
        print(f"worst relative deviation: {worst:.2e}")
    return rows, worst


def validate_symmetric_cancellation(verbose=True):
    """
    Physics check against the literature's central claim: at zero bias,
    two IDENTICAL contacts must give exactly zero net response, because
    F(x) = g_A(x) - g_B(x) is antisymmetric about x = L/2.

    Reported as a measured number, not asserted.
    """
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        N0 = net_response(W, W, E_bias=0.0)
        rows.append((metal, W, N0))
    worst = max(abs(n) for _, _, n in rows)
    if verbose:
        print("\nValidation: symmetric contacts at zero bias -> net response ~ 0")
        print("(Weiss & Duan 2013: 'equal positive and negative flow with a")
        print(" net zero photocurrent'; Shimomura et al. 2021: polarities 'canceled out')")
        for metal, W, N0 in rows:
            print(f"  {metal:<4} (W={W:.2f} eV) both contacts: N = {N0:+.3e}")
        print(f"worst |N| over symmetric pairs: {worst:.2e}")
    return rows, worst


# ---------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------
def symmetric_vs_single_contact_table(verbose=True):
    """
    The headline correction: for each metal, the single-contact model's
    collection efficiency at the working bias vs. the symmetric
    two-contact device's actual net response at the same bias.
    """
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        single = mean_collection_efficiency(W, include_doping=True)
        sym = net_response(W, W)
        rows.append((metal, W, single, sym, sym / single))
    if verbose:
        print(f"\nSymmetric two-contact device at V_bias = {V_bias} V "
              f"(E_bias = {E_BIAS:.3e} V/m)")
        print(f"{'metal':<6}{'W (eV)':>8}{'single-contact':>16}{'symmetric 2C':>15}{'ratio':>9}")
        for metal, W, single, sym, ratio in rows:
            print(f"{metal:<6}{W:>8.2f}{single:>16.6f}{sym:>15.6f}{ratio:>9.3f}")
    return rows


def asymmetric_pair_table(verbose=True, e_bias=0.0):
    """
    Zero-bias net response for every ordered metal pair (A, B). At zero
    bias a symmetric pair gives exactly zero, so any nonzero value here
    is purely the work-function asymmetry -- the Weiss & Duan mechanism.
    """
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    out = {}
    for mA, WA in metals:
        for mB, WB in metals:
            out[(mA, mB)] = net_response(WA, WB, E_bias=e_bias)
    if verbose:
        names = [m for m, _ in metals]
        print(f"\nZero-bias net response N(A,B), rows = contact A, cols = contact B")
        print("      " + "".join(f"{m:>9}" for m in names))
        for mA in names:
            print(f"{mA:<6}" + "".join(f"{out[(mA, mB)]:>9.4f}" for mB in names))
        best = max(out.items(), key=lambda kv: abs(kv[1]))
        print(f"largest |N| at zero bias: {best[0][0]}/{best[0][1]} -> {best[1]:+.4f}")
    return out


def make_plots(fname="photodetector_two_contact_net_response.png"):
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))

    # (a) signed field profiles: symmetric vs asymmetric, zero bias
    ax = axes[0]
    x = np.linspace(0.0, L_channel, N_POINTS)
    for (WA, WB, lab, st) in [
        (METAL_WORK_FUNCTIONS["Pt"], METAL_WORK_FUNCTIONS["Pt"], "Pt / Pt (symmetric)", "-"),
        (METAL_WORK_FUNCTIONS["Pt"], METAL_WORK_FUNCTIONS["Ti"], "Pt / Ti (asymmetric)", "--"),
        (METAL_WORK_FUNCTIONS["Pt"], METAL_WORK_FUNCTIONS["Cr"], "Pt / Cr (asymmetric)", ":"),
    ]:
        ax.plot(x * 1e9, signed_field(x, WA, WB, E_bias=0.0) / 1e6, st, label=lab)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("position x (nm), contact A at 0, contact B at L")
    ax.set_ylabel("signed field F(x)  (MV/m)\n[+ = sweeps toward contact A]")
    ax.set_title("(a) Two opposing doping fields, zero bias")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) single-contact vs symmetric two-contact, at working bias
    ax = axes[1]
    names = [m for m, _ in metals]
    single = [mean_collection_efficiency(W, include_doping=True) for _, W in metals]
    sym = [net_response(W, W) for _, W in metals]
    idx = np.arange(len(names))
    ax.bar(idx - 0.2, single, 0.4, label="single-contact model (2026-09-07)")
    ax.bar(idx + 0.2, sym, 0.4, label="symmetric two-contact (this model)")
    ax.set_xticks(idx); ax.set_xticklabels(names)
    ax.set_ylabel("net collection efficiency")
    ax.set_title(f"(b) The correction, at $V_{{bias}}$ = {V_bias} V")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    # (c) zero-bias net response vs work-function asymmetry
    ax = axes[2]
    pairs = asymmetric_pair_table(verbose=False, e_bias=0.0)
    dW, Nv = [], []
    for (mA, mB), N in pairs.items():
        dW.append(METAL_WORK_FUNCTIONS[mA] - METAL_WORK_FUNCTIONS[mB])
        Nv.append(N)
    ax.scatter(dW, Nv, s=22)
    ax.axhline(0.0, color="k", lw=0.8); ax.axvline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"$W_A - W_B$  (eV)")
    ax.set_ylabel("zero-bias net response N")
    ax.set_title("(c) Net response is driven by contact asymmetry")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"\nsaved {fname}")
    return fname


if __name__ == "__main__":
    print("=" * 72)
    print("Two-contact self-consistent photocarrier collection model")
    print("=" * 72)
    _, worst_v = validate_against_single_contact()
    _, worst_s = validate_symmetric_cancellation()
    symmetric_vs_single_contact_table()
    asymmetric_pair_table()
    make_plots()
    print("\nValidation summary:")
    print(f"  single-contact limit reproduced to {worst_v:.2e} relative")
    print(f"  symmetric zero-bias cancellation to |N| <= {worst_s:.2e}")
