"""
graphene_photodetector_collection_model.py

Spatially resolved, metal-dependent photocarrier collection-efficiency
model for the graphene photodetector, closing the Section 6.6 follow-on
item flagged in `thesis_draft/06-graphene-photodetectors.md` since
2026-08-24: "Replace the lumped EQE_bare = 0.15% parameter with a
spatially resolved diffusion-length collection model once the Chapter 4
spatially-resolved contact-doping-profile follow-on... is implemented."
That follow-on, `graphene_contact_doping_model.py`, has existed since
2026-08-26; this module is the integration that was always planned but
not yet done.

Physical picture (full literature review in
notes/2026-09-07-spatial-photocarrier-collection-model.md):

  Graphene's photocarrier lifetime is extremely short (~1 ps, Mueller,
  Xia & Avouris, Nature Photonics 4, 297 (2010)), so only photocarriers
  generated close enough to a contact to be swept out by a locally
  strong field survive to be collected. `graphene_photodetector_model.py`
  represents this "collection bottleneck" with a single literature
  midpoint, EQE_bare = 0.15%, applied identically regardless of contact
  metal. But the strength of the near-contact field is exactly what
  `graphene_contact_doping_model.py` already computes (from the same
  metal-work-function/quantum-capacitance physics used for Chapter 4's
  contact-resistance work) -- it has just never been fed into the
  photodetector side until now.

  This module superposes that contact-doping-induced field,
  E_doping(x) = |W_metal - W_graphene| / lambda_decay / (1+x/lambda_decay)^2,
  on the existing bare-device model's uniform bias field E_bias =
  V_bias/L_channel, and computes, for each candidate contact metal, the
  probability that a carrier generated at position x drifts to the
  contact before the ~1 ps lifetime elapses:

      v(x)    = mu * (E_bias + E_doping(x))
      t(x)    = integral_0^x dx'/v(x')          (transit time to contact)
      p(x)    = exp(-t(x)/tau_carrier)          (survival/collection probability)

  Averaging p(x) over the channel and comparing against a bias-field-only
  baseline (equivalent to a work-function-matched, zero-doping contact --
  already present in METAL_WORK_FUNCTIONS as Cr, W=4.50 eV vs. graphene's
  4.5 eV) gives a metal-dependent collection-efficiency ENHANCEMENT
  FACTOR, which is multiplied onto EQE_bare to get a metal-resolved
  modeled EQE.

  Scope limitation (stated explicitly, not hidden): a real two-terminal
  device has two contacts whose doping fields point in mirror-image
  directions (Weiss & Duan, NPG Asia Materials 5, e74 (2013): identical
  contacts give "equal and opposing" potentials and zero net zero-bias
  photocurrent). This module models only the single contact where the
  doping field REINFORCES the bias-driven sweep direction; the opposing
  contact's partially-cancelling contribution is not modeled here (see
  notes file Section 3 and thesis_draft/06-graphene-photodetectors.md
  Section 6.6 for this as an explicit open item, not a resolved one).
"""

import numpy as np
import matplotlib.pyplot as plt

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS, W_GRAPHENE, lambda_decay
from graphene_fet_model import mu
from graphene_photodetector_model import (
    L_channel, V_bias, TAU_TRANSIT, EQE_BARE, LAMBDA_NM, responsivity_bare,
)

# Graphene photocarrier lifetime (Mueller, Xia & Avouris, Nature Photonics
# 4, 297 (2010)) -- reused from graphene_photodetector_model.py's own
# docstring/notes citation rather than re-deriving.
TAU_CARRIER = 1.0e-12  # s

E_BIAS = V_bias / L_channel  # V/m, uniform bias-driven drift field

# Reference ("baseline") metal: Cr, whose work function (4.50 eV) is
# within 0.01 eV of graphene's own (4.5 eV) in this thesis's
# METAL_WORK_FUNCTIONS table -- i.e. already an approximately
# zero-doping, bias-field-only control case, not a separately invented one.
BASELINE_METAL = "Cr"

N_POINTS = 1000


def doping_field_magnitude(x, work_function_metal):
    """
    |E_doping(x)|, from differentiating the same saturating profile
    f(x) = 1/(1+x/lambda_decay) used for the *density* profile in
    graphene_contact_doping_model.doping_profile(). Reusing one profile
    shape for both the density and field is a stated simplification (see
    notes/2026-09-07-*.md, Section 3) -- a fully self-consistent
    treatment would derive both from one Poisson-like relation.
    """
    dV_mag = abs(work_function_metal - W_GRAPHENE)  # V (1 eV/e = 1 V)
    return dV_mag / lambda_decay / (1.0 + x / lambda_decay) ** 2


def collection_probability_profile(x, work_function_metal=None, include_doping=True):
    """
    Returns (p(x), t(x)): the collection (survival-to-contact) probability
    and cumulative transit time, for a carrier generated at each position
    in the array x (m, measured from the contact at x=0), under the
    superposed bias + (optionally) contact-doping field.
    """
    E_total = np.full_like(x, E_BIAS)
    if include_doping:
        if work_function_metal is None:
            raise ValueError("work_function_metal required when include_doping=True")
        E_total = E_total + doping_field_magnitude(x, work_function_metal)

    v = mu * E_total
    inv_v = 1.0 / v
    # Cumulative trapezoidal integral of 1/v(x') from 0 to each x.
    dt = 0.5 * (inv_v[1:] + inv_v[:-1]) * np.diff(x)
    t = np.concatenate([[0.0], np.cumsum(dt)])
    p = np.exp(-t / TAU_CARRIER)
    return p, t


def mean_collection_efficiency(work_function_metal=None, include_doping=True,
                                L=L_channel, n_points=N_POINTS):
    """Channel-averaged collection probability, <p(x)> over x in [0, L]."""
    x = np.linspace(0.0, L, n_points)
    p, _ = collection_probability_profile(x, work_function_metal, include_doping)
    return np.trapezoid(p, x) / L


def collection_enhancement_factor(work_function_metal, L=L_channel, n_points=N_POINTS):
    """
    eta_collect(metal) / eta_collect(bias-only baseline). >1 means the
    contact-doping field measurably speeds up collection relative to a
    work-function-matched (zero-doping) contact.
    """
    eta_metal = mean_collection_efficiency(work_function_metal, include_doping=True,
                                            L=L, n_points=n_points)
    eta_baseline = mean_collection_efficiency(include_doping=False, L=L, n_points=n_points)
    return eta_metal / eta_baseline


def modeled_eqe_by_metal(L=L_channel, n_points=N_POINTS):
    """
    Dict: metal -> (work_function, enhancement_factor, EQE_model).
    EQE_model = EQE_BARE * enhancement_factor -- see module docstring and
    notes/2026-09-07-*.md Section 3 for what this can/cannot claim
    (relative metal-to-metal trend is the load-bearing result; the
    absolute EQE_bare = 0.15% this is scaled from is not attributed to
    any specific metal in the literature reviewed for Section 6.2).
    """
    results = {}
    for metal, wf in METAL_WORK_FUNCTIONS.items():
        enh = collection_enhancement_factor(wf, L=L, n_points=n_points)
        results[metal] = (wf, enh, EQE_BARE * enh)
    return results


def plot_collection_profiles_and_eqe():
    """
    Two-panel figure:
      Left:  collection probability p(x) vs. distance from contact, for
             each metal plus the bias-only baseline.
      Right: modeled EQE by metal (bar chart), vs. the literature
             EQE_bare = 0.15% midpoint this thesis has used since Section 6.4.
    """
    x = np.linspace(0.0, L_channel, N_POINTS)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    p_baseline, _ = collection_probability_profile(x, include_doping=False)
    ax.plot(x * 1e9, p_baseline, 'k--', linewidth=2,
            label=f'bias-only baseline ($\\approx$ {BASELINE_METAL} contact)')
    for metal, wf in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        p, _ = collection_probability_profile(x, wf, include_doping=True)
        ax.plot(x * 1e9, p, linewidth=1.8, label=f'{metal} (W={wf:.2f} eV)')
    ax.set_xlabel('Distance from contact, x (nm)')
    ax.set_ylabel('Collection (survival) probability $p(x)$')
    ax.set_title('Photocarrier collection probability vs. contact metal\n'
                 f'($\\tau_{{carrier}}$={TAU_CARRIER*1e12:.1f} ps, '
                 f'$E_{{bias}}$={E_BIAS:.1e} V/m)')
    ax.legend(fontsize=7, loc='lower right')
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    results = modeled_eqe_by_metal()
    metals_sorted = sorted(results.keys(), key=lambda m: results[m][0])
    eqe_pct = [results[m][2] * 100 for m in metals_sorted]
    wfs = [results[m][0] for m in metals_sorted]
    colors = ['darkorange' if m == BASELINE_METAL else 'steelblue' for m in metals_sorted]
    x_pos = np.arange(len(metals_sorted))
    ax2.bar(x_pos, eqe_pct, color=colors)
    ax2.axhline(EQE_BARE * 100, color='gray', linestyle=':',
                label=f'literature EQE_bare midpoint ({EQE_BARE*100:.2f}%, Section 6.2)')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f'{m}\n(W={wfs[i]:.2f} eV)' for i, m in enumerate(metals_sorted)])
    ax2.set_ylabel('Modeled EQE (%)')
    ax2.set_title('Metal-resolved modeled EQE\n(this session\'s collection-enhancement factor)')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('photodetector_collection_efficiency_by_metal.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return results


def summary_numbers():
    print(f"E_bias = V_bias/L_channel = {E_BIAS:.3e} V/m "
          f"(V_bias={V_bias} V, L_channel={L_channel*1e9:.0f} nm, from graphene_fet_model.py)")
    print(f"tau_carrier (Mueller, Xia & Avouris 2010) = {TAU_CARRIER:.2e} s")
    print(f"TAU_TRANSIT (bias-field-only, graphene_photodetector_model.py) = {TAU_TRANSIT:.3e} s")
    eta_baseline = mean_collection_efficiency(include_doping=False)
    print(f"Baseline (bias-only, ~{BASELINE_METAL} contact) mean collection "
          f"efficiency = {eta_baseline:.4f}\n")

    print(f"{'Metal':6s} {'W (eV)':8s} {'|dV| (V)':10s} {'Enhancement':12s} {'EQE_model':10s}")
    results = modeled_eqe_by_metal()
    for metal, (wf, enh, eqe) in sorted(results.items(), key=lambda kv: kv[1][0]):
        print(f"{metal:6s} {wf:<8.2f} {abs(wf - W_GRAPHENE):<10.2f} {enh:<12.3f} {eqe*100:<10.4f}%")

    R_bare = responsivity_bare(LAMBDA_NM)
    best_metal = max(results, key=lambda m: results[m][1])
    wf_best, enh_best, eqe_best = results[best_metal]
    print(f"\nAt {LAMBDA_NM:.0f} nm: R_bare (literature EQE_bare) = {R_bare*1e3:.3f} mA/W; "
          f"best-modeled metal ({best_metal}, enhancement {enh_best:.2f}x) implies "
          f"R = {R_bare*1e3*enh_best:.3f} mA/W")
    print("\nCaveat: this is a single-contact (reinforcing-field-only) model; the "
          "opposing contact's partially-cancelling contribution is not modeled here "
          "(see notes/2026-09-07-*.md, Section 3, and thesis_draft/06-*.md Section 6.6).")


if __name__ == '__main__':
    print("Generating spatially resolved, metal-dependent photocarrier collection model...")
    summary_numbers()
    plot_collection_profiles_and_eqe()
    print("\nDone. Saved: photodetector_collection_efficiency_by_metal.png")
