"""
graphene_contact_doping_model.py

Spatially-resolved model of metal-work-function-dependent doping under and
near a graphene contact, closing the gap flagged repeatedly in
AUTOMATION_LOG.md since 2026-08-21: the existing contact-resistance
treatment (`graphene_fet_model.py`'s `Rc_total`) is a single lumped,
literature-calibrated series resistance and does not capture *why* contact
resistance is metal-specific, nor the fact that the contact locally dopes
the channel over a finite spatial extent (forming an in-plane p-n / p-p' /
n-n' junction), rather than acting as a resistance concentrated exactly at
the contact edge.

Physics summary (full literature review in
notes/2026-08-26-contact-induced-doping-profile.md):

  1. Work-function mismatch between the metal and graphene (~4.5 eV
     intrinsic) drives charge transfer at the interface. Because
     graphene's density of states is low near the Dirac point, the
     resulting Fermi-level shift is *not* the full ΔWF -- it is limited
     by graphene's own quantum capacitance, exactly as in the gated-
     channel case already modeled in graphene_fet_model.py. We reuse that
     module's quantum_capacitance() function here, replacing the ~90 nm
     SiO2 back-gate oxide capacitance with a much larger effective
     "interface capacitance" appropriate for a direct/weakly-bonded
     metal-graphene contact.
  2. The induced doping is not confined under the metal: Khomyakov et al.
     (Phys. Rev. B 82, 115437 (2010), arXiv:0911.2027) show the induced
     potential decays with distance x from the contact edge as x^-1 for
     doped graphene, extending hundreds of nm into the channel. We
     represent this with a saturating profile f(x) = 1/(1+x/lambda) that
     reduces to the literature x^-1 falloff for x >> lambda.
  3. Most common contact metals (Ti, Cr, Cu, Ag, Al) n-type dope
     graphene; only the highest-work-function metals (Au, Pt, and
     variably Pd) p-type dope it, because the literature n/p crossover
     work function (~5.4 eV) is higher than graphene's own work function
     (~4.5 eV) (Khomyakov et al. 2010).

This is a compact analytic model in the same spirit as graphene_fet_model.py
and rf_small_signal_model.py -- not a TCAD-grade electrostatic solver -- and
is intended to isolate the *additional* resistance contribution from the
contact-doping junction, on top of (not replacing) the lumped Rc_total
already used elsewhere in this thesis's device models.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, epsilon_0

from graphene_fet_model import quantum_capacitance, sheet_conductivity, C_ox, W, mu, n_puddle

# ---------------------------------------------------------------------------
# Metal work functions (eV) -- representative literature values for common
# graphene contact metals (see notes/2026-08-26-*.md, Section 1, for
# discussion of the spread and interface-chemistry caveats).
# ---------------------------------------------------------------------------
METAL_WORK_FUNCTIONS = {
    "Ti": 4.33,   # low-WF, n-type dopant, common adhesion-layer/contact metal
    "Cr": 4.50,   # ~graphene's own WF; weak/near-neutral doping
    "Cu": 4.65,   # common interconnect/contact metal; n-type dopant
    "Pd": 5.12,   # widely used for its low contact resistance (2026-08-21 notes)
    "Au": 5.10,   # near the reported n/p crossover; can go either way depending
                  # on interface bonding, commonly reported p-type in practice
    "Pt": 5.65,   # clearly above the ~5.4 eV crossover -> p-type
}

W_GRAPHENE = 4.5  # eV, intrinsic graphene work function (Khomyakov et al. 2010)

# Effective interface separation for a direct/weakly-bonded metal-graphene
# contact (sub-nm, vs. the 90 nm SiO2 back-gate dielectric in
# graphene_fet_model.py) -- physical vdW-type gap, ~3-4 Angstrom, consistent
# with reported equilibrium graphene-metal spacings for weakly-bonded metals
# (Au, Cu, Ag, Pt) in the DFT literature underlying Khomyakov et al. (2010).
d_interface = 3.5e-10  # m
eps_r_interface = 1.0  # vacuum-like gap; a simplification (no interface dipole layer)
C_interface = eps_r_interface * epsilon_0 / d_interface  # F/m^2

# Spatial decay length of the contact-induced doping profile (order "few
# hundred nm", per Khomyakov et al.'s reported doping extent).
lambda_decay = 250e-9  # m


def contact_edge_carrier_density(work_function_metal, T=300.0):
    """
    Self-consistent carrier density right at the contact edge (x=0), using
    the same quantum-capacitance series-combination approach as
    graphene_fet_model.carrier_density(), but with C_interface (large,
    direct-contact capacitance) in place of C_ox (back-gate oxide
    capacitance). Returns a signed density: positive = electron-doped
    (n-type), negative = hole-doped (p-type), consistent with
    "electrons flow from the lower-work-function material into the
    higher-work-function material."
    """
    dWF = work_function_metal - W_GRAPHENE  # eV
    # Sign convention: if the metal's work function is higher than graphene's,
    # electrons flow from graphene to the metal, leaving graphene hole-doped
    # (p-type, negative n). dWF > 0 -> p-type -> negative sign.
    dV = -dWF  # effective "gate voltage" equivalent, in volts (1 eV/e = 1 V)

    n_electrostatic = C_interface * dV / e
    C_q = quantum_capacitance(dV, T=T)
    series_factor = C_q / (C_q + C_interface)

    n_signed = n_electrostatic * series_factor
    return n_signed


def doping_profile(x, n_contact, n_bulk):
    """
    Spatial doping profile n(x) interpolating between the contact-edge
    density n_contact (at x=0) and the far-channel bulk density n_bulk
    (x -> infinity), using the saturating power-law decay
    f(x) = 1/(1+x/lambda) motivated by the Khomyakov et al. x^-1 asymptotic
    result for doped graphene.
    """
    f = 1.0 / (1.0 + x / lambda_decay)
    return n_bulk + (n_contact - n_bulk) * f


def _sheet_conductivity_puddle_regularized(n):
    """
    sheet_conductivity() as defined in graphene_fet_model.py goes to exactly
    zero when n=0, which is unphysical (real graphene has a residual
    disorder-induced "puddle" density n_puddle that regularizes the minimum
    conductivity point rather than letting it diverge to infinite sheet
    resistance -- the same regularization graphene_fet_model.carrier_density()
    already applies for the gate-swept case). Because our doping profile
    necessarily crosses through n=0 for any contact/bulk combination with
    opposite doping type (a genuine p-n junction), we apply the same
    regularization here so the local sheet resistance stays finite at the
    charge-neutrality crossing instead of producing a spurious 1/0 spike.
    """
    n_regularized = np.sqrt(np.asarray(n, dtype=float)**2 + n_puddle**2)
    return sheet_conductivity(n_regularized)


def junction_extra_resistance(work_function_metal, n_bulk, L_junction=1e-6, n_points=400):
    """
    Extra (above-bulk) sheet-resistance contribution of the contact-doping
    junction region, isolated from the lumped Rc_total already used in
    graphene_fet_model.py.

    Integrates 1/sigma_sheet(n(x)) - 1/sigma_sheet(n_bulk) over
    x in [0, L_junction], width-normalized (Ohm.um), so it can be compared
    directly against the Rc_per_width_ohm_um literature values used
    elsewhere in this thesis. Both terms use the puddle-regularized sheet
    conductivity so the integral stays finite across p-n crossings.
    """
    n_contact = contact_edge_carrier_density(work_function_metal)
    x = np.linspace(0, L_junction, n_points)
    n_x = doping_profile(x, n_contact, n_bulk)

    sigma_x = _sheet_conductivity_puddle_regularized(n_x)
    sigma_bulk = _sheet_conductivity_puddle_regularized(n_bulk)

    # Resistance per unit width of a segment dx of sheet resistance 1/sigma:
    # dR*W = dx / sigma(x). Extra resistance relative to a bulk-doped
    # channel of the same length.
    dR_extra_per_width = (1.0 / sigma_x) - (1.0 / sigma_bulk)
    R_extra_ohm_um = np.trapezoid(dR_extra_per_width, x) * 1e6  # Ohm (per um width)

    return R_extra_ohm_um, x, n_x, n_contact


def plot_doping_profiles_and_resistance():
    """Two-panel figure: (1) spatial doping profile n(x) for several
    contact metals against a representative n-type gated bulk channel, and
    (2) the resulting extra junction resistance vs. metal work function,
    compared against the lumped Rc range already used in this thesis."""
    n_bulk = 2.0e16  # 1/m^2, representative gate-induced n-type bulk channel
                      # density (order of the on-state density used in
                      # graphene_fet_model.py at moderate gate overdrive)

    x_plot = np.linspace(1e-9, 1.5e-6, 500)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    results = {}
    for metal, wf in METAL_WORK_FUNCTIONS.items():
        n_contact = contact_edge_carrier_density(wf)
        n_x = doping_profile(x_plot, n_contact, n_bulk)
        # Signed density for plotting (electron-doped positive convention
        # already baked into contact_edge_carrier_density / n_bulk).
        axes[0].plot(x_plot * 1e9, n_x * 1e-4, label=f"{metal} (W={wf:.2f} eV)",
                     linewidth=2)
        R_extra, _, _, _ = junction_extra_resistance(wf, n_bulk)
        results[metal] = (wf, R_extra, n_contact)

    axes[0].axhline(n_bulk * 1e-4, color='gray', linestyle=':',
                     label='bulk (gated) channel density')
    axes[0].axhline(0, color='black', linewidth=0.7)
    axes[0].set_xlabel('Distance from contact edge, x (nm)')
    axes[0].set_ylabel(r'Sheet carrier density $n(x)$ (cm$^{-2}$)')
    axes[0].set_title('Contact-induced doping profile vs. metal (n-type bulk channel)')
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)

    metals_sorted = sorted(results.keys(), key=lambda m: results[m][0])
    wfs = [results[m][0] for m in metals_sorted]
    R_extras = [results[m][1] for m in metals_sorted]
    colors = ['crimson' if r < 0 else 'steelblue' for r in R_extras]
    axes[1].bar(metals_sorted, R_extras, color=colors)
    axes[1].axhline(0, color='black', linewidth=0.8)
    axes[1].axhspan(110, 500, color='gray', alpha=0.15,
                     label='lumped Rc literature range\n(110-500 $\\Omega\\cdot\\mu$m,\n2026-08-21 notes)')
    axes[1].set_ylabel(r'Extra junction resistance ($\Omega\cdot\mu$m)')
    axes[1].set_title(f'Contact-doping junction resistance vs. metal\n'
                       f'(1 $\\mu$m junction length, $n_{{bulk}}$={n_bulk*1e-4:.1e} cm$^{{-2}}$)')
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('contact_doping_profile.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return results


def summary_numbers():
    """Print representative numbers for sanity-checking / logging."""
    n_bulk = 2.0e16
    print(f"Effective interface capacitance C_interface = {C_interface:.3e} F/m^2 "
          f"({C_interface / C_ox:.1f}x C_ox)")
    print(f"Bulk (gated) channel density used for comparison: {n_bulk*1e-4:.2e} cm^-2\n")
    print(f"{'Metal':6s} {'W (eV)':8s} {'n_contact (cm^-2)':20s} {'Type':6s} "
          f"{'R_extra (Ohm.um)':18s}")
    negative_R_metals = []
    for metal, wf in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        n_contact = contact_edge_carrier_density(wf)
        R_extra, _, _, _ = junction_extra_resistance(wf, n_bulk)
        dtype = "n-type" if n_contact >= 0 else "p-type"
        print(f"{metal:6s} {wf:<8.2f} {n_contact*1e-4:<20.3e} {dtype:6s} {R_extra:<18.1f}")
        if R_extra < 0:
            negative_R_metals.append(metal)

    if negative_R_metals:
        print(f"\nCaveat: {', '.join(negative_R_metals)} show *negative* R_extra "
              f"(i.e. the junction region integrates to lower resistance than the\n"
              f"bulk channel alone). This is not a bug: when |n_contact| >> n_bulk\n"
              f"(a very strongly doped contact against a lightly gated bulk channel),\n"
              f"most of the {lambda_decay*1e9:.0f} nm decay length is far better doped than the bulk\n"
              f"channel, so it locally conducts *better* than bulk over most of its\n"
              f"extent; only the narrow charge-neutrality crossing region is worse.\n"
              f"With the crossing spike now finite (puddle-regularized, see\n"
              f"_sheet_conductivity_puddle_regularized) rather than divergent, the\n"
              f"low-resistance well-doped stretch can dominate the net integral. This\n"
              f"is a real feature of the classical drift-sheet-conductance picture\n"
              f"used here, not a claim that the true device gets *easier* to contact\n"
              f"with such metals -- the model omits depletion-region/injection physics\n"
              f"right at the crossing itself, which would add resistance back in a way\n"
              f"this compact model cannot capture. See notes/2026-08-26-*.md, Section 4.")


if __name__ == '__main__':
    print("Generating contact-induced doping profile model...")
    summary_numbers()
    plot_doping_profiles_and_resistance()
    print("\nDone. Saved: contact_doping_profile.png")
