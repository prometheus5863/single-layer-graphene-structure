"""
graphene_edge_contact_model.py

Contact-*geometry* model: edge contacts vs. top (surface) contacts to
graphene, and a patterned (hole-array) contact that interpolates between
them. Closes the gap flagged in every "not yet covered" list since
2026-08-31: "Contact-geometry dependence (top vs. edge contact) is not
represented at all in graphene_contact_doping_model.py, despite Section
4.7's own finding that it is a larger lever (~11x) than metal choice (~3x)
at fixed geometry." Full literature review in
notes/2026-09-05-edge-vs-top-contact-geometry.md.

Physics summary:

  1. Top (surface) contacts inject current across a weak, van-der-Waals-
     coupled 2D overlap area -- graphene's carbon atoms are sp2-bonded with
     no out-of-plane dangling bonds, so there is no covalent path for a
     top contact to bond into (Wang et al., Science 342, 614 (2013)).
  2. At a graphene *edge*, carbon atoms CAN form sigma bonds directly with
     a contacting metal, giving stronger chemical coupling and higher
     carrier transmission. Passi et al. (arXiv:1807.04772) quantified this
     with DFT for Au: the metal-induced Fermi-level shift is 0.35 eV at an
     edge vs. only 0.14 eV at the flat surface -- a ~2.5x stronger
     doping-induced shift right at the edge.
  3. This module converts those two Fermi-level shifts directly into
     contact-edge carrier densities via graphene's linear dispersion
     (n(E_F) = sign(E_F) * E_F^2 / (pi*(hbar*v_F)^2)), then reuses
     graphene_contact_doping_model.py's existing doping-profile-integration
     machinery (factored out this session as
     junction_extra_resistance_from_ncontact()) to get an "extra junction
     resistance" for pure-edge-mode and pure-top-mode injection on an
     equal footing.
  4. A patterned (hole-array) contact -- Passi et al.'s actual device,
     metal deposited over graphene with an array of etched holes beneath
     it -- is modeled as an areal mix of edge-mode and top-mode regions,
     using a purely geometric relation between hole diameter/areal fill
     fraction and the fraction of the remaining graphene that lies close
     enough to a hole edge (within the existing lambda_decay from
     graphene_contact_doping_model.py) for edge-mode doping fronts from
     neighboring holes to overlap, plus a simple 1/(1-f) current-
     constriction penalty for the etched-away area (same area-dilution
     idiom as the liner-aware copper model in graphene_interconnect_model.py).

As with the rest of this thesis's device-physics modules, this is a compact
analytic model, not a TCAD-grade simulation. Section 4 of the notes file is
explicit about where it is only expected to match Passi et al.'s data
qualitatively (the large-hole, area-loss-dominated branch) and where it is
NOT expected to (the small-hole branch of their non-monotonic trend).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, e

from graphene_fet_model import v_F
from graphene_contact_doping_model import (
    junction_extra_resistance,
    junction_extra_resistance_from_ncontact,
    lambda_decay,
    METAL_WORK_FUNCTIONS,
)

# ---------------------------------------------------------------------------
# DFT-reported Fermi-level shifts at Au-graphene contacts, edge vs. surface
# (Passi et al., arXiv:1807.04772, Section 3 discussion / DFT calculations).
# Sign: Au is p-type (work function 5.10 eV, above the ~5.4 eV n/p crossover
# is NOT satisfied strictly, but Au is reported/treated as p-type in
# practice throughout this repo's existing METAL_WORK_FUNCTIONS comment and
# in the literature Rc source for Au already cited there) -- both shifts
# taken as downward (hole-doping) to match that existing convention.
# ---------------------------------------------------------------------------
AU_FERMI_SHIFT_EDGE_EV = 0.35
AU_FERMI_SHIFT_SURFACE_EV = 0.14

# Passi et al., Table (TLM measurements), on-state (V_BG = -40 V) and
# Dirac-point contact resistance vs. hole diameter, fixed 12 um x 5 um pad.
PASSI_HOLE_DIAMETER_NM = np.array([0, 50, 100, 200, 500, 1000])
PASSI_RC_ON_STATE_OHM_UM = np.array([519, 212, 352, 45, 410, 560])
PASSI_RC_DIRAC_OHM_UM = np.array([1372, 620, 732, 456, 1354, 1590])

# Same on-state bulk channel density used throughout graphene_contact_
# doping_model.py's recalibration section, for direct comparability.
N_BULK_ON_STATE = 2.0e16  # 1/m^2


def carrier_density_from_fermi_shift(E_F_eV, p_type=True):
    """
    Convert a metal-induced Fermi-level shift (eV, magnitude) into a
    graphene sheet carrier density via the standard linear-dispersion
    relation n(E_F) = E_F^2 / (pi * (hbar*v_F)^2), signed negative
    (p-type/hole-doped) or positive (n-type/electron-doped) to match the
    sign convention already used by
    graphene_contact_doping_model.contact_edge_carrier_density() (positive
    = electron-doped).
    """
    E_F_joules = abs(E_F_eV) * e
    n_magnitude = E_F_joules**2 / (np.pi * (hbar * v_F)**2)
    return -n_magnitude if p_type else n_magnitude


def edge_vs_top_extra_resistance_au(n_bulk=N_BULK_ON_STATE):
    """
    Pure-edge-mode and pure-top-mode extra junction resistance for Au,
    using the DFT Fermi-shift-derived carrier densities above, both run
    through graphene_contact_doping_model's existing doping-profile
    integration. Returns (R_extra_top, R_extra_edge, n_top, n_edge), all in
    Ohm.um / 1/m^2.
    """
    n_top = carrier_density_from_fermi_shift(AU_FERMI_SHIFT_SURFACE_EV, p_type=True)
    n_edge = carrier_density_from_fermi_shift(AU_FERMI_SHIFT_EDGE_EV, p_type=True)

    R_extra_top, _, _ = junction_extra_resistance_from_ncontact(n_top, n_bulk)
    R_extra_edge, _, _ = junction_extra_resistance_from_ncontact(n_edge, n_bulk)

    return R_extra_top, R_extra_edge, n_top, n_edge


def edge_influenced_area_fraction(D_nm, f, lambda_decay_nm=None):
    """
    Geometric model (Section 4 of this session's notes): for circular holes
    of diameter D_nm arranged on a square lattice at areal fill fraction f,
    return the fraction p_edge of the *remaining* graphene area close
    enough to a hole edge (within 2*lambda_decay of it) that neighboring
    holes' doping fronts overlap and the local area is treated as fully
    edge-dominated.

    f = (pi/4)(D/a)^2  =>  pitch a = D * sqrt(pi/(4f))
    ligament width = a - D  (remaining graphene between adjacent hole edges)
    p_edge = min(1, 2*lambda_decay / ligament)  when ligament > 0, else 1
    (f -> (pi/4) is the close-packing limit where ligament -> 0)
    """
    if lambda_decay_nm is None:
        lambda_decay_nm = lambda_decay * 1e9

    D_nm = np.asarray(D_nm, dtype=float)
    f = np.asarray(f, dtype=float)

    close_packing_limit = np.pi / 4.0
    f_clipped = np.clip(f, 1e-9, close_packing_limit - 1e-6)

    pitch = D_nm * np.sqrt(np.pi / (4.0 * f_clipped))
    ligament = pitch - D_nm
    ligament = np.where(ligament > 1e-9, ligament, 1e-9)

    p_edge = np.minimum(1.0, 2.0 * lambda_decay_nm / ligament)
    # D=0 (no holes at all) is a degenerate case of the formula (ligament
    # undefined/zero) that should just mean "no edges introduced" (f=0).
    p_edge = np.where((D_nm <= 0) | (f <= 0), 0.0, p_edge)
    return p_edge


def patterned_contact_extra_resistance(D_nm, f, R_extra_top, R_extra_edge,
                                        n_bulk=N_BULK_ON_STATE):
    """
    Effective extra junction resistance (Ohm.um) for a patterned contact
    with hole diameter D_nm and areal fill fraction f, mixing R_extra_top
    and R_extra_edge by edge-influenced area fraction (conductance-weighted
    parallel combination), then applying a 1/(1-f) current-constriction
    penalty for the etched-away area. Vectorized over D_nm and/or f.
    """
    p_edge = edge_influenced_area_fraction(D_nm, f)
    f = np.asarray(f, dtype=float)

    # Conductance-weighted mix: 1/R_eff = (1-p)/R_top + p/R_edge.
    # Both R_extra_top and R_extra_edge for Au are positive in this model
    # (checked numerically below in summary_numbers(); this mixing rule
    # assumes that -- see notes Section 4 for the caveat if a future metal
    # gives a negative R_extra in this framework).
    inv_R_eff = (1.0 - p_edge) / R_extra_top + p_edge / R_extra_edge
    R_eff = 1.0 / inv_R_eff

    constriction_penalty = 1.0 / (1.0 - np.clip(f, 0.0, 0.999))
    return R_eff * constriction_penalty, p_edge


def plot_edge_vs_top_and_patterned():
    """
    Two-panel figure: (1) pure edge-mode vs. pure top-mode extra junction
    resistance for Au (bar comparison, plus the Wang et al. ~100 Ohm.um
    literature edge-contact benchmark as a horizontal reference line), and
    (2) predicted patterned-contact resistance vs. hole diameter at a few
    illustrative (explicitly labeled as assumed) fill fractions, with
    Passi et al.'s five measured on-state data points overlaid.
    """
    R_extra_top, R_extra_edge, n_top, n_edge = edge_vs_top_extra_resistance_au()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # --- Panel 1: pure edge vs. pure top ---
    labels = ['Top contact\n(surface, DFT-shift-derived)',
              'Edge contact\n(DFT-shift-derived, this model)']
    values = [R_extra_top, R_extra_edge]
    colors = ['steelblue', 'darkorange']
    axes[0].bar(labels, values, color=colors)
    axes[0].axhline(100, color='crimson', linestyle='--', linewidth=1.5,
                     label='Wang et al. 2013 true 1D edge contact\n(~100 $\\Omega\\cdot\\mu$m, different metal/geometry)')
    axes[0].axhline(519, color='gray', linestyle=':', linewidth=1.2,
                     label='Passi et al. measured Au top contact\n(519 $\\Omega\\cdot\\mu$m, on-state)')
    axes[0].set_ylabel(r'Extra junction resistance ($\Omega\cdot\mu$m)')
    axes[0].set_title('Au: pure edge-mode vs. pure top-mode\ncontact-doping-junction resistance')
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3, axis='y')

    # --- Panel 2: patterned contact vs. hole diameter, swept fill fraction ---
    D_sweep = np.linspace(10, 1200, 300)
    assumed_fs = [0.05, 0.15, 0.30]
    for f_val, style in zip(assumed_fs, ['-', '--', '-.']):
        R_patterned, _ = patterned_contact_extra_resistance(
            D_sweep, f_val, R_extra_top, R_extra_edge)
        axes[1].plot(D_sweep, R_patterned, style,
                     label=f'model, assumed f={f_val:.2f}', linewidth=2)

    axes[1].scatter(PASSI_HOLE_DIAMETER_NM, PASSI_RC_ON_STATE_OHM_UM,
                    color='black', zorder=5, marker='o', s=50,
                    label='Passi et al. measured (on-state),\nactual f not reported')
    axes[1].axhline(R_extra_top, color='gray', linestyle=':', linewidth=1,
                     label='pure top-mode limit (D=0)')
    axes[1].set_xlabel('Hole diameter D (nm)')
    axes[1].set_ylabel(r'Extra junction resistance ($\Omega\cdot\mu$m)')
    axes[1].set_title('Patterned (hole-array) contact: model vs. Passi et al. data\n'
                       '(qualitative comparison only -- see notes Section 4)')
    axes[1].legend(fontsize=7.5)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim(0, 700)

    plt.tight_layout()
    plt.savefig('edge_vs_top_contact.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return R_extra_top, R_extra_edge


def summary_numbers():
    R_extra_top, R_extra_edge, n_top, n_edge = edge_vs_top_extra_resistance_au()

    print("Edge-vs-top contact geometry model (Au, DFT-Fermi-shift-derived)")
    print("=" * 78)
    print(f"Surface Fermi shift:  {AU_FERMI_SHIFT_SURFACE_EV:.2f} eV  ->  "
          f"n_top  = {n_top*1e-4:.3e} cm^-2 (p-type)")
    print(f"Edge Fermi shift:     {AU_FERMI_SHIFT_EDGE_EV:.2f} eV  ->  "
          f"n_edge = {n_edge*1e-4:.3e} cm^-2 (p-type)")
    print(f"Ratio n_edge/n_top = {n_edge/n_top:.2f} "
          f"(Fermi-shift ratio was {AU_FERMI_SHIFT_EDGE_EV/AU_FERMI_SHIFT_SURFACE_EV:.2f}; "
          f"n ~ E_F^2 so density ratio is the shift ratio squared)")
    print()
    print(f"R_extra, pure top-mode  = {R_extra_top:8.1f} Ohm.um")
    print(f"R_extra, pure edge-mode = {R_extra_edge:8.1f} Ohm.um "
          f"({'lower' if R_extra_edge < R_extra_top else 'HIGHER'} than top-mode)")
    print(f"Ratio R_top/R_edge = {R_extra_top/R_extra_edge:.2f}x "
          f"(cf. Passi et al.'s measured *device* ratio 519/45 = {519/45:.1f}x, "
          "which also includes their contact's patterning/geometry, not\n"
          "just the edge-vs-top doping-density effect isolated here)")
    print()

    if R_extra_edge <= 0 or R_extra_top <= 0:
        print("CAVEAT: at least one of R_extra_top/R_extra_edge came out "
              "<= 0 -- the conductance-weighted parallel-mixing rule used in\n"
              "patterned_contact_extra_resistance() assumes both are positive "
              "and would need revisiting (see graphene_contact_doping_model.py's\n"
              "own documented caveat about negative R_extra for strongly-doped "
              "metals).")

    print("\nPatterned-contact model vs. Passi et al. data, at illustrative "
          "assumed fill fractions (their actual f is not reported):")
    print(f"{'D (nm)':8s} {'f=0.05':>10s} {'f=0.15':>10s} {'f=0.30':>10s} "
          f"{'measured':>10s}")
    for D, meas in zip(PASSI_HOLE_DIAMETER_NM, PASSI_RC_ON_STATE_OHM_UM):
        row = [D]
        for f_val in (0.05, 0.15, 0.30):
            if D == 0:
                row.append(R_extra_top)
            else:
                r, _ = patterned_contact_extra_resistance(D, f_val, R_extra_top, R_extra_edge)
                row.append(float(r))
        print(f"{row[0]:<8d} {row[1]:>10.1f} {row[2]:>10.1f} {row[3]:>10.1f} {meas:>10.1f}")

    print("\nQualitative check: does the model produce a minimum at "
          "intermediate D (like the real device's D=200nm optimum), rather\n"
          "than a monotonic trend? ", end="")
    r_check, _ = patterned_contact_extra_resistance(
        np.array([50.0, 200.0, 800.0]), 0.15, R_extra_top, R_extra_edge)
    is_nonmonotonic_min_at_middle = r_check[1] < r_check[0] and r_check[1] < r_check[2]
    print("YES" if is_nonmonotonic_min_at_middle else "NO",
          f"(D=50/200/800nm at f=0.15 -> "
          f"{r_check[0]:.0f}/{r_check[1]:.0f}/{r_check[2]:.0f} Ohm.um)")
    print("As flagged in the notes, this model captures the large-D "
          "(area-loss-dominated) branch of the real non-monotonic trend but\n"
          "does NOT reproduce the small-D upturn Passi et al. observed at "
          "50-100nm -- it predicts monotonic improvement as D shrinks at\n"
          "fixed f, since smaller holes at fixed fill fraction geometrically "
          "always pack more edge length per unit area in this idealization.")


if __name__ == '__main__':
    print("Generating edge-vs-top contact geometry model...")
    summary_numbers()
    plot_edge_vs_top_and_patterned()
    print("\nDone. Saved: edge_vs_top_contact.png")
