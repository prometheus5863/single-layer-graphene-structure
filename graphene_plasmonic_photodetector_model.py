"""
graphene_plasmonic_photodetector_model.py

Extends graphene_photodetector_model.py with a plasmonic near-field
absorption-enhancement factor, closing the "Plasmonic-absorption-
enhancement factor for the photodetector model (Chapter 6, Section 6.5)"
item that has appeared in every "not yet covered" list since
2026-08-24's photodetector session, and referenced directly in
graphene_photodetector_model.py's own module docstring as a planned
follow-on. See notes/2026-09-06-plasmonic-enhancement-graphene-
photodetectors.md for full citations and the reasoning behind every
number pulled from the literature below.

Physical picture
-----------------
A resonant metal nanostructure (grating finger, bowtie antenna, etc.)
placed on or near a graphene sheet locally concentrates the incident
optical field. Since a 2D sheet's absorbed power density scales with the
local field *intensity* (|E_local|^2), not the far-field incident
intensity, this raises graphene's effective single-pass absorption --
and, assuming the collection efficiency downstream of photocarrier
generation is unchanged (an idealization, flagged below), the external
quantum efficiency (EQE) and responsivity by the same factor -- in a
narrow spectral window around the structure's plasmon resonance.

This module models that enhancement as a Lorentzian intensity-
enhancement factor

    F(lambda) = 1 + (F_max - 1) / (1 + (2*Q*(lambda-lambda0)/lambda0)^2)

(F -> F_max exactly on resonance, F -> 1 far off resonance) and applies
it multiplicatively to graphene_photodetector_model.EQE_BARE, i.e.
EQE_plasmonic(lambda) = EQE_BARE * F(lambda), then reuses
graphene_photodetector_model.responsivity_bare() unmodified with this
wavelength-dependent EQE. This is a deliberately compact model (one
Lorentzian per design, two free parameters F_max and Q), not a full
electromagnetic simulation -- the point is to make the *spectral shape*
and *scale* of the previously purely-qualitative "plasmonic enhancement"
follow-on item quantitative and checkable against literature, not to
reproduce any one paper's full simulated spectrum.

Two designs are calibrated directly against literature near-field/
absorption enhancement numbers (not device-level responsivity
comparisons that bundle in other effects -- see the notes file's
Section 3 for why a third literature design, Fang et al.'s bowtie
nanoantenna, is deliberately *not* fit here even though it is cited
descriptively):

  1. Echtermeyer et al., Nat. Commun. 2, 458 (2011) (Ti/Au finger
     grating at a graphene p-n junction, visible): FEM near-field
     simulations there report an amplitude enhancement of ~5 for the
     optimal geometry/polarization, i.e. an intensity (power)
     enhancement of ~5^2 = 25 -- this is what is used as F_max here,
     since graphene absorption is a |E|^2 (power) effect, not an
     amplitude effect. Two grating pitches were reported with resonances
     at 514 nm (110 nm finger width) and 633 nm (130 nm finger width).
  2. An arrayed bowtie-antenna-on-waveguide graphene photodetector,
     arXiv:1808.10823 / ACS Photonics (2018-2019 telecom-band device):
     FEM simulation there reports an absorption enhancement factor of
     8.5x for a single bowtie element vs. bare monolayer graphene of
     equal length, operating in the telecom C-band (~1550 nm). This is
     a *directly* usable absorption-enhancement number (not a
     device-to-device responsivity comparison), so F_max = 8.5 is used
     as-is.

Neither paper reports the plasmon resonance's spectral quality factor Q
(FWHM) in the material available this session, so Q is an *assumed*
representative value for lossy metal (Au/Ti) nanostructures at each
wavelength range (Q=8 visible, Q=6 telecom -- broader, consistent with
five weakly-coupled array elements rather than one sharp single
resonator) -- flagged explicitly as an assumption, not a fitted or
measured value, following this repo's established practice (e.g. the
2026-09-05 edge-vs-top-contact session's explicitly-assumed hole-array
fill fractions).
"""

import numpy as np
import matplotlib.pyplot as plt

from graphene_photodetector_model import (
    EQE_BARE, TAU_TRANSIT, responsivity_bare, photoconductive_gain,
)

# ---------------------------------------------------------------------------
# Literature-calibrated plasmonic designs
# ---------------------------------------------------------------------------

PLASMONIC_DESIGNS = [
    {
        "label": "Echtermeyer Ti/Au grating, 110 nm finger (514 nm)",
        "lambda0_nm": 514.0,
        "F_max": 25.0,   # (near-field amplitude ~5)^2, Echtermeyer 2011
        "Q_assumed": 8.0,
        "source": "Echtermeyer et al., Nat. Commun. 2, 458 (2011)",
    },
    {
        "label": "Echtermeyer Ti/Au grating, 130 nm finger (633 nm)",
        "lambda0_nm": 633.0,
        "F_max": 25.0,
        "Q_assumed": 8.0,
        "source": "Echtermeyer et al., Nat. Commun. 2, 458 (2011)",
    },
    {
        "label": "Bowtie array on Si waveguide, telecom C-band (1550 nm)",
        "lambda0_nm": 1550.0,
        "F_max": 8.5,    # simulated single-element absorption enhancement
        "Q_assumed": 6.0,
        "source": "arXiv:1808.10823 (ACS Photonics)",
    },
]

# Descriptive-only literature point: Fang et al., Appl. Phys. Lett. 105,
# 241114 (2014), Au nano-antenna, LSPR at 580 nm, tested at 635 nm.
# Deliberately NOT given an F_max/fit here -- see notes file Section 3
# for why its reported "~17 nA/uW, four orders of magnitude higher than
# previously reported single-antenna graphene photodetectors" figure is
# a device-to-device responsivity comparison (bundling contact geometry,
# bias, and collection differences along with the antenna's near-field
# gain), not an isolated absorption/near-field enhancement factor
# comparable to the two designs above.
FANG_ANTENNA_RESONANCE_NM = 580.0
FANG_ANTENNA_TEST_WAVELENGTH_NM = 635.0


def enhancement_factor(wavelength_nm, lambda0_nm, F_max, Q):
    """
    Lorentzian intensity (power) enhancement factor vs. wavelength,
    F(lambda0) = F_max, F -> 1 far from resonance. `Q` sets the
    resonance linewidth: FWHM = lambda0 / Q.
    """
    detuning = 2.0 * Q * (wavelength_nm - lambda0_nm) / lambda0_nm
    return 1.0 + (F_max - 1.0) / (1.0 + detuning ** 2)


def eqe_plasmonic(wavelength_nm, design, eqe_bare=EQE_BARE):
    """Plasmonically-enhanced EQE at a given wavelength for one design
    dict from PLASMONIC_DESIGNS. Capped at 1 (100%) as a physical bound
    -- checked in summary_numbers() below; not reached by either design
    here, but the cap is applied so the model fails safely if a much
    larger F_max were substituted."""
    F = enhancement_factor(wavelength_nm, design["lambda0_nm"],
                            design["F_max"], design["Q_assumed"])
    return min(1.0, eqe_bare * F)


def responsivity_plasmonic(wavelength_nm, design, eqe_bare=EQE_BARE):
    """Bare (gain-free) responsivity with the plasmonic EQE enhancement
    applied, by reusing graphene_photodetector_model.responsivity_bare()
    with the wavelength-dependent EQE computed above -- i.e. the
    plasmonic effect is modeled purely as an EQE multiplier at fixed
    photon energy, not a separate formula."""
    return responsivity_bare(wavelength_nm, eqe=eqe_plasmonic(wavelength_nm, design, eqe_bare))


def responsivity_plasmonic_with_gain(wavelength_nm, tau_trap, design,
                                      eqe_bare=EQE_BARE, tau_transit=TAU_TRANSIT):
    """Combines this module's plasmonic EQE enhancement with
    graphene_photodetector_model's photoconductive gain -- the two
    mechanisms are independent in this compact picture (one boosts
    photocarrier generation, the other boosts collection via gain), so
    they simply multiply."""
    return responsivity_plasmonic(wavelength_nm, design, eqe_bare) * \
        photoconductive_gain(tau_trap, tau_transit)


def plot_plasmonic_enhancement():
    """Two-panel figure: (1) visible-range comparison (bare graphene
    detector vs. the two Echtermeyer-grating designs) and (2) telecom
    C-band comparison (bare vs. the bowtie/waveguide design), each
    showing bare-device responsivity_bare(lambda) alongside the
    plasmonically-enhanced curve, with the design's resonance
    wavelength(s) marked."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

    # --- Panel 1: visible-range finger-grating designs ---
    ax = axes[0]
    lam_vis = np.linspace(420.0, 750.0, 600)
    R_bare_vis = np.array([responsivity_bare(l) for l in lam_vis]) * 1e3  # mA/W
    ax.plot(lam_vis, R_bare_vis, color="0.4", lw=1.8, ls="--",
            label="Bare graphene MGM detector (no plasmonics)")
    colors = ["#1f77b4", "#d62728"]
    for design, color in zip(PLASMONIC_DESIGNS[:2], colors):
        R_plasm = np.array([responsivity_plasmonic(l, design) for l in lam_vis]) * 1e3
        ax.plot(lam_vis, R_plasm, color=color, lw=2.0, label=design["label"])
        ax.axvline(design["lambda0_nm"], color=color, lw=0.8, ls=":")
    ax.axvline(FANG_ANTENNA_RESONANCE_NM, color="0.2", lw=0.8, ls=":")
    ax.annotate("Fang et al. antenna\nLSPR (580 nm) --\nresonance shown,\nnot fit here",
                xy=(FANG_ANTENNA_RESONANCE_NM, ax.get_ylim()[1] if False else 0),
                xytext=(FANG_ANTENNA_RESONANCE_NM + 8, 0.55),
                textcoords=("data", "axes fraction"), fontsize=7.5, color="0.3")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Responsivity (mA/W)")
    ax.set_title("Visible-range finger-grating designs\n(Echtermeyer et al. 2011)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    # --- Panel 2: telecom C-band bowtie/waveguide design ---
    ax = axes[1]
    lam_tc = np.linspace(1400.0, 1700.0, 600)
    R_bare_tc = np.array([responsivity_bare(l) for l in lam_tc]) * 1e3
    design = PLASMONIC_DESIGNS[2]
    R_plasm_tc = np.array([responsivity_plasmonic(l, design) for l in lam_tc]) * 1e3
    ax.plot(lam_tc, R_bare_tc, color="0.4", lw=1.8, ls="--",
            label="Bare graphene MGM detector (no plasmonics)")
    ax.plot(lam_tc, R_plasm_tc, color="#2ca02c", lw=2.0, label=design["label"])
    ax.axvline(design["lambda0_nm"], color="#2ca02c", lw=0.8, ls=":")
    ax.axvspan(1530, 1565, color="0.85", zorder=0, label="Telecom C-band")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Responsivity (mA/W)")
    ax.set_title("Telecom-band bowtie/waveguide design\n(arXiv:1808.10823)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    fig.suptitle(
        "Plasmonic near-field absorption enhancement: bare vs. resonant "
        "graphene photodetector designs\n(EQE_plasmonic = EQE_bare * "
        "Lorentzian F(lambda); F_max literature-calibrated, Q assumed -- "
        "see module docstring)",
        fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig("plasmonic_photodetector_enhancement.png", dpi=150)
    print("Saved plasmonic_photodetector_enhancement.png")


def summary_numbers():
    print("=" * 78)
    print("Plasmonic photodetector enhancement -- summary")
    print("=" * 78)
    for design in PLASMONIC_DESIGNS:
        lam0 = design["lambda0_nm"]
        eqe0 = eqe_plasmonic(lam0, design)
        R0_bare = responsivity_bare(lam0) * 1e3
        R0_plasm = responsivity_plasmonic(lam0, design) * 1e3
        print(f"\n{design['label']}  [{design['source']}]")
        print(f"  On-resonance F_max            = {design['F_max']:.1f}x "
              f"(assumed Q = {design['Q_assumed']:.1f}, FWHM = "
              f"{lam0/design['Q_assumed']:.0f} nm)")
        print(f"  EQE_bare -> EQE_plasmonic(res) = {EQE_BARE*100:.3f}% -> "
              f"{eqe0*100:.3f}%  (capped at 100%: "
              f"{'YES -- check model!' if eqe0 >= 1.0 else 'no, well under cap'})")
        print(f"  R_bare(res) -> R_plasmonic(res) = {R0_bare:.2f} mA/W -> "
              f"{R0_plasm:.2f} mA/W  ({R0_plasm/R0_bare:.1f}x, matches F_max "
              f"by construction on-resonance)")
        # Combined with the strongest literature photogating point
        # (extended trap-lifetime device, tau=1s, from
        # graphene_photodetector_model.py's LITERATURE_POINTS) as an
        # illustrative (not literature-matched) upper-bound combination.
        tau_trap_demo = 1.0
        R_combined = responsivity_plasmonic_with_gain(lam0, tau_trap_demo, design)
        print(f"  Illustrative combination with tau_trap=1s photogating "
              f"gain: {R_combined:.3e} A/W (not a literature device -- "
              f"multiplicative combination of two independent mechanisms "
              f"in this compact model, not verified against a real "
              f"gain+plasmonics device)")

    print("\n" + "-" * 78)
    print("Fang et al. antenna (580 nm LSPR, tested at 635 nm) is cited "
          "descriptively only -- its reported enhancement is a device-to-"
          "device responsivity comparison, not an isolated near-field "
          "absorption factor, so it is not fit into this model. See notes "
          "file Section 3.")
    print("=" * 78)


if __name__ == "__main__":
    summary_numbers()
    plot_plasmonic_enhancement()
