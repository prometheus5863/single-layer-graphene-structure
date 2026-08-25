"""
graphene_photodetector_model.py

Compact quantitative model of the graphene photodetector responsivity /
response-time / gain tradeoff discussed qualitatively in
thesis_draft/06-graphene-photodetectors.md (Section 6.3) and reviewed in
notes/2026-08-24-photodetector-responsivity.md. This closes out the
"Planned follow-on work" item flagged in both of those documents:
implement a compact responsivity model reusing the existing 2.3%
universal-absorption result (graphene_transport_properties.py) combined
with a parametrized photoconductive gain factor, to reproduce the
qualitative responsivity-vs-response-time tradeoff curve quantitatively.

Physical picture
-----------------
Bare, ungated graphene metal-graphene-metal (MGM) photodetectors are
absorption- and lifetime-limited: only the fixed ~2.3% single-pass
absorption fraction (pi*alpha, wavelength-independent -- the same result
already derived in graphene_transport_properties.py::calculate_optical_
conductivity and plotted in optical_absorption.png) generates
photocarriers, and only those generated within the ~100-200 nm built-in-
field region near a contact are collected before the ~1 ps graphene
carrier lifetime lets them recombine (Mueller, Xia & Avouris, Nature
Photonics 2010). This model uses the literature-reported ~0.1-0.2% net
external quantum efficiency of bare MGM devices directly (rather than
re-deriving the collection efficiency from a diffusion equation, which
would need doping-profile inputs beyond this thesis's current scope --
see the Chapter 4 note on this being a planned spatially-resolved
follow-on).

Photoconductive gain in graphene has no intrinsic (avalanche/bandgap)
origin, so essentially all reported gain is extrinsic photogating: a
trap state captures one carrier species for a lifetime tau_trap, while
the other recirculates through the channel with transit time
tau_transit, giving a classic photoconductor gain
    G = tau_trap / tau_transit
which trades off directly against response bandwidth
    f_3dB = 1 / (2*pi*tau_trap)
so that the gain-bandwidth product G * f_3dB = 1 / (2*pi*tau_transit) is,
in this simple picture, a device constant set purely by the (fast,
transit-time-limited) bare-device physics -- independent of tau_trap.
Section "Model vs. literature" below checks this invariant against the
three literature anchor points from notes/2026-08-24-*.md and finds the
literature devices exceed it by a roughly constant factor of ~4-15x
across six decades of tau_trap, i.e. the simple transit-time picture
gets the *scaling* right (same tradeoff slope) but underestimates the
*absolute* gain-bandwidth product by a device-specific factor this model
does not capture (trap capture cross-section, channel geometry, contact
placement). That gap, not the slope, is what more detailed follow-on
modeling (Section 6.4's plasmonic-enhancement and spatial-doping
extensions) would need to close.

References: see notes/2026-08-24-photodetector-responsivity.md for full
citations. Anchor data points used below:
  - Interfacial photogating: tau ~= 400 ns,  R ~= 1e3   A/W
  - Alternating-channel (2025): tau ~= 3.5 us, R ~= 1.7e4 A/W (1.7e7 mA/W)
  - Extended trap-lifetime photogating: tau ~= 1 s,   R ~= 1e10  A/W
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, h, c

# ---------------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------------

# Universal single-pass absorption fraction, pi*alpha (Chapters 2-3 /
# graphene_transport_properties.py), wavelength-independent in the
# visible/NIR.
ALPHA_ABS = np.pi / 137.036  # ~= 0.0229

# Net external quantum efficiency reported for bare (no-gain) MGM
# graphene photodetectors, ~0.1-0.2% (notes Section 1); take the
# midpoint. This already folds in the collection-bottleneck loss on top
# of the absorption bottleneck.
EQE_BARE = 0.0015

# Bare-device carrier transit time: reuse the 200 nm channel length and
# 0.4 m^2/(V.s) mobility already established in graphene_fet_model.py,
# at the same 0.1 V bias used for that model's transfer-characteristic
# sweep, so tau_transit here is consistent with the rest of this thesis's
# device parameters rather than a new free parameter.
L_channel = 200e-9      # m (graphene_fet_model.py channel length)
mu = 0.4                # m^2/(V.s) (graphene_fet_model.py mobility)
V_bias = 0.1             # V (graphene_fet_model.py Vds used for GFET plots)
E_field = V_bias / L_channel
TAU_TRANSIT = L_channel / (mu * E_field)  # = L^2 / (mu * V_bias)

# Representative visible wavelength for numeric examples
LAMBDA_NM = 550.0


def responsivity_bare(wavelength_nm, eqe=EQE_BARE):
    """
    Bare (gain-free) responsivity R = EQE * e * lambda / (h*c) [A/W].
    Implemented directly from R = EQE * e / E_photon, with E_photon in
    joules from hc/lambda (equivalent to the common R[A/W] = EQE *
    lambda[nm]/1240 shortcut, 1240 = hc/e in eV.nm).
    """
    E_photon_J = h * c / (wavelength_nm * 1e-9)
    return eqe * e / E_photon_J


def photoconductive_gain(tau_trap, tau_transit=TAU_TRANSIT):
    """Classic photoconductor gain G = tau_trap / tau_transit."""
    return tau_trap / tau_transit


def bandwidth_3db(tau_trap):
    """3 dB electrical bandwidth for a single-pole response of time constant tau_trap."""
    return 1.0 / (2 * np.pi * tau_trap)


def responsivity_with_gain(tau_trap, wavelength_nm=LAMBDA_NM, eqe=EQE_BARE,
                            tau_transit=TAU_TRANSIT):
    """Gain-boosted responsivity R(tau_trap) = R_bare(lambda) * G(tau_trap)."""
    return responsivity_bare(wavelength_nm, eqe) * photoconductive_gain(tau_trap, tau_transit)


def gain_bandwidth_invariant(tau_transit=TAU_TRANSIT):
    """
    Model-predicted gain-bandwidth product, independent of tau_trap in
    this simple picture: GBP = 1 / (2*pi*tau_transit).
    """
    return 1.0 / (2 * np.pi * tau_transit)


# ---------------------------------------------------------------------------
# Literature anchor points (notes/2026-08-24-photodetector-responsivity.md)
# ---------------------------------------------------------------------------

LITERATURE_POINTS = [
    {"label": "Interfacial photogating", "tau_s": 400e-9, "R_AW": 1.0e3},
    {"label": "Alternating-channel (2025)", "tau_s": 3.5e-6, "R_AW": 1.7e4},
    {"label": "Extended-trap-lifetime photogating", "tau_s": 1.0, "R_AW": 1.0e10},
]


def plot_responsivity_tradeoff():
    """
    Two-panel figure:
      Left:  model responsivity-vs-response-time curve (log-log) with the
             three literature anchor points overlaid.
      Right: gain-bandwidth product for the model (constant line) vs. the
             effective gain-bandwidth product implied by each literature
             point, showing the model captures the correct tradeoff slope
             but underestimates absolute GBP by a roughly constant factor.
    """
    tau_range = np.logspace(-12, 1, 500)  # 1 ps to 10 s
    R_model = responsivity_with_gain(tau_range)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # --- Left: responsivity vs. response time ---
    ax = axes[0]
    ax.loglog(tau_range, R_model, 'b-', linewidth=2,
              label=r'Model: $R = R_{bare}\times(\tau_{trap}/\tau_{transit})$')
    ax.axhline(responsivity_bare(LAMBDA_NM), color='gray', linestyle=':',
               label=r'$R_{bare}$ (no gain, $\tau=\tau_{transit}$)')
    for pt in LITERATURE_POINTS:
        ax.plot(pt["tau_s"], pt["R_AW"], 'r*', markersize=16, zorder=5)
        ax.annotate(pt["label"], (pt["tau_s"], pt["R_AW"]),
                    textcoords="offset points", xytext=(8, -4), fontsize=8)
    ax.set_xlabel(r'Response time $\tau_{trap}$ (s)')
    ax.set_ylabel('Responsivity (A/W)')
    ax.set_title('Graphene Photodetector Responsivity vs. Response Time')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=9, loc='upper left')

    # --- Right: gain-bandwidth product, model (constant) vs. literature ---
    ax2 = axes[1]
    gbp_model = gain_bandwidth_invariant()
    ax2.axhline(gbp_model, color='blue', linewidth=2,
                label=f'Model invariant GBP = {gbp_model:.2e} Hz')
    labels, gbp_lit = [], []
    for pt in LITERATURE_POINTS:
        G_lit = pt["R_AW"] / responsivity_bare(LAMBDA_NM)
        BW_lit = bandwidth_3db(pt["tau_s"])
        gbp = G_lit * BW_lit
        labels.append(pt["label"])
        gbp_lit.append(gbp)
    x_pos = np.arange(len(labels))
    ax2.bar(x_pos, gbp_lit, color='crimson', alpha=0.7, label='Literature-implied GBP')
    ax2.set_yscale('log')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(labels, rotation=20, ha='right', fontsize=8)
    ax2.set_ylabel('Gain x Bandwidth product (Hz)')
    ax2.set_title('Gain-Bandwidth Invariant: Model vs. Literature')
    ax2.legend(fontsize=9)
    ax2.grid(True, which='both', axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('photodetector_responsivity_gain_tradeoff.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def summary_numbers():
    """Print model predictions vs. literature anchor points for logging/sanity-checking."""
    print(f"tau_transit (200 nm channel, mu=0.4 m^2/Vs, Vds=0.1V) = {TAU_TRANSIT:.3e} s")
    print(f"R_bare at {LAMBDA_NM:.0f} nm (EQE={EQE_BARE*100:.2f}%)        = "
          f"{responsivity_bare(LAMBDA_NM)*1e3:.3f} mA/W")
    print(f"Model gain-bandwidth invariant GBP           = {gain_bandwidth_invariant():.3e} Hz "
          f"({gain_bandwidth_invariant()/1e9:.1f} GHz)")
    print()
    print(f"{'Literature point':<32}{'tau (s)':>12}{'R_lit (A/W)':>14}{'R_model (A/W)':>16}{'ratio':>10}")
    for pt in LITERATURE_POINTS:
        R_model = responsivity_with_gain(pt["tau_s"])
        ratio = pt["R_AW"] / R_model
        print(f"{pt['label']:<32}{pt['tau_s']:>12.2e}{pt['R_AW']:>14.2e}"
              f"{R_model:>16.2e}{ratio:>10.1f}x")


if __name__ == '__main__':
    print("Generating graphene photodetector responsivity/gain tradeoff model...")
    plot_responsivity_tradeoff()
    summary_numbers()
    print("Done. Saved: photodetector_responsivity_gain_tradeoff.png")
