"""
graphene_fet_model.py

Device-oriented model of a back-gated single-layer graphene field-effect
transistor (GFET), built on top of the fundamental band-structure and
transport results already established elsewhere in this repo
(graphene_band_structure.py, graphene_transport_properties.py).

This script adds three device-physics ingredients that are absent from a
pure band-structure treatment but are what semiconductor/device engineers
actually evaluate a 2D-material transistor on:

  1. Quantum capacitance C_q(V_g) of the graphene channel (finite DOS effect)
  2. Series combination with oxide capacitance C_ox to get the true gate
     capacitance C_total(V_g), and the resulting self-consistent channel
     carrier density n(V_g)
  3. A simple long-channel drift model for the ambipolar Id-Vg transfer
     characteristic, including a back-of-envelope contact-resistance term
     calibrated against literature values (see notes/2026-08-21-*.md)

This is a compact analytic model (not a TCAD-grade solver) intended to
reproduce the qualitative, literature-known shape of a GFET transfer curve:
a V-shaped ambipolar characteristic with a minimum conductance point (the
"Dirac point") that shifts with any built-in doping, and current asymmetry
introduced by drain-bias-dependent local charge neutrality points along the
channel.

References (see notes/2026-08-21-contact-resistance-and-quantum-capacitance.md
for full citations and discussion):
  - Phenomenological quantum capacitance model, arXiv:1105.5827
  - Quantum-capacitance-limited vertical scaling, ResearchGate 49838119
  - Contact resistance literature range ~110-500 Ohm.um (PubMed 21297624,
    Nature Sci. Rep. 2024 s41598-024-58360-9)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, hbar, epsilon_0, Boltzmann as kB

# ---------------------------------------------------------------------------
# Physical / device parameters
# ---------------------------------------------------------------------------

T = 300.0                    # Temperature (K)
v_F = 1.0e6                  # Graphene Fermi velocity (m/s)

# Back-gate oxide stack: SiO2, thickness t_ox
eps_r_ox = 3.9                # relative permittivity of SiO2
t_ox = 90e-9                  # oxide thickness (m) -- common back-gate value
C_ox = eps_r_ox * epsilon_0 / t_ox   # geometric oxide capacitance (F/m^2)

# Channel geometry
W = 1e-6                      # channel width (m) -> 1 um, results normalized per um
L = 200e-9                    # channel length (m) -> 200 nm, short-channel regime

# Mobility (diffusive, literature-typical for SiO2-supported exfoliated/CVD graphene)
mu = 0.4                      # m^2/(V.s)  == 4000 cm^2/(V.s)

# Residual carrier density at the Dirac point (disorder-induced puddles)
n_puddle = 5e15               # 1/m^2  (~5e11 cm^-2, typical for SiO2-supported graphene)

# Dirac point gate voltage (built-in doping offset, e.g. from substrate/contacts)
V_dirac = 0.8                 # V

# Contact resistance, width-normalized (Ohm.um), literature range ~110-500 Ohm.um
Rc_per_width_ohm_um = 300.0
Rc_total = 2 * (Rc_per_width_ohm_um * 1e-6) / W   # two contacts, converted to Ohms


def quantum_capacitance(V_g_minus_Vdirac, T=300.0):
    """
    Graphene quantum capacitance as a function of (V_g - V_Dirac), including
    finite-temperature broadening near the Dirac point so C_q does not
    diverge to zero exactly at the neutrality point.

    C_q = (2 e^2 / (pi * hbar^2 * v_F^2)) * kB*T * F(eta)

    where F(eta) is evaluated numerically from the graphene DOS integrated
    against the derivative of the Fermi function; here we use the standard
    closed-form high-|eta| limit patched with a thermal rounding term near
    eta = 0, following the phenomenological approach of arXiv:1105.5827.
    """
    kT = kB * T
    eta = e * V_g_minus_Vdirac / kT  # normalized energy relative to thermal energy

    # Prefactor: e^2 * DOS-scale for graphene
    prefactor = (2 * e**3) / (np.pi * hbar**2 * v_F**2)

    # Thermally-broadened |V_g - V_dirac|-like term:
    # smooth approximation that -> |V| for |eta| >> 1 and -> ~ (kT/e)*ln(2) at eta=0
    kT_over_e = kT / e
    smoothed_V = kT_over_e * np.log(2 * (1 + np.cosh(eta))) 

    return prefactor * smoothed_V


def carrier_density(V_g, V_ch=0.0):
    """
    Self-consistent channel sheet carrier density n(V_g, V_ch) [1/m^2],
    combining the oxide capacitance (electrostatic charge from V_g - V_ch)
    with the quantum-capacitance-limited channel response, plus a residual
    puddle density that regularizes n at the Dirac point.

    We solve the series-capacitance charge relation self-consistently:
        Q = C_ox * (V_g - V_ch - V_dirac - V_channel_shift)
    where V_channel_shift is determined by requiring Q/e also equals the
    charge implied by the quantum capacitance branch. For a closed analytic
    form we use the standard local (no-cross-coupling) approximation:
        n(V_g) = C_ox * (V_g - V_ch - V_dirac) / e   [electrostatic estimate]
    corrected by the quantum-capacitance series factor
        n_eff = n_electrostatic * C_q / (C_q + C_ox)
    which correctly suppresses induced charge exactly where C_q is small
    (near the Dirac point), reproducing the vertical-scaling-limited
    behavior discussed in the notes.
    """
    dV = V_g - V_ch - V_dirac
    n_electrostatic = C_ox * dV / e

    C_q = quantum_capacitance(dV, T=T)
    series_factor = C_q / (C_q + C_ox)

    n_eff = n_electrostatic * series_factor

    # Add residual puddle density in quadrature (regularizes conductivity
    # minimum at the Dirac point, matching experimentally observed minimum
    # conductivity plateau rather than a true zero). Only the magnitude is
    # used downstream (sheet_conductivity takes abs(n)), so np.sign() would
    # incorrectly zero the result exactly at n_eff = 0; return the magnitude
    # directly instead.
    n_total = np.sqrt(n_eff**2 + n_puddle**2)
    return n_total


def sheet_conductivity(n):
    """
    Drude sheet conductivity sigma = n * e * mu  [S] (per square, i.e. S/sq)
    """
    return np.abs(n) * e * mu


def channel_resistance(V_g, V_ch=0.0):
    """
    Channel resistance for given gate/channel bias point, R = (L/W) / sigma_sheet
    """
    n = carrier_density(V_g, V_ch)
    sigma_sheet = sheet_conductivity(n)
    return (L / W) / sigma_sheet


def transfer_characteristic(Vg_range, Vds=0.05):
    """
    Compute Id vs Vg at fixed (small) Vds, using a simple long-channel
    drift approximation with the channel divided into segments to capture
    the drain-bias-dependent local charge-neutrality-point shift, plus a
    lumped series contact resistance Rc_total.

    This reproduces the well-known ambipolar V-shaped GFET transfer curve.
    """
    Id = np.zeros_like(Vg_range)
    n_segments = 50
    V_channel_profile = np.linspace(0, Vds, n_segments)

    for i, V_g in enumerate(Vg_range):
        # Average the local channel resistance along the channel to account
        # for the fact that the local carrier density (and hence local
        # resistivity) varies with position due to the Vds drop -- this is
        # what produces the characteristic asymmetric/kinked GFET Id-Vg
        # curve rather than a symmetric V shape once Vds is non-negligible.
        R_channel_local = np.mean([
            channel_resistance(V_g, V_ch) for V_ch in V_channel_profile
        ])
        R_total = R_channel_local + Rc_total
        Id[i] = Vds / R_total

    return Id


def plot_quantum_capacitance():
    """Plot C_q and C_total vs (Vg - V_Dirac), illustrating the
    quantum-capacitance-limited scaling discussed in the notes."""
    dV = np.linspace(-2, 2, 500)
    Cq = quantum_capacitance(dV, T=T)
    C_total = 1.0 / (1.0 / Cq + 1.0 / C_ox)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(dV, Cq * 1e4, label=r'$C_q$ (quantum capacitance)', linewidth=2)
    ax.axhline(C_ox * 1e4, color='gray', linestyle='--',
               label=r'$C_{ox}$ (90 nm SiO$_2$)')
    ax.plot(dV, C_total * 1e4, label=r'$C_{total} = (C_{ox}^{-1}+C_q^{-1})^{-1}$',
            linewidth=2, color='crimson')
    ax.set_xlabel(r'$V_g - V_{Dirac}$ (V)')
    ax.set_ylabel(r'Capacitance ($\mu$F/cm$^2$)')
    ax.set_title('Graphene Quantum Capacitance vs. Gate Overdrive (T = 300 K)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('quantum_capacitance.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_transfer_characteristics():
    """Plot Id-Vg transfer characteristics at several Vds, and separately
    show the effect of contact resistance on the curve."""
    Vg_range = np.linspace(-2.0, 3.5, 300)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: transfer curves at multiple Vds
    for Vds in [0.05, 0.2, 0.5]:
        Id = transfer_characteristic(Vg_range, Vds=Vds)
        axes[0].plot(Vg_range, Id * 1e3, label=f'$V_{{ds}}$ = {Vds} V', linewidth=2)
    axes[0].axvline(V_dirac, color='gray', linestyle=':', label='Dirac point')
    axes[0].set_xlabel(r'$V_g$ (V)')
    axes[0].set_ylabel(r'$I_d$ (mA/$\mu$m width)')
    axes[0].set_title('GFET Transfer Characteristics (200 nm channel)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Right panel: effect of contact resistance at fixed Vds
    global Rc_total
    Rc_saved = Rc_total
    Vds = 0.1
    for Rc_um in [0, 110, 300, 500]:
        Rc_total = 2 * (Rc_um * 1e-6) / W
        Id = transfer_characteristic(Vg_range, Vds=Vds)
        axes[1].plot(Vg_range, Id * 1e3, label=f'$R_c$ = {Rc_um} $\\Omega\\cdot\\mu$m', linewidth=2)
    Rc_total = Rc_saved
    axes[1].axvline(V_dirac, color='gray', linestyle=':', label='Dirac point')
    axes[1].set_xlabel(r'$V_g$ (V)')
    axes[1].set_ylabel(r'$I_d$ (mA/$\mu$m width)')
    axes[1].set_title(f'Effect of Contact Resistance ($V_{{ds}}$ = {Vds} V)')
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('gfet_transfer_characteristics.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def summary_numbers():
    """Print a few representative numbers for sanity-checking / logging."""
    n_max = carrier_density(3.5)
    R_ch_on = channel_resistance(3.5)
    R_ch_dirac = channel_resistance(V_dirac)
    print(f"C_ox (90 nm SiO2)          = {C_ox*1e4:.3f} uF/cm^2")
    print(f"n at Vg=3.5V (on-state)    = {n_max:.3e} 1/m^2  ({n_max*1e-4:.3e} 1/cm^2)")
    print(f"R_channel at Vg=3.5V       = {R_ch_on:.1f} Ohm (L=200nm, W=1um)")
    print(f"R_channel at Dirac point   = {R_ch_dirac:.1f} Ohm")
    print(f"R_contact (total, 2x{Rc_per_width_ohm_um} Ohm.um / 1um width) = {Rc_total:.1f} Ohm")
    print(f"Contact resistance fraction of total R at Dirac point = "
          f"{Rc_total/(Rc_total+R_ch_dirac)*100:.1f}%")


if __name__ == '__main__':
    print("Generating graphene FET quantum capacitance and transfer characteristic plots...")
    plot_quantum_capacitance()
    plot_transfer_characteristics()
    summary_numbers()
    print("Done. Saved: quantum_capacitance.png, gfet_transfer_characteristics.png")
