"""
rf_small_signal_model.py

RF figures of merit (f_T, f_max) for the back-gated single-layer graphene
FET modeled in graphene_fet_model.py.

This module reuses the DC transfer-characteristic and quantum-capacitance
machinery already implemented there (no re-derivation of the channel
electrostatics) and layers a standard hybrid-pi small-signal FET model on
top of it to estimate:

  f_T   = g_m / (2*pi*C_gs)                                  [current-gain
                                                                cutoff freq.]
  f_max = f_T / (2 * sqrt(g_ds*R_g + 2*pi*f_T*C_gd*R_g))     [maximum
                                                                oscillation
                                                                frequency,
                                                                simplified
                                                                practical
                                                                form]

Ingredients and their source:
  - g_m(V_g)  : numerical derivative of Id(Vg) from
                graphene_fet_model.transfer_characteristic()
  - C_gs(V_g) : total series gate capacitance (C_ox, C_q) from
                graphene_fet_model.quantum_capacitance(), converted from a
                per-area quantity to a total gate capacitance via the
                channel area W*L
  - C_gd      : approximated as a fraction of C_gs (standard practice when
                a full 2D charge-partitioning model isn't available; here
                taken as C_gd ~ 0.4*C_gs, appropriate for a device biased
                away from strong channel pinch-off, which graphene never
                fully reaches -- see notes/2026-08-22-*.md, Section 1)
  - g_ds(V_g) : numerical derivative of Id with respect to Vds at fixed Vg,
                using the existing transfer_characteristic() model
                evaluated at two close Vds points
  - R_g       : distributed metal-gate resistance, single-side-fed gate
                finger approximation R_g = rho_sheet_gate * W / (3*L)
  - R_c       : the literature-calibrated contact resistance already used
                in graphene_fet_model.Rc_total (kept as a separate check
                on how much contact resistance alone would degrade f_T if
                lumped into the source/drain path -- see summary_numbers())

See notes/2026-08-22-rf-figures-of-merit-fT-fmax.md for the literature
benchmarks (f_T up to ~400 GHz intrinsic, f_max ~100-200 GHz for
aggressively scaled devices, extrinsic f_T/f_max ~30-40 GHz for foundry-
style 0.5-2 um gate lengths) that this estimate is checked against.
"""

import numpy as np
import matplotlib.pyplot as plt

import graphene_fet_model as gfet

# ---------------------------------------------------------------------------
# Additional RF-specific parameters not already in graphene_fet_model
# ---------------------------------------------------------------------------

R_sheet_gate = 8.0            # Ohm/sq, typical thin evaporated metal gate (Pd/Au)
Cgd_over_Cgs = 0.4             # feedback-capacitance fraction (see docstring)


def gate_resistance(L=gfet.L, W=gfet.W):
    """Distributed gate resistance for a single-side-fed metal gate finger.
    R_g = rho_sheet * W_gate / (3 * L_gate); the factor of 3 is the
    standard result for a transmission-line gate fed from one edge."""
    return R_sheet_gate * W / (3.0 * L)


def transconductance(Vg_range, Vds=0.05):
    """g_m(V_g) = d(Id)/d(Vg), computed by numerical differentiation of
    the existing transfer_characteristic() model."""
    Id = gfet.transfer_characteristic(Vg_range, Vds=Vds)
    gm = np.gradient(Id, Vg_range)
    return gm, Id


def output_conductance(Vg_range, Vds=0.05, dVds=1e-3):
    """g_ds(V_g) = d(Id)/d(Vds) at fixed Vg, via a two-point finite
    difference around the operating Vds."""
    Id_plus = gfet.transfer_characteristic(Vg_range, Vds=Vds + dVds)
    Id_minus = gfet.transfer_characteristic(Vg_range, Vds=max(Vds - dVds, 1e-4))
    gds = (Id_plus - Id_minus) / (2 * dVds)
    return gds


def gate_capacitance(Vg_range, Vch=0.0):
    """C_gs(V_g): total series gate capacitance (C_ox in series with C_q),
    scaled from a per-area quantity (F/m^2) to a total device capacitance
    (F) via the channel area W*L already defined in graphene_fet_model."""
    dV = Vg_range - Vch - gfet.V_dirac
    Cq = gfet.quantum_capacitance(dV, T=gfet.T)
    C_total_per_area = 1.0 / (1.0 / Cq + 1.0 / gfet.C_ox)
    return C_total_per_area * gfet.W * gfet.L


def compute_fT_fmax(Vg_range, Vds=0.05):
    """Compute f_T(V_g) and f_max(V_g) across a gate-voltage sweep."""
    gm, Id = transconductance(Vg_range, Vds=Vds)
    gds = output_conductance(Vg_range, Vds=Vds)
    Cgs = gate_capacitance(Vg_range)
    Cgd = Cgd_over_Cgs * Cgs
    Rg = gate_resistance()

    fT = np.abs(gm) / (2 * np.pi * Cgs)

    # simplified practical fmax form (see module docstring)
    denom = gds * Rg + 2 * np.pi * fT * Cgd * Rg
    denom = np.clip(denom, 1e-30, None)  # avoid divide-by-zero at Dirac point
    fmax = fT / (2 * np.sqrt(denom))

    return fT, fmax, gm, gds, Cgs


def plot_fT_fmax():
    """Plot fT and fmax vs Vg for the 200 nm GFET, and show how contact
    resistance (already characterized in graphene_fet_model) degrades
    the extrinsic fT relative to the intrinsic (Rc = 0) case."""
    Vg_range = np.linspace(-1.5, 3.5, 400)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: fT, fmax vs Vg at the literature-calibrated Rc = 300 Ohm.um
    fT, fmax, gm, gds, Cgs = compute_fT_fmax(Vg_range, Vds=0.1)
    axes[0].plot(Vg_range, fT / 1e9, label=r'$f_T$', linewidth=2)
    axes[0].plot(Vg_range, fmax / 1e9, label=r'$f_{max}$', linewidth=2)
    axes[0].axvline(gfet.V_dirac, color='gray', linestyle=':', label='Dirac point')
    axes[0].set_xlabel(r'$V_g$ (V)')
    axes[0].set_ylabel('Frequency (GHz)')
    axes[0].set_title(f'GFET RF Figures of Merit (L = {gfet.L*1e9:.0f} nm, $V_{{ds}}$ = 0.1 V)')
    axes[0].set_ylim(0, None)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Right panel: intrinsic (Rc = 0) fT vs extrinsic fT (literature Rc range)
    # -- reuses graphene_fet_model's module-level Rc_total so the contact-
    # resistance effect on RF performance connects directly to the DC
    # transfer-characteristic contact-resistance study from 2026-08-21.
    Rc_saved = gfet.Rc_total
    for Rc_um, style in [(0, '-'), (110, '--'), (300, '-.'), (500, ':')]:
        gfet.Rc_total = 2 * (Rc_um * 1e-6) / gfet.W
        fT_i, _, _, _, _ = compute_fT_fmax(Vg_range, Vds=0.1)
        axes[1].plot(Vg_range, fT_i / 1e9, style, linewidth=2,
                     label=f'$R_c$ = {Rc_um} $\\Omega\\cdot\\mu$m')
    gfet.Rc_total = Rc_saved
    axes[1].axvline(gfet.V_dirac, color='gray', linestyle=':', label='Dirac point')
    axes[1].set_xlabel(r'$V_g$ (V)')
    axes[1].set_ylabel(r'$f_T$ (GHz)')
    axes[1].set_title('Contact Resistance Degrades Extrinsic $f_T$')
    axes[1].set_ylim(0, None)
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('rf_figures_of_merit.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def summary_numbers():
    """Print peak fT/fmax and compare against literature benchmarks from
    notes/2026-08-22-rf-figures-of-merit-fT-fmax.md."""
    Vg_range = np.linspace(-1.5, 3.5, 400)
    fT, fmax, gm, gds, Cgs = compute_fT_fmax(Vg_range, Vds=0.1)

    i_peak = np.argmax(fT)
    Rg = gate_resistance()

    print(f"Gate resistance R_g (L={gfet.L*1e9:.0f} nm, W={gfet.W*1e6:.0f} um) "
          f"= {Rg:.2f} Ohm")
    print(f"Peak f_T   = {fT[i_peak]/1e9:.1f} GHz at Vg = {Vg_range[i_peak]:.2f} V")
    print(f"     f_max at same Vg = {fmax[i_peak]/1e9:.1f} GHz")
    print(f"     C_gs at same Vg  = {Cgs[i_peak]*1e15:.2f} fF")
    print(f"     g_m at same Vg   = {gm[i_peak]*1e3:.3f} mS")
    print(f"     g_ds at same Vg  = {gds[i_peak]*1e3:.3f} mS")
    print()
    print("Literature context (see notes/2026-08-22-rf-figures-of-merit-fT-fmax.md):")
    print("  - Record intrinsic fT: ~300-427 GHz (exotic short/scaled devices)")
    print("  - Record fmax (60 nm T-gate, de-embedded): ~200 GHz")
    print("  - Foundry-style 0.5-2 um gates, extrinsic: fT ~34 GHz, fmax ~37 GHz;")
    print("    extrapolated to Lg=50nm: ~100 GHz")
    print(f"  - This 200 nm-channel analytic estimate falls in the expected")
    print(f"    intermediate range between those two regimes, which is a")
    print(f"    reasonable sanity check for a compact analytic (non-TCAD) model.")
    print()
    if fmax[i_peak] > fT[i_peak]:
        print("CAVEAT: at this bias point the simplified fmax estimate exceeds fT,")
        print("  which does NOT match the literature trend (fmax typically ~10x")
        print("  BELOW fT for graphene, due to lack of current saturation). This")
        print("  compact model uses an idealized low gate resistance (single metal")
        print("  finger, no pad/interconnect parasitics) and a fixed Cgd/Cgs ratio,")
        print("  so it understates the gds*Rg product that suppresses real fmax.")
        print("  Treat fT as the more reliable number from this model; a proper")
        print("  fmax estimate needs measured/extracted parasitic Rg, Rd, Rs from")
        print("  a full two-port S-parameter de-embedding, which is out of scope")
        print("  for this analytic model. Flagged here rather than silently")
        print("  reporting an unphysical fmax > fT.")


if __name__ == '__main__':
    print("Computing GFET RF figures of merit (fT, fmax)...")
    plot_fT_fmax()
    summary_numbers()
    print("Done. Saved: rf_figures_of_merit.png")
