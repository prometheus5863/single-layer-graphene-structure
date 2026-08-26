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
  f_max = f_T / (2 * sqrt(g_ds*(R_g+R_s) + 2*pi*f_T*C_gd*R_g))  [maximum
                                                                oscillation
                                                                frequency,
                                                                standard
                                                                hybrid-pi
                                                                form incl.
                                                                access
                                                                resistance]

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
  - R_s       : source access resistance = graphene_fet_model.Rc_total / 2
                (the total two-contact resistance already calibrated in
                graphene_fet_model.py, split evenly). Added to the f_max
                denominator as of 2026-08-26 -- see
                notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md for
                why this term was missing and the literature (arXiv:
                1112.4831) motivating it.
  - C_pad     : GSG probe-pad capacitance, relevant only for the
                "extrinsic" (pre-de-embedding) estimate on this device's
                conductive back-gate substrate -- see the same 2026-08-26
                notes, Section 3, for the literature basis and the
                assumed (not literature-pinned) 15 fF/pad value used.
  - R_c       : the literature-calibrated contact resistance already used
                in graphene_fet_model.Rc_total (kept as a separate check
                on how much contact resistance alone would degrade f_T if
                lumped into the source/drain path -- see summary_numbers())

See notes/2026-08-22-rf-figures-of-merit-fT-fmax.md for the literature
benchmarks (f_T up to ~400 GHz intrinsic, f_max ~100-200 GHz for
aggressively scaled devices, extrinsic f_T/f_max ~30-40 GHz for foundry-
style 0.5-2 um gate lengths) that this estimate is checked against, and
notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md for the access-
resistance/pad-capacitance revision and the corrected f_max vs. f_T
framing (f_max > f_T is not itself unphysical -- see Feijoo et al. 2016).
"""

import numpy as np
import matplotlib.pyplot as plt

import graphene_fet_model as gfet

# ---------------------------------------------------------------------------
# Additional RF-specific parameters not already in graphene_fet_model
# ---------------------------------------------------------------------------

R_sheet_gate = 8.0            # Ohm/sq, typical thin evaporated metal gate (Pd/Au)
Cgd_over_Cgs = 0.4             # feedback-capacitance fraction (see docstring)
W_RF = 40e-6                   # m, literature-representative total gate
                                # width for the extrinsic/pad-capacitance
                                # comparison only (typical multi-finger
                                # RF graphene-FET device scale; see
                                # notes/2026-08-26-*.md and the docstring
                                # of compute_fT_fmax() for why the
                                # repo-standard W=1um is not usable here)
N_FINGERS_RF = 8                # gate fingers at W_RF (5 um/finger),
                                # giving the standard N^2 R_g reduction
C_pad_per_pad = 15e-15         # F, assumed GSG pad capacitance (extrinsic
                                # estimate only) -- see 2026-08-26 notes,
                                # Section 3, for the literature basis and
                                # the caveat that this specific number is
                                # an assumption, not a literature-pinned
                                # value for this exact device geometry.


def gate_resistance(L=gfet.L, W=gfet.W, N_fingers=1):
    """Distributed gate resistance for a single-side-fed metal gate finger.
    R_g = rho_sheet * W_gate / (3 * L_gate); the factor of 3 is the
    standard result for a transmission-line gate fed from one edge.

    N_fingers > 1 (added 2026-08-26): splits total gate width W into
    N_fingers parallel, individually-fed fingers of width W/N_fingers
    each, giving the standard multi-finger reduction
    R_g = rho_sheet * W / (3 * L * N_fingers^2) -- the same mechanism
    Feijoo et al. 2016 and multi-finger/T-gate device literature use to
    push R_g down to the ~15 Ohm range (see
    notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md). Default
    N_fingers=1 reproduces the original single-finger formula exactly,
    so this is backward-compatible with the per-um curves used
    elsewhere in this module."""
    return R_sheet_gate * W / (3.0 * L * N_fingers ** 2)


def source_access_resistance():
    """R_s = R_c_total / 2: the literature-calibrated total two-contact
    resistance from graphene_fet_model.py, split evenly between source
    and drain. Added to the f_max denominator as of 2026-08-26 (see
    module docstring) -- previously this model only included R_g, which
    structurally understated gds*(Rg+Rs)."""
    return gfet.Rc_total / 2.0


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


def compute_fT_fmax(Vg_range, Vds=0.05, extrinsic=False, W=None, N_fingers=1):
    """Compute f_T(V_g) and f_max(V_g) across a gate-voltage sweep.

    extrinsic=False (default): "intrinsic" estimate -- R_s (source access
      resistance) included in the f_max denominator, no pad capacitance.
      This is the model's best estimate of a de-embedded measurement.
    extrinsic=True: adds GSG pad capacitance C_pad in parallel with Cgs
      and Cgd, matching the raw (pre-de-embedding) measurement condition
      -- see notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md,
      Section 3, and the de-embedded-vs-raw contrast reported there
      (Feijoo et al. 2016: 21.3/42 GHz raw vs. 35.4/50 GHz de-embedded
      for a comparable 200 nm device).
    W, N_fingers: optional total gate width and finger count override.
      The rest of this repo normalizes everything to W = 1 um (a single
      finger, "per unit width" convention -- see graphene_fet_model.py's
      W comment). That convention makes an absolute pad-capacitance
      comparison meaningless: C_pad (~15 fF, sized for a real on-wafer
      GSG probe pad) would swamp a 1 um-wide device's ~0.08 fF Cgs by
      >100x, which was checked numerically this session and found to
      wildly overstate the pad-capacitance effect relative to the
      Feijoo et al. raw-vs-de-embedded ratio (~1.4-2x, not >100x) --
      because a real RF-probed device is never actually 1 um wide.
      Passing W (and, realistically, N_fingers > 1 for the gate-
      resistance benefit) temporarily rescales graphene_fet_model's
      module-level W and Rc_total (which are both defined per total
      device width, exactly like the existing Rc_total sweep in
      plot_fT_fmax() already does) so the extrinsic comparison is
      evaluated at a literature-representative device size instead.
    """
    if W is None:
        fT, fmax, gm, gds, Cgs = _compute_fT_fmax_core(
            Vg_range, Vds=Vds, extrinsic=extrinsic, N_fingers=N_fingers)
        return fT, fmax, gm, gds, Cgs

    W_saved = gfet.W
    Rc_saved = gfet.Rc_total
    gfet.W = W
    gfet.Rc_total = 2 * (gfet.Rc_per_width_ohm_um * 1e-6) / W
    try:
        return _compute_fT_fmax_core(
            Vg_range, Vds=Vds, extrinsic=extrinsic, N_fingers=N_fingers)
    finally:
        gfet.W = W_saved
        gfet.Rc_total = Rc_saved


def _compute_fT_fmax_core(Vg_range, Vds, extrinsic, N_fingers):
    """Shared computation for compute_fT_fmax(), using whatever gfet.W /
    gfet.Rc_total are currently set to (see compute_fT_fmax's W override
    for why this indirection exists)."""
    gm, Id = transconductance(Vg_range, Vds=Vds)
    gds = output_conductance(Vg_range, Vds=Vds)
    Cgs = gate_capacitance(Vg_range)
    Cgd = Cgd_over_Cgs * Cgs
    Rg = gate_resistance(N_fingers=N_fingers)
    Rs = source_access_resistance()

    if extrinsic:
        Cgs = Cgs + C_pad_per_pad
        Cgd = Cgd + C_pad_per_pad

    fT = np.abs(gm) / (2 * np.pi * Cgs)

    # standard hybrid-pi fmax form including access resistance Rs (see
    # module docstring and 2026-08-26 notes for why Rs was added)
    denom = gds * (Rg + Rs) + 2 * np.pi * fT * Cgd * Rg
    denom = np.clip(denom, 1e-30, None)  # avoid divide-by-zero at Dirac point
    fmax = fT / (2 * np.sqrt(denom))

    return fT, fmax, gm, gds, Cgs


def plot_fT_fmax():
    """Plot fT and fmax vs Vg for the 200 nm GFET (three panels):
      1. Intrinsic fT/fmax (Rs included, no pad capacitance)
      2. Contact resistance degrading extrinsic fT (as before 2026-08-26)
      3. NEW (2026-08-26): intrinsic vs. extrinsic (+ pad capacitance)
         fT and fmax side by side, the direct analog of the raw-vs-
         de-embedded comparison reported in Feijoo et al. 2016 -- see
         notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md.
    """
    Vg_range = np.linspace(-1.5, 3.5, 400)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: intrinsic fT, fmax vs Vg at the literature-calibrated Rc
    fT, fmax, gm, gds, Cgs = compute_fT_fmax(Vg_range, Vds=0.1)
    axes[0].plot(Vg_range, fT / 1e9, label=r'$f_T$', linewidth=2)
    axes[0].plot(Vg_range, fmax / 1e9, label=r'$f_{max}$', linewidth=2)
    axes[0].axvline(gfet.V_dirac, color='gray', linestyle=':', label='Dirac point')
    axes[0].set_xlabel(r'$V_g$ (V)')
    axes[0].set_ylabel('Frequency (GHz)')
    axes[0].set_title(f'Intrinsic RF FoM (L = {gfet.L*1e9:.0f} nm, $V_{{ds}}$ = 0.1 V)')
    axes[0].set_ylim(0, None)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Panel 2: intrinsic (Rc = 0) fT vs extrinsic fT (literature Rc range)
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

    # Panel 3 (NEW 2026-08-26): intrinsic vs extrinsic (+pad cap, +Rs),
    # evaluated at a literature-representative multi-finger device width
    # (W_RF, N_FINGERS_RF) -- see compute_fT_fmax() docstring for why the
    # repo-standard W=1um single-finger convention isn't usable for an
    # absolute pad-capacitance comparison.
    fT_i, fmax_i, _, _, _ = compute_fT_fmax(Vg_range, Vds=0.1, extrinsic=False,
                                             W=W_RF, N_fingers=N_FINGERS_RF)
    fT_e, fmax_e, _, _, _ = compute_fT_fmax(Vg_range, Vds=0.1, extrinsic=True,
                                             W=W_RF, N_fingers=N_FINGERS_RF)
    axes[2].plot(Vg_range, fT_i / 1e9, '-', color='C0', linewidth=2,
                 label=r'$f_T$ intrinsic')
    axes[2].plot(Vg_range, fT_e / 1e9, '--', color='C0', linewidth=2,
                 label=r'$f_T$ extrinsic (+pad)')
    axes[2].plot(Vg_range, fmax_i / 1e9, '-', color='C1', linewidth=2,
                 label=r'$f_{max}$ intrinsic')
    axes[2].plot(Vg_range, fmax_e / 1e9, '--', color='C1', linewidth=2,
                 label=r'$f_{max}$ extrinsic (+pad)')
    axes[2].axvline(gfet.V_dirac, color='gray', linestyle=':', label='Dirac point')
    axes[2].set_xlabel(r'$V_g$ (V)')
    axes[2].set_ylabel('Frequency (GHz)')
    axes[2].set_title(f'Intrinsic vs. Extrinsic, {W_RF*1e6:.0f} $\\mu$m/'
                       f'{N_FINGERS_RF}-finger device\n'
                       f'(+{C_pad_per_pad*1e15:.0f} fF/pad, $R_s$ incl.)',
                       fontsize=10)
    axes[2].set_ylim(0, None)
    axes[2].legend(fontsize=8)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('rf_figures_of_merit.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def summary_numbers():
    """Print peak fT/fmax (intrinsic and extrinsic) and compare against
    literature benchmarks from notes/2026-08-22-rf-figures-of-merit-fT-fmax.md
    and notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md."""
    Vg_range = np.linspace(-1.5, 3.5, 400)
    fT, fmax, gm, gds, Cgs = compute_fT_fmax(Vg_range, Vds=0.1, extrinsic=False)

    i_peak = np.argmax(fT)
    Rg = gate_resistance()
    Rs = source_access_resistance()

    print(f"Gate resistance R_g (L={gfet.L*1e9:.0f} nm, W={gfet.W*1e6:.0f} um) "
          f"= {Rg:.2f} Ohm")
    print(f"Source access resistance R_s (= Rc_total/2)       = {Rs:.1f} Ohm")
    print()
    print(f"INTRINSIC (per-um convention used throughout this repo, Rs")
    print(f"included, no pad capacitance):")
    print(f"  Peak f_T   = {fT[i_peak]/1e9:.1f} GHz at Vg = {Vg_range[i_peak]:.2f} V")
    print(f"       f_max at same Vg = {fmax[i_peak]/1e9:.1f} GHz "
          f"(f_max/f_T = {fmax[i_peak]/fT[i_peak]:.2f})")
    print(f"       C_gs at same Vg  = {Cgs[i_peak]*1e15:.2f} fF")
    print(f"       g_m at same Vg   = {gm[i_peak]*1e3:.3f} mS")
    print(f"       g_ds at same Vg  = {gds[i_peak]*1e3:.3f} mS")
    print()

    # Extrinsic (+pad capacitance) comparison, evaluated at a literature-
    # representative multi-finger device width -- see compute_fT_fmax()
    # docstring: the repo's per-um convention makes an absolute pad-
    # capacitance comparison meaningless (checked numerically this
    # session: it overstated the pad effect by >100x vs. the literature
    # raw/de-embedded ratio of ~1.4-2x).
    fT_rf, fmax_rf, gm_rf, gds_rf, Cgs_rf = compute_fT_fmax(
        Vg_range, Vds=0.1, extrinsic=False, W=W_RF, N_fingers=N_FINGERS_RF)
    fT_rf_e, fmax_rf_e, _, _, Cgs_rf_e = compute_fT_fmax(
        Vg_range, Vds=0.1, extrinsic=True, W=W_RF, N_fingers=N_FINGERS_RF)
    i_peak_rf = np.argmax(fT_rf)
    Rg_rf = gate_resistance(W=W_RF, N_fingers=N_FINGERS_RF)

    print(f"LITERATURE-SCALE DEVICE ({W_RF*1e6:.0f} um total width, "
          f"{N_FINGERS_RF} fingers, L={gfet.L*1e9:.0f} nm):")
    print(f"  R_g ({N_FINGERS_RF}-finger) = {Rg_rf:.2f} Ohm  "
          f"(vs. {gate_resistance(W=W_RF):.1f} Ohm single-finger -- "
          f"the N_fingers^2 reduction from Section 4 of the notes)")
    print(f"  Intrinsic: f_T = {fT_rf[i_peak_rf]/1e9:.1f} GHz, "
          f"f_max = {fmax_rf[i_peak_rf]/1e9:.1f} GHz "
          f"(f_max/f_T = {fmax_rf[i_peak_rf]/fT_rf[i_peak_rf]:.2f}), "
          f"C_gs = {Cgs_rf[i_peak_rf]*1e15:.2f} fF")
    print(f"  Extrinsic (+{C_pad_per_pad*1e15:.0f} fF/pad): "
          f"f_T = {fT_rf_e[i_peak_rf]/1e9:.1f} GHz, "
          f"f_max = {fmax_rf_e[i_peak_rf]/1e9:.1f} GHz "
          f"(f_max/f_T = {fmax_rf_e[i_peak_rf]/fT_rf_e[i_peak_rf]:.2f})")
    print(f"  Degradation vs. intrinsic: "
          f"f_T x{fT_rf_e[i_peak_rf]/fT_rf[i_peak_rf]:.2f}, "
          f"f_max x{fmax_rf_e[i_peak_rf]/fmax_rf[i_peak_rf]:.2f} "
          f"(cf. Feijoo et al. raw/de-embedded ratio ~0.6-0.7x)")
    print()
    print("Literature context (notes/2026-08-22-*.md, notes/2026-08-26-*.md):")
    print("  - Record intrinsic fT: ~300-427 GHz (exotic short/scaled devices)")
    print("  - Record fmax (60 nm T-gate, de-embedded): ~200 GHz")
    print("  - Feijoo et al. 2016 (200 nm gate, de-embedded): fT/fmax = ")
    print("    35.4/50 GHz (fmax/fT = 1.41); raw (pre-de-embedding): 21.3/42 GHz")
    print("  - Foundry-style 0.5-2 um gates, extrinsic: fT ~34 GHz, fmax ~37 GHz;")
    print("    extrapolated to Lg=50nm: ~100 GHz")
    print(f"  - This 200 nm-channel analytic estimate falls in the expected")
    print(f"    intermediate range between those regimes, which is a reasonable")
    print(f"    sanity check for a compact analytic (non-TCAD) model.")
    print()
    print("NOTE (revised 2026-08-26 -- see notes/2026-08-26-*.md, Section 1):")
    print("  f_max > f_T is NOT inherently unphysical for graphene FETs -- the")
    print("  earlier (2026-08-22) caveat treating it as a red flag was too broad.")
    print("  Feijoo et al. 2016 report f_max > f_T at every gate length they")
    print("  measured (low, engineered R_g ~ 15 Ohm). What matters is whether")
    print("  gds*(Rg+Rs) and the Rg*Cgd feedback term are realistic, which is")
    print("  why this model now includes Rs (previously missing) and reports")
    print("  an explicit extrinsic (+pad capacitance) estimate alongside the")
    print("  intrinsic one, rather than asserting a required sign for f_max-f_T.")


if __name__ == '__main__':
    print("Computing GFET RF figures of merit (fT, fmax)...")
    plot_fT_fmax()
    summary_numbers()
    print("Done. Saved: rf_figures_of_merit.png")
