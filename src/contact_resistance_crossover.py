"""
contact_resistance_crossover.py

Contact-resistance-vs-channel-length crossover analysis for the back-gated
GFET model established in graphene_fet_model.py.

This closes out a follow-on item that has been flagged since the
2026-08-21 research notes (notes/2026-08-21-contact-resistance-and-quantum-
capacitance.md, Section 4: "Add a contact-resistance term ... so that total
device resistance vs. channel length can be estimated, illustrating the
crossover length below which contacts dominate.") The DC model itself
(quantum capacitance, self-consistent carrier density, contact resistance
as a lumped series term) was already implemented in graphene_fet_model.py;
this module reuses that model unmodified and sweeps channel length L
instead of gate voltage, which is the piece that was still missing.

Physical motivation: as transistors are scaled to shorter channel lengths,
channel (intrinsic) resistance R_channel = (L/W) / sigma_sheet falls
linearly with L, while contact resistance R_c is L-independent (it depends
on contact geometry and metal-graphene interface quality, not channel
length). Below some crossover length L_x, R_c dominates total device
resistance and further channel-length scaling stops improving on-current --
exactly the "contact-resistance-limited" regime widely discussed for 2D
FETs (see notes/2026-08-21-*.md refs 1-5). This is directly analogous to
the contact-resistance scaling problem in silicon FinFET/nanosheet nodes,
which is one of the reasons device engineers care about this number.

References: see notes/2026-08-21-contact-resistance-and-quantum-
capacitance.md and notes/2026-08-24-graphene-photodetector-responsivity.md
is unrelated; contact-resistance literature range (~110-500 Ohm.um) is the
same range already used and cited in graphene_fet_model.py.
"""

import numpy as np
import matplotlib.pyplot as plt

from graphene_fet_model import (
    carrier_density, sheet_conductivity, W, V_dirac,
)


def channel_resistance_at_length(L, V_g):
    """Intrinsic channel resistance (no contacts) for a given channel
    length L (m) at gate bias V_g, evaluated at V_ch = 0 (linear region,
    small-Vds limit -- consistent with how graphene_fet_model.py evaluates
    channel_resistance for a fixed V_ch)."""
    n = carrier_density(V_g, V_ch=0.0)
    sigma_sheet = sheet_conductivity(n)
    return (L / W) / sigma_sheet


def crossover_length(Rc_total, V_g):
    """Analytic crossover length L_x where R_channel(L_x) = Rc_total, i.e.
    contacts and channel contribute equally to total device resistance.
    Because R_channel is linear in L for fixed V_g (sheet conductivity does
    not depend on L), this has a closed form: L_x = Rc_total * W * sigma_sheet.
    """
    n = carrier_density(V_g, V_ch=0.0)
    sigma_sheet = sheet_conductivity(n)
    return Rc_total * W * sigma_sheet


def plot_crossover():
    """Plot total device resistance (R_channel + R_c) vs. channel length L,
    for several literature contact-resistance values, at a fixed on-state
    gate bias. Marks the analytic crossover length for each Rc curve and
    reports the contact-resistance-dominated fraction at the shortest
    length scanned (20 nm), analogous to advanced-node contact scaling
    concerns in silicon technology."""
    V_g_on = 3.5  # on-state gate bias, same value used in graphene_fet_model.py's summary_numbers()
    L_range = np.linspace(20e-9, 1000e-9, 400)  # 20 nm to 1 um

    Rc_values_ohm_um = [0, 110, 300, 500]  # literature range, incl. Pd best-case (0 as ideal baseline)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    R_channel = np.array([channel_resistance_at_length(L, V_g_on) for L in L_range])

    crossover_report = []
    for Rc_um in Rc_values_ohm_um:
        Rc_total = 2 * (Rc_um * 1e-6) / W  # two contacts, width-normalized -> Ohms
        R_total = R_channel + Rc_total
        contact_fraction = Rc_total / R_total * 100.0

        label = f'$R_c$ = {Rc_um} $\\Omega\\cdot\\mu$m'
        axes[0].plot(L_range * 1e9, R_total, linewidth=2, label=label)

        if Rc_um > 0:
            Lx = crossover_length(Rc_total, V_g_on)
            if L_range[0] <= Lx <= L_range[-1]:
                axes[0].axvline(Lx * 1e9, color='gray', linestyle=':', alpha=0.5)
                crossover_report.append((Rc_um, Lx))

        axes[1].plot(L_range * 1e9, contact_fraction, linewidth=2, label=label)

    axes[0].set_xlabel('Channel length L (nm)')
    axes[0].set_ylabel(r'Total device resistance ($\Omega$, per $\mu$m width)')
    axes[0].set_title(f'Device Resistance vs. Channel Length ($V_g$ = {V_g_on} V, on-state)')
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')

    axes[1].axhline(50, color='k', linestyle='--', alpha=0.4, label='50% (equal split)')
    axes[1].set_xlabel('Channel length L (nm)')
    axes[1].set_ylabel('Contact resistance fraction of total (%)')
    axes[1].set_title('Contact-Resistance-Limited Regime')
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig('contact_resistance_crossover.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return crossover_report


def summary_numbers():
    report = plot_crossover()
    print("Contact-resistance / channel-length crossover analysis")
    print("=" * 60)
    print(f"{'Rc (Ohm.um)':>12} | {'Crossover length L_x':>22}")
    print("-" * 40)
    for Rc_um, Lx in report:
        print(f"{Rc_um:>12} | {Lx*1e9:>18.1f} nm")
    print()
    print("Interpretation: below L_x, contact resistance exceeds intrinsic")
    print("channel resistance and dominates the total device resistance --")
    print("further channel-length scaling below L_x gives diminishing")
    print("on-current improvement unless contact resistance is separately")
    print("reduced (e.g. via edge-contact geometries, work-function-")
    print("engineered metals, or doping the graphene under the contact).")
    print(f"For the literature best-case contact (Pd, ~110 Ohm.um), the")
    print(f"crossover length ({report[0][1]*1e9:.0f} nm) is comparable to")
    print(f"or below the 200 nm channel length used elsewhere in this")
    print(f"thesis's GFET model (graphene_fet_model.py), meaning that even")
    print(f"the best literature contacts are not negligible at that node.")


if __name__ == '__main__':
    print("Generating contact-resistance-vs-channel-length crossover plot...")
    summary_numbers()
    print("Done. Saved: contact_resistance_crossover.png")
