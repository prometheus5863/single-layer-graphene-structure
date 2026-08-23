"""
graphene_interconnect_model.py

Device-oriented model of graphene nanoribbon (GNR) interconnect
resistivity as a function of linewidth W, built to complement the GFET
device-physics work already in this repo (graphene_fet_model.py,
rf_small_signal_model.py). Where those scripts cover transistor-level
device physics (contacts, gate capacitance, RF performance), this script
covers the other classic "device application" pillar semiconductor
companies evaluate 2D materials on: back-end-of-line (BEOL) interconnect
resistivity, and specifically how it scales as linewidth is shrunk.

Physical model
---------------
Graphene's linewidth-dependent resistivity is modeled via Matthiessen's
rule combination of a bulk (phonon/substrate-limited) mean free path and
a width-dependent edge-scattering mean free path set by a phenomenological
specularity parameter p (0 = fully diffuse edge scattering, 1 = fully
specular / non-resistive), following the approach used for graphene
nanoribbon interconnects in the literature (see
notes/2026-08-23-interconnect-resistivity-vs-linewidth.md for full
citations and discussion):

    1 / lambda_eff(W) = 1 / lambda_bulk + [(1 - p) / (1 + p)] * (1 / W)

Sheet/3D resistivity is then obtained from a simple Drude-type relation
using this width-dependent mean free path in place of a fixed bulk mean
free path, calibrated so that the W -> infinity (wide-ribbon) limit
reproduces the literature phonon-scattering-limited intrinsic graphene
resistivity of ~1.2 uOhm.cm at n = 5e12 cm^-2 (Murali et al.,
arXiv:0906.0924), and the finite-width behavior is calibrated against the
same paper's measured 15-25 uOhm.cm cluster for 18 nm < W < 52 nm GNRs.

The edge-scattering + residual-impurity curves are calibrated once (at
p = 0.15, W = 22 nm) against the *best individual GNR* result in Murali
et al. (resistivity ~3x the phonon-scattering-limited intrinsic value),
representing a "best-achievable, edge-quality-limited" family of curves.
The *average* measured cluster from the same paper (15-25 uOhm.cm at
W = 18-52 nm) -- reflecting additional real-world process-induced
LER/impurity scattering beyond a clean specularity model -- is overlaid
separately on the plot as an empirical reference band, not used as a fit
target. This is a compact, literature-calibrated analytic model (not an
ab-initio transport calculation) intended to reproduce the qualitative
and order-of-magnitude quantitative trend: graphene resistivity rises
sharply as W is scaled down due to edge scattering, and only becomes
competitive with (eventually, in the idealized single-layer limit,
better than) Cu at aggressively scaled linewidths.

References: see notes/2026-08-23-interconnect-resistivity-vs-linewidth.md
  - Murali et al., "Resistivity of Graphene Nanoribbon Interconnects,"
    arXiv:0906.0924
  - Naeemi & Meindl, "Conductance modeling for graphene nanoribbon (GNR)
    interconnects," IEEE EDL 28, 428 (2007)
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Model parameters (literature-calibrated, see notes/2026-08-23-*.md)
# ---------------------------------------------------------------------------

# Phonon/substrate-scattering-limited intrinsic resistivity of wide-ribbon
# (effectively 2D) graphene on SiO2 at n = 5e12 cm^-2, 300 K.
rho_bulk_uohm_cm = 1.2

# Bulk (width-independent) mean free path implied by rho_bulk. We do not
# need the literal value of lambda_bulk to reproduce rho(W); instead we
# work directly with the mean-free-path *ratio* lambda_eff(W)/lambda_bulk,
# since resistivity ~ 1/lambda for a fixed-density Drude conductor. This
# avoids introducing an extra, poorly-constrained absolute length scale.
#
# lambda_bulk is nonetheless a physically meaningful quantity (room-
# temperature literature values for supported single-layer graphene are
# ~50-60 nm; see e.g. the mean-free-path discussion in
# notes/2026-08-23-*.md). We use 55 nm as a representative value purely
# for reporting/plotting the implied mean free path -- it cancels out of
# the resistivity ratio calculation.
lambda_bulk_nm = 55.0

# Specularity parameter p: 0 = fully diffuse edge scattering (rough,
# plasma/lithographically-etched edges), 1 = fully specular (idealized,
# atomically smooth armchair/zigzag edges). Real lithographically
# patterned GNR edges are diffuse; we use p = 0.15 as a representative
# "realistic etched edge" value and also sweep a range for comparison.
p_default = 0.15

# Additional, width-independent impurity/residual-disorder mean free path,
# calibrated ONCE (at p = p_default, not re-solved per p -- doing so would
# defeat the purpose of a specularity sweep by silently absorbing all of
# the p-dependence into a compensating impurity term) against the "best
# GNR" result in Murali et al.: the best individual GNR at a given width
# had resistivity ~3x the phonon-scattering-limited intrinsic value, i.e.
# best-case, cleanly-patterned samples are edge-scattering-dominated with
# only a small residual impurity contribution. This makes the calibrated
# curves below a "best-achievable, edge-quality-limited" family, distinct
# from the *typical* (average, more disorder-affected) measured cluster of
# 15-25 uOhm.cm at W = 18-52 nm, which is overlaid separately on the plot
# as an empirical reference band rather than used as the fit target -- see
# notes/2026-08-23-*.md for why average real devices sit well above a
# clean specularity-only model (additional process-induced LER/impurity
# scattering not captured by this compact model).
W_calibration_nm = 22.0
rho_calibration_uohm_cm = 3.6  # "best GNR" ~ 3x rho_bulk (Murali et al.)


def _edge_inverse_mfp_ratio(W_nm, p):
    """
    Edge-scattering contribution to the inverse mean free path, expressed
    as a multiple of 1/lambda_bulk: lambda_bulk/lambda_edge.
    """
    return ((1 - p) / (1 + p)) * (lambda_bulk_nm / W_nm)


# Impurity mean free path solved once, at p_default and W_calibration_nm,
# and then held fixed (width- and p-independent) for all curves below.
def _solve_lambda_impurity_nm():
    target_ratio = rho_calibration_uohm_cm / rho_bulk_uohm_cm  # = lambda_bulk/lambda_eff at calibration point
    edge_term = _edge_inverse_mfp_ratio(W_calibration_nm, p_default)
    residual = target_ratio - 1.0 - edge_term
    if residual <= 0:
        raise ValueError(
            "Calibration point implies negative impurity scattering "
            "contribution -- check p_default / calibration target."
        )
    return lambda_bulk_nm / residual


lambda_impurity_nm = _solve_lambda_impurity_nm()


def resistivity_vs_width(W_nm, p=p_default):
    """
    3D resistivity (uOhm.cm) of a graphene nanoribbon of width W_nm (nm),
    for edge specularity p, including the fixed (p-independent) residual
    impurity term calibrated above.
    """
    edge_term = _edge_inverse_mfp_ratio(W_nm, p)
    impurity_term = lambda_bulk_nm / lambda_impurity_nm
    ratio = 1.0 + edge_term + impurity_term
    return rho_bulk_uohm_cm * ratio


def cu_resistivity_vs_width(W_nm):
    """
    Simple projected Cu resistivity-vs-linewidth trend, in the spirit of
    the ITRS-2007-style projection used as the comparison baseline in
    Murali et al. (arXiv:0906.0924): bulk Cu resistivity plus a
    Fuchs-Sondheimer-like surface/grain-boundary scattering rise as W
    shrinks, using a Cu bulk mean free path of ~40 nm and a Fuchs
    specularity parameter of 0.6 (representative literature value; see
    notes/2026-08-23-*.md). This is a simplified single-parameter
    surface-scattering model, not a full Mayadas-Shatzkes grain-boundary
    treatment, and is intended only as an order-of-magnitude comparison
    baseline -- real scaled Cu resistivity is also strongly affected by
    liner/barrier thickness, which does not scale down with W and is not
    modeled here.
    """
    rho_cu_bulk_uohm_cm = 1.68  # bulk Cu resistivity at 300 K
    lambda_cu_bulk_nm = 40.0
    p_cu = 0.6
    # Fuchs-Sondheimer thin-film-style approximate correction factor
    kappa = lambda_cu_bulk_nm / W_nm
    F = 1.0 + (3.0 / 8.0) * (1 - p_cu) * kappa
    return rho_cu_bulk_uohm_cm * F


def plot_resistivity_vs_linewidth():
    """Plot graphene (several specularity values) and Cu resistivity vs.
    linewidth, reproducing the qualitative crossover behavior discussed
    in the notes: graphene is less conductive than Cu at the
    lithographically-relevant widths measured to date (tens of nm), and
    the model's idealized specular-edge (p -> 1) case is required to
    approach Cu-competitive resistivity, consistent with the Naeemi &
    Meindl theoretical crossover being pushed to very narrow widths in
    the idealized single-layer limit. Note this compact model's p=0.9
    (near-specular) curve is already below the Cu model across the whole
    plotted range -- see summary_numbers() / AUTOMATION_LOG.md for the
    realistic (p=0.15, diffuse-edge) crossover width, which is the more
    process-relevant number."""
    W_nm = np.linspace(5, 300, 600)

    fig, ax = plt.subplots(figsize=(9, 6.5))

    for p, style in [(0.0, '--'), (0.15, '-'), (0.5, '-.'), (0.9, ':')]:
        rho = resistivity_vs_width(W_nm, p=p)
        ax.plot(W_nm, rho, style, linewidth=2, label=f'Graphene GNR, p = {p}')

    rho_cu = cu_resistivity_vs_width(W_nm)
    ax.plot(W_nm, rho_cu, color='black', linewidth=2.5, label='Cu (Fuchs-Sondheimer-style projection)')

    # Overlay the literature measured cluster (Murali et al.) as a shaded band
    ax.axvspan(18, 52, color='gray', alpha=0.12, label='Murali et al. measured range (18-52 nm)')
    ax.axhspan(15, 25, color='gray', alpha=0.12)

    ax.set_xlabel('Linewidth W (nm)')
    ax.set_ylabel(r'Resistivity ($\mu\Omega\cdot$cm)')
    ax.set_title('Graphene Nanoribbon Interconnect Resistivity vs. Linewidth')
    ax.set_ylim(0, 30)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('interconnect_resistivity_vs_linewidth.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def summary_numbers():
    """Print representative numbers for sanity-checking / logging."""
    print(f"Fixed residual impurity mean free path (calibrated once at "
          f"p={p_default}, W={W_calibration_nm}nm against the 'best GNR "
          f"~3x bulk limit' result): lambda_impurity = {lambda_impurity_nm:.1f} nm")
    for p in [0.0, p_default, 0.5, 0.9]:
        rho_22 = resistivity_vs_width(22.0, p=p)
        rho_8 = resistivity_vs_width(8.0, p=p)
        print(f"p = {p:.2f}: rho(W=22nm) = {rho_22:.2f} uOhm.cm, "
              f"rho(W=8nm) = {rho_8:.2f} uOhm.cm")

    rho_cu_22 = cu_resistivity_vs_width(22.0)
    rho_cu_8 = cu_resistivity_vs_width(8.0)
    print(f"Cu model: rho(W=22nm) = {rho_cu_22:.2f} uOhm.cm, "
          f"rho(W=8nm) = {rho_cu_8:.2f} uOhm.cm")

    # Find the actual crossover width (if any, within a wide scan range)
    # where each graphene curve drops below the Cu model, via a genuine
    # sign change of (rho_graphene - rho_cu) rather than a naive
    # closest-point search (which can misleadingly report a "crossover"
    # at the edge of the scan range even when the curves never cross).
    W_scan = np.linspace(2, 500, 20000)
    for p, label in [(p_default, "realistic, diffuse edges"), (0.9, "idealized, near-specular edges")]:
        rho_g = resistivity_vs_width(W_scan, p=p)
        rho_cu = cu_resistivity_vs_width(W_scan)
        diff = rho_g - rho_cu
        sign_changes = np.where(np.diff(np.sign(diff)) != 0)[0]
        if len(sign_changes) == 0:
            state = "below" if diff[0] < 0 else "above"
            print(f"p = {p:.2f} ({label}): graphene stays {state} the Cu "
                  f"model across the full W = 2-500 nm scan range (no crossover).")
        else:
            W_cross = W_scan[sign_changes[0]]
            print(f"p = {p:.2f} ({label}): graphene crosses below Cu at "
                  f"approx. W = {W_cross:.0f} nm.")


if __name__ == '__main__':
    print("Generating graphene interconnect resistivity-vs-linewidth model and plot...")
    plot_resistivity_vs_linewidth()
    summary_numbers()
    print("Done. Saved: interconnect_resistivity_vs_linewidth.png")
