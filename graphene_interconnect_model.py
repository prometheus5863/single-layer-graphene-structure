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

2026-08-28 update: the Cu comparison baseline below originally modeled
surface/grain-boundary scattering only, and explicitly noted that this
made it a conservative (upper-bound) estimate of scaled-Cu performance
because it omitted the diffusion-barrier/adhesion-liner layer that a
real damascene Cu interconnect requires -- a layer whose thickness does
not scale down with linewidth and so consumes a growing cross-section
fraction as W shrinks. `cu_resistivity_with_liner()` adds that effect as
a separate, explicit correction (see
notes/2026-08-28-copper-liner-barrier-thickness-effect.md for the
literature basis and derivation), so both the original surface-
scattering-only curve and the new liner-aware curve are available for
comparison.

References: see notes/2026-08-23-interconnect-resistivity-vs-linewidth.md
  - Murali et al., "Resistivity of Graphene Nanoribbon Interconnects,"
    arXiv:0906.0924
  - Naeemi & Meindl, "Conductance modeling for graphene nanoribbon (GNR)
    interconnects," IEEE EDL 28, 428 (2007)

References for the 2026-08-28 liner/barrier addition: see
notes/2026-08-28-copper-liner-barrier-thickness-effect.md
  - Domenichini et al., "Selecting alternative metals for advanced
    interconnects," arXiv:2406.09106
  - "Mechanisms of Scaling Effect for Emerging Nanoscale Interconnect
    Materials," Nanomaterials 12(10), 1760 (2022)

2026-08-29 update: `cu_resistivity_with_liner()` treated the barrier/liner
as strictly non-conducting -- flagged as a simplification in its own
docstring and in notes/2026-08-28-*.md, Section 6.
`cu_resistivity_with_liner_parallel()` replaces that assumption with an
explicit core+liner parallel-conduction model, and is used to evaluate
thinner Ru- and Co-liner Cu scenarios (also left open in the 2026-08-28
notes) alongside the original 3 nm TaN/Co case. See
notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md.
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


# ---------------------------------------------------------------------------
# Cu liner/barrier-thickness effect (added 2026-08-28; see
# notes/2026-08-28-copper-liner-barrier-thickness-effect.md for the
# literature basis and full derivation). The Cu model above is
# surface-scattering-only and was explicitly flagged in the Chapter 5
# draft and this file's original docstring as omitting the diffusion-
# barrier/adhesion-liner layer that a real damascene Cu interconnect
# requires -- a layer whose thickness does not scale down with linewidth
# and therefore consumes a growing fraction of the wire's cross-section
# as W shrinks. This section adds that effect as a separate, explicit
# correction rather than folding it into cu_resistivity_vs_width(), so
# the "surface-scattering-only" and "surface-scattering + liner" curves
# can still be compared directly.
# ---------------------------------------------------------------------------

# Combined barrier (TaN) + adhesion-liner (Co) thickness for a
# conventional Cu damascene stack, per sidewall/dimension. Literature
# value from the Nanomaterials 12(10), 1760 (2022) resistance
# calculations (3 nm for a Cu/TaN-Co stack), consistent with the ~2-3 nm
# functional floor reported by Domenichini et al., arXiv:2406.09106.
t_liner_nm_default = 3.0


def cu_resistivity_with_liner(W_nm, t_liner_nm=t_liner_nm_default):
    """
    Effective Cu resistivity (uOhm.cm) *including* the non-conducting
    barrier/liner's consumption of the drawn cross-section, evaluated
    across a drawn linewidth W_nm (nm) that includes the liner on both
    sides -- i.e. W_nm is the "linewidth (w + 2t)" quantity used in the
    Nanomaterials review this is calibrated against, not the bare Cu
    core width.

    Model (see notes/2026-08-28-*.md, Section 4, for the full
    derivation): treating the liner/barrier as non-conducting and
    assuming the wire's two in-plane dimensions scale together, the
    remaining Cu core has effective width/height W_eff = W - 2*t_liner.
    The reported resistivity combines two compounding effects:

      1. The remaining Cu core is itself narrower, so it is *more*
         surface-scattering-limited than a same-drawn-width, no-liner
         wire would be -- captured by evaluating the existing
         Fuchs-Sondheimer-style cu_resistivity_vs_width() at W_eff
         rather than at W.
      2. That core-limited resistivity is measured across the full
         drawn cross-section (the quantity a process/design actually
         allocates), which dilutes the effective conductance by the
         ratio of drawn area to conducting area, (W / W_eff)^2, for a
         cross-section scaling isotropically in both in-plane
         dimensions.

    Below W = 2*t_liner_nm, the liner/barrier consumes the entire drawn
    cross-section and there is no remaining Cu core -- this is flagged
    explicitly by returning NaN (not a large-but-finite number), since a
    real process simply cannot form a conducting Cu wire in that regime.
    """
    W_arr = np.atleast_1d(np.asarray(W_nm, dtype=float))
    W_eff = W_arr - 2.0 * t_liner_nm

    with np.errstate(divide='ignore', invalid='ignore'):
        rho_core = cu_resistivity_vs_width(np.where(W_eff > 0, W_eff, np.nan))
        area_dilution = (W_arr / np.where(W_eff > 0, W_eff, np.nan)) ** 2
        rho_eff = rho_core * area_dilution

    rho_eff = np.where(W_eff > 0, rho_eff, np.nan)

    if np.isscalar(W_nm) or (hasattr(W_nm, 'ndim') and W_nm.ndim == 0):
        return float(rho_eff[0])
    return rho_eff


# ---------------------------------------------------------------------------
# Parallel-conduction refinement + thinner-liner scenarios (added 2026-08-29;
# see notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md
# for the literature basis and derivation). `cu_resistivity_with_liner()`
# above (2026-08-28) treats the barrier/liner as strictly non-conducting --
# explicitly flagged in that function's own docstring and in the 2026-08-28
# notes (Section 6, first open item) as a simplification, since TaN and
# especially thin Co/Ru liners carry some non-zero current in reality. This
# section replaces that assumption with an explicit parallel-conduction
# (core + liner) model, and uses it to evaluate two liner material/thickness
# scenarios left unexplored on 2026-08-28: thinner Ru- or Co-liner-enabled Cu
# stacks (using the same `t_liner_nm` override mechanism, but now paired
# with each liner's own, much lower resistivity than TaN's -- a materially
# different scenario from just shrinking t_liner_nm while keeping TaN's
# resistivity, since the whole point of switching liner material is that Ru/
# Co liners are far more conductive than TaN).
# ---------------------------------------------------------------------------

# Representative *effective* (thin-film, few-nm-scale) resistivities of
# barrier/liner materials, in uOhm.cm. These are order-of-magnitude
# literature-representative values, not read off a single traceable
# thickness-resolved data table: WebSearch/WebFetch were used this session
# to look for a quantitative Ru/Co/TaN thin-film resistivity-vs-thickness
# source, and found strong *qualitative* confirmation of the ordering used
# here (TaN's barrier resistivity is far higher than Ru's or Co's -- this is
# the standard reason Ru/Co are being pursued as thinner-liner-compatible
# alternatives; see e.g. imec's public Ru/Co interconnect liner work and
# Domenichini et al., arXiv:2406.09106 / pubs.aip.org/aip/jap/article/136/17/171101,
# already cited in notes/2026-08-28-*.md) but no single open-access
# quantitative table for the specific few-nm-thickness values (several
# candidate sources were paywalled or blocked by robots.txt during this
# session's fetch attempts, recorded in today's notes rather than silently
# skipped). The values below are therefore representative figures assembled
# from general semiconductor-interconnect knowledge, exposed as adjustable
# function parameters rather than hard-coded, so a future session can
# tighten them against a specific source without changing the model
# structure.
rho_liner_tan_uohm_cm_default = 400.0   # TaN barrier, ~1-3 nm: strongly
                                         # size-effect-elevated above its
                                         # already-high bulk value (bulk TaN
                                         # is itself only a modest conductor,
                                         # ~125-250 uOhm.cm depending on phase)
rho_liner_ru_uohm_cm_default = 30.0     # Ru liner, ~1-3 nm (bulk Ru ~7.1 uOhm.cm)
rho_liner_co_uohm_cm_default = 20.0     # Co liner, ~1-3 nm (bulk Co ~6.2 uOhm.cm)

# Per-liner default thicknesses, from the Nanomaterials 12(10), 1760 (2022)
# figures already cited in notes/2026-08-28-*.md: Cu/TaN-Co conventional
# stack 3 nm combined, Ru 0.3 nm, Co(-only, thin-liner scheme) 1 nm.
t_liner_nm_ru_default = 0.3
t_liner_nm_co_default = 1.0


def cu_resistivity_with_liner_parallel(W_nm, t_liner_nm=t_liner_nm_default,
                                        rho_liner_uohm_cm=rho_liner_tan_uohm_cm_default):
    """
    Effective Cu resistivity (uOhm.cm) including the liner/barrier's
    cross-section consumption *and* its own (finite, non-zero) parallel
    conduction path -- the refinement flagged as an open item in
    notes/2026-08-28-*.md, Section 6. Same drawn-linewidth convention and
    isotropic-scaling assumption as cu_resistivity_with_liner() (W_nm
    includes the liner on both sides; W_eff = W - 2*t_liner is the
    remaining Cu core).

    Model: treat the drawn cross-section as two conductors in parallel --
    the Cu core (area ~ W_eff^2, resistivity from the existing
    surface-scattering model evaluated at W_eff) and the liner "frame"
    (area ~ W^2 - W_eff^2, resistivity rho_liner_uohm_cm, assumed uniform
    across the whole frame rather than spatially resolved). For a fixed
    length L, conductances add:

        G_total = A_core/(rho_core * L) + A_liner/(rho_liner * L)
        rho_eff = A_drawn / (G_total * L)
                = W^2 / (W_eff^2/rho_core + (W^2 - W_eff^2)/rho_liner)

    This strictly generalizes cu_resistivity_with_liner(): taking the
    rho_liner -> infinity limit recovers exactly the non-conducting-liner
    formula, rho_core * (W/W_eff)^2 (verified numerically in
    notes/2026-08-29-*.md and in this module's __main__ self-check below).
    Below W = 2*t_liner_nm there is still no remaining Cu core, so the
    function returns NaN there for consistency with
    cu_resistivity_with_liner() -- even a perfectly conducting liner alone,
    with no Cu core, is a materially different (and not modeled here)
    all-liner-conductor regime, not a continuation of this Cu-core-centric
    model.
    """
    W_arr = np.atleast_1d(np.asarray(W_nm, dtype=float))
    W_eff = W_arr - 2.0 * t_liner_nm

    with np.errstate(divide='ignore', invalid='ignore'):
        rho_core = cu_resistivity_vs_width(np.where(W_eff > 0, W_eff, np.nan))
        A_core = np.where(W_eff > 0, W_eff, np.nan) ** 2
        A_liner = W_arr ** 2 - np.where(W_eff > 0, W_eff, 0.0) ** 2
        G_total = A_core / rho_core + A_liner / rho_liner_uohm_cm
        rho_eff = (W_arr ** 2) / G_total

    rho_eff = np.where(W_eff > 0, rho_eff, np.nan)

    if np.isscalar(W_nm) or (hasattr(W_nm, 'ndim') and W_nm.ndim == 0):
        return float(rho_eff[0])
    return rho_eff


def plot_resistivity_vs_linewidth():
    """Plot graphene (several specularity values) against *two* Cu
    resistivity models vs. linewidth: the original surface-scattering-
    only baseline, and the surface-scattering + liner/barrier model added
    2026-08-28 (see notes/2026-08-28-copper-liner-barrier-thickness-
    effect.md). Against the surface-scattering-only Cu baseline, realistic
    (p=0.15, diffuse-edge) graphene only wins at wide W (see
    summary_numbers() for the crossover width); against the liner-aware
    Cu model, realistic-edge graphene wins across the entire valid W
    range plotted here -- see summary_numbers() / AUTOMATION_LOG.md for
    both crossover results side by side, which is the more complete,
    process-relevant comparison than either Cu curve alone."""
    W_nm = np.linspace(5, 300, 600)

    fig, ax = plt.subplots(figsize=(9, 6.5))

    for p, style in [(0.0, '--'), (0.15, '-'), (0.5, '-.'), (0.9, ':')]:
        rho = resistivity_vs_width(W_nm, p=p)
        ax.plot(W_nm, rho, style, linewidth=2, label=f'Graphene GNR, p = {p}')

    rho_cu = cu_resistivity_vs_width(W_nm)
    ax.plot(W_nm, rho_cu, color='black', linewidth=2.5, linestyle='--',
            label='Cu, surface scattering only (no liner effect)')

    rho_cu_liner = cu_resistivity_with_liner(W_nm)
    ax.plot(W_nm, rho_cu_liner, color='black', linewidth=2.5,
            label=f'Cu, {t_liner_nm_default:.0f}nm TaN/Co liner (non-conducting approx., 2026-08-28)')

    rho_cu_liner_parallel = cu_resistivity_with_liner_parallel(W_nm)
    ax.plot(W_nm, rho_cu_liner_parallel, color='dimgray', linewidth=2.0, linestyle='-.',
            label=f'Cu, {t_liner_nm_default:.0f}nm TaN liner (parallel conduction, 2026-08-29)')

    rho_cu_liner_ru = cu_resistivity_with_liner_parallel(
        W_nm, t_liner_nm=t_liner_nm_ru_default, rho_liner_uohm_cm=rho_liner_ru_uohm_cm_default)
    ax.plot(W_nm, rho_cu_liner_ru, color='firebrick', linewidth=2.0, linestyle='-.',
            label=f'Cu, {t_liner_nm_ru_default:.1f}nm Ru liner (parallel conduction, 2026-08-29)')

    rho_cu_liner_co = cu_resistivity_with_liner_parallel(
        W_nm, t_liner_nm=t_liner_nm_co_default, rho_liner_uohm_cm=rho_liner_co_uohm_cm_default)
    ax.plot(W_nm, rho_cu_liner_co, color='darkorange', linewidth=2.0, linestyle='-.',
            label=f'Cu, {t_liner_nm_co_default:.1f}nm Co liner (parallel conduction, 2026-08-29)')

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
    print(f"Cu model (surface scattering only): rho(W=22nm) = {rho_cu_22:.2f} uOhm.cm, "
          f"rho(W=8nm) = {rho_cu_8:.2f} uOhm.cm")

    rho_cu_liner_22 = cu_resistivity_with_liner(22.0)
    rho_cu_liner_12 = cu_resistivity_with_liner(12.0)
    print(f"Cu model (surface scattering + {t_liner_nm_default:.0f}nm liner/barrier, "
          f"added 2026-08-28): rho(W=22nm) = {rho_cu_liner_22:.2f} uOhm.cm, "
          f"rho(W=12nm) = {rho_cu_liner_12:.2f} uOhm.cm (W=8nm is below the "
          f"W = 2*t_liner = {2*t_liner_nm_default:.0f}nm floor -- no Cu core remains, "
          f"model returns NaN as flagged in cu_resistivity_with_liner()'s docstring)")

    # Find the actual crossover width (if any, within a wide scan range)
    # where each graphene curve drops below a given Cu model, via a
    # genuine sign change of (rho_graphene - rho_cu) rather than a naive
    # closest-point search (which can misleadingly report a "crossover"
    # at the edge of the scan range even when the curves never cross).
    # Repeated for both the original surface-scattering-only Cu model and
    # the 2026-08-28 liner-aware Cu model, so the shift in crossover width
    # caused specifically by the liner effect is visible directly.
    def _report_crossovers(cu_model_fn, cu_model_label):
        W_scan = np.linspace(2, 500, 20000)
        rho_cu_scan = cu_model_fn(W_scan)
        valid = ~np.isnan(rho_cu_scan)
        for p, label in [(p_default, "realistic, diffuse edges"), (0.9, "idealized, near-specular edges")]:
            rho_g = resistivity_vs_width(W_scan, p=p)
            diff = np.where(valid, rho_g - rho_cu_scan, np.nan)
            finite = np.isfinite(diff)
            sign_changes = np.where(np.diff(np.sign(diff[finite])) != 0)[0]
            W_valid = W_scan[finite]
            if len(sign_changes) == 0:
                state = "below" if diff[finite][0] < 0 else "above"
                print(f"  p = {p:.2f} ({label}): graphene stays {state} {cu_model_label} "
                      f"across the full valid scan range (no crossover).")
            else:
                W_cross = W_valid[sign_changes[0]]
                print(f"  p = {p:.2f} ({label}): graphene crosses below {cu_model_label} at "
                      f"approx. W = {W_cross:.0f} nm.")

    print("Crossover vs. Cu, surface scattering only:")
    _report_crossovers(cu_resistivity_vs_width, "the surface-scattering-only Cu model")
    print("Crossover vs. Cu, surface scattering + liner/barrier (2026-08-28):")
    _report_crossovers(cu_resistivity_with_liner, "the liner-aware Cu model")

    print("Crossover vs. Cu, parallel-conduction TaN liner (2026-08-29):")
    _report_crossovers(cu_resistivity_with_liner_parallel, "the parallel-conduction TaN-liner Cu model")

    def _ru_model(W):
        return cu_resistivity_with_liner_parallel(
            W, t_liner_nm=t_liner_nm_ru_default, rho_liner_uohm_cm=rho_liner_ru_uohm_cm_default)

    def _co_model(W):
        return cu_resistivity_with_liner_parallel(
            W, t_liner_nm=t_liner_nm_co_default, rho_liner_uohm_cm=rho_liner_co_uohm_cm_default)

    print(f"Crossover vs. Cu, {t_liner_nm_ru_default:.1f}nm Ru liner, parallel conduction (2026-08-29):")
    _report_crossovers(_ru_model, "the Ru-liner Cu model")
    print(f"Crossover vs. Cu, {t_liner_nm_co_default:.1f}nm Co liner, parallel conduction (2026-08-29):")
    _report_crossovers(_co_model, "the Co-liner Cu model")


def self_check_parallel_limit():
    """
    Sanity check (2026-08-29): cu_resistivity_with_liner_parallel() should
    reduce to cu_resistivity_with_liner() in the rho_liner -> infinity
    (non-conducting) limit. Verified numerically here (not just claimed in
    the docstring) at a representative width, using a very large but finite
    liner resistivity as a stand-in for infinity.
    """
    W_test = 20.0
    rho_nonconducting = cu_resistivity_with_liner(W_test)
    rho_parallel_huge_rho = cu_resistivity_with_liner_parallel(
        W_test, rho_liner_uohm_cm=1.0e12)
    rel_diff = abs(rho_parallel_huge_rho - rho_nonconducting) / rho_nonconducting
    status = "OK" if rel_diff < 1e-6 else "FAILED"
    print(f"Self-check: parallel-conduction model at rho_liner->inf matches "
          f"non-conducting model at W={W_test:.0f}nm: "
          f"{rho_nonconducting:.4f} vs {rho_parallel_huge_rho:.4f} uOhm.cm "
          f"(rel. diff {rel_diff:.2e}) -- {status}")
    if status == "FAILED":
        raise AssertionError("Parallel-conduction limiting-case self-check failed")

    # A second check in the opposite direction: with a liner resistivity
    # *equal to* the Cu core's own resistivity at that width, the parallel
    # model should reduce to the plain area-averaged (liner indistinguishable
    # from core) resistivity -- i.e. rho_eff should equal rho_core exactly,
    # independent of t_liner_nm, since both "conductors" are then identical.
    rho_core_at_W = cu_resistivity_vs_width(W_test - 2.0 * t_liner_nm_default)
    rho_parallel_matched = cu_resistivity_with_liner_parallel(
        W_test, rho_liner_uohm_cm=rho_core_at_W)
    rel_diff2 = abs(rho_parallel_matched - rho_core_at_W) / rho_core_at_W
    status2 = "OK" if rel_diff2 < 1e-6 else "FAILED"
    print(f"Self-check: parallel-conduction model at rho_liner == rho_core "
          f"matches rho_core exactly: {rho_core_at_W:.4f} vs "
          f"{rho_parallel_matched:.4f} uOhm.cm (rel. diff {rel_diff2:.2e}) -- {status2}")
    if status2 == "FAILED":
        raise AssertionError("Parallel-conduction matched-resistivity self-check failed")


if __name__ == '__main__':
    print("Generating graphene interconnect resistivity-vs-linewidth model and plot...")
    self_check_parallel_limit()
    plot_resistivity_vs_linewidth()
    summary_numbers()
    print("Done. Saved: interconnect_resistivity_vs_linewidth.png")
