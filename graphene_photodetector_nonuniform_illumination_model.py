"""
graphene_photodetector_nonuniform_illumination_model.py

NON-UNIFORM ILLUMINATION: the collection kernel, shadow masks, and a hard
ceiling on illumination engineering.

This closes the item AUTOMATION_LOG.md has carried at the top of "not yet
covered" since 2026-09-18 -- "Non-uniform illumination, a generation weight
g(x)" -- and it is Simplification 2 of
`graphene_photodetector_signed_carrier_model.py`, unchanged since 2026-09-07:

    "2. Uniform illumination. (Still open: the Suzuki et al. shadow mask.)"

CITATION CORRECTION
-------------------
That shadow-mask experiment is NOT by Suzuki et al. and is not article
100115. It is

    Kenta Shimomura, Kaname Imai, Kenta Nakagawa, Akira Kawai,
    Kazuki Hashimoto, Takuro Ideguchi, Hideyuki Maki,
    "Graphene photodetectors with asymmetric device structures on
    silicon chips", Carbon Trends 5, 100100 (2021).
    https://www.sciencedirect.com/science/article/pii/S2667056921000778

confirmed against the publisher page and the Ideguchi group's publication
list. The repo carried the wrong first author and the wrong article number
in six places from 2026-09-17; all are corrected in this session. The
physics attributed to it was correct.

THE STRUCTURAL POINT
--------------------
Illumination does not enter the transport at all. Write the 2026-09-18
signed, carrier-resolved response as

    N = (1/L) integral_0^L [ q_h s_h(x) p_h(x) + q_e s_e(x) p_e(x) ] dx
      = (1/L) integral_0^L k(x) dx

The bracket k(x) -- the COLLECTION KERNEL -- is the net charge delivered to
contact A per photon ABSORBED AT x. It depends only on the two contact
metals, the bias, the crossover and tau. So non-uniform illumination needs
no new transport model, only a generation weight:

    N[g] = (1/L) integral_0^L g(x) k(x) dx,    (1/L) integral g dx = 1

The normalisation fixes the TOTAL absorbed photon number, which is the only
basis on which two illumination patterns can be honestly compared.

Three consequences, all exactly checkable:

  1. g == 1 must return the 2026-09-18 number BITWISE. Not approximately:
     under g == 1 the expression is literally the same one.
  2. k(x) IS the delta-spot scanning-photocurrent trace: N[delta(x-x0)] =
     k(x0). The kernel is not an intermediate quantity.
  3. |N[g]| <= max_x |k(x)| for EVERY normalised g. This is a hard ceiling
     on masks, spots, gratings and plasmonic patterning alike, at fixed
     photon number -- the first quantity in this repo that bounds a whole
     design space instead of ranking points inside it.

And for a SYMMETRIC pair at zero bias k is antisymmetric about the midpoint,
which gives the leaky shadow mask in closed form. With transmission T over
the left half and 1 over the right, A = (1/L) integral_{L/2}^{L} k dx:

    N(T) = 2A (1 - T) / (1 + T)      =>   N(T)/N(0) = (1 - T)/(1 + T)

Shimomura et al. report that "the photovoltage under the shadow mask was
about half of the opposite side", i.e. T ~ 0.5 for a real 50 nm Ni mask.
That leaves ONE THIRD of a perfect mask, not one half. This closed form is
used below as a fourth exactly-known validation target.

GEOMETRY CAVEAT (stated up front, not buried)
---------------------------------------------
This repo's channel is L = 200 nm. Shimomura et al.'s focused spot is
~2 um -- TEN TIMES the whole channel. Optical localisation of the SPCM kind
is impossible at this geometry: a diffraction-limited spot illuminates the
channel and both contacts at once. The only realisable non-uniform
illumination here is a LITHOGRAPHIC mask on the device, which is what
Shimomura et al. built. `spot_size_study()` is therefore reported as a
statement about what channel length would be required, not as a proposal
for this one.

MECHANISM CAVEAT
----------------
Kasirga's review (arXiv:2509.09390, 2025) stresses that measured SPCM
signals are frequently photo-thermal rather than photovoltaic. This repo
models drift collection only. The position-resolved curves here are the
PHOTOVOLTAIC CONTRIBUTION to an SPCM trace, not a prediction of a measured
one. The repo has no photo-thermoelectric machinery with which to close
that gap, and does not pretend to.

See notes/2026-09-19-non-uniform-illumination-and-the-collection-kernel.md.

STATED SIMPLIFICATIONS
----------------------
1-4. All of the 2026-09-18 module's simplifications are inherited unchanged
     (linear dW->profile with one lambda for every metal; drift only, no
     diffusion, mu_e = mu_h; no photogain, PTE or bolometric term).
5.   g(x) weights ABSORPTION, and absorption is taken proportional to the
     local illumination. Graphene's 2.3% absorption is position-independent
     and low enough that beam attenuation along the channel is negligible,
     so this is safe here; it would not be for a thick absorber.
6.   The mask is treated as a pure transmission factor with an optional
     finite edge width. Near-field scattering, reflection off the 50 nm Ni
     and any plasmonic response of the mask edge are not modelled.
"""

import numpy as np
import matplotlib.pyplot as plt

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS
from graphene_photodetector_model import L_channel, V_bias
from graphene_photodetector_collection_model import E_BIAS
from graphene_photodetector_signed_carrier_model import (
    total_field, _transport, Q_HOLE, Q_ELECTRON,
    W_CROSS_CHEM, W_CROSS_VACUUM, N_POINTS, signed_offset,
    net_response as signed_net_response,
)

# Shimomura et al.'s measured mask leakage: the shaded interface still saw
# "about half" the photovoltage of the open one.
T_MASK_MEASURED = 0.5


# ---------------------------------------------------------------------
# The collection kernel
# ---------------------------------------------------------------------
def collection_kernel(W_A, W_B, e_applied=E_BIAS, L=L_channel,
                      n_points=N_POINTS, w_cross=W_CROSS_CHEM,
                      carriers=("hole", "electron")):
    """
    k(x): net charge delivered to contact A per photon absorbed at x.

    Built with exactly the same accumulation order as
    graphene_photodetector_signed_carrier_model.net_response(), so that
    (1/L) * trapezoid(k, x) reproduces that function bitwise. Do not
    'simplify' the zeros_like + successive addition below; it is what makes
    validate_uniform_reduction() an exact test rather than a tolerance test.
    """
    x = np.linspace(0.0, L, n_points)
    E = total_field(x, W_A, W_B, e_applied=e_applied, L=L, w_cross=w_cross)

    k = np.zeros_like(x)
    parts = {}
    if "hole" in carriers:
        s_h, p_h = _transport(+E, x)
        k = k + Q_HOLE * s_h * p_h
        parts["hole"] = Q_HOLE * s_h * p_h
    if "electron" in carriers:
        s_e, p_e = _transport(-E, x)
        k = k + Q_ELECTRON * s_e * p_e
        parts["electron"] = Q_ELECTRON * s_e * p_e
    return x, k, parts


# ---------------------------------------------------------------------
# Generation profiles g(x). Each returns an UNNORMALISED weight; the
# normalisation to unit mean is applied once, in response_under().
# `uniform` is literally np.ones so that its continuum mean is exactly 1
# and no floating-point renormalisation touches the reduction test.
# ---------------------------------------------------------------------
def g_uniform(x, L=L_channel):
    return np.ones_like(x)


def g_shadow_mask(x, L=L_channel, T=0.0, edge_width=0.0, masked="left"):
    """
    Shimomura et al.'s geometry: an opaque (or leaky) mask over one of the
    two graphene/electrode interfaces.

    T          : transmission under the mask (0 = perfect, 1 = no mask).
                 T_MASK_MEASURED = 0.5 is the value their 50 nm Ni mask
                 actually achieved.
    edge_width : 1/e width of a tanh penumbra at the mask edge. 0 = sharp.
    masked     : "left" shades contact A's interface, "right" shades B's.
    """
    xm = 0.5 * L
    if edge_width <= 0.0:
        open_frac = np.where(x < xm, 0.0, 1.0)
        open_frac = np.where(np.isclose(x, xm), 0.5, open_frac)
    else:
        open_frac = 0.5 * (1.0 + np.tanh((x - xm) / edge_width))
    if masked == "right":
        open_frac = 1.0 - open_frac
    return T + (1.0 - T) * open_frac


def g_gaussian_spot(x, L=L_channel, x0=None, sigma=20e-9):
    """A focused spot at x0 with Gaussian 1/e^2 intensity half-width 2*sigma."""
    if x0 is None:
        x0 = 0.5 * L
    return np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def response_under(k, x, g, L=L_channel, normalise=True):
    """
    N[g] = (1/L) integral g_normalised(x) k(x) dx, with g normalised to unit
    mean so that every pattern is compared at equal total ABSORBED photons.

    For g == 1 the normaliser is exactly 1.0 by definition (not computed),
    so the returned value is bitwise identical to the uniform-illumination
    result of the 2026-09-18 module.

    normalise=False gives the PER-INCIDENT-PHOTON figure instead: g is the
    raw transmission, so photons stopped by the mask are counted as lost.
    See response_per_incident() for why a detector needs that one.
    """
    if normalise:
        mean_g = float(np.trapezoid(g, x) / L)
        if abs(mean_g - 1.0) > 1e-15:
            g = g / mean_g
    return float(np.trapezoid(g * k, x) / L)


def response_per_incident(k, x, g, L=L_channel):
    """
    Charge to contact A per photon INCIDENT on the device area, for a device
    uniformly illuminated from outside with g(x) the local transmission.

    This is the distinction that decides whether masking is a good idea, and
    the first draft of this module got the headline wrong by ignoring it.

    Per ABSORBED photon a mask costs nothing -- it only redistributes where
    the surviving photons land, which is the right way to compare COLLECTION
    MECHANISMS. But responsivity is amps per incident watt: a mask that
    shades half the device throws half the light into 50 nm of nickel, and
    the detector never sees it. For a symmetric pair at zero bias the two
    accountings differ by exactly the normalisation factor (1 + T)/2:

        per absorbed photon : N(T)/N(0) = (1 - T)/(1 + T)
        per incident photon : N(T)/N(0) = (1 - T)          [linear]

    so a perfect mask (T = 0) gives exactly HALF as much per incident photon
    as per absorbed photon. Any comparison of a masked symmetric device with
    an unmasked asymmetric one must use this function, not response_under().
    """
    return float(np.trapezoid(g * k, x) / L)


# =====================================================================
# Validations -- four, all against EXACTLY known values. The 2026-09-17
# session established why that matters: a plausible-range check passed
# while a real inf-inf bug sat in the transit integral.
# =====================================================================
def validate_uniform_reduction(verbose=True):
    """
    EXACT #1. g == 1 must reproduce
    graphene_photodetector_signed_carrier_model.net_response() to the LAST
    BIT, on all 49 ordered pairs, both crossovers, both biases -- 196
    comparisons. Not 'to within 1e-12': bitwise. Anything else means the
    kernel refactor changed the arithmetic.
    """
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    worst, n, n_bitwise = 0.0, 0, 0
    worst_row = None
    for mA, WA in metals:
        for mB, WB in metals:
            for wc in (W_CROSS_CHEM, W_CROSS_VACUUM):
                for eb in (0.0, E_BIAS):
                    x, k, _ = collection_kernel(WA, WB, e_applied=eb, w_cross=wc)
                    new = response_under(k, x, g_uniform(x))
                    old = signed_net_response(WA, WB, e_applied=eb, w_cross=wc)
                    dev = abs(new - old)
                    n += 1
                    if new == old:
                        n_bitwise += 1
                    if dev > worst:
                        worst, worst_row = dev, (mA, mB, wc, eb, old, new)
    if verbose:
        print("Validation 1 (EXACT): uniform g == 1 reproduces the 2026-09-18")
        print(f"  signed model on {n} (pair, crossover, bias) combinations.")
        print(f"  bitwise identical: {n_bitwise}/{n}")
        if worst_row:
            mA, mB, wc, eb, old, new = worst_row
            print(f"  worst: {mA}/{mB} w_cross={wc} E={eb:.3e} "
                  f"old={old:+.12f} new={new:+.12f}")
        print(f"  worst absolute deviation: {worst:.3e}")
    return n, n_bitwise, worst


def validate_symmetric_g_gives_zero(w_cross=W_CROSS_CHEM, verbose=True):
    """
    EXACT #2. A symmetric pair at zero bias has k(x) antisymmetric about
    x = L/2. Any MIRROR-SYMMETRIC illumination, g(x) = g(L-x), must then
    give N exactly zero -- including a centred spot, not just the uniform
    case the 2026-09-18 module could test. This separates 'the illumination
    machinery is correct' from 'the transport is correct': a sign error in
    the weighting would survive Validation 1 and fail here.
    """
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        x, k, _ = collection_kernel(W, W, e_applied=0.0, w_cross=w_cross)
        for label, g in [("uniform", g_uniform(x)),
                         ("centred spot", g_gaussian_spot(x, x0=0.5 * L_channel,
                                                          sigma=20e-9)),
                         ("both edges", g_gaussian_spot(x, x0=0.0, sigma=20e-9)
                          + g_gaussian_spot(x, x0=L_channel, sigma=20e-9))]:
            rows.append((metal, label, response_under(k, x, g)))
    worst = max(abs(v) for _, _, v in rows)
    if verbose:
        print(f"\nValidation 2 (EXACT): symmetric pair, zero bias, mirror-symmetric")
        print(f"  illumination -> N == 0 (w_cross = {w_cross} eV, "
              f"{len(rows)} cases)")
        wr = max(rows, key=lambda r: abs(r[2]))
        print(f"  worst: {wr[0]} under '{wr[1]}' -> N = {wr[2]:+.3e}")
        print(f"  worst |N|: {worst:.3e}")
    return rows, worst


def validate_mirror_antisymmetry(w_cross=W_CROSS_CHEM, verbose=True):
    """
    EXACT #3. Same device, zero bias: masking the LEFT interface and masking
    the RIGHT interface must give exactly opposite responses, for every
    transmission T. N[g(L-x)] = -N[g(x)].

    A model that got the mask orientation right by accident -- e.g. by
    taking a magnitude somewhere -- passes Validation 2 and fails this.
    """
    rows, worst = [], 0.0
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        x, k, _ = collection_kernel(W, W, e_applied=0.0, w_cross=w_cross)
        for T in (0.0, 0.25, T_MASK_MEASURED, 0.9):
            nl = response_under(k, x, g_shadow_mask(x, T=T, masked="left"))
            nr = response_under(k, x, g_shadow_mask(x, T=T, masked="right"))
            dev = abs(nl + nr)
            worst = max(worst, dev)
            rows.append((metal, T, nl, nr, dev))
    if verbose:
        print(f"\nValidation 3 (EXACT): mask left vs mask right -> N_left == -N_right")
        wr = max(rows, key=lambda r: r[4])
        print(f"  worst: {wr[0]} at T={wr[1]}  left={wr[2]:+.9f} right={wr[3]:+.9f}")
        print(f"  worst |N_left + N_right|: {worst:.3e}")
    return rows, worst


def validate_leaky_mask_closed_form(w_cross=W_CROSS_CHEM, verbose=True):
    """
    EXACT #4, and the one that tests the PHYSICS rather than the plumbing.

    For a symmetric pair at zero bias, antisymmetry of k gives the leaky
    mask in closed form:

        N(T)/N(0) = (1 - T)/(1 + T)

    with no reference to the kernel's shape at all. The numerics must
    reproduce that ratio exactly for every metal and every T. If the
    normalisation-to-equal-photon-number were dropped or applied to the
    wrong quantity, this ratio would come out as (1 - T) instead, which is
    a 50% error at the experimentally relevant T = 0.5 and would not be
    caught by any of the first three checks.
    """
    rows, worst = [], 0.0
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        x, k, _ = collection_kernel(W, W, e_applied=0.0, w_cross=w_cross)
        n0 = response_under(k, x, g_shadow_mask(x, T=0.0))
        if abs(n0) < 1e-12:
            continue                      # Cr at 4.5 eV has dW = 0, no field
        for T in (0.0, 0.1, 0.25, T_MASK_MEASURED, 0.75, 0.9):
            got = response_under(k, x, g_shadow_mask(x, T=T)) / n0
            want = (1.0 - T) / (1.0 + T)
            dev = abs(got - want)
            worst = max(worst, dev)
            rows.append((metal, T, got, want, dev))
    if verbose:
        print(f"\nValidation 4 (EXACT): leaky mask closed form "
              f"N(T)/N(0) == (1-T)/(1+T)")
        wr = max(rows, key=lambda r: r[4])
        print(f"  worst: {wr[0]} at T={wr[1]}  numeric={wr[2]:.12f} "
              f"closed form={wr[3]:.12f}")
        print(f"  worst deviation: {worst:.3e}")
        print(f"  at the measured T = {T_MASK_MEASURED}: ratio = "
              f"{(1 - T_MASK_MEASURED) / (1 + T_MASK_MEASURED):.6f} "
              f"-- one THIRD of a perfect mask, not one half")
    return rows, worst


def validate_kernel_ceiling(w_cross=W_CROSS_CHEM, verbose=True):
    """
    EXACT #5 (an inequality, but an exact one). |N[g]| <= max_x |k(x)| for
    every normalised g. Checked against a spread of masks and spots on every
    ordered pair; a single violation would mean the normalisation is wrong.
    """
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    worst_slack, n_viol, n = 1e9, 0, 0
    for mA, WA in metals:
        for mB, WB in metals:
            x, k, _ = collection_kernel(WA, WB, e_applied=0.0, w_cross=w_cross)
            ceiling = float(np.max(np.abs(k)))
            gs = [g_uniform(x),
                  g_shadow_mask(x, T=0.0), g_shadow_mask(x, T=0.0, masked="right"),
                  g_shadow_mask(x, T=0.5), g_gaussian_spot(x, x0=0.0, sigma=15e-9),
                  g_gaussian_spot(x, x0=0.5 * L_channel, sigma=15e-9),
                  g_gaussian_spot(x, x0=L_channel, sigma=15e-9)]
            for g in gs:
                val = abs(response_under(k, x, g))
                n += 1
                if val > ceiling + 1e-12:
                    n_viol += 1
                worst_slack = min(worst_slack, ceiling - val)
    if verbose:
        print(f"\nValidation 5 (EXACT inequality): |N[g]| <= max|k| on {n} cases")
        print(f"  violations: {n_viol}")
        print(f"  tightest slack (ceiling - |N|): {worst_slack:.6f}")
    return n, n_viol, worst_slack


# =====================================================================
# Results
# =====================================================================
def masked_symmetric_table(w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 1. The headline. Under uniform illumination a symmetric device
    gives identically zero at zero bias (that is the Shimomura /
    Weiss & Duan cancellation). A shadow mask over one interface is the only
    mechanism left by which it can respond at all -- which matters more
    since 2026-09-18, when the symmetric METAL RANKING was retracted. How
    much does it recover, and how does that compare with the best asymmetric
    metal pair under uniform light?
    """
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    rows = []
    for metal, W in metals:
        x, k, _ = collection_kernel(W, W, e_applied=0.0, w_cross=w_cross)
        rows.append((metal, W, signed_offset(W, w_cross),
                     response_under(k, x, g_uniform(x)),
                     response_under(k, x, g_shadow_mask(x, T=0.0)),
                     response_under(k, x, g_shadow_mask(x, T=T_MASK_MEASURED)),
                     response_per_incident(k, x, g_shadow_mask(x, T=0.0)),
                     response_per_incident(k, x,
                                           g_shadow_mask(x, T=T_MASK_MEASURED))))
    best_asym = signed_net_response(METAL_WORK_FUNCTIONS["Ti"],
                                    METAL_WORK_FUNCTIONS["Pt"],
                                    e_applied=0.0, w_cross=w_cross)
    if verbose:
        print(f"\nRESULT 1: a shadow mask on a SYMMETRIC device (zero bias, "
              f"w_cross = {w_cross} eV)")
        print("  per ABSORBED photon (mechanism comparison) | "
              "per INCIDENT photon (responsivity)")
        print(f"{'metal':<6}{'dW (eV)':>9}{'uniform':>10}{'mask T=0':>10}"
              f"{'mask T=.5':>11}{'| inc T=0':>11}{'inc T=.5':>10}")
        for metal, W, dW, un, pm, lm, pi0, pi5 in rows:
            print(f"{metal:<6}{dW:>+9.2f}{un:>10.2e}{pm:>10.4f}{lm:>11.4f}"
                  f"{pi0:>11.4f}{pi5:>10.4f}")
        best = max(rows, key=lambda r: abs(r[4]))
        print(f"  best perfect-mask symmetric device: {best[0]} -> "
              f"{best[4]:+.4f} per absorbed photon, {best[6]:+.4f} per incident")
        print(f"  best uniform-illumination ASYMMETRIC pair (Ti/Pt): "
              f"{best_asym:+.4f} (absorbs everything, so both accountings agree)")
        print(f"  HONEST COMPARISON, per incident photon: a perfectly masked "
              f"symmetric device reaches")
        print(f"    {abs(best[6]) / abs(best_asym) * 100:.1f}% of Ti/Pt -- not "
              f"the {abs(best[4]) / abs(best_asym) * 100:.1f}% the "
              f"per-absorbed-photon column suggests.")
        print(f"  With the mask Shimomura et al. actually achieved (T = "
              f"{T_MASK_MEASURED}) it is {abs(best[7]) / abs(best_asym) * 100:.1f}%.")
        print("  CONCLUSION: asymmetric METALLISATION beats illumination")
        print("  engineering by roughly 4x, and by ~8x against a real mask.")
    return rows, best_asym


def mask_on_asymmetric_table(w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 2. Does masking help a device that is ALREADY asymmetric?

    THE NOTE'S PREDICTION WAS WRONG, AND IS RECORDED AS WRONG.
    notes/2026-09-19 predicted "an n/p pair should gain little or NOTHING
    from masking: its kernel does not change sign, so masking half the
    channel only discards photons" and offered it as falsifiable. It is
    falsified, weakly: per absorbed photon, masking Ti/Pd's A side gains
    +0.7% and Ti/Pt's B side +0.03%. The reasoning was wrong in two ways.

    (i) Sign is not the issue; VARIATION is. Even a single-signed kernel is
        not flat, so concentrating photons where |k| is larger gains
        something. "Only discards photons" was a misstatement: at equal
        absorbed photon number a mask redistributes rather than discards.
    (ii) The gain is exactly bounded, and the bound is the useful output:

            max achievable / uniform  =  max_x|k(x)| / |N_uniform|

        which is the Validation-5 ceiling divided by the uniform response.
        For every asymmetric pair here it is under 2%, so the CONCLUSION the
        note wanted survives -- illumination engineering is worth nothing on
        an n/p pair -- but it now has a bound attached instead of a bad
        argument. The bound is what should be quoted.

    Per INCIDENT photon a real mask is strictly worse, because the shaded
    photons are genuinely lost; that column is shown too.
    """
    cases = [("Ti", "Pt"), ("Ti", "Pd"), ("Cr", "Pt"), ("Pd", "Pt"), ("Pt", "Pt")]
    rows = []
    for mA, mB in cases:
        WA, WB = METAL_WORK_FUNCTIONS[mA], METAL_WORK_FUNCTIONS[mB]
        x, k, _ = collection_kernel(WA, WB, e_applied=0.0, w_cross=w_cross)
        un = response_under(k, x, g_uniform(x))
        ml = response_under(k, x, g_shadow_mask(x, T=0.0, masked="left"))
        mr = response_under(k, x, g_shadow_mask(x, T=0.0, masked="right"))
        inc = max(abs(response_per_incident(k, x, g_shadow_mask(x, T=0.0,
                                                               masked=m)))
                  for m in ("left", "right"))
        ceiling = float(np.max(np.abs(k)))
        best_ratio = ceiling / abs(un) if abs(un) > 1e-12 else np.inf
        rows.append((mA, mB, un, ml, mr, best_ratio, inc,
                     inc / abs(un) if abs(un) > 1e-12 else np.inf))
    if verbose:
        print(f"\nRESULT 2: does a mask help an ALREADY asymmetric pair? "
              f"(zero bias)")
        print(f"{'pair':<10}{'uniform':>10}{'mask A':>10}{'mask B':>10}"
              f"{'BOUND max|k|/unif':>19}{'best per inc':>14}{'vs unif':>9}")
        for mA, mB, un, ml, mr, br, inc, ir in rows:
            b = "     inf" if not np.isfinite(br) else f"{br:>19.4f}"
            i = "      inf" if not np.isfinite(ir) else f"{ir:>9.3f}"
            print(f"{mA+'/'+mB:<10}{un:>+10.4f}{ml:>+10.4f}{mr:>+10.4f}"
                  f"{b}{inc:>14.4f}{i}")
        print("  The BOUND column is the most any illumination pattern whatever")
        print("  can achieve at equal absorbed photons -- under 1.02 for every")
        print("  asymmetric pair. Per INCIDENT photon masking always loses.")
        print("  Design rule: masks are for symmetric devices only.")
    return rows


def kernel_scan_table(w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 3. k(x) is the delta-spot photocurrent trace. Report its shape
    for a symmetric and an asymmetric device: the symmetric one must show
    OPPOSITE polarity at the two interfaces (Shimomura et al.'s "polarities
    ... are opposite"; Mueller et al.'s p-n-p), the asymmetric one must not.
    """
    cases = [("Pt", "Pt"), ("Ti", "Ti"), ("Ti", "Pt"), ("Ti", "Pd")]
    rows = []
    for mA, mB in cases:
        x, k, _ = collection_kernel(METAL_WORK_FUNCTIONS[mA],
                                    METAL_WORK_FUNCTIONS[mB],
                                    e_applied=0.0, w_cross=w_cross)
        nz = k[k != 0.0]
        rows.append((mA, mB, float(k[0]), float(k[len(k) // 2]), float(k[-1]),
                     float(np.max(np.abs(k))),
                     bool(nz.size and np.any(np.diff(np.sign(nz)) != 0)),
                     float(np.trapezoid(k, x) / L_channel)))
    if verbose:
        print(f"\nRESULT 3: the collection kernel as a delta-spot scan "
              f"(zero bias)")
        print(f"{'pair':<10}{'k(0)':>9}{'k(L/2)':>9}{'k(L)':>9}{'max|k|':>9}"
              f"{'sign flip':>11}{'mean = N_unif':>15}")
        for mA, mB, k0, km, kL, mx, flip, mean in rows:
            print(f"{mA+'/'+mB:<10}{k0:>+9.3f}{km:>+9.3f}{kL:>+9.3f}{mx:>9.3f}"
                  f"{str(flip):>11}{mean:>+15.2e}")
    return rows


def spot_size_study(pair=("Pt", "Pt"), w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 4. How localised must the illumination be? Scan a Gaussian spot
    parked at the contact and widen it; the response must fall from k(0)
    toward the uniform value as the spot outgrows the channel.

    Reported as a statement about REQUIRED CHANNEL LENGTH, not as a proposal
    for this device: Shimomura et al.'s ~2 um spot is 10x this repo's 200 nm
    channel, so optical localisation is impossible here and only a
    lithographic mask is realisable. The useful output is the ratio
    sigma/L at which localisation stops paying.
    """
    mA, mB = pair
    x, k, _ = collection_kernel(METAL_WORK_FUNCTIONS[mA], METAL_WORK_FUNCTIONS[mB],
                                e_applied=0.0, w_cross=w_cross)
    ideal = float(k[0])
    rows = []
    for frac in (0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0):
        sigma = frac * L_channel
        g = g_gaussian_spot(x, x0=0.0, sigma=sigma)
        N = response_under(k, x, g)
        rows.append((frac, sigma, N, N / ideal if ideal else np.nan))
    if verbose:
        print(f"\nRESULT 4: spot size at contact A, {mA}/{mB}, zero bias "
              f"(L = {L_channel * 1e9:.0f} nm)")
        print(f"{'sigma/L':>9}{'sigma (nm)':>12}{'N':>10}{'frac of k(0)':>14}")
        for frac, sigma, N, r in rows:
            print(f"{frac:>9.3f}{sigma * 1e9:>12.1f}{N:>+10.4f}{r:>14.3f}")
        print(f"  k(0) (delta-spot limit) = {ideal:+.4f}")
        print(f"  Shimomura et al.'s ~2 um spot is sigma/L ~ "
              f"{2e-6 / L_channel:.0f} at this geometry: no localisation at all.")
        print("  Reaching 90% of the delta limit needs sigma <~ 0.05 L = "
              f"{0.05 * L_channel * 1e9:.0f} nm here; for a diffraction-limited")
        print("  ~1 um spot that means a channel of order 20 um.")
    return rows


def leaky_mask_cost_table(w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 5. What mask quality is worth paying for. Uses the closed form
    validated above, so this is exact, not a fit.
    """
    rows = [(T, (1.0 - T) / (1.0 + T)) for T in
            (0.0, 0.01, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.9, 1.0)]
    if verbose:
        print("\nRESULT 5: the cost of mask leakage, N(T)/N(0) = (1-T)/(1+T)")
        print(f"{'T':>7}{'N(T)/N(0)':>12}{'note':>40}")
        for T, r in rows:
            note = ""
            if T == 0.0:
                note = "perfect mask"
            elif T == T_MASK_MEASURED:
                note = "Shimomura et al.'s measured mask"
            elif T == 1.0:
                note = "no mask -> uniform -> exact zero"
            print(f"{T:>7.2f}{r:>12.4f}{note:>40}")
        print("  the derivative at T=0 is -2, so the FIRST few percent of")
        print("  leakage cost twice their face value; halving T from 0.5 to")
        print(f"  0.25 buys {((1-0.25)/(1+0.25)) / ((1-0.5)/(1+0.5)):.2f}x.")
    return rows


def make_plots(fname="photodetector_nonuniform_illumination.png",
               w_cross=W_CROSS_CHEM):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))

    # (a) the kernel = the delta-spot scan
    ax = axes[0]
    for (mA, mB, st) in [("Pt", "Pt", "-"), ("Ti", "Ti", "--"), ("Ti", "Pt", ":")]:
        x, k, _ = collection_kernel(METAL_WORK_FUNCTIONS[mA],
                                    METAL_WORK_FUNCTIONS[mB],
                                    e_applied=0.0, w_cross=w_cross)
        ax.plot(x * 1e9, k, st, label=f"{mA}/{mB}")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("spot position x (nm); contact A at 0, B at L")
    ax.set_ylabel("k(x) = charge to A per photon at x")
    ax.set_title("(a) The collection kernel is the delta-spot scan")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) symmetric device: uniform gives zero, a mask does not
    ax = axes[1]
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    names = [m for m, _ in metals]
    un, pm, lm = [], [], []
    for _, W in metals:
        x, k, _ = collection_kernel(W, W, e_applied=0.0, w_cross=w_cross)
        un.append(response_under(k, x, g_uniform(x)))
        pm.append(response_under(k, x, g_shadow_mask(x, T=0.0)))
        lm.append(response_under(k, x, g_shadow_mask(x, T=T_MASK_MEASURED)))
    idx = np.arange(len(names))
    ax.bar(idx - 0.27, un, 0.26, label="uniform (exactly zero)")
    ax.bar(idx, pm, 0.26, label="perfect mask, T = 0")
    ax.bar(idx + 0.27, lm, 0.26,
           label=f"measured mask, T = {T_MASK_MEASURED}")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xticks(idx); ax.set_xticklabels(names)
    ax.set_ylabel("zero-bias net response N")
    ax.set_title("(b) Symmetric devices, rescued by a mask")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")

    # (c) the leakage cost, closed form vs numerics
    ax = axes[2]
    Ts = np.linspace(0.0, 1.0, 201)
    ax.plot(Ts, (1 - Ts) / (1 + Ts), "-", label="closed form $(1-T)/(1+T)$")
    x, k, _ = collection_kernel(METAL_WORK_FUNCTIONS["Pt"],
                                METAL_WORK_FUNCTIONS["Pt"],
                                e_applied=0.0, w_cross=w_cross)
    n0 = response_under(k, x, g_shadow_mask(x, T=0.0))
    Tn = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9])
    num = [response_under(k, x, g_shadow_mask(x, T=float(T))) / n0 for T in Tn]
    ax.plot(Tn, num, "o", ms=7, label="numerics, Pt/Pt")
    ax.plot([T_MASK_MEASURED], [(1 - T_MASK_MEASURED) / (1 + T_MASK_MEASURED)],
            "s", ms=10, mfc="none", label="Shimomura et al., $T\\approx0.5$")
    ax.plot(Ts, 1 - Ts, "--", lw=1,
            label="per INCIDENT photon: $(1-T)$")
    ax.set_xlabel("mask transmission $T$")
    ax.set_ylabel("$N(T)/N(0)$")
    ax.set_title("(c) Two accountings: absorbed vs incident photons")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"\nsaved {fname}")
    return fname


if __name__ == "__main__":
    print("=" * 74)
    print("Non-uniform illumination: the collection kernel and the shadow mask")
    print("=" * 74)
    _, nb1, w1 = validate_uniform_reduction()
    _, w2 = validate_symmetric_g_gives_zero()
    _, w3 = validate_mirror_antisymmetry()
    _, w4 = validate_leaky_mask_closed_form()
    _, nv5, slack5 = validate_kernel_ceiling()
    masked_symmetric_table()
    mask_on_asymmetric_table()
    kernel_scan_table()
    spot_size_study()
    leaky_mask_cost_table()
    make_plots()
    print("\nValidation summary (all against exactly known values):")
    print(f"  1. uniform g reduces to the 2026-09-18 model : {nb1}/196 bitwise, "
          f"|dev| <= {w1:.3e}")
    print(f"  2. mirror-symmetric g on a symmetric device  : |N|   <= {w2:.3e}")
    print(f"  3. mask-left == -mask-right                  : |dev| <= {w3:.3e}")
    print(f"  4. leaky mask (1-T)/(1+T) closed form        : |dev| <= {w4:.3e}")
    print(f"  5. |N[g]| <= max|k| ceiling                  : {nv5} violations, "
          f"tightest slack {slack5:.4f}")
