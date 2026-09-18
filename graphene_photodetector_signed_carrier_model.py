"""
graphene_photodetector_signed_carrier_model.py

SIGNED, CARRIER-RESOLVED two-contact photocarrier collection model.

This closes the item that AUTOMATION_LOG.md has carried at the top of
"not yet covered" since 2026-09-17, and that
`graphene_photodetector_two_contact_model.py` itself lists as its own
Simplification 1:

    "Both doping fields are taken as sweeping carriers toward their own
     contact, via |W_metal - W_gr| ... A fully signed treatment would
     distinguish n-type (W_metal < W_gr, e.g. Ti, Cu) from p-type
     (W_metal > W_gr, e.g. Pt, Pd) contacts and track electrons and
     holes separately; a p-n pair (e.g. Ti/Pt) would then ADD rather
     than partially cancel for one carrier species."

That is exactly what this module does. It supersedes the magnitude
convention of the 2026-09-17 module for every ASYMMETRIC pair. For
IDENTICAL contacts the two agree on the sign structure, but not on the
magnitude, because this model also collects the second carrier species
(see Result 1 below).

TWO THINGS WERE WRONG WITH THE MAGNITUDE CONVENTION, NOT ONE
------------------------------------------------------------
Writing |W_metal - W_graphene| hid two separate physical facts.

(i) The SIGN. A contact that n-dopes graphene and one that p-dopes it
    build electrostatic fields that point the same way along the
    channel, not opposite ways. The magnitude convention forced every
    contact to behave like a p-type one, so every pair partially
    cancelled and every pair developed an interior stagnation point.

(ii) The CROSSOVER. The sign of the contact doping is NOT
     sign(W_metal - W_graphene). Graphene's own work function is
     4.5 eV, but the n/p crossover for a metal ON graphene sits at

         W_c ~ 5.4 eV

     because the short-range metal-graphene chemical interaction adds
     a ~0.9 eV interface dipole on top of simple vacuum-level
     alignment:

       Giovannetti, Khomyakov, Brocks, Karpan, van den Brink & Kelly,
       "Doping Graphene with Metal Contacts", Phys. Rev. Lett. 101,
       026803 (2008). arXiv:0802.2267.
       https://arxiv.org/abs/0802.2267
       Crossover at "a metal work function of ~5.4 eV"; metals below
       it (Al, Ag, Cu, Au, and by extension Ti, Cr, Ni, Pd here)
       n-dope graphene, metals above it (Pt) p-dope it.

     `graphene_contact_doping_model.py`'s own METAL_WORK_FUNCTIONS
     table already knows this -- its comments call Cu "n-type dopant"
     (despite 4.65 > 4.5) and Pt "clearly above the ~5.4 eV crossover
     -> p-type". The magnitude convention simply could not express it,
     because |W - 4.5| discards the crossover entirely.

     Consequence: under the physical crossover, SIX of the seven metals
     in the table n-dope graphene and only Pt p-dopes it. Pt is the
     only available partner for a true n/p pair.

To keep the comparison honest this module carries BOTH crossovers as an
explicit parameter, `w_cross`:

    W_CROSS_VACUUM = 4.5 eV  -- pure vacuum-level alignment; the
                               convention the 2026-09-17 magnitudes
                               implicitly assume. Used for the exact
                               reduction test, NOT recommended physics.
    W_CROSS_CHEM   = 5.4 eV  -- Giovannetti et al. Default.

Results are reported under both, and where they disagree that is said
out loud rather than averaged away.

THE MODEL
---------
Contact A at x=0, contact B at x=L. Signed offset for a contact metal:

    dW = W_metal - w_cross      (volts; >0 p-dopes, <0 n-dopes)

The Dirac point in the channel, measured from the Fermi level, is
rigidly shifted near each contact with the SAME saturating profile the
rest of the repo uses (Khomyakov et al. 2010, lambda_decay = 250 nm):

    E_D(x) = dW_A / (1 + x/lam)  +  dW_B / (1 + (L-x)/lam)

A rigid band shift is an electron potential energy U(x) = E_D(x), so
the electrostatic potential is phi = -E_D/e and the field (positive =
pointing along +x) is E = -dphi/dx = (1/e) dE_D/dx:

    E(x) = -E_applied
           - dW_A / lam / (1 + x/lam)**2
           + dW_B / lam / (1 + (L-x)/lam)**2

with E_applied the uniform bias field, signed so that a POSITIVE
E_applied sweeps HOLES toward contact A -- the same convention as the
two previous modules, so the numbers stay comparable.

Both carriers then drift in that one field, in opposite directions:

    v_hole(x)     = +mu * E(x)
    v_electron(x) = -mu * E(x)

Each is transported independently by the same stagnation-aware
collection logic as the 2026-09-17 module: a carrier is collected only
if the velocity keeps its sign over the whole path to the boundary it
started toward, with survival exp(-t/tau), tau = 1 ps.

The figure of merit is the net charge delivered to contact A per
absorbed photon, averaged over uniform illumination:

    N = (1/L) * integral [ (+1) s_h(x) p_h(x) + (-1) s_e(x) p_e(x) ] dx

    s = +1 collected at A, -1 collected at B, 0 stalled at a null.

The (-1) on the electron term is its charge, so an electron reaching B
(s_e = -1) contributes +1, the same sign as a hole reaching A: that is
ordinary photovoltaic charge separation. N therefore runs over
[-2, +2], not [-1, +1] -- BOTH carriers can now be collected, which the
magnitude convention structurally could not represent.

EXPERIMENTAL ANCHOR
-------------------
The n/p pairing this model rewards is not hypothetical. Mueller, Xia &
Avouris, Nature Photonics 4, 297 (2010) -- already cited in this repo
for tau = 1 ps -- built exactly this device: interdigitated fingers,
"One set of fingers was made of palladium/gold (20/25 nm in thickness),
and the other of titanium/gold (20/25 nm)", because "if both electrodes
... consist of the same metal, the built-in electric field profile in
the channel between two neighbouring fingers is symmetric, and the
total photocurrent is zero." Asymmetric metallization gave 6.1 mA/W at
1.55 um, a 15-fold improvement.

See notes/2026-09-18-signed-carrier-resolved-contact-fields.md.

STATED SIMPLIFICATIONS (unchanged from 2026-09-17 unless noted)
---------------------------------------------------------------
1. The doping-profile magnitude is taken as linear in dW with the same
   lambda_decay for every metal. Giovannetti et al. find the relation
   between W_metal and the graphene Fermi-level shift is only roughly
   linear and is not linear at all for the strongly bonded metals
   (Ti, Ni, Pd chemisorb). This model does not attempt that; it changes
   only the sign structure and the crossover, and says so.
2. Uniform illumination. (Still open: the Suzuki et al. shadow mask.)
3. Drift only, no diffusion, single mobility for both carriers. Taking
   mu_e = mu_h is standard for graphene's symmetric bands but means the
   electron and hole terms differ only through the field geometry.
4. No photogain, photo-thermoelectric or bolometric contribution.
"""

import numpy as np
import matplotlib.pyplot as plt

from graphene_contact_doping_model import METAL_WORK_FUNCTIONS, W_GRAPHENE, lambda_decay
from graphene_fet_model import mu
from graphene_photodetector_model import L_channel, V_bias
from graphene_photodetector_collection_model import TAU_CARRIER, E_BIAS
from graphene_photodetector_two_contact_model import net_response as magnitude_net_response

N_POINTS = 4001  # odd, so x = L/2 is on the grid (matters for symmetry tests)

W_CROSS_VACUUM = W_GRAPHENE   # 4.5 eV -- vacuum-level alignment only
W_CROSS_CHEM = 5.4            # eV -- Giovannetti et al. PRL 101, 026803 (2008)

Q_HOLE = +1.0
Q_ELECTRON = -1.0


def signed_offset(W_metal, w_cross=W_CROSS_CHEM):
    """dW = W_metal - w_cross, in volts. >0 p-dopes graphene, <0 n-dopes."""
    return W_metal - w_cross


def total_field(x, W_A, W_B, e_applied=E_BIAS, L=L_channel,
                w_cross=W_CROSS_CHEM, lam=lambda_decay, use_magnitude=False):
    """
    Electrostatic field E(x) in V/m, positive = pointing along +x.
    `e_applied` is signed so that positive sweeps HOLES toward contact A,
    matching the convention of the two earlier modules.

    use_magnitude=True forces |dW| for both contacts, reproducing the
    2026-09-17 convention; used only by the reduction test.
    """
    dW_A = signed_offset(W_A, w_cross)
    dW_B = signed_offset(W_B, w_cross)
    if use_magnitude:
        dW_A, dW_B = abs(dW_A), abs(dW_B)
    return (-e_applied
            - dW_A / lam / (1.0 + x / lam) ** 2
            + dW_B / lam / (1.0 + (L - x) / lam) ** 2)


def _transport(v, x, tau=TAU_CARRIER):
    """
    Stagnation-aware drift collection for ONE carrier species.

    `v` is the CARRIER-SIGNED FIELD, i.e. the drift velocity divided by the
    mobility: pass +E for holes and -E for electrons. Its sign is the
    direction of motion (v > 0 moves toward x = L = contact B) and mu*|v|
    is the drift speed, which is how `speed` is formed below.

    Returns (s, p) arrays: s = +1 collected at contact A, -1 collected at
    contact B, 0 stalled at an interior null; p = exp(-t/tau) survival.

    The two cumulative transit integrals are built independently in each
    direction. They must NOT be derived from one another: 1/|v| diverges
    at a stagnation point, so a cumulative sum that has crossed a null is
    infinite from there on and the difference of two such values is
    inf - inf. That exact bug (2026-09-17) silently marked every carrier
    beyond a null uncollectable and was caught only by the exactly-known
    zero of the symmetric test.
    """
    speed = mu * np.abs(v)
    with np.errstate(divide="ignore"):
        inv_v = np.where(speed > 0.0, 1.0 / np.maximum(speed, 1e-300), np.inf)

    dx = np.diff(x)
    seg = 0.5 * (inv_v[1:] + inv_v[:-1]) * dx
    cum_fwd = np.concatenate([[0.0], np.cumsum(seg)])          # 0   -> x_i
    cum_rev = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])  # x_i -> L

    min_from_left = np.minimum.accumulate(v)      # min of v on [0, x_i]
    max_from_left = np.maximum.accumulate(v)      # max of v on [0, x_i]
    min_from_right = np.minimum.accumulate(v[::-1])[::-1]   # min on [x_i, L]
    max_from_right = np.maximum.accumulate(v[::-1])[::-1]   # max on [x_i, L]

    s = np.zeros_like(x)
    t = np.full_like(x, np.inf)

    # v < 0 -> heads toward contact A at x = 0; needs v < 0 on all [0, x_i]
    ok_A = (v < 0.0) & (max_from_left < 0.0)
    s[ok_A] = +1.0
    t[ok_A] = cum_fwd[ok_A]

    # v > 0 -> heads toward contact B at x = L; needs v > 0 on all [x_i, L]
    ok_B = (v > 0.0) & (min_from_right > 0.0)
    s[ok_B] = -1.0
    t[ok_B] = cum_rev[ok_B]

    # min_from_left / max_from_right are computed for symmetry of the
    # expression but only the two used above are needed; reference them so
    # a future edit does not silently drop the wrong one.
    assert min_from_left.shape == max_from_right.shape == x.shape

    p = np.where(np.isfinite(t), np.exp(-t / tau), 0.0)
    return s, p


def net_response(W_A, W_B, e_applied=E_BIAS, L=L_channel, n_points=N_POINTS,
                 w_cross=W_CROSS_CHEM, use_magnitude=False,
                 carriers=("hole", "electron"), return_profiles=False):
    """
    Signed net charge delivered to contact A per absorbed photon, averaged
    over uniform illumination. Positive = net positive charge to contact A.

    carriers: restrict to ("hole",) to reproduce the single-species
    magnitude convention exactly (see validate_reduction_to_magnitude).
    """
    x = np.linspace(0.0, L, n_points)
    E = total_field(x, W_A, W_B, e_applied=e_applied, L=L,
                    w_cross=w_cross, use_magnitude=use_magnitude)

    terms = {}
    integrand = np.zeros_like(x)
    if "hole" in carriers:
        s_h, p_h = _transport(+E, x)   # v_hole = +mu E
        terms["hole"] = (s_h, p_h)
        integrand = integrand + Q_HOLE * s_h * p_h
    if "electron" in carriers:
        s_e, p_e = _transport(-E, x)   # v_electron = -mu E
        terms["electron"] = (s_e, p_e)
        integrand = integrand + Q_ELECTRON * s_e * p_e

    N = float(np.trapezoid(integrand, x) / L)
    if return_profiles:
        prof = {"x": x, "E": E, "N": N}
        for name, (s, p) in terms.items():
            prof[f"s_{name}"] = s
            prof[f"p_{name}"] = p
            prof[f"stalled_frac_{name}"] = float(np.count_nonzero(s == 0.0)) / len(x)
            prof[f"to_A_{name}"] = float(np.trapezoid(np.where(s > 0, p, 0.0), x) / L)
            prof[f"to_B_{name}"] = float(np.trapezoid(np.where(s < 0, p, 0.0), x) / L)
        return N, prof
    return N


# ---------------------------------------------------------------------
# Validations -- all three against EXACTLY known values, not plausible
# ranges. The 2026-09-17 session recorded why that distinction matters:
# a range check passed while a real inf-inf bug sat in the transit
# integral, and only an exact zero exposed it.
# ---------------------------------------------------------------------
def validate_reduction_to_magnitude(verbose=True):
    """
    EXACT check #1. Restricting this model to holes only, forcing both
    contacts p-type via |dW|, and putting the crossover back at graphene's
    own work function must reproduce
    `graphene_photodetector_two_contact_model.net_response()` to machine
    precision -- because under those three restrictions the two modules
    are algebraically the same expression:

        E(x) = -(E_bias + g_A - g_B) = -F_old,  v_hole = mu E,

    so 'hole moves toward A' (v<0) is exactly 'F_old > 0', and the transit
    integrals are over identical intervals of an identical 1/|v|.

    Any nonzero deviation here means one of the two modules is wrong.
    """
    rows, worst = [], 0.0
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    for mA, WA in metals:
        for mB, WB in metals:
            for eb in (0.0, E_BIAS):
                new = net_response(WA, WB, e_applied=eb,
                                   w_cross=W_CROSS_VACUUM,
                                   use_magnitude=True, carriers=("hole",))
                old = magnitude_net_response(WA, WB, E_bias=eb)
                dev = abs(new - old)
                worst = max(worst, dev)
                rows.append((mA, mB, eb, old, new, dev))
    if verbose:
        print("Validation 1 (EXACT): hole-only + |dW| + w_cross=4.5 eV must")
        print("reproduce the 2026-09-17 magnitude model on all 49 ordered")
        print(f"pairs at both biases ({len(rows)} comparisons).")
        worst_row = max(rows, key=lambda r: r[5])
        print(f"  worst case: {worst_row[0]}/{worst_row[1]} at E={worst_row[2]:.3e} V/m, "
              f"old={worst_row[3]:+.9f} new={worst_row[4]:+.9f}")
        print(f"  worst absolute deviation: {worst:.3e}")
    return rows, worst


def validate_symmetric_cancellation(w_cross=W_CROSS_CHEM, verbose=True):
    """
    EXACT check #2. Identical contacts at zero bias: E(x) is antisymmetric
    about x = L/2 for any dW, so holes and electrons each split evenly and
    N must be exactly zero. This is Weiss & Duan's 'net zero photocurrent'
    for symmetric metal-graphene-metal devices, and it must survive the
    move to a signed, two-carrier treatment.
    """
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        rows.append((metal, W, net_response(W, W, e_applied=0.0, w_cross=w_cross)))
    worst = max(abs(n) for _, _, n in rows)
    if verbose:
        print(f"\nValidation 2 (EXACT): symmetric pair at zero bias -> N == 0 "
              f"(w_cross = {w_cross} eV)")
        for metal, W, N0 in rows:
            print(f"  {metal:<4} (W={W:.2f} eV) both contacts: N = {N0:+.3e}")
        print(f"  worst |N|: {worst:.3e}")
    return rows, worst


def validate_charge_conjugation(w_cross=W_CROSS_CHEM, verbose=True):
    """
    EXACT check #3, and the one that is only meaningful once the model is
    signed -- the magnitude convention could not even state it.

    At zero bias, flipping the sign of BOTH contact offsets (every p-type
    contact made n-type and vice versa) flips E -> -E, which exchanges the
    hole and electron trajectories exactly. Since the two carriers enter N
    with opposite charges,

        N(-dW_A, -dW_B) = -N(dW_A, dW_B)      at zero bias, exactly.

    Implemented by reflecting the work functions through the crossover,
    W -> 2*w_cross - W, which negates dW by construction. A model that
    mixed carriers incorrectly, or that leaked a magnitude anywhere, would
    fail this even though it passes Validation 2.
    """
    rows, worst = [], 0.0
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    for mA, WA in metals:
        for mB, WB in metals:
            N = net_response(WA, WB, e_applied=0.0, w_cross=w_cross)
            Nflip = net_response(2 * w_cross - WA, 2 * w_cross - WB,
                                 e_applied=0.0, w_cross=w_cross)
            dev = abs(N + Nflip)
            worst = max(worst, dev)
            rows.append((mA, mB, N, Nflip, dev))
    if verbose:
        print(f"\nValidation 3 (EXACT): charge conjugation at zero bias, "
              f"N(-dW) == -N(dW) (w_cross = {w_cross} eV)")
        wr = max(rows, key=lambda r: r[4])
        print(f"  worst case: {wr[0]}/{wr[1]}  N={wr[2]:+.9f}  N_flip={wr[3]:+.9f}")
        print(f"  worst |N + N_flip|: {worst:.3e}")
    return rows, worst


# ---------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------
def doping_type_table(verbose=True):
    """Which metals n-dope and which p-dope graphene, under each crossover."""
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        rows.append((metal, W,
                     signed_offset(W, W_CROSS_VACUUM),
                     signed_offset(W, W_CROSS_CHEM)))
    if verbose:
        print("\nContact doping type under the two crossovers")
        print(f"{'metal':<6}{'W (eV)':>8}{'dW vs 4.5':>11}{'type':>7}"
              f"{'dW vs 5.4':>12}{'type':>7}")
        for metal, W, dv, dc in rows:
            print(f"{metal:<6}{W:>8.2f}{dv:>+11.2f}{('p' if dv > 0 else 'n'):>7}"
                  f"{dc:>+12.2f}{('p' if dc > 0 else 'n'):>7}")
        n_p = sum(1 for _, _, _, dc in rows if dc > 0)
        print(f"  under the physical crossover only {n_p} of {len(rows)} metals "
              f"p-dope graphene")
    return rows


def symmetric_pair_table(w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 1. Symmetric device at the working bias: magnitude convention
    (holes only) vs. this model (both carriers).
    """
    rows = []
    for metal, W in sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1]):
        mag = magnitude_net_response(W, W)
        sgn = net_response(W, W, w_cross=w_cross)
        holes = net_response(W, W, w_cross=w_cross, carriers=("hole",))
        rows.append((metal, W, mag, holes, sgn))
    if verbose:
        print(f"\nRESULT 1: symmetric pairs at V_bias = {V_bias} V "
              f"(w_cross = {w_cross} eV)")
        print(f"{'metal':<6}{'W (eV)':>8}{'magnitude':>12}{'signed,h only':>15}"
              f"{'signed, h+e':>13}")
        for metal, W, mag, holes, sgn in rows:
            print(f"{metal:<6}{W:>8.2f}{mag:>12.6f}{holes:>15.6f}{sgn:>13.6f}")
    return rows


def asymmetric_pair_table(w_cross=W_CROSS_CHEM, e_applied=0.0, verbose=True):
    """
    RESULT 2. Zero-bias net response for every ordered pair, signed model.
    """
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    out = {}
    for mA, WA in metals:
        for mB, WB in metals:
            out[(mA, mB)] = net_response(WA, WB, e_applied=e_applied,
                                         w_cross=w_cross)
    if verbose:
        names = [m for m, _ in metals]
        print(f"\nRESULT 2: zero-bias net response N(A,B), signed model "
              f"(w_cross = {w_cross} eV)")
        print("rows = contact A (x=0), cols = contact B (x=L)")
        print("      " + "".join(f"{m:>9}" for m in names))
        for mA in names:
            print(f"{mA:<6}" + "".join(f"{out[(mA, mB)]:>9.4f}" for mB in names))
        best = max(out.items(), key=lambda kv: abs(kv[1]))
        print(f"  largest |N|: {best[0][0]}/{best[0][1]} -> {best[1]:+.4f}")
        mag = {}
        for mA, WA in metals:
            for mB, WB in metals:
                mag[(mA, mB)] = magnitude_net_response(WA, WB, E_bias=0.0)
        bmag = max(mag.items(), key=lambda kv: abs(kv[1]))
        print(f"  magnitude convention said: {bmag[0][0]}/{bmag[0][1]} "
              f"-> {bmag[1]:+.4f}")
    return out


def stagnation_audit(w_cross=W_CROSS_CHEM, verbose=True):
    """
    RESULT 3. The mechanism behind Result 2: an n/p pair builds two fields
    that point the SAME way, so the interior stagnation point disappears
    and both carrier species are collected, at opposite ends.
    """
    cases = [("Pt", "Pt"), ("Pd", "Pd"), ("Ti", "Pd"), ("Ti", "Pt"), ("Cr", "Pt")]
    rows = []
    for mA, mB in cases:
        N, prof = net_response(METAL_WORK_FUNCTIONS[mA], METAL_WORK_FUNCTIONS[mB],
                               e_applied=0.0, w_cross=w_cross,
                               return_profiles=True)
        E = prof["E"]
        crossings = int(np.count_nonzero(np.diff(np.sign(E)) != 0))
        rows.append((mA, mB, N, crossings,
                     prof["stalled_frac_hole"], prof["stalled_frac_electron"],
                     prof["to_A_hole"], prof["to_B_hole"],
                     prof["to_A_electron"], prof["to_B_electron"]))
    if verbose:
        print(f"\nRESULT 3: stagnation audit at zero bias (w_cross = {w_cross} eV)")
        print(f"{'pair':<10}{'N':>9}{'E nulls':>9}{'stall h':>9}{'stall e':>9}"
              f"{'h->A':>8}{'h->B':>8}{'e->A':>8}{'e->B':>8}")
        for mA, mB, N, cr, sh, se, hA, hB, eA, eB in rows:
            print(f"{mA+'/'+mB:<10}{N:>+9.4f}{cr:>9d}{sh:>9.2f}{se:>9.2f}"
                  f"{hA:>8.3f}{hB:>8.3f}{eA:>8.3f}{eB:>8.3f}")
    return rows


def make_plots(fname="photodetector_signed_carrier_response.png",
               w_cross=W_CROSS_CHEM):
    metals = sorted(METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])
    names = [m for m, _ in metals]
    x = np.linspace(0.0, L_channel, N_POINTS)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))

    # (a) the field itself: same-type pair keeps a null, n/p pair does not
    ax = axes[0]
    for (mA, mB, st) in [("Pt", "Pt", "-"), ("Pd", "Pd", "--"), ("Ti", "Pt", ":")]:
        E = total_field(x, METAL_WORK_FUNCTIONS[mA], METAL_WORK_FUNCTIONS[mB],
                        e_applied=0.0, w_cross=w_cross)
        dA = signed_offset(METAL_WORK_FUNCTIONS[mA], w_cross)
        dB = signed_offset(METAL_WORK_FUNCTIONS[mB], w_cross)
        lab = (f"{mA}({'p' if dA > 0 else 'n'}) / {mB}({'p' if dB > 0 else 'n'})")
        ax.plot(x * 1e9, E / 1e6, st, label=lab)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("position x (nm); contact A at 0, contact B at L")
    ax.set_ylabel("field E(x)  (MV/m)")
    ax.set_title("(a) An n/p pair has no interior null")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) zero-bias response against a Pt counter-electrode: the two
    #     conventions, side by side
    ax = axes[1]
    W_pt = METAL_WORK_FUNCTIONS["Pt"]
    magv = [magnitude_net_response(W, W_pt, E_bias=0.0) for _, W in metals]
    sgnv = [net_response(W, W_pt, e_applied=0.0, w_cross=w_cross) for _, W in metals]
    idx = np.arange(len(names))
    ax.bar(idx - 0.2, magv, 0.4, label="magnitude convention (2026-09-17)")
    ax.bar(idx + 0.2, sgnv, 0.4, label="signed, carrier-resolved (this model)")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xticks(idx); ax.set_xticklabels(names)
    ax.set_ylabel("zero-bias net response N")
    ax.set_title("(b) Contact A = metal shown, contact B = Pt")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    # (c) N vs the signed offset product: cancellation vs. reinforcement
    ax = axes[2]
    for w_c, mk, lab in [(W_CROSS_CHEM, "o", f"w_cross = {W_CROSS_CHEM} eV (physical)"),
                         (W_CROSS_VACUUM, "^", f"w_cross = {W_CROSS_VACUUM} eV (vacuum only)")]:
        xs, ys = [], []
        for mA, WA in metals:
            for mB, WB in metals:
                xs.append(signed_offset(WA, w_c) * signed_offset(WB, w_c))
                ys.append(net_response(WA, WB, e_applied=0.0, w_cross=w_c))
        ax.scatter(xs, ys, s=20, marker=mk, alpha=0.75, label=lab)
    ax.axhline(0.0, color="k", lw=0.8); ax.axvline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"$\Delta W_A \cdot \Delta W_B$  (eV$^2$);  $<0$ = n/p pair")
    ax.set_ylabel("zero-bias net response N")
    ax.set_title("(c) Only opposite-type pairs reinforce")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"\nsaved {fname}")
    return fname


if __name__ == "__main__":
    print("=" * 74)
    print("Signed, carrier-resolved two-contact photocarrier collection model")
    print("=" * 74)
    _, w1 = validate_reduction_to_magnitude()
    _, w2 = validate_symmetric_cancellation()
    _, w3 = validate_charge_conjugation()
    doping_type_table()
    symmetric_pair_table()
    asymmetric_pair_table()
    stagnation_audit()
    make_plots()
    print("\nValidation summary (all three against exactly known values):")
    print(f"  1. reduction to the magnitude model : |dev| <= {w1:.3e}")
    print(f"  2. symmetric zero-bias cancellation : |N|   <= {w2:.3e}")
    print(f"  3. charge conjugation at zero bias  : |dev| <= {w3:.3e}")
