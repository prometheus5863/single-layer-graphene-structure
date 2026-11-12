"""
graphene_self_heating_model.py

Self-consistent electro-thermal model of the back-gated single-layer graphene
FET of graphene_fet_model.py: the channel dissipates power, the power raises
the channel temperature through a series thermal network, and the temperature
degrades both the mobility and the saturation velocity, which changes the
current, which changes the power.  Solved to a fixed point.

This discharges the standing-list item "A SELF-HEATING TERM", created
2026-10-06 and carried unchanged for 34 sessions.  See
notes/2026-11-10-self-heating-the-thermal-network-and-two-exact-benchmarks.md
for the literature reading, the two exact benchmarks below, and the two
extraction ambiguities (one settled by arithmetic, one carried as a spread).

PROVENANCE OF EVERY THERMAL AND HIGH-FIELD NUMBER
-------------------------------------------------
All of the following are from Dorgan, Bae and Pop, "Mobility and saturation
velocity in graphene on SiO2", Appl. Phys. Lett. 97, 082112 (2010), read as
arXiv:1005.2711 (preprint + supplement):

  Eq. (5)  dT = P * (R_B + R_ox + R_Si),  R_B = 1/(h A),
           R_ox = t_ox/(kappa_ox A),  R_Si = 1/(2 kappa_Si sqrt(A)),  A = L W
  supp. 6  kappa_ox = 0.52*ln(T_ox) - 1.687  [W/m.K],  T_ox = (T0 + T)/2
           kappa_Si = 2.4e4 / T0             [W/m.K]
  Eq. (4)  mu(n,T) -- see MU0/N_REF/ALPHA/BETA below
  Eq. (6)  v(F) = mu F / [1 + (mu F/v_sat)^gamma]^(1/gamma),  gamma ~ 2
  Eq. (7)  v_sat(n,T) -- see V_SAT_PHONON_MEV below

h ~ 1e8 W m^-2 K^-1 is quoted by Dorgan FROM Chen et al., APL 95, 161910
(2009), which was NOT fetched here.  It is therefore a SECONDARY CITATION and
is labelled as one, per the standing 10-01 rule.

TWO EXACT EXTERNAL BENCHMARKS (verify_benchmarks(), run by __main__)
-------------------------------------------------------------------
A. On Dorgan's OWN geometry (W=7um, L=4um, t_ox=300nm) the network must
   reproduce the paper's stated term shares 4 / 84 / 12 %, total ~1e4 K/W,
   and -- added 2026-11-12 -- the ~2.8e-7 m^2.K/W per unit of device area
   printed in the same sentence as the total.  The areal statement carries
   two figures where the total carries one, so it is the tighter half of the
   SAME constraint and not a further one; A2c asserts that rather than
   leaving the reader to assume independence (which the 11-11 standing list
   did assume, wrongly -- see graphene_areal_thermal_resistance.py).  That is agreement on three independent ratios, so no single
   fudged constant can buy it.
B. The Eq. (7) branch-cut threshold n_min = (w_OP/v_F)^2/(2 pi) must
   reproduce the paper's footnote-18 values 9.4e11 cm^-2 at 160 meV and
   1.1e11 cm^-2 at 55 meV.  This pins v_F, hbar and the factor 2 pi
   together: a factor error in the v_sat prefactor survives A and dies here.

These are external and exact.  On 2026-09-17 a plausible-range check in this
repository passed a real bug; a published three-ratio identity does not have
that failure mode, which is the whole reason the item waited for a source
that carried its own thermal budget.

WHAT THIS MODEL IS NOT
----------------------
  - The channel is LUMPED and the lateral field UNIFORM (F = V_ch/L).  The
    repo's graphene_fet_model.transfer_characteristic() instead averages the
    local resistance along a 50-segment channel.  The two agree in the
    low-field limit and overlap_with_isothermal_model() MEASURES the residual
    rather than asserting it.
  - Eq. (5) is stated valid for W, L >> t_ox.  This device has L/t_ox = 2.22,
    so lateral spreading in the oxide is neglected and R_ox is therefore
    OVER-estimated.  The error is DIRECTIONAL: this network is an upper bound
    on R_th, hence on dT, hence a worst-case bound on the degradation.  A
    bound of known sign, not a centred estimate of unknown sign.
  - No diffusion current (the RF thread's binding question since 10-06), no
    spreading-corrected R_ox, no lateral heat flow into the contacts (Dorgan
    supp. 5 measures the probe thermal resistance at ~1e6 K/W, two orders
    above the network, so this omission is small for a 4um device and is
    NOT re-justified for a 200nm one -- see the standing list).
"""

import numpy as np

import graphene_fet_model as gfet

# ---------------------------------------------------------------------------
# Constants (SI unless a name says otherwise)
# ---------------------------------------------------------------------------
e = 1.602176634e-19
k_B = 1.380649e-23
hbar = 1.054571817e-34

# --- Dorgan Eq. (4) mobility fit -------------------------------------------
MU0 = 4650.0 * 1e-4          # cm^2/V.s -> m^2/V.s
N_REF = 1.1e13 * 1e4         # cm^-2 -> m^-2
T_REF = 300.0                # K
ALPHA = 2.2
BETA = 3.0

# --- Dorgan Eq. (7) saturation velocity -----------------------------------
V_SAT_PHONON_MEV = 82.0      # the paper's OWN suggestion for device
                             # simulation; 55 meV (SiO2) and 160 meV
                             # (graphene OP) bracket the data and are swept
V_F = 1.0e6                  # m/s -- the paper's v_F ~ 1e8 cm/s
GAMMA_VF = 2.0               # Eq. (6) exponent

# --- Thermal network -------------------------------------------------------
H_BOUNDARY = 1.0e8           # W/m^2.K, graphene/SiO2 -- SECONDARY CITATION
T_AMBIENT = 300.0            # K

# --- Dorgan's own fit-exclusion ceiling -----------------------------------
# Data with dT > ~200 K are excluded from the paper's fits because heating
# past that point causes irreproducible behaviour (a V_0 shift).  Carried as
# a REPORTED regime boundary, never as a clamp on the solve.
DT_REPRODUCIBILITY_CEILING = 200.0

# Fixed-point solver controls
_T_TOL = 1e-9    # K; the solved dT is O(1..1e2) K, so this is far
                 # tighter than any physics in the model
_T_MAXIT = 400
_VCH_TOL = 1e-14


# ---------------------------------------------------------------------------
# Thermal network
# ---------------------------------------------------------------------------
def kappa_ox(T_ox):
    """SiO2 thermal conductivity, Dorgan supp. Sec. 6.

    The extracted text reads ln(T_ox^0.52) - 1.687, which is IDENTICALLY
    0.52*ln(T_ox) - 1.687; see the note's Section 4.  Returns W/m.K.
    """
    return 0.52 * np.log(T_ox) - 1.687


def kappa_Si(T0):
    """Doped-Si thermal conductivity, Dorgan supp. Sec. 6.  W/m.K."""
    return 2.4e4 / T0


def thermal_resistances(W, L, t_ox, T=T_AMBIENT, T0=T_AMBIENT,
                        h=H_BOUNDARY):
    """The three terms of Dorgan Eq. (5), in K/W, as a dict.

    kappa_ox is evaluated at the paper's mean oxide temperature (T0+T)/2.
    """
    A = L * W
    k_ox = kappa_ox(0.5 * (T0 + T))
    R_B = 1.0 / (h * A)
    R_ox = t_ox / (k_ox * A)
    R_Si = 1.0 / (2.0 * kappa_Si(T0) * np.sqrt(A))
    return {"R_B": R_B, "R_ox": R_ox, "R_Si": R_Si,
            "R_th": R_B + R_ox + R_Si, "A": A, "kappa_ox": k_ox}


# ---------------------------------------------------------------------------
# Temperature- and density-dependent transport parameters
# ---------------------------------------------------------------------------
def mobility(n, T, reading="product"):
    """Dorgan Eq. (4).  n in m^-2, T in K, returns m^2/V.s.

    TWO READINGS of the extracted temperature factor, deliberately both
    implemented (note Section 4).  They are bitwise identical at T = T_REF,
    so no room-temperature number can discriminate them.

      "product" : mu0 / (1 + (n/n_ref)^a) * 1/(1 + (T/T_ref - 1)^b)
      "sum"     : mu0 / (1 + (n/n_ref)^a + (T/T_ref - 1)^b)
    """
    n_term = (np.abs(n) / N_REF) ** ALPHA
    t_term = (T / T_REF - 1.0) ** BETA
    if reading == "product":
        return MU0 / (1.0 + n_term) / (1.0 + t_term)
    if reading == "sum":
        return MU0 / (1.0 + n_term + t_term)
    raise ValueError("reading must be 'product' or 'sum', got %r" % (reading,))


def n_min_branch_cut(phonon_meV):
    """Dorgan footnote 18: Eq. (7)'s square root is real only above this
    density.  Returns m^-2.  This is Benchmark B."""
    w_OP = phonon_meV * 1e-3 * e / hbar
    return (w_OP / V_F) ** 2 / (2.0 * np.pi)


def saturation_velocity(n, T, phonon_meV=V_SAT_PHONON_MEV):
    """Dorgan Eq. (7).  n in m^-2, T in K, returns m/s.

    Below the Eq. (7) branch cut the square root is imaginary and the
    expression has no meaning; this REFUSES rather than clamping, because a
    clamp would silently manufacture a velocity in a regime the paper
    excludes.
    """
    n = np.abs(n)
    n_cut = n_min_branch_cut(phonon_meV)
    if np.any(n <= n_cut):
        raise ValueError(
            "saturation_velocity called below the Eq. (7) branch cut: "
            "n_min = %.4e m^-2 at %.1f meV, got min(n) = %.4e m^-2"
            % (n_cut, phonon_meV, float(np.min(n))))
    w_OP = phonon_meV * 1e-3 * e / hbar
    N_OP = 1.0 / (np.expm1(hbar * w_OP / (k_B * T)))
    root = np.sqrt(1.0 - w_OP ** 2 / (4.0 * np.pi * n * V_F ** 2))
    return (2.0 / np.pi) * (w_OP / np.sqrt(np.pi * n)) * root / (N_OP + 1.0)


def drift_velocity(F, mu, v_sat, gamma=GAMMA_VF):
    """Dorgan Eq. (6).  F in V/m."""
    if mu * F == 0.0:
        return 0.0 * F
    return mu * F / (1.0 + (mu * F / v_sat) ** gamma) ** (1.0 / gamma)


# ---------------------------------------------------------------------------
# The self-consistent electro-thermal solve
# ---------------------------------------------------------------------------
def _mean_density(V_g, V_ch, n_segments=50):
    """Segment-averaged |n| along the channel, matching the averaging spirit
    of graphene_fet_model.transfer_characteristic()."""
    if V_ch == 0.0:
        return abs(float(gfet.carrier_density(V_g, 0.0)))
    Vs = np.linspace(0.0, V_ch, n_segments)
    # Vectorised over the segment array: gfet.carrier_density is pure numpy,
    # and the thermal fixed point re-enters this function ~1e4 times per
    # sweep, so the scalar loop this replaces dominated the runtime.
    return float(np.mean(np.abs(gfet.carrier_density(V_g, Vs))))


def _current_at_T(V_g, V_ds, T, reading, phonon_meV, W, L, Rc):
    """Solve the contact-resistance voltage division at FIXED temperature.

    Returns (Id, V_ch, n_bar, mu, v_sat).  Bisection on V_ch in [0, V_ds]
    of g(V_ch) = V_ch + Id(V_ch)*Rc - V_ds, which is monotone increasing.
    """
    if V_ds == 0.0:
        n_bar = _mean_density(V_g, 0.0)
        return (0.0, 0.0, n_bar, mobility(n_bar, T, reading),
                saturation_velocity(n_bar, T, phonon_meV))

    def Id_of(V_ch):
        n_bar = _mean_density(V_g, V_ch)
        mu = mobility(n_bar, T, reading)
        vs = saturation_velocity(n_bar, T, phonon_meV)
        v = drift_velocity(V_ch / L, mu, vs)
        return W * n_bar * e * v, n_bar, mu, vs

    lo, hi = 0.0, V_ds
    for _ in range(60):   # bisection on [0, V_ds] to double precision
        mid = 0.5 * (lo + hi)
        Id, n_bar, mu, vs = Id_of(mid)
        if mid + Id * Rc - V_ds > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < _VCH_TOL * max(1.0, V_ds):
            break
    V_ch = 0.5 * (lo + hi)
    Id, n_bar, mu, vs = Id_of(V_ch)
    return Id, V_ch, n_bar, mu, vs


def solve_bias_point(V_g, V_ds, self_heating=True, reading="product",
                     phonon_meV=V_SAT_PHONON_MEV, W=None, L=None, t_ox=None,
                     Rc=None, T0=T_AMBIENT, h=H_BOUNDARY,
                     R_th_override=None):
    """Self-consistent electro-thermal bias point.

    With self_heating=False (or R_th_override=0.0) the temperature is pinned
    at T0 and the result is the isothermal model; isothermal_recovery() uses
    that for an EXACT recovery check.

    P is taken as Id*V_ds -- the TOTAL dissipation including the contacts,
    not just the channel -- because the contacts sit on the same oxide.  Both
    are reported; using the total is the conservative choice and is stated
    rather than silently assumed.
    """
    W = gfet.W if W is None else W
    L = gfet.L if L is None else L
    t_ox = gfet.t_ox if t_ox is None else t_ox
    Rc = gfet.Rc_total if Rc is None else Rc

    T = T0
    iters = 0
    converged = True
    for iters in range(1, _T_MAXIT + 1):
        Id, V_ch, n_bar, mu, vs = _current_at_T(V_g, V_ds, T, reading,
                                                phonon_meV, W, L, Rc)
        P_total = Id * V_ds
        P_channel = Id * V_ch
        if not self_heating:
            R_th = 0.0
        elif R_th_override is not None:
            R_th = R_th_override
        else:
            R_th = thermal_resistances(W, L, t_ox, T=T, T0=T0, h=h)["R_th"]
        T_new = T0 + P_total * R_th
        if abs(T_new - T) < _T_TOL:
            T = T_new
            break
        # light damping: the kappa_ox(T) feedback is weak but the P(T)
        # feedback is not, and undamped iteration rings at high V_ds
        T = 0.5 * (T + T_new)
    else:
        converged = False

    Id, V_ch, n_bar, mu, vs = _current_at_T(V_g, V_ds, T, reading,
                                            phonon_meV, W, L, Rc)
    P_total = Id * V_ds
    terms = thermal_resistances(W, L, t_ox, T=T, T0=T0, h=h)
    return {
        "V_g": V_g, "V_ds": V_ds, "Id": Id, "V_ch": V_ch, "n_bar": n_bar,
        "mu": mu, "v_sat": vs, "F": (V_ch / L) if V_ds else 0.0,
        "T": T, "dT": T - T0, "P_total": P_total, "P_channel": Id * V_ch,
        "R_th": 0.0 if not self_heating else (
            R_th_override if R_th_override is not None else terms["R_th"]),
        "R_B": terms["R_B"], "R_ox": terms["R_ox"], "R_Si": terms["R_Si"],
        "iters": iters, "converged": converged,
        "past_reproducibility_ceiling": (T - T0) > DT_REPRODUCIBILITY_CEILING,
    }


# ---------------------------------------------------------------------------
# CHECKS.  Every one is an exact identity, an exact zero, a bitwise recovery
# or a published two-significant-figure value -- never a plausible range.
# ---------------------------------------------------------------------------
_CHECKS = []


def _check(name, ok, detail):
    _CHECKS.append((name, bool(ok), detail))
    return bool(ok)


def verify_benchmarks():
    """BENCHMARK A and BENCHMARK B: the two exact external checks."""
    # --- A: Dorgan's own device, his own stated term shares --------------
    t = thermal_resistances(W=7e-6, L=4e-6, t_ox=300e-9,
                            T=T_AMBIENT, T0=T_AMBIENT)
    shares = {k: 100.0 * t[k] / t["R_th"] for k in ("R_B", "R_ox", "R_Si")}
    # A1: TWICE CORRECTED ON 2026-11-11, AND THE SECOND CORRECTION RETRACTS
    # PART OF THE FIRST.  Neither number is deleted.
    #
    # STEP 1, from graphene_tolerance_census_audit.py.  A1 read
    # "kappa_ox(300K) == accepted thermal-SiO2 value 1.2794 W/m.K" and tested
    # it with an absolute tolerance of 5e-4 where five quoted figures license
    # 5e-5.  The transcribed formula gives
    #
    #     kappa_ox(300) = 0.52*ln(300) - 1.687 = 1.2789668868 W/m.K
    #
    # which rounds to 1.2790 at five figures, not 1.2794.  Residual 4.3311e-4,
    # 0.0339 %, between the licensed half-width and the coded one -- so A1 was
    # green while the equality in its own label was false.  That commit left
    # the tolerance standing and called which of the two numbers was wrong an
    # open question, listing an intercept of 1.686567 or a coefficient of
    # 0.52007593 as candidate roundings of "either published coefficient".
    #
    # STEP 2, from reading the paper.  THE QUESTION WAS WRONGLY POSED AND
    # 1.2794 IS NOT A PUBLISHED NUMBER AT ALL.  Dorgan supp. Sec. 6 prints
    #
    #     kappa_ox = ln(T_ox^0.52) - 1.687        (identically 0.52*ln T_ox
    #                                              - 1.687; the module's
    #                                              existing docstring is
    #                                              right about this)
    #
    # with 0.52 and 1.687 as the only digits given, and the only numerical
    # conductivity it prints is in supp. Sec. 5:
    #
    #     "kappa_ox ~ 1.3 Wm-1K-1 is the oxide thermal conductivity (at 300 K)"
    #
    # 1.2794 appears nowhere in the paper, to any precision.  So the candidate
    # roundings offered in step 1 were roundings of coefficients the paper
    # never quotes beyond the two digits this module already uses, and
    # "accepted thermal-SiO2 value" named an authority that does not exist.
    # The likeliest provenance of 1.2794 is an earlier evaluation of this
    # module's OWN formula written down to five figures with a slip in the
    # fourth decimal -- i.e. a self-comparison wearing a citation.
    #
    # This is the worse of the two defects the census could have found here.
    # A vacuous tolerance passes a wrong value; a check mislabelled as
    # external validation reports agreement with the literature that was never
    # tested.  This repository's only check against the "accepted" SiO2
    # conductivity had no external referent.
    #
    # WHAT THE CHECKS NOW DO.  A1 is the comparison the paper actually
    # licenses, against the 1.3 it prints, by the significant-figure test --
    # and it PASSES, so the transcription was never in doubt; only its
    # warrant was.  A1b keeps the five-figure arithmetic that exposed the
    # defect.  A1c keeps 1.2794 itself, with its tolerance untouched and its
    # label no longer claiming publication, so the number survives for
    # whoever wants to trace where it came from.  A1c is deliberately NOT
    # worded as a literature comparison any more, which is why the tolerance
    # census no longer lists it: that is the repair, not a blinding.
    _check("A1 kappa_ox(300K) rounds to the 1.3 W/m.K Dorgan supp. Sec. 5 "
           "prints for the oxide at 300 K, at the 2 s.f. that value is "
           "printed to",
           _rounds_to(kappa_ox(300.0), 1.3, 2),
           "kappa_ox(300) = %.10f -> %r at 2 s.f. (paper: 1.3)"
           % (kappa_ox(300.0), _round_sig(kappa_ox(300.0), 2)))
    _check("A1b kappa_ox(300K) = 1.2789668868 does NOT round to 1.2794 at "
           "five figures (it gives 1.2790), so A1's pre-2026-11-11 wording "
           "'== accepted value 1.2794' was false while A1 was green",
           _round_sig(kappa_ox(300.0), 5) == 1.2790
           and not _rounds_to(kappa_ox(300.0), 1.2794, 5),
           "%.10f -> %r at 5 s.f." % (kappa_ox(300.0),
                                      _round_sig(kappa_ox(300.0), 5)))
    _check("A1c the 1.2794 of the old label is retained and still agrees to "
           "5e-4 W/m.K, but it is NOT a value this paper prints at any "
           "precision and carries no citation -- provenance unresolved",
           abs(kappa_ox(300.0) - 1.2794) < 5e-4,
           "residual %.4e W/m.K; the paper prints only ~1.3 (supp. Sec. 5) "
           "and the fit ln(T^0.52) - 1.687 (supp. Sec. 6)"
           % abs(kappa_ox(300.0) - 1.2794))
    # A2 REPAIRED ON 2026-11-12, AND THE REPAIR GOES THE OPPOSITE WAY TO THE
    # ONE THE 11-11 LIST PROPOSED.  Neither number is deleted.
    #
    # The 11-11 standing list ranked third "REPAIR A2's OVER-TIGHT TOLERANCE,
    # OR DECIDE NOT TO", on the reasoning that +-5 % against a value the paper
    # prints to one figure ("~1e4 K/W") is a tighter claim than the source
    # supports, with widening to the licensed +-50 % the alternative and
    # "nearly vacuous" the objection to it.  THE DILEMMA WAS FALSE, because
    # the paper states the same fact twice in the SAME SENTENCE:
    #
    #   "At 300 K for our geometry R_th ~ 10^4 K/W, or ~2.8x10^-7 m^2 K/W
    #    per unit of device area"
    #
    # and the second half carries TWO figures.  Since A = L*W = 2.8e-11 m^2 is
    # given by the stated geometry and not fitted, 2.8e-7/A is 1.0e4 to within
    # one ulp, so the two halves are one constraint and the tighter half
    # licenses +-0.05/2.8 = +-1.7857 % on it.  See
    # graphene_areal_thermal_resistance.py, checks I1/I2/F4, which assert the
    # two residuals agree to 14 significant figures.
    #
    # So +-5 % was never over-tight.  It was 2.8x LOOSER than the source
    # supports, and A2's defect was its CITATION rather than its width: it
    # quoted the one-figure half of a two-figure sentence and then used a
    # width that neither half licenses.  A2 is therefore TIGHTENED to the
    # licensed +-1.7857 % and NOT widened to +-50 %.  It passes with 0.94 %
    # of the licensed width to spare.  The old +-0.05 is preserved in A2a,
    # which records what was being asserted before today.
    _AREAL_PAPER, _AREAL_SIGFIGS = 2.8e-7, 2
    _AREAL_HALFWIDTH = 0.05 / 2.8          # what 2 printed figures license
    _check("A2 Dorgan R_th total agrees with the paper's ~1e4 K/W to the "
           "+-1.7857 %% that the SAME SENTENCE's two-figure areal statement "
           "licenses (tightened from +-5 %% on 2026-11-12; the 11-11 list "
           "had this backwards and proposed widening to +-50 %%)",
           abs(t["R_th"] - 1.0e4) / 1.0e4 < _AREAL_HALFWIDTH,
           "R_th = %.4f K/W, residual %.4f %% of 1e4, licensed +-%.4f %%"
           % (t["R_th"], 100.0 * (t["R_th"] - 1.0e4) / 1.0e4,
              100.0 * _AREAL_HALFWIDTH))
    _check("A2a the pre-2026-11-12 form of A2 -- the same comparison at the "
           "+-5 %% it coded against the one-figure ~1e4 -- also passed, so "
           "today's repair changed the WARRANT and not the verdict",
           abs(t["R_th"] - 1.0e4) / 1.0e4 < 0.05,
           "residual %.4f %% against the old +-5 %%; the one figure in ~1e4 "
           "would itself have licensed +-50 %%"
           % (100.0 * (t["R_th"] - 1.0e4) / 1.0e4))
    _check("A2b R_th * A rounds to the paper's ~2.8e-7 m^2.K/W at the 2 s.f. "
           "it is printed to -- the per-unit-area half of the sentence, "
           "which this module had never checked before 2026-11-12.  NOTE: "
           "the quantity is R_th TIMES area, not divided by it; the quotient "
           "is 21 decades away",
           _rounds_to(t["R_th"] * t["A"], _AREAL_PAPER, _AREAL_SIGFIGS),
           "R_th*A = %.6e -> %.1e at 2 s.f. (paper: 2.8e-7); R_th/A = %.3e"
           % (t["R_th"] * t["A"],
              float(_round_sig(t["R_th"] * t["A"], _AREAL_SIGFIGS)),
              t["R_th"] / t["A"]))
    _check("A2c A2 and A2b are ONE constraint, not two: their relative "
           "residuals agree to within one ulp of 1.0, so neither can fail "
           "while the other passes at equal width.  A2b is kept because the "
           "paper prints it to more figures, NOT because it is independent",
           abs((t["R_th"] * t["A"] - _AREAL_PAPER) / _AREAL_PAPER
               - (t["R_th"] - 1.0e4) / 1.0e4) <= 4.0 * np.spacing(1.0),
           "areal residual %.14e vs total residual %.14e, difference %.3e"
           % ((t["R_th"] * t["A"] - _AREAL_PAPER) / _AREAL_PAPER,
              (t["R_th"] - 1.0e4) / 1.0e4,
              (t["R_th"] * t["A"] - _AREAL_PAPER) / _AREAL_PAPER
              - (t["R_th"] - 1.0e4) / 1.0e4))
    _check("A3 R_B share rounds to the paper's 4 %%",
           round(shares["R_B"]) == 4,
           "R_B = %.1f K/W, %.2f %%" % (t["R_B"], shares["R_B"]))
    _check("A4 R_ox share rounds to the paper's 84 %%",
           round(shares["R_ox"]) == 84,
           "R_ox = %.1f K/W, %.2f %%" % (t["R_ox"], shares["R_ox"]))
    _check("A5 R_Si share rounds to the paper's 12 %%",
           round(shares["R_Si"]) == 12,
           "R_Si = %.1f K/W, %.2f %%" % (t["R_Si"], shares["R_Si"]))
    _check("A6 the three shares sum to exactly 100 %% (no fourth term hiding)",
           abs(sum(shares.values()) - 100.0) < 1e-9,
           "sum = %.12f %%" % sum(shares.values()))
    # --- B: Eq. (7) branch cut at the paper's two phonon energies -------
    #
    # TOLERANCE DEFECT, FOUND BY THIS CHECK ON ITSELF (2026-11-10).
    # The first form of B2 compared the computed 1.1113e11 cm^-2 against the
    # paper's quoted 1.1e11 with a 1 % RELATIVE tolerance, and failed at
    # 1.03 %.  The model was right and the CHECK was wrong: the paper quotes
    # two significant figures, and 1.1113e11 rounds to 1.1e11 exactly as
    # required.  A relative tolerance is the wrong instrument for comparing
    # against a rounded published number -- it tests a precision the source
    # never claimed.  Replaced with a round-to-the-quoted-figures test, which
    # is what "matches the paper" actually means.
    #
    # This is the first live instance of the standing item AUDIT EVERY
    # TOLERANCE (created 2026-10-09, never exercised until now), and it
    # arrived from the direction that item did not anticipate: not a
    # tolerance too loose to catch a defect, but one too tight against a
    # number quoted to fewer figures than the comparison assumed.
    for meV, paper_cm2, sigfigs, label in ((160.0, 9.4e11, 2, "B1"),
                                           (55.0, 1.1e11, 2, "B2")):
        got = n_min_branch_cut(meV) * 1e-4   # m^-2 -> cm^-2
        _check("%s n_min(%.0f meV) rounds to the paper's footnote-18 "
               "%.1e cm^-2 at %d s.f."
               % (label, meV, paper_cm2, sigfigs),
               _rounds_to(got, paper_cm2, sigfigs),
               "n_min = %.4e cm^-2 -> %.1e at %d s.f. (paper: %.1e); "
               "relative difference %.3f %%"
               % (got, _round_sig(got, sigfigs), sigfigs, paper_cm2,
                  100.0 * abs(got - paper_cm2) / paper_cm2))
    _check("B3 n_min scales exactly as w_OP^2 between the two energies",
           abs(n_min_branch_cut(160.0) / n_min_branch_cut(55.0)
               - (160.0 / 55.0) ** 2) < 1e-9,
           "ratio = %.10f, (160/55)^2 = %.10f"
           % (n_min_branch_cut(160.0) / n_min_branch_cut(55.0),
              (160.0 / 55.0) ** 2))
    return t, shares


def verify_exact_zeros_and_limits():
    """Exact zeros, exact identities and bitwise recoveries."""
    # --- Z: zero dissipated power must give EXACTLY zero heating --------
    z = solve_bias_point(V_g=3.0, V_ds=0.0)
    _check("Z1 V_ds = 0 gives Id EXACTLY 0.0", z["Id"] == 0.0,
           "Id = %r" % z["Id"])
    _check("Z2 V_ds = 0 gives P EXACTLY 0.0", z["P_total"] == 0.0,
           "P = %r" % z["P_total"])
    _check("Z3 V_ds = 0 gives dT EXACTLY 0.0 (not merely small)",
           z["dT"] == 0.0, "dT = %r" % z["dT"])
    _check("Z4 V_ds = 0 gives T bitwise equal to T0",
           z["T"] == T_AMBIENT, "T = %r" % z["T"])

    # --- R: R_th -> 0 must recover the isothermal model BITWISE ---------
    for vg, vds in ((3.0, 0.05), (3.0, 0.5), (1.5, 1.0)):
        hot0 = solve_bias_point(vg, vds, self_heating=True, R_th_override=0.0)
        iso = solve_bias_point(vg, vds, self_heating=False)
        _check("R1 R_th=0 recovers isothermal T bitwise at "
               "(Vg=%.1f,Vds=%.2f)" % (vg, vds),
               hot0["T"] == iso["T"] == T_AMBIENT,
               "T_hot0 = %r, T_iso = %r" % (hot0["T"], iso["T"]))
        _check("R2 R_th=0 recovers isothermal Id bitwise at "
               "(Vg=%.1f,Vds=%.2f)" % (vg, vds),
               hot0["Id"] == iso["Id"],
               "Id_hot0 = %r, Id_iso = %r" % (hot0["Id"], iso["Id"]))

    # --- M: the Eq. (4) temperature factor at T_REF ---------------------
    for n in (1e15, N_REF, 1e17):
        mp = mobility(n, T_REF, "product")
        ms = mobility(n, T_REF, "sum")
        _check("M1 the two Eq. (4) readings are BITWISE equal at T=T_REF "
               "(n=%.3e)" % n, mp == ms,
               "product = %r, sum = %r" % (mp, ms))
    _check("M2 T_REF temperature factor is EXACTLY 1 "
           "(mu(n,T_REF) == mu0/(1+(n/n_ref)^a))",
           mobility(N_REF, T_REF, "product") == MU0 / 2.0,
           "mu(n_ref,300) = %r, mu0/2 = %r"
           % (mobility(N_REF, T_REF, "product"), MU0 / 2.0))
    # M4/M5 exist because the first draft of this suite was BLIND to both
    # Eq. (4) exponents: every other exact check is anchored at T = T_REF,
    # where the temperature factor is identically 1 for any beta, and at
    # n = N_REF, where the density factor is identically 2 for any alpha.
    # A suite that cannot see a parameter cannot defend it.  Both constants
    # below are hand-computed literals, NOT expressions in ALPHA/BETA, or
    # the check would be circular and pass for any value.
    _check("M4 mu(2*n_ref, T_REF) == mu0/(1 + 2^2.2) exactly "
           "(this is what pins ALPHA; nothing else does)",
           abs(mobility(2.0 * N_REF, T_REF, "product")
               / (MU0 / (1.0 + 4.5947934199881175)) - 1.0) < 1e-12,
           "mu = %.12e, mu0/(1+2^2.2) = %.12e"
           % (mobility(2.0 * N_REF, T_REF, "product"),
              MU0 / (1.0 + 4.5947934199881175)))
    _check("M5 mu(n->0, 450 K) == mu0/1.125 exactly "
           "(this is what pins BETA; T = 2*T_REF would not, since "
           "1^beta = 1 for any beta)",
           abs(mobility(0.0, 450.0, "product") / (MU0 / 1.125) - 1.0)
           < 1e-12,
           "mu = %.12e, mu0/1.125 = %.12e"
           % (mobility(0.0, 450.0, "product"), MU0 / 1.125))
    _check("M3 the two readings DIVERGE off T_REF "
           "(so M1 is not vacuous)",
           mobility(N_REF, 500.0, "product")
           != mobility(N_REF, 500.0, "sum"),
           "500K product = %.6e, sum = %.6e, ratio = %.6f"
           % (mobility(N_REF, 500.0, "product"),
              mobility(N_REF, 500.0, "sum"),
              mobility(N_REF, 500.0, "product")
              / mobility(N_REF, 500.0, "sum")))

    # --- V: Eq. (7) limits ----------------------------------------------
    _check("V1 v_sat REFUSES below the Eq. (7) branch cut rather than "
           "clamping", _raises(lambda: saturation_velocity(
               0.5 * n_min_branch_cut(V_SAT_PHONON_MEV), 300.0)),
           "ValueError raised as required")
    _check("V2 the T->0 Bose limit N_OP->0 gives v_sat EXACTLY the "
           "unoccupied form",
           abs(saturation_velocity(1e17, 1.0)
               / _v_sat_unoccupied(1e17) - 1.0) < 1e-12,
           "ratio = %.16f" % (saturation_velocity(1e17, 1.0)
                              / _v_sat_unoccupied(1e17)))
    _check("V3 v_sat is monotonically DECREASING in T (Dorgan's ~20%% / "
           "~2%% statement has this sign)",
           saturation_velocity(1e17, 500.0) < saturation_velocity(1e17, 300.0),
           "v_sat(300K) = %.4e, v_sat(500K) = %.4e m/s"
           % (saturation_velocity(1e17, 300.0),
              saturation_velocity(1e17, 500.0)))
    _check("V4 dT >= 0 at every bias swept (a cooling channel would be a "
           "sign error)", True, "checked in sweep_vds(), see table")

    # --- D: the drift-velocity form's own low-field limit ---------------
    mu_t, vs_t = 0.4, 5e5
    F_small = 1.0
    _check("D1 Eq. (6) -> mu*F exactly in the low-field limit",
           abs(drift_velocity(F_small, mu_t, vs_t) / (mu_t * F_small) - 1.0)
           < 1e-12,
           "v/muF - 1 = %.3e"
           % (drift_velocity(F_small, mu_t, vs_t) / (mu_t * F_small) - 1.0))
    _check("D2 Eq. (6) -> v_sat from below as F -> inf",
           drift_velocity(1e12, mu_t, vs_t) < vs_t
           and drift_velocity(1e12, mu_t, vs_t) / vs_t > 0.999999,
           "v(F=1e12)/v_sat = %.12f"
           % (drift_velocity(1e12, mu_t, vs_t) / vs_t))


def _v_sat_unoccupied(n):
    """Eq. (7) with the Bose factor dropped (N_OP = 0), for the V2 check."""
    w_OP = V_SAT_PHONON_MEV * 1e-3 * e / hbar
    root = np.sqrt(1.0 - w_OP ** 2 / (4.0 * np.pi * n * V_F ** 2))
    return (2.0 / np.pi) * (w_OP / np.sqrt(np.pi * n)) * root


def _round_sig(x, sigfigs):
    """Round x to the given number of significant figures."""
    if x == 0.0:
        return 0.0
    from math import floor, log10
    return round(x, -int(floor(log10(abs(x)))) + (sigfigs - 1))


def _rounds_to(got, quoted, sigfigs):
    """Does `got` agree with a value `quoted` to `sigfigs` significant
    figures?  This is the correct comparison against a published rounded
    number; a relative tolerance is not (see the B2 note)."""
    return _round_sig(got, sigfigs) == _round_sig(quoted, sigfigs)


def _raises(fn):
    try:
        fn()
    except ValueError:
        return True
    except Exception:
        return False
    return False


def overlap_with_isothermal_model(V_g=3.0, V_ds=1e-4):
    """Validate the NEW lumped model against the EXISTING repo model in the
    limit where both are supposed to agree: low field, isothermal, so that
    Eq. (6) -> mu*F and the channel is ohmic.

    The residual is MEASURED AND REPORTED, not asserted away.  The two
    models differ in two declared ways -- this one uses Dorgan's
    mu(n,T) where gfet uses a constant mobility, and this one is
    lumped where gfet averages 50 segments -- so a nonzero residual is
    EXPECTED and its size is the point.
    """
    r = solve_bias_point(V_g, V_ds, self_heating=False)
    G_new = r["Id"] / r["V_ch"] if r["V_ch"] else float("nan")
    R_old = float(gfet.channel_resistance(V_g, 0.0))
    G_old = 1.0 / R_old
    mu_dorgan = r["mu"]
    # with mu scaled to the repo's constant mu, the two must agree to the
    # segment-averaging difference alone -- which at V_ds = 1e-4 is tiny
    G_new_mu_matched = G_new * (gfet.mu / mu_dorgan)
    rel = abs(G_new_mu_matched - G_old) / G_old
    _check("O1 low-field conductance matches the existing repo model to "
           "<0.5 %% once the mobility CONSTANT is matched",
           rel < 5e-3,
           "G_new(mu-matched) = %.6e S, G_old = %.6e S, rel = %.3e"
           % (G_new_mu_matched, G_old, rel))
    return {"G_new": G_new, "G_old": G_old, "rel_mu_matched": rel,
            "mu_dorgan": mu_dorgan, "mu_repo": gfet.mu,
            "mu_ratio": mu_dorgan / gfet.mu, "n_bar": r["n_bar"]}


def sweep_vds(V_g=3.0, V_ds_max=1.5, N=61, reading="product",
              phonon_meV=V_SAT_PHONON_MEV):
    """Self-consistent dT, Id and parameter degradation vs V_ds."""
    vds = np.linspace(0.0, V_ds_max, N)
    rows = [solve_bias_point(V_g, float(v), reading=reading,
                             phonon_meV=phonon_meV) for v in vds]
    iso = [solve_bias_point(V_g, float(v), self_heating=False,
                            reading=reading, phonon_meV=phonon_meV)
           for v in vds]
    if not all(r["dT"] >= 0.0 for r in rows):
        raise AssertionError("dT < 0 somewhere: sign error in the network")
    if not all(r["converged"] for r in rows):
        bad = [r["V_ds"] for r in rows if not r["converged"]]
        raise AssertionError("thermal fixed point did not converge at "
                             "V_ds = %r" % bad)
    return vds, rows, iso


def vg_for_branch_cut(phonon_meV, V_dirac=None):
    """Gate voltage at which this device's carrier density first reaches the
    Eq. (7) branch cut for a given phonon energy.

    Found by bisection on gfet.carrier_density, so it uses the repository's
    own electrostatics (quantum capacitance and puddle floor included)
    rather than the plain C_ox*V/e estimate.
    """
    V_dirac = gfet.V_dirac if V_dirac is None else V_dirac
    n_cut = n_min_branch_cut(phonon_meV)
    lo, hi = V_dirac, V_dirac + 200.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if abs(float(gfet.carrier_density(mid, 0.0))) < n_cut:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def verify_branch_cut_reachability(V_g=3.0):
    """The 160 meV branch's inadmissibility at this device's bias is a
    RESULT, so it gets a check rather than a sentence.

    Both-sided on purpose: the branch must be refused at the operating
    point AND admitted above vg_for_branch_cut(), or the guard would be
    indistinguishable from a function that always refuses.
    """
    n_op = _mean_density(V_g, 0.0)
    _check("C1 160 meV branch is INADMISSIBLE at the operating point "
           "(this is a physics result, not a bug)",
           n_op < n_min_branch_cut(160.0),
           "n(V_g=%.1f) = %.4e m^-2 < n_min(160 meV) = %.4e m^-2"
           % (V_g, n_op, n_min_branch_cut(160.0)))
    _check("C2 82 meV branch IS admissible at the operating point "
           "(so C1 is a statement about 160 meV, not about the guard)",
           n_op > n_min_branch_cut(82.0),
           "n(V_g=%.1f) = %.4e m^-2 > n_min(82 meV) = %.4e m^-2"
           % (V_g, n_op, n_min_branch_cut(82.0)))
    vg_cut = vg_for_branch_cut(160.0)
    _check("C3 above vg_for_branch_cut(160 meV) the SAME call succeeds "
           "(the guard is not a constant refusal)",
           saturation_velocity(_mean_density(vg_cut + 0.5, 0.0), 300.0,
                               160.0) > 0.0,
           "V_g = %.4f V admits 160 meV; v_sat = %.4e m/s"
           % (vg_cut + 0.5,
              saturation_velocity(_mean_density(vg_cut + 0.5, 0.0),
                                  300.0, 160.0)))
    _check("C4 and it is refused just BELOW that gate voltage "
           "(two-sided, so C3 is not vacuous)",
           _raises(lambda: saturation_velocity(
               _mean_density(vg_cut - 0.5, 0.0), 300.0, 160.0)),
           "V_g = %.4f V refuses 160 meV as required" % (vg_cut - 0.5))
    return vg_cut


def thermal_current_peak(V_g=3.0, V_ds_max=1.5, N=151):
    """Locate the self-heating-induced current maximum.

    The isothermal model's Id is monotone in V_ds; the self-consistent one
    is not, because mobility loss outruns the added bias.  That turning
    point is a thermal negative-differential-resistance onset and it is a
    prediction of the coupled model that neither half makes alone.
    """
    vds, rows, iso = sweep_vds(V_g=V_g, V_ds_max=V_ds_max, N=N)
    Id = np.array([r["Id"] for r in rows])
    Id_iso = np.array([r["Id"] for r in iso])
    i_pk = int(np.argmax(Id))
    i_pk_iso = int(np.argmax(Id_iso))
    #
    # OVERCLAIM, CAUGHT BY THIS CHECK BEFORE THE COMMIT (2026-11-10).
    # N1 was first written as "the isothermal current is MONOTONE in V_ds,
    # so N2's turning point is thermal in origin".  It FAILED:
    # min diff(Id_iso) = -2.0e-07 A.  The isothermal model turns over too,
    # near V_ds ~ 1.1 V, because the drain bias lowers the channel carrier
    # density faster than it raises the field -- the ordinary GFET
    # quasi-saturation kink, nothing to do with temperature.
    #
    # So self-heating does NOT create the turning point.  It moves it
    # earlier and deepens it, which is a smaller claim and the true one.
    # The checks below state the true claim and are two-sided about it.
    _check("N1 the ISOTHERMAL current ALSO has an interior maximum "
           "(self-heating does NOT create the turning point)",
           0 < i_pk_iso < len(Id_iso) - 1,
           "isothermal peak at V_ds = %.4f V, Id = %.4f uA; "
           "min diff(Id_iso) = %.3e A"
           % (vds[i_pk_iso], Id_iso[i_pk_iso] * 1e6,
              float(np.min(np.diff(Id_iso)))))
    _check("N2 self-heating moves the peak to STRICTLY LOWER V_ds",
           vds[i_pk] < vds[i_pk_iso],
           "thermal peak %.4f V vs isothermal peak %.4f V "
           "(shift %.4f V earlier)"
           % (vds[i_pk], vds[i_pk_iso], vds[i_pk_iso] - vds[i_pk]))
    _check("N3 and to STRICTLY LOWER peak current",
           Id[i_pk] < Id_iso[i_pk_iso],
           "thermal peak %.4f uA vs isothermal peak %.4f uA "
           "(%.2f %% lower)"
           % (Id[i_pk] * 1e6, Id_iso[i_pk_iso] * 1e6,
              100.0 * (1.0 - Id[i_pk] / Id_iso[i_pk_iso])))
    _check("N4 the two peaks are DISTINGUISHABLE on this sweep grid "
           "(so N2/N3 are not reporting grid noise)",
           abs(i_pk - i_pk_iso) >= 3,
           "peak indices %d vs %d on an N=%d grid (dV_ds = %.5f V)"
           % (i_pk, i_pk_iso, len(vds), vds[1] - vds[0]))
    return {"V_ds_peak": float(vds[i_pk]), "Id_peak": float(Id[i_pk]),
            "V_ds_peak_iso": float(vds[i_pk_iso]),
            "Id_peak_iso": float(Id_iso[i_pk_iso]),
            "Id_end": float(Id[-1]), "V_ds_end": float(vds[-1])}


# ---------------------------------------------------------------------------
# Figure and transcript
# ---------------------------------------------------------------------------
def plot_self_heating(fname="self_heating_electrothermal.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    V_g = 3.0
    vds, rows, iso = sweep_vds(V_g=V_g)
    dT = np.array([r["dT"] for r in rows])
    Id = np.array([r["Id"] for r in rows])
    Id_iso = np.array([r["Id"] for r in iso])
    mu_h = np.array([r["mu"] for r in rows])
    mu_i = np.array([r["mu"] for r in iso])
    vs_h = np.array([r["v_sat"] for r in rows])
    vs_i = np.array([r["v_sat"] for r in iso])

    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.4))

    a = ax[0, 0]
    a.plot(vds, dT, "r-", lw=2)
    a.axhline(DT_REPRODUCIBILITY_CEILING, color="k", ls="--", lw=1.2)
    a.text(0.03, DT_REPRODUCIBILITY_CEILING * 1.05,
           "Dorgan fit-exclusion ceiling  $\\Delta T \\approx 200$ K",
           fontsize=8)
    a.set_xlabel("$V_{ds}$ (V)")
    a.set_ylabel("$\\Delta T$ (K)")
    a.set_title("Self-consistent channel temperature rise\n"
                "$W$=1 $\\mu$m, $L$=200 nm, $t_{ox}$=90 nm, "
                "$R_{th}\\approx$ %.2e K/W" % rows[-1]["R_th"], fontsize=9)
    a.grid(alpha=0.3)

    a = ax[0, 1]
    a.plot(vds, Id_iso * 1e6, "b--", lw=2, label="isothermal (300 K)")
    a.plot(vds, Id * 1e6, "r-", lw=2, label="self-consistent")
    a.set_xlabel("$V_{ds}$ (V)")
    a.set_ylabel("$I_d$ ($\\mu$A)")
    a.set_title("Drain current: the cost of self-heating", fontsize=9)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[1, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        a.plot(vds, 100.0 * (1.0 - Id / Id_iso), "r-", lw=2, label="$I_d$")
        a.plot(vds, 100.0 * (1.0 - mu_h / mu_i), "g-", lw=2, label="$\\mu$")
        a.plot(vds, 100.0 * (1.0 - vs_h / vs_i), "m-", lw=2,
               label="$v_{sat}$")
    a.set_xlabel("$V_{ds}$ (V)")
    a.set_ylabel("degradation vs isothermal (%)")
    a.set_title("Where the current loss comes from", fontsize=9)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[1, 1]
    for meV, c in ((55.0, "tab:blue"), (82.0, "k"), (160.0, "tab:orange")):
        try:
            _, rr, _ = sweep_vds(V_g=V_g, phonon_meV=meV)
        except ValueError:
            # Eq. (7) has no meaning below its branch cut at this density;
            # the curve is ABSENT and said to be absent, never clamped in.
            a.plot([], [], color=c, lw=1.8, ls="-.",
                   label="$\\hbar\\omega_{OP}$ = %.0f meV: inadmissible\n"
                         "(below Eq. (7) branch cut)" % meV)
            continue
        a.plot(vds, [r["dT"] for r in rr], color=c, lw=1.8,
               label="$\\hbar\\omega_{OP}$ = %.0f meV" % meV)
    for rd, ls in (("product", "-"), ("sum", ":")):
        _, rr, _ = sweep_vds(V_g=V_g, reading=rd)
        a.plot(vds, [r["dT"] for r in rr], "r", ls=ls, lw=1.2,
               label="Eq.(4) %s reading" % rd)
    a.set_xlabel("$V_{ds}$ (V)")
    a.set_ylabel("$\\Delta T$ (K)")
    a.set_title("Spread from the two UNSETTLED model choices", fontsize=9)
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    return fname


def summary_numbers():
    print("=" * 74)
    print("SELF-CONSISTENT ELECTRO-THERMAL MODEL -- graphene FET on SiO2")
    print("Dorgan, Bae & Pop, APL 97, 082112 (2010) / arXiv:1005.2711")
    print("=" * 74)

    t, shares = verify_benchmarks()
    print("\nBENCHMARK A -- the paper's OWN device (W=7um, L=4um, "
          "t_ox=300nm)")
    print("  %-10s %14s %9s   %s" % ("term", "value (K/W)", "share",
                                     "paper"))
    for k, paper in (("R_B", "4 %"), ("R_ox", "84 %"), ("R_Si", "12 %")):
        print("  %-10s %14.1f %8.2f %%   %s" % (k, t[k], shares[k], paper))
    print("  %-10s %14.1f %8s     %s" % ("TOTAL", t["R_th"], "",
                                         "~1e4 K/W"))
    print("  kappa_ox(300 K) = %.6f W/m.K   kappa_Si(300 K) = %.3f W/m.K"
          % (kappa_ox(300.0), kappa_Si(300.0)))

    print("\nBENCHMARK B -- Eq. (7) branch cut n_min = (w_OP/v_F)^2/(2 pi)")
    for meV, paper in ((160.0, "9.4e11"), (82.0, "(not quoted)"),
                       (55.0, "1.1e11")):
        print("  hbar*w_OP = %5.1f meV  ->  n_min = %.4e cm^-2   "
              "paper: %s" % (meV, n_min_branch_cut(meV) * 1e-4, paper))

    print("\nTHIS repository's device (W=%.2g m, L=%.2g m, t_ox=%.2g m)"
          % (gfet.W, gfet.L, gfet.t_ox))
    td = thermal_resistances(gfet.W, gfet.L, gfet.t_ox)
    for k in ("R_B", "R_ox", "R_Si"):
        print("  %-6s %12.4e K/W  (%5.2f %%)"
              % (k, td[k], 100.0 * td[k] / td["R_th"]))
    print("  %-6s %12.4e K/W   -- %.1fx Dorgan's device"
          % ("R_th", td["R_th"], td["R_th"] / t["R_th"]))
    print("  A = %.3e m^2, %.1fx SMALLER in area than Dorgan's"
          % (td["A"], t["A"] / td["A"]))
    print("  CAVEAT: Eq. (5) assumes W, L >> t_ox; here L/t_ox = %.2f, so"
          % (gfet.L / gfet.t_ox))
    print("          lateral spreading is neglected and R_ox is an "
          "OVER-estimate.")
    print("          R_th above is therefore an UPPER BOUND, and every dT "
          "below")
    print("          is a WORST CASE, not a centred estimate.")

    print("\nOVERLAP with the existing isothermal repo model "
          "(low field, 300 K)")
    o = overlap_with_isothermal_model()
    print("  G_new = %.6e S   G_old = %.6e S" % (o["G_new"], o["G_old"]))
    print("  mu(Dorgan, n=%.3e) = %.4f  vs  gfet.mu = %.4f  "
          "(ratio %.4f)" % (o["n_bar"], o["mu_dorgan"], o["mu_repo"],
                            o["mu_ratio"]))
    print("  residual once the mobility constant is matched: %.3e"
          % o["rel_mu_matched"])
    print("  NOTE: the mobility ratio is NOT a defect -- Dorgan's Eq. (4)")
    print("        at this density gives a different number from the repo's")
    print("        hand-set mu, and that disagreement is itself a result.")

    print("\nBIAS SWEEP at V_g = 3.0 V (self-consistent)")
    print("  %-7s %9s %9s %9s %9s %8s %7s %s"
          % ("V_ds", "Id(uA)", "Id_iso", "dT(K)", "P(uW)", "mu", "v_sat",
             "iters"))
    vds, rows, iso = sweep_vds(V_g=3.0)
    for i in range(0, len(vds), 6):
        r, ri = rows[i], iso[i]
        print("  %-7.3f %9.3f %9.3f %9.3f %9.3f %8.4f %7.3e %5d%s"
              % (r["V_ds"], r["Id"] * 1e6, ri["Id"] * 1e6, r["dT"],
                 r["P_total"] * 1e6, r["mu"], r["v_sat"], r["iters"],
                 "  <-- PAST CEILING" if r["past_reproducibility_ceiling"]
                 else ""))

    rf = solve_bias_point(3.0, 0.05)
    rf_iso = solve_bias_point(3.0, 0.05, self_heating=False)
    print("\nTHE RF CHAPTER'S OWN BIAS (V_ds = 50 mV), which is what the")
    print("10-06 item was created to bound:")
    print("  dT = %.4f K, P = %.4f uW, Id = %.4f uA (isothermal %.4f uA)"
          % (rf["dT"], rf["P_total"] * 1e6, rf["Id"] * 1e6,
             rf_iso["Id"] * 1e6))
    print("  relative current loss from self-heating: %.3e"
          % (1.0 - rf["Id"] / rf_iso["Id"]))
    first_past = [r for r in rows if r["past_reproducibility_ceiling"]]
    print("  lowest swept V_ds past Dorgan's ~200 K ceiling: %s"
          % ("%.3f V" % first_past[0]["V_ds"] if first_past
             else "none in sweep"))

    print("\nUNSETTLED CHOICES, reported as a SPREAD rather than picked")
    for meV in (55.0, 82.0, 160.0):
        try:
            r = solve_bias_point(3.0, 1.0, phonon_meV=meV)
        except ValueError as exc:
            print("  hbar*w_OP = %5.1f meV: REFUSED -- %s"
                  % (meV, str(exc).split(": ", 1)[1]))
            print("           this is the Eq. (7) validity guard working, "
                  "not a defect:")
            print("           at V_g = 3.0 V this device's density is BELOW "
                  "the branch cut,")
            print("           so the graphene-OP branch cannot bracket "
                  "v_sat(T) here at all.")
            print("           V_g required to reach it: %.3f V "
                  "(see vg_for_branch_cut())" % vg_for_branch_cut(meV))
            continue
        print("  hbar*w_OP = %5.1f meV: dT(V_ds=1V) = %8.3f K, "
              "Id = %8.3f uA" % (meV, r["dT"], r["Id"] * 1e6))
    for rd in ("product", "sum"):
        r = solve_bias_point(3.0, 1.0, reading=rd)
        print("  Eq. (4) %-8s reading: dT(V_ds=1V) = %8.3f K, "
              "Id = %8.3f uA" % (rd, r["dT"], r["Id"] * 1e6))

    vg_cut = verify_branch_cut_reachability()
    print("\n  V_g needed to reach the 160 meV branch cut on THIS device: "
          "%.4f V" % vg_cut)
    print("  (operating point is V_g = 3.0 V, so the graphene-OP branch is")
    print("   out of reach and the v_sat(T) spread cannot be bracketed by")
    print("   the paper's two endpoints at this bias.)")

    pk = thermal_current_peak()
    print("\nWHERE THE CURRENT TURNS OVER -- and what self-heating does "
          "and does NOT do")
    print("  isothermal      peak: V_ds = %.4f V, Id = %.3f uA"
          % (pk["V_ds_peak_iso"], pk["Id_peak_iso"] * 1e6))
    print("  self-consistent peak: V_ds = %.4f V, Id = %.3f uA"
          % (pk["V_ds_peak"], pk["Id_peak"] * 1e6))
    print("  CORRECTED CLAIM: self-heating does NOT create the turning")
    print("  point -- the isothermal model turns over too (check N1), from")
    print("  drain-induced density loss, the ordinary GFET kink.  What")
    print("  self-heating does is move it %.4f V EARLIER and %.2f %% LOWER."
          % (pk["V_ds_peak_iso"] - pk["V_ds_peak"],
             100.0 * (1.0 - pk["Id_peak"] / pk["Id_peak_iso"])))
    print("  The first draft of check N1 asserted monotonicity and FAILED,")
    print("  which is how the overclaim was caught before being committed;")
    print("  see the comment at thermal_current_peak().")

    verify_exact_zeros_and_limits()
    print("\n" + "=" * 74)
    print("CHECKS")
    print("=" * 74)
    npass = 0
    for name, ok, detail in _CHECKS:
        print("  [%s] %s" % ("PASS" if ok else "FAIL", name))
        print("         %s" % detail)
        npass += ok
    print("\n  %d/%d passed, %d failed"
          % (npass, len(_CHECKS), len(_CHECKS) - npass))
    return npass, len(_CHECKS)


if __name__ == "__main__":
    import sys as _sys
    # --checks-only skips the figure, which costs five extra bias sweeps.
    # graphene_self_heating_mutation.py uses it so that eleven subprocess
    # runs fit in one session step; it changes NO check.
    _checks_only = "--checks-only" in _sys.argv
    npass, ntot = summary_numbers()
    if not _checks_only:
        f = plot_self_heating()
        print("\nFigure written: %s" % f)
    raise SystemExit(0 if npass == ntot else 1)
