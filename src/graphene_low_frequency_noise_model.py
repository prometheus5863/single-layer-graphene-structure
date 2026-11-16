"""
graphene_low_frequency_noise_model.py

Low-frequency (1/f) noise in a back-gated single-layer graphene FET:
charge-number fluctuation (CNF, dN), correlated mobility fluctuation (CMF)
and Hooge mobility fluctuation (dmu), evaluated on THIS repository's own
GFET transfer characteristic.

WHY THIS IS A DEVICE-PHYSICS TOPIC AND NOT AN AFTERTHOUGHT
----------------------------------------------------------
1/f noise is the quantity that decides whether a transistor technology is
usable in analog and RF front ends: it up-converts to oscillator phase
noise (as (1/f)^3 near carrier), sets the 1/f corner, and -- unlike thermal
noise -- cannot be averaged away by measuring longer.  For graphene the
question is sharper than for Si, because graphene has no gap, so the
device is never "off" and the carrier number N passes through a minimum
(not through zero) at the charge-neutrality point.  Both noise mechanisms
below are functions of N, so the gate dependence of the noise is a direct
probe of which mechanism is operating -- which is why the experimental
literature classifies graphene noise by the SHAPE of S_I/I^2 versus V_g
(V-shape, M-shape, Lambda-shape) rather than by a single number.

MODELS IMPLEMENTED
------------------
(i) CNF / number fluctuation, simplified homogeneous-channel form:

        S_I/I^2 |_dN = (g_m/I_d)^2 * S_Vfb,
        S_Vfb        = q^2 * k_B * T * lambda * N_T / (W * L * C_ox^2 * f)

    N_T = dielectric volumetric trap density per unit energy [1/(eV cm^3)],
    lambda = trap tunnelling attenuation length (~0.1 nm).

(ii) CNF+CMF, the correlated form, in which a trapped charge changes the
     mobility as well as the carrier number:

        S_I/I^2 |_dN+dmu = (g_m/I_d + alpha_c * mu * C_ox)^2
                           * q^2 k_B T lambda N_T / (W L C_ox^2 f)

     alpha_c = Coulomb scattering coefficient [V s / C].  Setting
     alpha_c = 0 must return (i) EXACTLY; that is check V3 below.

     *** A RESULT OF THIS SESSION, AND THE REASON TWO VARIANTS OF (ii)
     *** ARE IMPLEMENTED RATHER THAN ONE.
     That form is written for a unipolar MOSFET, where g_m > 0 at every
     usable bias, so the sign of the correlated term is a convention and
     the literature writes it as "+-".  Graphene is AMBIPOLAR: g_m changes
     sign at the charge-neutrality point, so on the hole branch g_m/I_d is
     negative and the sum can pass through zero.  Transcribed verbatim, the
     form therefore predicts a bias at which a trap-limited device has NO
     1/f noise from its dominant mechanism -- here at g_m/I_d
     = -alpha_c*mu*C_ox = -6.906e-2 1/V, i.e. Vg ~ -1.3 V, on the hole side
     of a device whose electron-side mirror bias is one of the two noise
     MAXIMA.  No measurement reports such a notch, and nothing physical
     selects that particular gate voltage: alpha_c*mu*C_ox is a property of
     the oxide and the mobility, while g_m/I_d is a property of the bias
     point, and their crossing is an artefact of the imported sign.
     Both are kept:
        s_cnf_cmf_signed()  -- the verbatim transcription, retained so the
                               artefact is reproducible and checkable (V7
                               pins its exact zero);
        s_cnf_cmf()         -- branch-corrected, (|g_m/I_d| + alpha_c mu
                               C_ox)^2, which is what every result below
                               uses.  Coulomb scattering from a trapped
                               charge raises the scattering rate on both
                               branches, so the correlated term adds in
                               magnitude; taking |g_m/I_d| is the minimal
                               change that makes the unipolar form
                               branch-agnostic and leaves the electron
                               branch bitwise unchanged (check V8).

(iii) Hooge / uncorrelated mobility fluctuation:

        S_I/I^2 |_Hooge = alpha_H / (f * N),   N = n_sheet * W * L

PROVENANCE OF EVERY IMPORTED NUMBER, AND ITS STRENGTH
-----------------------------------------------------
See notes/2026-11-15-low-frequency-noise-in-gfets.md for the full citation
discussion.  Two strengths are distinguished there and here:

  VERBATIM  (read in the source's own text)
    - alpha_H = 1.3e-3, N_T = 1.1e20 1/(eV cm^3): Mavredakis et al.,
      arXiv:1806.06589, Table 1, L = 13 um electrolyte-top-gated SLG row.
    - lambda ~ 0.1 nm: same paper, main text.
    - S_I/I^2 = 1e-9 .. 1e-7 Hz^-1 at f = 10 Hz, and area-normalised
      1e-8 .. 1e-7 um^2/Hz: Balandin, arXiv:1307.4797, Sec. II.
    - 1/f corner 1 .. 100 kHz; (1/f)^3 phase-noise up-conversion: same
      review, Sec. II and Introduction.
  SUMMARY-ONLY  (reported from the paper's SUPPLEMENT, which was not read
  verbatim -- so it is used but flagged, not asserted)
    - the (g_m/I_d + alpha_c mu C_t)^2 form (Supp. Eq. A25) and
      alpha_c = 450 V s / C.
  THIS REPOSITORY (not from any paper)
    - C_ox, W, L, mu, V_dirac, n_puddle and the I_d(V_g) curve itself all
      come from graphene_fet_model.py, imported unchanged.

The N_T = 1.1e20 value is NOT a SiO2 number -- it is fitted to an
electrolyte-gated device and the same paper says it is far above the Si
range of 1e16 .. 1e18.  So both are evaluated below, and the difference
between them turns out to decide the SHAPE of the gate dependence, which
is this script's main result.
"""

import numpy as np; from graphene.paths import results_path  # noqa: E702 -- same line on purpose: no line numbers move (docs/INVENTORY.md §6)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.constants import e as q, Boltzmann as kB

import graphene_fet_model as fet

# ---------------------------------------------------------------- parameters
T = 300.0                 # K
F_MEAS = 10.0             # Hz -- the frequency the experimental literature
                          # quotes S_I/I^2 at (Balandin Sec. II)
LAMBDA_TUNNEL = 0.1e-9    # m
ALPHA_H = 1.3e-3          # dimensionless (arXiv:1806.06589 Table 1, L=13um)
ALPHA_C = 450.0           # V s / C  (SUMMARY-ONLY provenance, see above)
NT_PAPER_CM = 1.1e20      # 1/(eV cm^3), fitted, electrolyte-gated SLG
NT_SIO2_CM = 1.0e17       # 1/(eV cm^3), mid Si/SiO2 range quoted by the
                          # same paper as the conventional interval

VDS = 0.05                # V, the same small drain bias the repo's
                          # transfer_characteristic() is exercised at
VDS_SYMMETRIC = 1.0e-6    # V, the Vds -> 0 limit used by check V2


def nt_si(nt_per_eV_cm3):
    """Convert a trap density from 1/(eV cm^3) to SI 1/(J m^3)."""
    return nt_per_eV_cm3 * (1.0 / q) * 1.0e6


def s_vfb(nt_per_eV_cm3, f=F_MEAS, W=None, L=None, C_ox=None):
    """Flat-band-voltage noise PSD [V^2/Hz]."""
    W = fet.W if W is None else W
    L = fet.L if L is None else L
    C_ox = fet.C_ox if C_ox is None else C_ox
    return (q**2 * kB * T * LAMBDA_TUNNEL * nt_si(nt_per_eV_cm3)
            / (W * L * C_ox**2 * f))


def gm_over_id(Vg, Id):
    """g_m/I_d on a uniform Vg grid, by central differences.

    np.gradient is exact for a linear I(Vg) up to floating-point rounding,
    which is what check V1 exploits.
    """
    return np.gradient(Id, Vg) / Id


def s_cnf(gm_id, nt_per_eV_cm3, f=F_MEAS, **kw):
    """CNF only: (g_m/I_d)^2 * S_Vfb  [1/Hz]."""
    return gm_id**2 * s_vfb(nt_per_eV_cm3, f=f, **kw)


def s_cnf_cmf_signed(gm_id, nt_per_eV_cm3, f=F_MEAS, alpha_c=ALPHA_C,
                     mu=None, C_ox=None, **kw):
    """CNF+CMF exactly as transcribed from the unipolar literature form:
    (g_m/I_d + alpha_c mu C_ox)^2 * S_Vfb.

    Retained ONLY so the ambipolar artefact it produces stays reproducible
    (see the module docstring and check V7).  Not used for any result.
    """
    mu = fet.mu if mu is None else mu
    C_ox = fet.C_ox if C_ox is None else C_ox
    return (gm_id + alpha_c * mu * C_ox)**2 * s_vfb(
        nt_per_eV_cm3, f=f, C_ox=C_ox, **kw)


def s_cnf_cmf(gm_id, nt_per_eV_cm3, f=F_MEAS, alpha_c=ALPHA_C, mu=None,
              C_ox=None, **kw):
    """Branch-corrected CNF+CMF: (|g_m/I_d| + alpha_c mu C_ox)^2 * S_Vfb.

    Identical to s_cnf_cmf_signed() wherever g_m >= 0 (check V8), and free
    of the hole-branch cancellation notch (check V7b).
    """
    mu = fet.mu if mu is None else mu
    C_ox = fet.C_ox if C_ox is None else C_ox
    return (np.abs(gm_id) + alpha_c * mu * C_ox)**2 * s_vfb(
        nt_per_eV_cm3, f=f, C_ox=C_ox, **kw)


def cmf_term():
    """alpha_c * mu * C_ox  [1/V] -- the bias-independent correlated term."""
    return ALPHA_C * fet.mu * fet.C_ox


def id_zero_vds_limit(Vg, Vds=1.0):
    """I_d in the EXACT Vds -> 0 limit: Vds/(R_channel(Vg, V_ch=0) + R_c).

    transfer_characteristic() averages the local channel resistance over a
    0..Vds quasi-Fermi profile, so its output is only asymptotically even in
    (Vg - V_dirac); this is the limit itself, which is exactly even because
    carrier_density() depends on (Vg - V_dirac) through its magnitude only.
    Used by check V2, which needs an exact symmetry rather than a converging
    one.  The Vds prefactor cancels in S_I/I^2.
    """
    return Vds / (np.array([fet.channel_resistance(v, 0.0) for v in Vg])
                  + fet.Rc_total)


def s_hooge(n_sheet, f=F_MEAS, alpha_H=ALPHA_H, W=None, L=None):
    """Hooge: alpha_H/(f N), N = n_sheet * W * L  [1/Hz]."""
    W = fet.W if W is None else W
    L = fet.L if L is None else L
    return alpha_H / (f * n_sheet * W * L)


# ------------------------------------------------------------------- sweeps
def noise_sweep(Vg, nt_per_eV_cm3, Vds=VDS, alpha_H=ALPHA_H,
                alpha_c=ALPHA_C, f=F_MEAS):
    """Everything on one gate sweep.  Returns a dict of arrays."""
    Id = fet.transfer_characteristic(Vg, Vds=Vds)
    n = fet.carrier_density(Vg)
    gm_id = gm_over_id(Vg, Id)
    cnf = s_cnf(gm_id, nt_per_eV_cm3, f=f)
    cnfcmf = s_cnf_cmf(gm_id, nt_per_eV_cm3, f=f, alpha_c=alpha_c)
    signed = s_cnf_cmf_signed(gm_id, nt_per_eV_cm3, f=f, alpha_c=alpha_c)
    hooge = s_hooge(n, f=f, alpha_H=alpha_H)
    return dict(Vg=Vg, Id=Id, n=n, gm_id=gm_id, cnf=cnf,
                cnf_cmf=cnfcmf, cnf_cmf_signed=signed, hooge=hooge,
                total=cnfcmf + hooge, total_signed=signed + hooge)


def dip_depth(Vg, S, V_cnp):
    """Depth of the central dip in S(Vg), the M-shape discriminator.

    Positive  -> S has a local minimum at the CNP with a higher maximum on
                 BOTH sides: an M (or V) shape.
    Negative  -> the CNP is at or near the global maximum: a Lambda shape.

    Returned in decades, so the number is scale-free:
        log10( min(max_left, max_right) / S(CNP) )
    """
    i_cnp = int(np.argmin(np.abs(Vg - V_cnp)))
    left, right = S[:i_cnp + 1], S[i_cnp:]
    if len(left) < 2 or len(right) < 2:
        raise ValueError("CNP not interior to the sweep")
    return np.log10(min(left.max(), right.max()) / S[i_cnp])


def critical_nt(Vg, Vds=VDS, alpha_H=ALPHA_H, lo=1e15, hi=1e23, tol=1e-4,
                key="total"):
    """Trap density at which the central dip appears, by bisection on
    log10(N_T).

    NOT a root find.  dip_depth() is identically ZERO over the whole
    Lambda-shaped region (the maximum IS the CNP there, so the ratio is
    exactly 1), so the threshold is the boundary of a region rather than a
    sign change, and a sign-bracketing bisection cannot start -- the first
    draft of this module tried and raised.  Bisection is therefore on the
    PREDICATE dip_depth > 0.

    Returns (N_T*, iterations, dip depth just above, dip depth just below).
    """
    def has_dip(nt):
        return dip_depth(Vg, noise_sweep(Vg, nt, Vds=Vds,
                                         alpha_H=alpha_H)[key],
                         fet.V_dirac) > 0.0
    if has_dip(lo) or not has_dip(hi):
        raise ValueError("predicate does not change over [%g, %g]" % (lo, hi))
    a, b, it = np.log10(lo), np.log10(hi), 0
    while b - a > tol:
        m = 0.5 * (a + b)
        if has_dip(10.0**m):
            b = m
        else:
            a = m
        it += 1
    d_above = dip_depth(Vg, noise_sweep(Vg, 10.0**b, Vds=Vds,
                                        alpha_H=alpha_H)[key], fet.V_dirac)
    d_below = dip_depth(Vg, noise_sweep(Vg, 10.0**a, Vds=Vds,
                                        alpha_H=alpha_H)[key], fet.V_dirac)
    return 10.0**b, it, d_above, d_below


# -------------------------------------------------------------- validation
def validate():
    """V1-V5 and V7-V8 are against EXACTLY known values (an analytic
    derivative, an exact symmetry, exact algebraic identities, an exact
    power law, an exact area scaling, an exact algebraic zero).  V6 is a
    RANGE comparison against published measurements and V2b a convergence
    rate; both are weaker and are labelled so they cannot be read as exact.
    """
    out, ok = [], True

    def check(name, detail, passed, exact=True):
        nonlocal ok
        ok &= bool(passed)
        out.append((name, detail, bool(passed), exact))

    # V1  g_m/I_d on an analytic linear I(Vg) = A (Vg - V0) must return
    #     1/(Vg - V0) exactly (central differences are exact on a line).
    V0, A = 0.8, 3.0e-5
    Vg = np.linspace(2.0, 12.0, 401)
    gm_id = gm_over_id(Vg, A * (Vg - V0))
    err = np.max(np.abs(gm_id - 1.0 / (Vg - V0)) * (Vg - V0))
    check("V1 EXACT  g_m/I_d of a linear I(Vg) equals 1/(Vg-V0)",
          "max rel err %.3e" % err, err < 1e-12)

    # V2  EXACT ZERO BY SYMMETRY, in the Vds -> 0 LIMIT ITSELF (not
    #     approached numerically).  carrier_density() depends on
    #     (Vg - V_dirac) through its magnitude, so at V_ch = 0 the channel
    #     resistance is exactly even, I_d is exactly even, g_m is exactly
    #     odd, and (g_m/I_d)^2 -- hence the CNF noise -- must be EXACTLY
    #     even about the CNP.  A sign error anywhere in carrier_density or
    #     in the differentiation breaks this and nothing else here would.
    d = np.linspace(0.5, 25.0, 60)
    Vg_sym = np.concatenate([fet.V_dirac - d[::-1], fet.V_dirac + d])
    Id_lim = id_zero_vds_limit(Vg_sym)
    gm_id_lim = gm_over_id(Vg_sym, Id_lim)
    cnf_lim = s_cnf(gm_id_lim, NT_PAPER_CM)
    half = len(d)
    lo, hi = cnf_lim[:half][::-1], cnf_lim[half:]
    asym = float(np.max(np.abs(lo - hi) / np.maximum(lo, hi)))
    check("V2 EXACT  CNF noise is even about the CNP in the Vds->0 limit",
          "max asymmetry %.3e" % asym, asym < 1e-12)

    # V2b CONVERGENCE (not exact): at finite Vds the asymmetry is the drain
    #     bias itself, so it must fall LINEARLY with Vds.  This is what
    #     establishes that the 1e-6-scale asymmetry of
    #     transfer_characteristic() is the gradual-channel average and not a
    #     defect -- the first draft of this module asserted an exact
    #     symmetry on it and failed at 3.4e-6.
    def asym_at(vds):
        sw = noise_sweep(Vg_sym, NT_PAPER_CM, Vds=vds)
        a, b = sw["cnf"][:half][::-1], sw["cnf"][half:]
        return float(np.max(np.abs(a - b) / np.maximum(a, b)))
    a1, a2 = asym_at(2.0e-6), asym_at(1.0e-6)
    ratio = a1 / a2
    check("V2b RANGE  finite-Vds asymmetry halves when Vds halves "
          "(ratio 2 +-2%)",
          "a(2e-6)=%.3e a(1e-6)=%.3e ratio %.4f" % (a1, a2, ratio),
          abs(ratio - 2.0) < 0.04, exact=False)

    # V3  EXACT ALGEBRAIC IDENTITY.  alpha_c = 0 must reduce BOTH CMF
    #     variants to CNF with no residue at all.
    Vg = np.linspace(-30.0, 30.0, 241)
    sw = noise_sweep(Vg, NT_PAPER_CM)
    r_signed = float(np.max(np.abs(
        s_cnf_cmf_signed(sw["gm_id"], NT_PAPER_CM, alpha_c=0.0) - sw["cnf"])))
    r_abs = float(np.max(np.abs(
        s_cnf_cmf(sw["gm_id"], NT_PAPER_CM, alpha_c=0.0) - sw["cnf"])))
    check("V3 EXACT  both CMF variants at alpha_c=0 are bitwise CNF",
          "max |residue| signed %.1e, branch-corrected %.1e"
          % (r_signed, r_abs), r_signed == 0.0 and r_abs == 0.0)

    # V4  EXACT POWER LAW.  Every term is 1/f by construction.
    r_cnf = s_cnf(sw["gm_id"], NT_PAPER_CM, f=1.0) / s_cnf(
        sw["gm_id"], NT_PAPER_CM, f=100.0)
    r_h = s_hooge(sw["n"], f=1.0) / s_hooge(sw["n"], f=100.0)
    err = max(np.max(np.abs(r_cnf - 100.0)),
              np.max(np.abs(r_h - 100.0))) / 100.0
    check("V4 EXACT  both terms scale as 1/f (S(1Hz)/S(100Hz)=100)",
          "max rel err %.3e" % err, err < 1e-13)

    # V5  EXACT AREA SCALING.  alpha_H/(f n W L) at fixed SHEET density must
    #     divide by exactly 4 when both W and L double.
    n0 = 1.0e16
    r = s_hooge(n0) / s_hooge(n0, W=2 * fet.W, L=2 * fet.L)
    err = abs(r - 4.0) / 4.0
    check("V5 EXACT  Hooge noise scales as 1/(W L): ratio is 4",
          "ratio %.15f, rel err %.1e" % (r, err), err < 1e-14)

    # V6  RANGE comparison against published measurements (NOT exact).
    #     Balandin Sec. II: area-normalised S_I/I^2 at 10 Hz is
    #     1e-8 .. 1e-7 um^2/Hz for graphene FETs.
    Vg = np.linspace(-30.0, 30.0, 241)
    sw6 = noise_sweep(Vg, NT_SIO2_CM)
    area_um2 = fet.W * fet.L * 1e12
    norm = sw6["total"] * area_um2
    check("V6 RANGE  area-normalised S_I/I^2 overlaps the published "
          "1e-8..1e-7 um^2/Hz band",
          "model spans %.2e .. %.2e um^2/Hz" % (norm.min(), norm.max()),
          (norm.min() <= 1e-7) and (norm.max() >= 1e-8), exact=False)

    # V7  EXACT ALGEBRAIC ZERO of the verbatim signed form.  Feeding it
    #     g_m/I_d = -alpha_c mu C_ox must return EXACTLY zero: that is the
    #     artefact, pinned as an identity rather than described in prose.
    gz = np.array([-cmf_term()])
    z = float(s_cnf_cmf_signed(gz, NT_PAPER_CM)[0])
    check("V7 EXACT  signed CMF form is exactly 0 at g_m/I_d = "
          "-alpha_c*mu*C_ox", "S = %.1e 1/Hz (alpha_c mu C_ox = %.5e 1/V)"
          % (z, cmf_term()), z == 0.0)

    # V7b EXACT NON-ZERO of the branch-corrected form at the same point:
    #     (|-x| + x)^2 = 4x^2, so the ratio to S_Vfb must be exactly
    #     4 (alpha_c mu C_ox)^2.
    v = float(s_cnf_cmf(gz, NT_PAPER_CM)[0]) / s_vfb(NT_PAPER_CM)
    want = 4.0 * cmf_term()**2
    err = abs(v - want) / want
    check("V7b EXACT branch-corrected form at that point is exactly "
          "4(alpha_c mu C_ox)^2 S_Vfb",
          "ratio %.6e vs %.6e, rel err %.1e" % (v, want, err), err < 1e-14)

    # V8  EXACT AGREEMENT ON THE ELECTRON BRANCH.  The correction must be a
    #     no-op wherever g_m >= 0, or it is not a minimal change.
    pos = sw["gm_id"] >= 0.0
    dif = float(np.max(np.abs(
        s_cnf_cmf(sw["gm_id"][pos], NT_PAPER_CM)
        - s_cnf_cmf_signed(sw["gm_id"][pos], NT_PAPER_CM))))
    check("V8 EXACT  the branch correction is bitwise a no-op where g_m>=0",
          "%d points, max |difference| %.1e" % (int(pos.sum()), dif),
          dif == 0.0)
    return out, ok


# ------------------------------------------------------------------ figure
def make_figure(Vg, sw_paper, sw_sio2, nt_star, nt_scan, dips):
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

    for a, sw, title, nt in (
            (ax[0], sw_paper, "(a) N_T = 1.1e20 /eV/cm3 (fitted, "
             "electrolyte-gated)", NT_PAPER_CM),
            (ax[1], sw_sio2, "(b) N_T = 1e17 /eV/cm3 (mid Si/SiO2 range)",
             NT_SIO2_CM)):
        a.semilogy(sw["Vg"], sw["cnf"], lw=1.3, ls="--",
                   label=r"CNF  $(g_m/I_d)^2 S_{Vfb}$")
        a.semilogy(sw["Vg"], sw["cnf_cmf"], lw=1.6,
                   label=r"CNF+CMF  $(g_m/I_d+\alpha_c\mu C_{ox})^2 S_{Vfb}$")
        a.semilogy(sw["Vg"], sw["hooge"], lw=1.6, ls=":",
                   label=r"Hooge  $\alpha_H/(fN)$")
        a.semilogy(sw["Vg"], sw["cnf_cmf_signed"], lw=1.1, ls="-.",
                   color="C3",
                   label=r"CNF+CMF, verbatim signed form (artefact)")
        a.semilogy(sw["Vg"], sw["total"], lw=2.4, color="k",
                   label="total (branch-corrected + Hooge)")
        a.axvline(fet.V_dirac, color="0.6", lw=0.9)
        a.set_xlabel(r"$V_g$  (V)")
        a.set_ylabel(r"$S_I/I^2$  at 10 Hz  (1/Hz)")
        a.set_title(title, fontsize=10)
        a.grid(alpha=0.3, which="both")
        a.legend(fontsize=7.5, loc="lower left")
        a.set_ylim(1e-14, 3e-6)

    ax[2].semilogx(nt_scan, dips, lw=1.8, color="C3")
    ax[2].axhline(0.0, color="0.4", lw=1.0)
    ax[2].axvline(nt_star, color="C0", lw=1.2, ls="--",
                  label=r"$N_T^{*}=%.2e$" % nt_star)
    ax[2].set_xlabel(r"$N_T$  (1/eV/cm$^3$)")
    ax[2].set_ylabel("central dip depth  (decades)")
    ax[2].set_title("(c) M-shape onset: the dip exists only above $N_T^{*}$",
                    fontsize=10)
    ax[2].grid(alpha=0.3, which="both")
    ax[2].legend(fontsize=9)
    ax[2].annotate(r"$\Lambda$-shape" "\n(peak AT the CNP)", (3e15, 0.06),
                   fontsize=9, color="0.3")
    ax[2].annotate("M-shape\n(dip at the CNP)", (3e21, 0.35), fontsize=9,
                   color="0.3")

    fig.suptitle("Low-frequency 1/f noise of a back-gated single-layer "
                 "graphene FET: which mechanism the gate dependence reveals",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(results_path('ch4', "low_frequency_noise_gate_dependence.png"), dpi=150)
    plt.close(fig)


# -------------------------------------------------------------------- main
def main():
    print("=" * 78)
    print("LOW-FREQUENCY 1/f NOISE IN A SINGLE-LAYER GRAPHENE FET")
    print("CNF / CNF+CMF / Hooge, on this repository's own I_d(V_g)")
    print("=" * 78)
    print("Device (imported unchanged from graphene_fet_model.py):")
    print("  W = %.3g um, L = %.3g nm, area = %.4g um^2" %
          (fet.W * 1e6, fet.L * 1e9, fet.W * fet.L * 1e12))
    print("  C_ox = %.4e F/m^2 (t_ox = %.3g nm), mu = %.4g cm^2/Vs" %
          (fet.C_ox, fet.t_ox * 1e9, fet.mu * 1e4))
    print("  V_dirac = %.3g V, n_puddle = %.3e m^-2, Vds = %.3g V" %
          (fet.V_dirac, fet.n_puddle, VDS))
    print("Noise parameters:")
    print("  f = %.3g Hz, lambda = %.3g nm, alpha_H = %.3g, "
          "alpha_c = %.4g V s/C" %
          (F_MEAS, LAMBDA_TUNNEL * 1e9, ALPHA_H, ALPHA_C))
    print("  S_Vfb(N_T=1.1e20) = %.4e V^2/Hz" % s_vfb(NT_PAPER_CM))
    print("  S_Vfb(N_T=1e17)   = %.4e V^2/Hz" % s_vfb(NT_SIO2_CM))

    print("\n" + "-" * 78)
    print("SECTION 1.  VALIDATION")
    print("-" * 78)
    rows, ok = validate()
    n_exact = sum(1 for r in rows if r[3])
    for name, detail, passed, exact in rows:
        print("  [%s] %s\n          %s" %
              ("PASS" if passed else "FAIL", name, detail))
    print("  verdict: %s over %d checks (%d exact, %d range/convergence)"
          % ("PASS" if ok else "FAIL", len(rows), n_exact,
             len(rows) - n_exact))

    Vg = np.linspace(-30.0, 30.0, 601)
    sw_paper = noise_sweep(Vg, NT_PAPER_CM)
    sw_sio2 = noise_sweep(Vg, NT_SIO2_CM)

    print("\n" + "-" * 78)
    print("SECTION 2.  GATE DEPENDENCE AT THE TWO TRAP DENSITIES")
    print("-" * 78)
    for label, sw, nt in (("N_T = 1.1e20 (fitted, electrolyte-gated)",
                           sw_paper, NT_PAPER_CM),
                          ("N_T = 1.0e17 (mid Si/SiO2 range)",
                           sw_sio2, NT_SIO2_CM)):
        i_cnp = int(np.argmin(np.abs(Vg - fet.V_dirac)))
        i_max = int(np.argmax(sw["total"]))
        dd = dip_depth(Vg, sw["total"], fet.V_dirac)
        print("  %s" % label)
        print("    S_total at the CNP          : %.4e 1/Hz" %
              sw["total"][i_cnp])
        print("    S_total global max           : %.4e 1/Hz at Vg = %+.2f V"
              % (sw["total"][i_max], Vg[i_max]))
        print("    CNF/Hooge ratio at the max   : %.4g" %
              (sw["cnf_cmf"][i_max] / sw["hooge"][i_max]))
        print("    central dip depth            : %+.4f decades  -> %s" %
              (dd, "M-shape" if dd > 0 else "Lambda-shape"))
        print("    area-normalised S at the CNP : %.4e um^2/Hz"
              " (published band 1e-8..1e-7)" %
              (sw["total"][i_cnp] * fet.W * fet.L * 1e12))

    print("\n" + "-" * 78)
    print("SECTION 3.  THE IMPORTED SIGN CONVENTION NOTCHES THE HOLE BRANCH")
    print("-" * 78)
    x = cmf_term()
    print("  alpha_c*mu*C_ox = %.5e 1/V (bias-INDEPENDENT: oxide and"
          " mobility only)" % x)
    i_notch = int(np.argmin(np.abs(sw_paper["gm_id"] + x)))
    print("  g_m/I_d crosses -alpha_c*mu*C_ox at Vg = %+.3f V"
          " (g_m/I_d = %+.5e 1/V)" % (Vg[i_notch], sw_paper["gm_id"][i_notch]))
    i_mirror = int(np.argmin(np.abs(Vg - (2 * fet.V_dirac - Vg[i_notch]))))
    print("  at that bias the two forms disagree by a factor %.4g:"
          % (sw_paper["cnf_cmf"][i_notch] /
             sw_paper["cnf_cmf_signed"][i_notch]))
    print("    verbatim signed form   : %.4e 1/Hz" %
          sw_paper["cnf_cmf_signed"][i_notch])
    print("    branch-corrected form  : %.4e 1/Hz" %
          sw_paper["cnf_cmf"][i_notch])
    print("  its mirror bias Vg = %+.3f V (same |Vg-V_dirac|, electron side)"
          % Vg[i_mirror])
    print("    both forms there       : %.4e 1/Hz (V8: bitwise equal)" %
          sw_paper["cnf_cmf"][i_mirror])
    print("  So the verbatim form makes two mirror-image bias points differ")
    print("  by %.4g x in trap-limited 1/f noise, with nothing physical"
          % (sw_paper["cnf_cmf_signed"][i_mirror] /
             sw_paper["cnf_cmf_signed"][i_notch]))
    print("  distinguishing them.  Totals (CMF+Hooge), which is what a")
    print("  measurement would see: signed %.4e vs corrected %.4e 1/Hz."
          % (sw_paper["total_signed"][i_notch], sw_paper["total"][i_notch]))
    print("  The Hooge floor hides part of it -- the notch is a factor %.4g"
          % (sw_paper["total"][i_notch] / sw_paper["total_signed"][i_notch]))
    print("  in the TOTAL, which is why it would be missed by eye on a log")
    print("  plot but not by a fit.  Checks V7/V7b/V8 pin all three claims.")

    print("\n" + "-" * 78)
    print("SECTION 4.  THE M-SHAPE ONSET IS A TRAP-DENSITY THRESHOLD")
    print("-" * 78)
    nt_star, iters, d_above, d_below = critical_nt(Vg)
    print("  bisection on the PREDICATE dip>0 (dip_depth is identically 0")
    print("  over the whole Lambda region, so this is a region boundary and")
    print("  not a root), tol = 1e-4 decades, %d iterations" % iters)
    print("  N_T* = %.4e 1/(eV cm^3)" % nt_star)
    print("  dip depth just above N_T*: %+.3e decades" % d_above)
    print("  dip depth just below N_T*: %+.3e decades (exactly 0 => no dip)"
          % d_below)
    print("  Below N_T* the predicted gate dependence is Lambda-shaped (the")
    print("  noise PEAKS at the charge-neutrality point); above it the CNF")
    print("  term wins on the flanks and opens the central dip, giving the")
    print("  M-shape the experimental literature reports.")
    print("  N_T* / 1e18 (top of the conventional Si/SiO2 interval) = %.3g"
          % (nt_star / 1e18))
    print("  N_T* / 1.1e20 (the fitted electrolyte-gated value)     = %.3g"
          % (nt_star / NT_PAPER_CM))
    print("  => N_T* lies INSIDE the conventional 1e16..1e18 interval, at")
    print("     %.0f%% of its top end." % (100 * nt_star / 1e18))
    print("  SENSITIVITY TO THE HOOGE FLOOR.  N_T* is where CNF overtakes")
    print("  the Hooge term at the flanks, so it moves with alpha_H.  The")
    print("  three alpha_H values the same Table 1 lists for other devices:")
    for ah, src in ((3.0e-4, "its Ref. 33, lowest"),
                    (ALPHA_H, "this model, L=13um row"),
                    (7.0e-3, "its Ref. 32, highest")):
        nt_a, _, _, _ = critical_nt(Vg, alpha_H=ah)
        print("    alpha_H = %-8.1e -> N_T* = %.3e  (%s)" % (ah, nt_a, src))
    nt_signed, _, _, _ = critical_nt(Vg, key="total_signed")
    print("  DOES THE SECTION-3 ARTEFACT MOVE THIS THRESHOLD?  Recomputing")
    print("  N_T* on the verbatim signed total instead of the corrected one:")
    print("    N_T*(branch-corrected) = %.4e" % nt_star)
    print("    N_T*(verbatim signed)  = %.4e   ratio %.4f"
          % (nt_signed, nt_signed / nt_star))
    print("  The notch sits on one flank, and dip depth takes the MINIMUM of")
    print("  the two flank maxima, so the notched flank is the one that")
    print("  decides -- the artefact %s the threshold."
          % ("does not move" if abs(nt_signed / nt_star - 1) < 1e-3
             else "DOES move"))
    print("  Over that whole factor-23 span in alpha_H, N_T* stays below")
    print("  the fitted 1.1e20 and spans the top of the conventional")
    print("  interval rather than sitting far above it.")

    nt_scan = np.logspace(15.0, 23.0, 65)
    dips = [dip_depth(Vg, noise_sweep(Vg, nt)["total"], fet.V_dirac)
            for nt in nt_scan]

    print("\n" + "-" * 78)
    print("SECTION 5.  WHAT THIS DOES AND DOES NOT ESTABLISH")
    print("-" * 78)
    print("  DOES, AND IT IS THE OPPOSITE OF WHAT THIS MODULE WAS DRAFTED")
    print("  EXPECTING.  The draft prose asserted that the M-shape onset")
    print("  would sit far ABOVE the conventional SiO2 trap-density range,")
    print("  and therefore that an M-shape on a SiO2-gated GFET could not")
    print("  be homogeneous CNF.  The computed threshold refutes that:")
    print("  N_T* = %.3e 1/(eV cm^3) is INSIDE the 1e16..1e18 interval."
          % nt_star)
    print("  So on this device, homogeneous charge-number fluctuation at a")
    print("  perfectly ordinary oxide trap density already produces the")
    print("  M-shape, and the M-shape is NOT by itself evidence of the")
    print("  electron/hole-puddle inhomogeneity that Balandin's review")
    print("  attributes it to.  The shape discriminates mechanism only")
    print("  together with an independent trap-density measurement.")
    print("  DOES NOT: the puddle mechanism is not in this model at all --")
    print("  n_puddle here REGULARISES the density, it does not fluctuate.")
    print("  This result weakens the inference FROM shape TO mechanism; it")
    print("  says nothing about whether puddles also contribute.")
    print("  DOES NOT: contact noise (negligible near the CNP per")
    print("  arXiv:1806.06589, dominant at high |Vg|) is absent, so the")
    print("  flanks are a lower bound and the dip depth an upper one.")
    print("  DOES NOT: the CMF form and alpha_c are SUMMARY-ONLY provenance")
    print("  (that paper's supplement was not read verbatim).  V3/V7b pin")
    print("  the form's algebra, not the form.")
    print("  DOES NOT: N_T* is a threshold in THIS geometry at THIS mu and")
    print("  t_ox.  S_Vfb carries 1/(W L C_ox^2), so a thinner oxide or a")
    print("  larger area moves it; no scaling study is done here.")

    make_figure(Vg, sw_paper, sw_sio2, nt_star, nt_scan, dips)
    print("\n  figure written: low_frequency_noise_gate_dependence.png")
    print("=" * 78)
    print("VERDICT: %s over V1..V8" % ("PASS" if ok else "FAIL"))
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
