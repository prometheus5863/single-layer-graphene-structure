"""
graphene_branch_asymmetric_mobility.py

ELECTRON-HOLE MOBILITY ASYMMETRY IN THE GFET, AND HOW TO TELL IT APART FROM
CONTACT ASYMMETRY.

This closes the follow-on item created 2026-10-09 ("ASYMMETRIC ELECTRON/HOLE
MOBILITY"), which has been carried untouched for twenty sessions, and it is
the item `graphene_fet_signed_carrier_model.py` deferred in writing:

    "the electron-hole mobility asymmetry that would separate them is a real
     effect in contact-doped devices and is NOT introduced here, because this
     module's contract is to add a branch without moving a number.  Asymmetric
     mobility is logged as a follow-on item instead."
            -- sheet_conductivity_two_carrier.__doc__

That module supplied the branch (n_e, n_h signed decomposition) and left the
mobility hooks `mu_e`, `mu_h` deliberately defaulted to the single shipped
`gfet.mu`.  This module is what the hooks were for.

------------------------------------------------------------------------------
1. WHY THE TWO BRANCHES ARE NOT EQUIVALENT (the physics, and the conflict)
------------------------------------------------------------------------------
Within the Born approximation, Coulomb scattering of Dirac fermions from a
charged impurity is electron-hole SYMMETRIC: the Born cross-section depends on
the square of the impurity charge, so the sign drops out.  The asymmetry seen
in experiment therefore has to come from somewhere outside leading-order Born,
and the literature offers two mechanisms that are still not reconciled:

  (A) CHANNEL ORIGIN.  Beyond Born (Novikov's relativistic treatment), the
      sign of the impurity potential matters: carriers of charge opposite to
      the impurity are scattered more strongly.  Zion, Butenko, Sharoni,
      Kaveh and Shlimak, "Charge carrier transport asymmetry in monolayer
      graphene", arXiv:1801.03565, measure a Hall bar with non-invasive
      voltage probes 30 um from the contacts, find mu_h < mu_e, and attribute
      it to negatively charged acceptors in the SiO2 (about 1.5e12 cm^-2,
      with a puddle density 0.25e12 cm^-2 and potential amplitude ~58 meV).
      They report no numerical mu_e/mu_h; they report that the asymmetry
      SHRINKS with increasing |n| as the long-range potential is screened and
      short-range neutral-defect scattering takes over (hole crossover at
      p ~ 3e12 cm^-2).

  (B) CONTACT ORIGIN.  Huard, Stander, Sulpizio and Goldhaber-Gordon,
      "Evidence of the role of contacts on the observed electron-hole
      asymmetry in graphene", arXiv:0804.2040, show the asymmetry appears with
      INVASIVE probes and vanishes with external ones: rho_e/rho_h =
      1 +/- 0.03 measured away from the metal, against ~1.20 inferred when the
      contacts are in the current path.  The metal pins the graphene beneath it
      p-type (Pd ~5.1 eV vs graphene ~4.5 eV), so the junction is p-n for an
      n-type channel and p-p' for a p-type one, and the two junction
      resistances differ.  Their odd resistance collapses when scaled by the
      INTERFACE WIDTH w, not by the channel length, and their specific contact
      resistance is lambda = 135 Ohm.um -- which sits at the bottom of the
      110-500 Ohm.um band this repository already calibrates `Rc_total` on.

These two papers disagree about sign as well as about mechanism: (A) finds the
holes worse, (B) finds the electrons more resistive for p-doping contacts.
This module does not adjudicate the physics.  It does something cheaper and
decidable: it shows that the two mechanisms have DIFFERENT CHANNEL-LENGTH
SIGNATURES in this repository's own model, and that the signatures separate
exactly rather than approximately (Section 5).

------------------------------------------------------------------------------
2. THE MODEL, AND THE ONE CLOSED FORM IT TURNS ON
------------------------------------------------------------------------------
Keep the electrostatics untouched -- mobility is a transport parameter and does
not enter the charge relation -- and put the asymmetry only in the two-carrier
Drude sum:

    sigma(n_net) = e*(mu_e*n_e + mu_h*n_h),   n_e, n_h = (S +/- n_net)/2,
                   S = sqrt(n_net^2 + n_p^2)                             (A1)

Substituting the split and writing

    mu_bar = (mu_e + mu_h)/2,   dmu = mu_e - mu_h,   r = dmu/(2*mu_bar)  (A2)

gives the form this module works in:

    sigma(n_net) = e*(mu_bar*S + (dmu/2)*n_net)                          (A3)

S is EVEN in n_net and n_net is ODD in dV, so (A3) splits sigma into an even
part the shipped model already had and an odd part that is new and is linear
in the mobility difference.  Minimising (A3):

    d(sigma)/d(n_net) = e*(mu_bar*n_net/S + dmu/2) = 0
      =>  n_net/S = -r
      =>  n_net* = -n_p * r / sqrt(1 - r^2)                              (A4)
      =>  sigma_min = e*n_p*mu_bar*sqrt(1 - r^2) = e*n_p*sqrt(mu_e*mu_h) (A5)

(A5) is the result worth having.  Because 1 - r^2 = 4*mu_e*mu_h/(mu_e+mu_h)^2
identically, the minimum sheet conductivity of a disordered graphene channel is
set by the GEOMETRIC MEAN of the two mobilities, not the arithmetic mean --
exactly, with no expansion.

  SCOPE NARROWED 2026-10-30.  "Exactly, with no expansion" is correct and
  stands, but its premise is a CONSTANT pair, which is the premise of this
  module and is no longer the premise of the repository.  Once mu = mu(S) the
  derivation above loses that premise -- mu_bar and dmu can no longer be
  pulled through d/dn -- and D10 of graphene_density_dependent_mobility.py
  measures what survives: e*n_p*sqrt(mu_e(S*)*mu_h(S*)) is accurate to
  1.17e-3 at worst, with the error SECOND ORDER in n*/S*
  (dev/(n*/S*)^2 constant to a factor 1.044 across a 48.6x range in
  (n*/S*)^2).  The geometric-mean reading therefore survives as a
  second-order-accurate statement rather than an exact one.  No number below
  moves: they are the r = const case, and for that case (A5) is exact.

Three consequences, all checked below:

  (i)  CORRECTION OF RECORD.  `graphene_fet_signed_carrier_model.sigma_min_ch4`
       returns `e*mu*n_puddle` and its docstring calls that "EXACTLY".  It is
       exact only on the mu_e == mu_h diagonal; (A5) is the general statement
       and reduces to it.  The committed number does not move (the shipped
       model is on the diagonal), so this is a scope correction, not a numeric
       one -- see Section 6 and the annotation added in that file.
  (ii) THE FLOOR IS SECOND-ORDER, THE DIRAC POINT IS FIRST-ORDER.  Since
       sqrt(1-r^2) = 1 - r^2/2 + O(r^4) while (A4) is linear in r, a mobility
       asymmetry moves the apparent Dirac point at FIRST order in r and the
       conductivity floor only at SECOND order.  That asymmetry of orders is
       the reason an experiment can show a large Dirac-point offset and a
       floor that still looks universal, and it is why the Chapter 3/Chapter 4
       floor adjudication is robust to an asymmetry that would nonetheless
       wreck a doping extraction.
  (iii) A MOBILITY ASYMMETRY IMPERSONATES A DOPING OFFSET.  n_net* != 0 means
       the resistance maximum no longer sits at V_g = V_dirac.  Section 4
       converts (A4) into the apparent gate offset through this repository's
       own series-factor charge relation and compares it with the V_dirac =
       0.8 V the model attributes to built-in doping.

------------------------------------------------------------------------------
3. WHAT IS CHECKED, AND WHICH CHECKS ARE EXACT
------------------------------------------------------------------------------
Five of the checks below are exact-value checks, not plausible-range checks,
which is the standing requirement in this repository since 2026-09-17:

  C1  mu_e == mu_h reproduces the shipped sheet conductivity over a sweep, to
      a MEASURED ulp distance (not asserted bitwise: (A1) and the shipped
      `abs(n)*e*mu` differ in association order, so the honest claim is a
      bound, and the bound is printed).
  C2  mu_e == mu_h makes sigma(+dV) - sigma(-dV) EXACTLY 0.0, bitwise, in the
      n_net_of_dV parameterisation -- the parameterisation caveat is the one
      `graphene_fet_signed_carrier_model.n_net_of_dV` exists to enforce.
  C3  n_e*n_h == n_p^2/4 exactly, for every mobility pair, because mobility is
      absent from the electrostatics.  A mobility that leaked into the charge
      relation would show up here as a nonzero residual.
  C4  d(sigma)/d(n_net) vanishes at the analytic n_net* of (A4) to machine
      precision, evaluated on the analytic derivative -- the zero that must be
      zero.
  C5  Exchanging mu_e <-> mu_h maps n_net* -> -n_net* bitwise and leaves
      sigma_min bitwise unchanged, since (A5) is symmetric in the two
      mobilities.
  C6  AM-GM: sqrt(mu_e*mu_h) <= mu_bar with equality iff mu_e == mu_h, so the
      asymmetric floor is never above the symmetric one.
  C7  The Section 5 extraction recovers PLANTED mu_e, mu_h, Rc_n, Rc_p from a
      length sweep to ~1e-12 relative.  This is the strong check: the odd part
      of the total resistance is exactly linear in L, so a planted truth must
      come back out of a straight-line fit.

------------------------------------------------------------------------------
4. WHAT THIS MODULE DOES NOT CLAIM
------------------------------------------------------------------------------
mu_e and mu_h are constants here.  Zion et al.'s central observation is that
the asymmetry is DENSITY-DEPENDENT (it shrinks as |n| grows), which a constant
pair cannot represent; a density-dependent mobility pair is logged as a
follow-on item rather than guessed at.

    SUPERSEDED 2026-10-30, and the numbers it affects are annotated where they
    are made rather than changed.  graphene_density_dependent_mobility.py
    closes that item with 1/mu_b(n) = 1/mu_L,b + |n|/(n_c*mu_L_bar) and
    derives, exactly,

        r(n) = r0 / (1 + (1 - r0^2) * |n| / n_c)

    so the constant pair here is the n_c -> infinity limit of that model, and
    the deviation from it falls as 1/n_c (measured: asymptotic log-log slope
    -0.99996865, with the departure from -1 accounted for by the subleading
    term to 0.01 %).  What makes the simplification harmless FOR THIS DEVICE
    rather than merely convenient is a device fact and not a modelling one:
    the shipped 90 nm SiO2 back gate needs 12.53 V of overdrive to reach
    n_c = 3e12 cm^-2, and every sweep in this repository runs inside +/-5 V.
    The same simplification is NOT harmless for a 10 nm HfO2 top-gated
    device, which crosses the same density at 0.217 V.

The contact asymmetry is likewise a
constant resistance difference, not Huard et al.'s n^(-1/2) junction law
(r_nc(n) > 2.06*(h/4e^2)*n^(-1/2)); Section 5's discriminator only needs the
contact term to be L-INDEPENDENT, which the n^(-1/2) law also is, so the
separation argument survives the simplification even though the fitted
intercept would then be bias-dependent.
"""

import os
import sys

import numpy as np

import graphene_fet_model as gfet
import graphene_fet_signed_carrier_model as sgn

E_CHARGE = 1.602176634e-19

passed = 0
failed = 0


def check(label, ok, detail=''):
    global passed, failed
    if ok:
        passed += 1
        print('  [PASS] %s' % label)
    else:
        failed += 1
        print('  [FAIL] %s' % label)
    if detail:
        print('         %s' % detail)


def ulp_distance(a, b):
    """Integer-representable distance between two finite float arrays."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ia = a.view(np.int64).copy()
    ib = b.view(np.int64).copy()
    ia[ia < 0] = np.int64(np.iinfo(np.int64).min) - ia[ia < 0]
    ib[ib < 0] = np.int64(np.iinfo(np.int64).min) - ib[ib < 0]
    return np.abs(ia - ib)


# ---------------------------------------------------------------------------
# Section 1: the asymmetric two-carrier conductivity
# ---------------------------------------------------------------------------

def mu_pair(ratio, mu_bar=None):
    """
    (mu_e, mu_h) with mu_e/mu_h == ratio and arithmetic mean mu_bar.

    The mean is held fixed rather than mu_h, so that a sweep over `ratio`
    varies ONLY the asymmetry and not the overall mobility scale; otherwise
    every effect below would be contaminated by a change in mu_bar.
    """
    if mu_bar is None:
        mu_bar = gfet.mu
    mu_h = 2.0 * mu_bar / (1.0 + ratio)
    return ratio * mu_h, mu_h


def r_of(mu_e, mu_h):
    """The asymmetry parameter r = (mu_e - mu_h)/(mu_e + mu_h) of (A2)."""
    return (mu_e - mu_h) / (mu_e + mu_h)


def sigma_asym(n_net_val, mu_e, mu_h, n_puddle=None):
    """
    Sheet conductivity [S/sq] from the species split, eq. (A1).  Routed through
    the signed-carrier module's own split and two-carrier sum so that this
    module cannot drift from the decomposition it builds on.
    """
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    n_e, n_h = sgn.split_species(n_net_val, n_puddle=n_puddle)
    return sgn.sheet_conductivity_two_carrier(n_e, n_h, mu_e=mu_e, mu_h=mu_h)


def sigma_asym_closed(n_net_val, mu_e, mu_h, n_puddle=None):
    """The same conductivity in the even/odd form (A3)."""
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    n_net_val = np.asarray(n_net_val, dtype=float)
    mu_bar = 0.5 * (mu_e + mu_h)
    dmu = mu_e - mu_h
    S = np.sqrt(n_net_val ** 2 + n_puddle ** 2)
    return E_CHARGE * (mu_bar * S + 0.5 * dmu * n_net_val)


def dsigma_dn(n_net_val, mu_e, mu_h, n_puddle=None):
    """Analytic d(sigma)/d(n_net) of (A3), used for the C4 zero."""
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    n_net_val = np.asarray(n_net_val, dtype=float)
    mu_bar = 0.5 * (mu_e + mu_h)
    dmu = mu_e - mu_h
    S = np.sqrt(n_net_val ** 2 + n_puddle ** 2)
    return E_CHARGE * (mu_bar * n_net_val / S + 0.5 * dmu)


def n_net_at_minimum(mu_e, mu_h, n_puddle=None):
    """n_net* of (A4): the net density at which sigma is minimal."""
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    r = r_of(mu_e, mu_h)
    return -n_puddle * r / np.sqrt(1.0 - r * r)


def sigma_min_closed(mu_e, mu_h, n_puddle=None):
    """sigma_min of (A5): e*n_p*sqrt(mu_e*mu_h), the geometric mean."""
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    return E_CHARGE * n_puddle * np.sqrt(mu_e * mu_h)


# ---------------------------------------------------------------------------
# Section 2: from n_net* back to a gate voltage
# ---------------------------------------------------------------------------

def dV_of_n_net(target, lo=-5.0, hi=5.0, iters=200):
    """
    Invert this repository's own series-factor charge relation n_net(dV) for
    the overdrive dV that produces `target`.

    Bisection rather than a library root finder, deliberately: the bracket is
    known from the sign of `target` and bisection's convergence does not depend
    on the derivative of the quantum-capacitance factor, whose thermal rounding
    is the subject of three open items in this repository.  Monotonicity of
    n_net(dV) is checked in Section 2 of main() before this is trusted.
    """
    target = float(target)
    f = lambda x: float(sgn.n_net_of_dV(x)) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0.0:
        raise ValueError('bracket does not straddle the target')
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if flo * fmid <= 0.0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def apparent_dirac_offset(mu_e, mu_h, n_puddle=None):
    """
    The gate-voltage offset, relative to V_dirac, at which the channel
    resistance peaks when the two branches have different mobilities.

    Positive means the measured Dirac point sits ABOVE the built-in V_dirac.
    """
    return dV_of_n_net(n_net_at_minimum(mu_e, mu_h, n_puddle=n_puddle))


# ---------------------------------------------------------------------------
# Section 3: total two-terminal resistance with branch-dependent contacts
# ---------------------------------------------------------------------------

def R_total(dV, mu_e, mu_h, L, Rc_n, Rc_p, W=None):
    """
    Two-terminal resistance [Ohm] at local overdrive dV: channel sheet
    resistance times L/W, plus the contact term of the branch the channel is
    in.  The contact term is selected by sign(dV) -- i.e. by the branch -- which
    is the whole content of mechanism (B).
    """
    if W is None:
        W = gfet.W
    dV = np.asarray(dV, dtype=float)
    nn = sgn.n_net_of_dV(dV)
    sigma = sigma_asym_closed(nn, mu_e, mu_h)
    Rc = np.where(nn >= 0.0, Rc_n, Rc_p)
    return (L / W) / sigma + Rc


def R_odd_vs_length(dV, mu_e, mu_h, lengths, Rc_n, Rc_p, W=None):
    """
    The odd part of the resistance, R(-dV) - R(+dV), over a length sweep.

    Exactly linear in L:
        R_odd(L) = (L/W)*(1/sigma(-dV) - 1/sigma(+dV)) + (Rc_p - Rc_n)
    so the SLOPE carries the channel-mobility asymmetry and the INTERCEPT
    carries the contact asymmetry.  This is the discriminator: Huard et al.'s
    contact term collapses on interface width and is L-independent, while a
    channel asymmetry must grow with L.
    """
    lengths = np.asarray(lengths, dtype=float)
    return np.array([
        float(R_total(-abs(dV), mu_e, mu_h, L, Rc_n, Rc_p, W=W)
              - R_total(+abs(dV), mu_e, mu_h, L, Rc_n, Rc_p, W=W))
        for L in lengths
    ])


def extract_from_length_sweep(dV, lengths, R_odd, W=None):
    """
    Recover (channel-asymmetry slope, contact asymmetry) from a length sweep by
    exact two-point line solve on the first and last length, then report the
    worst residual over the interior points.  A two-point solve is used instead
    of a least-squares fit so that the recovery is an identity when the
    underlying relation is linear, which is what C7 tests.
    """
    if W is None:
        W = gfet.W
    lengths = np.asarray(lengths, dtype=float)
    R_odd = np.asarray(R_odd, dtype=float)
    slope = (R_odd[-1] - R_odd[0]) / (lengths[-1] - lengths[0])
    intercept = R_odd[0] - slope * lengths[0]
    resid = R_odd - (slope * lengths + intercept)
    return slope, intercept, float(np.max(np.abs(resid)))


# ---------------------------------------------------------------------------
# Section 4: figure
# ---------------------------------------------------------------------------

def make_figure(path='branch_asymmetric_mobility.png'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ratios = [1.0, 1.1, 1.2, 1.3, 1.5]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))

    # (a) sheet conductivity vs overdrive, with the closed-form minima marked
    ax = axes[0][0]
    dV = np.linspace(-3.0, 3.0, 1201)
    for ratio in ratios:
        mu_e, mu_h = mu_pair(ratio)
        s = sigma_asym_closed(sgn.n_net_of_dV(dV), mu_e, mu_h)
        ax.plot(dV, s * 1e6, label=r'$\mu_e/\mu_h$ = %.1f' % ratio)
        dv_star = apparent_dirac_offset(mu_e, mu_h)
        ax.plot([dv_star], [sigma_min_closed(mu_e, mu_h) * 1e6], 'k.', ms=7)
    ax.set_xlabel(r'$V_g - V_{dirac}$  (V)')
    ax.set_ylabel(r'$\sigma_{sheet}$  ($\mu$S/sq)')
    ax.set_title('(a) Asymmetry tilts the V, and moves its vertex off zero\n'
                 r'dots: closed form $\sigma_{min}=e n_p\sqrt{\mu_e\mu_h}$, '
                 r'$n^*=-n_p r/\sqrt{1-r^2}$')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) the two orders: Dirac shift linear in r, floor quadratic in r
    ax = axes[0][1]
    rr = np.linspace(0.0, 0.25, 81)
    shift = []
    floor_rel = []
    for r in rr:
        mu_e = gfet.mu * (1.0 + r)
        mu_h = gfet.mu * (1.0 - r)
        shift.append(apparent_dirac_offset(mu_e, mu_h))
        floor_rel.append(sigma_min_closed(mu_e, mu_h)
                         / sigma_min_closed(gfet.mu, gfet.mu) - 1.0)
    ax.plot(rr, np.abs(shift), 'C0-', lw=2,
            label=r'$|\Delta V_{dirac}^{app}|$  (V, left)')
    ax.axhline(gfet.V_dirac, color='C3', ls='--', lw=1.2,
               label=r'built-in $V_{dirac}$ = %.2f V' % gfet.V_dirac)
    ax.set_xlabel(r'asymmetry  $r=(\mu_e-\mu_h)/(\mu_e+\mu_h)$')
    ax.set_ylabel('apparent Dirac-point offset  (V)')
    ax2 = ax.twinx()
    ax2.plot(rr, 100.0 * np.array(floor_rel), 'C2-', lw=2)
    ax2.plot(rr, -100.0 * 0.5 * rr ** 2, 'C2:', lw=1.4)
    ax2.set_ylabel(r'$\sigma_{min}$ change (%), dotted: $-r^2/2$', color='C2')
    ax.set_title('(b) First order in the Dirac point, second order in the floor')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)

    # (c) transfer characteristics, equal-mean mobility
    ax = axes[1][0]
    Vg = np.linspace(-2.0, 4.0, 601)
    for ratio in ratios:
        mu_e, mu_h = mu_pair(ratio)
        R = R_total(Vg - gfet.V_dirac, mu_e, mu_h, gfet.L,
                    gfet.Rc_total / 2.0, gfet.Rc_total / 2.0)
        ax.plot(Vg, 1e6 * 0.05 / R, label=r'$\mu_e/\mu_h$ = %.1f' % ratio)
    ax.axvline(gfet.V_dirac, color='k', ls=':', lw=1)
    ax.set_xlabel(r'$V_g$  (V)')
    ax.set_ylabel(r'$I_d$  ($\mu$A/$\mu$m)  at $V_{ds}$ = 50 mV')
    ax.set_title('(c) Transfer characteristic: the branch slopes split and the\n'
                 'current minimum walks away from the dotted built-in $V_{dirac}$')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (d) the length discriminator
    ax = axes[1][1]
    lengths = np.linspace(100e-9, 5e-6, 60)
    dV_probe = 2.5
    mu_e, mu_h = mu_pair(1.3)
    Rc_half = gfet.Rc_total / 2.0
    cases = [
        ('channel only  ($\\mu_e/\\mu_h$=1.3, $R_c$ equal)',
         (mu_e, mu_h, Rc_half, Rc_half), 'C0'),
        ('contact only  ($\\mu$ equal, $\\Delta R_c$=120 $\\Omega$)',
         (gfet.mu, gfet.mu, Rc_half, Rc_half + 120.0), 'C1'),
        ('both', (mu_e, mu_h, Rc_half, Rc_half + 120.0), 'C2'),
    ]
    for label, (me, mh, rn, rp), c in cases:
        R_odd = R_odd_vs_length(dV_probe, me, mh, lengths, rn, rp)
        slope, intercept, _ = extract_from_length_sweep(dV_probe, lengths, R_odd)
        ax.plot(lengths * 1e6, R_odd, c + '-', lw=2, label=label)
        ax.plot(lengths * 1e6, slope * lengths + intercept, 'k:', lw=0.9)
    ax.axhline(0.0, color='k', lw=0.8)
    ax.set_xlabel(r'channel length $L$  ($\mu$m)')
    ax.set_ylabel(r'$R(-\Delta V) - R(+\Delta V)$  ($\Omega$)  at $|\Delta V|$=2.5 V')
    ax.set_title('(d) The discriminator: slope is the channel, intercept is the\n'
                 'contact.  Exactly linear, so a planted truth is recoverable.')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle('Electron-hole mobility asymmetry in the GFET channel '
                 '(closes the 2026-10-09 follow-on item)', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('ELECTRON-HOLE MOBILITY ASYMMETRY IN THE GFET CHANNEL')
    print('closes the follow-on item created 2026-10-09 (twenty sessions carried)')
    print('=' * 78)
    print('Shipped parameters: mu = %.4g m^2/Vs, n_puddle = %.4g m^-2, '
          'V_dirac = %.4g V' % (gfet.mu, gfet.n_puddle, gfet.V_dirac))
    print('                    C_ox = %.6g F/m^2, L = %.4g m, W = %.4g m'
          % (gfet.C_ox, gfet.L, gfet.W))

    # -- Section 1 ----------------------------------------------------------
    print('\nSECTION 1.  The symmetric limit must reproduce the shipped model')
    Vg = np.linspace(-2.0, 4.0, 601)
    shipped = gfet.sheet_conductivity(gfet.carrier_density(Vg))
    mine = sigma_asym(sgn.n_net(Vg), gfet.mu, gfet.mu)
    mine_closed = sigma_asym_closed(sgn.n_net(Vg), gfet.mu, gfet.mu)
    u_split = int(np.max(ulp_distance(shipped, mine)))
    u_closed = int(np.max(ulp_distance(shipped, mine_closed)))
    rel_split = float(np.max(np.abs(mine / shipped - 1.0)))
    rel_closed = float(np.max(np.abs(mine_closed / shipped - 1.0)))

    # The admissible bound is COUNTED, not guessed.  The shipped expression
    # abs(n)*e*mu rounds twice after S.  The species-split form rounds five
    # times after S (the two halves, the two mobility products, their sum, and
    # the charge factor), and the even/odd form (A3) rounds four times (mu_bar*S,
    # the odd term, the sum, the charge factor).  Seven roundings of relative
    # size eps/2 each bound the disagreement at 7*eps/2 relative, i.e. ~8 ulp
    # for mantissas in [1, 2).  A FIRST attempt at this check used 1 ulp, which
    # is not a bound derived from anything; it failed at 3 ulp and the threshold
    # is corrected here rather than the measurement being talked around.
    ULP_BUDGET = 8
    EPS = float(np.finfo(np.float64).eps)
    check('C1 species-split form matches shipped sigma within the counted '
          '%d-ulp budget' % ULP_BUDGET,
          u_split <= ULP_BUDGET,
          'max ulp distance %d, max relative deviation %.3e (= %.2f eps)'
          % (u_split, rel_split, rel_split / EPS))
    check('C1 even/odd closed form (A3) matches shipped sigma within the '
          'counted %d-ulp budget' % ULP_BUDGET,
          u_closed <= ULP_BUDGET,
          'max ulp distance %d, max relative deviation %.3e (= %.2f eps)'
          % (u_closed, rel_closed, rel_closed / EPS))
    print('         Exact bitwise equality is NOT claimed: the forms differ in')
    print('         association order.  The claim is the counted bound above.')

    # A bound alone cannot tell a rounding floor from a small MODEL difference,
    # because both can hide under 8 ulp on one grid.  The discriminator is
    # refinement: a rounding floor does not care how finely the sweep is
    # sampled, while a model difference evaluated at new bias points would
    # generally move.  Refine 10x and require the deviation not to grow.
    Vg_fine = np.linspace(-2.0, 4.0, 6001)
    rel_fine = float(np.max(np.abs(
        sigma_asym(sgn.n_net(Vg_fine), gfet.mu, gfet.mu)
        / gfet.sheet_conductivity(gfet.carrier_density(Vg_fine)) - 1.0)))
    check('C1 deviation is a rounding floor, not a model difference '
          '(10x refinement does not grow it)',
          rel_fine <= max(rel_split, 8.0 * EPS),
          '601-point sweep %.3e, 6001-point sweep %.3e, budget %.3e'
          % (rel_split, rel_fine, 8.0 * EPS))

    # -- Section 2 ----------------------------------------------------------
    print('\nSECTION 2.  The zeros that must be exactly zero')
    dV = np.linspace(0.05, 3.0, 400)
    s_plus = sigma_asym_closed(sgn.n_net_of_dV(+dV), gfet.mu, gfet.mu)
    s_minus = sigma_asym_closed(sgn.n_net_of_dV(-dV), gfet.mu, gfet.mu)
    worst = float(np.max(np.abs(s_plus - s_minus)))
    check('C2 mu_e == mu_h gives sigma(+dV) - sigma(-dV) == 0.0 bitwise',
          worst == 0.0, 'max |difference| = %r (must be exactly 0.0)' % worst)

    nn = sgn.n_net_of_dV(np.linspace(-3.0, 3.0, 401))
    for ratio in (1.0, 1.3, 2.0):
        mu_e, mu_h = mu_pair(ratio)
        n_e, n_h = sgn.split_species(nn)
        prod_resid = float(np.max(np.abs(n_e * n_h - gfet.n_puddle ** 2 / 4.0)))
        rel = prod_resid / (gfet.n_puddle ** 2 / 4.0)
        check('C3 mass action n_e*n_h == n_p^2/4 at mu_e/mu_h = %.1f' % ratio,
              rel < 1e-15,
              'max absolute residual %.4e (relative %.2e) -- mobility must not '
              'enter the electrostatics' % (prod_resid, rel))

    print('\n  Monotonicity of n_net(dV), required before inverting it:')
    nn_mono = sgn.n_net_of_dV(np.linspace(-5.0, 5.0, 2001))
    strictly_up = bool(np.all(np.diff(nn_mono) > 0.0))
    check('n_net(dV) strictly increasing on [-5, 5] V', strictly_up,
          'min forward difference %.4e m^-2' % float(np.min(np.diff(nn_mono))))

    # -- Section 3 ----------------------------------------------------------
    print('\nSECTION 3.  The closed form (A4)/(A5)')
    print('  %-9s %-11s %-11s %-14s %-14s %-11s' %
          ('mu_e/mu_h', 'r', 'n*/n_p', 'sigma_min num', 'sigma_min (A5)', 'rel. dev.'))
    for ratio in (1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0):
        mu_e, mu_h = mu_pair(ratio)
        r = r_of(mu_e, mu_h)
        n_star = n_net_at_minimum(mu_e, mu_h)
        # numerical minimum on a fine grid, refined by golden section
        grid = np.linspace(-20.0 * gfet.n_puddle, 20.0 * gfet.n_puddle, 400001)
        s = sigma_asym_closed(grid, mu_e, mu_h)
        s_num = float(np.min(s))
        s_cf = float(sigma_min_closed(mu_e, mu_h))
        dev = abs(s_num / s_cf - 1.0)
        print('  %-9.2f %-11.6f %-11.6f %-14.7e %-14.7e %-11.2e'
              % (ratio, r, n_star / gfet.n_puddle, s_num, s_cf, dev))
        check('C4 d(sigma)/d(n) == 0 at n* for mu_e/mu_h = %.2f' % ratio,
              abs(float(dsigma_dn(n_star, mu_e, mu_h)))
              < 1e-16 * E_CHARGE * gfet.mu * gfet.n_puddle * 1e3,
              'derivative at n* = %.4e  (scale e*mu*n_p = %.4e)'
              % (float(dsigma_dn(n_star, mu_e, mu_h)),
                 E_CHARGE * gfet.mu * gfet.n_puddle))
        check('    grid minimum agrees with (A5) for mu_e/mu_h = %.2f' % ratio,
              dev < 1e-9, 'relative deviation %.3e' % dev)

    print('\n  C5 exchange symmetry mu_e <-> mu_h:')
    for ratio in (1.3, 2.0):
        mu_e, mu_h = mu_pair(ratio)
        n1 = n_net_at_minimum(mu_e, mu_h)
        n2 = n_net_at_minimum(mu_h, mu_e)
        s1 = float(sigma_min_closed(mu_e, mu_h))
        s2 = float(sigma_min_closed(mu_h, mu_e))
        check('    n*(mu_h, mu_e) == -n*(mu_e, mu_h) bitwise at ratio %.1f' % ratio,
              float(n2) == -float(n1), 'n* = %r vs %r' % (float(n1), float(n2)))
        check('    sigma_min symmetric in the pair, bitwise, at ratio %.1f' % ratio,
              s1 == s2, 'sigma_min = %r vs %r' % (s1, s2))

    print('\n  C6 AM-GM, i.e. the asymmetric floor never exceeds the symmetric one:')
    worst_gap = None
    for ratio in (1.0, 1.05, 1.2, 1.5, 2.0, 4.0):
        mu_e, mu_h = mu_pair(ratio)
        gm = float(np.sqrt(mu_e * mu_h))
        am = 0.5 * (mu_e + mu_h)
        ok = (gm < am) if ratio != 1.0 else (gm == am)
        if worst_gap is None or (am - gm) > worst_gap:
            worst_gap = am - gm
        check('    sqrt(mu_e*mu_h) %s mu_bar at ratio %.2f'
              % ('==' if ratio == 1.0 else '<', ratio), ok,
              'GM = %.9g, AM = %.9g, GM/AM = %.9f' % (gm, am, gm / am))

    print('\n  CORRECTION OF RECORD.  sigma_min_ch4() in')
    print('  graphene_fet_signed_carrier_model.py returns e*mu*n_puddle and its')
    print('  docstring calls that EXACT.  Measured here:')
    s_ch4 = float(sgn.sigma_min_ch4())
    s_diag = float(sigma_min_closed(gfet.mu, gfet.mu))
    check('    (A5) reduces to the shipped sigma_min on the diagonal within 1 ulp',
          int(ulp_distance(np.array([s_ch4]), np.array([s_diag]))[0]) <= 1,
          'sigma_min_ch4 = %r, (A5) at mu_e=mu_h = %r, ulp distance %d'
          % (s_ch4, s_diag,
             int(ulp_distance(np.array([s_ch4]), np.array([s_diag]))[0])))
    print('         So no committed number moves; what was wrong was the SCOPE of')
    print('         the word "exactly".  e*mu*n_p is the mu_e == mu_h case of')
    print('         e*n_p*sqrt(mu_e*mu_h), and off the diagonal it is an')
    print('         OVERESTIMATE of the floor by mu_bar/sqrt(mu_e*mu_h).')

    # -- Section 4 ----------------------------------------------------------
    print('\nSECTION 4.  A mobility asymmetry impersonates a doping offset')
    print('  %-9s %-9s %-16s %-16s %-14s' %
          ('mu_e/mu_h', 'r', 'dV_app (V)', 'as frac of V_dirac',
           'n_mis (cm^-2)'))
    for ratio in (1.05, 1.1, 1.2, 1.3, 1.5):
        mu_e, mu_h = mu_pair(ratio)
        dv = apparent_dirac_offset(mu_e, mu_h)
        n_mis = gfet.C_ox * abs(dv) / E_CHARGE
        print('  %-9.2f %-9.6f %-16.6f %-16.4f %-14.4e'
              % (ratio, r_of(mu_e, mu_h), dv, abs(dv) / gfet.V_dirac,
                 n_mis * 1e-4))
    print('\n  Reading: at the shipped 90 nm SiO2 back gate (C_ox = %.4g F/m^2)'
          % gfet.C_ox)
    print('  a 30 % electron-hole mobility asymmetry moves the resistance')
    print('  maximum by %.3f V, which is %.1f %% of the %.2f V this repository'
          % (abs(apparent_dirac_offset(*mu_pair(1.3))),
             100.0 * abs(apparent_dirac_offset(*mu_pair(1.3))) / gfet.V_dirac,
             gfet.V_dirac))
    print('  attributes entirely to built-in doping.')
    print('  [SUPERSEDED IN SCOPE 2026-10-30: this offset is an OVERESTIMATE')
    print('   by 0.07 %% at mu_e/mu_h = 1.3 once the mobility pair is allowed')
    print('   to depend on density.  It is only 0.07 %% because n_p/n_c =')
    print('   0.167 for this device -- the DIRAC POINT is unscreened even')
    print('   though the high-bias branches are not.  The same figure for a')
    print('   thin-oxide top-gated stack is not small.  See')
    print('   graphene_density_dependent_mobility.py Section 9.]')
    print('  A doping extraction that')
    print('  reads V_dirac off a transfer curve without an independent mobility')
    print('  measurement on each branch cannot separate the two.')
    print('  The sign is informative: mu_e > mu_h puts the minimum on the HOLE')
    print('  side (n* < 0), i.e. the apparent Dirac point moves to LOWER V_g.')

    sign_ok = apparent_dirac_offset(*mu_pair(1.3)) < 0.0
    check('direction: mu_e > mu_h shifts the apparent Dirac point to lower V_g',
          sign_ok, 'dV_app = %.6f V' % apparent_dirac_offset(*mu_pair(1.3)))

    # -- Section 5 ----------------------------------------------------------
    print('\nSECTION 5.  The length discriminator, and C7')
    print('  R_odd(L) = (L/W)*(1/sigma(-dV) - 1/sigma(+dV)) + (Rc_p - Rc_n)')
    print('  exactly linear in L: slope = channel asymmetry, intercept = contact.')
    lengths = np.linspace(100e-9, 5e-6, 25)
    dV_probe = 2.5
    Rc_half = gfet.Rc_total / 2.0
    dRc_planted = 120.0
    scenarios = [
        ('channel only', mu_pair(1.3), Rc_half, Rc_half),
        ('contact only', (gfet.mu, gfet.mu), Rc_half, Rc_half + dRc_planted),
        ('both', mu_pair(1.3), Rc_half, Rc_half + dRc_planted),
    ]
    print('\n  %-14s %-16s %-16s %-14s %-12s' %
          ('scenario', 'slope (Ohm/m)', 'intercept (Ohm)', 'max resid', 'verdict'))
    for name, (mu_e, mu_h), rn, rp in scenarios:
        R_odd = R_odd_vs_length(dV_probe, mu_e, mu_h, lengths, rn, rp)
        slope, intercept, resid = extract_from_length_sweep(dV_probe, lengths, R_odd)
        s_plus = float(sigma_asym_closed(sgn.n_net_of_dV(+dV_probe), mu_e, mu_h))
        s_minus = float(sigma_asym_closed(sgn.n_net_of_dV(-dV_probe), mu_e, mu_h))
        slope_truth = (1.0 / s_minus - 1.0 / s_plus) / gfet.W
        inter_truth = rp - rn
        rel_slope = (abs(slope / slope_truth - 1.0) if slope_truth != 0.0
                     else abs(slope))
        abs_inter = abs(intercept - inter_truth)
        verdict = ('CHANNEL' if abs_inter < 1e-6 and abs(slope) > 1.0
                   else 'CONTACT' if abs(slope) < 1e-6
                   else 'BOTH')
        print('  %-14s %-16.6e %-16.6f %-14.3e %-12s'
              % (name, slope, intercept, resid, verdict))
        check('C7 slope recovers the planted channel asymmetry (%s)' % name,
              rel_slope < 1e-12,
              'fitted %.9e vs planted %.9e, relative %.2e'
              % (slope, slope_truth, rel_slope))
        check('C7 intercept recovers the planted contact asymmetry (%s)' % name,
              abs_inter < 1e-9,
              'fitted %.9f vs planted %.9f Ohm' % (intercept, inter_truth))
        check('C7 residual to the straight line is at the rounding floor (%s)'
              % name, resid < 1e-8,
              'max residual %.3e Ohm over %d lengths' % (resid, len(lengths)))

    print('\n  Why this is the discriminator rather than a fit convenience:')
    print('  the contact term enters R_odd with NO factor of L, so it cannot be')
    print('  absorbed into the slope at any length, and the channel term enters')
    print('  with L and cannot be absorbed into the intercept.  Huard et al.')
    print('  (arXiv:0804.2040) report exactly this: their odd resistance')
    print('  collapses when scaled by the metal/graphene INTERFACE WIDTH, and the')
    print('  bulk-only reading of the same data gives a resistivity asymmetry')
    print('  that drifts with probe spacing -- the signature of a nonzero')
    print('  intercept being misread as a slope.')
    print('\n  Cost of the simplification: a single-length measurement cannot do')
    print('  this at all.  At the shipped L = %.0f nm the contact-only and' % (gfet.L * 1e9))
    print('  channel-only scenarios differ in R_odd by only:')
    ro_ch = float(R_odd_vs_length(dV_probe, *mu_pair(1.3), [gfet.L], Rc_half, Rc_half)[0])
    ro_ct = float(R_odd_vs_length(dV_probe, gfet.mu, gfet.mu, [gfet.L],
                                  Rc_half, Rc_half + dRc_planted)[0])
    print('    channel-only R_odd = %10.4f Ohm' % ro_ch)
    print('    contact-only R_odd = %10.4f Ohm' % ro_ct)
    print('    ratio              = %10.4f' % (ro_ch / ro_ct))
    print('  i.e. at 200 nm the two mechanisms are within a factor %.2f of each'
          % (ro_ct / ro_ch if ro_ct > ro_ch else ro_ch / ro_ct))
    print('  other and a single-length experiment is blind between them.')

    # -- Section 6 ----------------------------------------------------------
    print('\nSECTION 6.  What the shipped model is, in the new vocabulary')
    print('  The committed GFET model sits exactly on the mu_e == mu_h diagonal,')
    print('  i.e. r = 0, so every number already in this repository is the r = 0')
    print('  case of the forms above and NONE of them moves.  What changes is')
    print('  what they are claims ABOUT: the floor e*mu*n_p is a geometric mean')
    print('  evaluated on a diagonal, and the V_dirac = %.2f V offset is a' % gfet.V_dirac)
    print('  doping offset only under an assumption (r = 0) that the model makes')
    print('  silently and that no committed check tests.')

    path = make_figure()
    print('\n  Figure written: %s' % path)

    print('\n' + '=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
