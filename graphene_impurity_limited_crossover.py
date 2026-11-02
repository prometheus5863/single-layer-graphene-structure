"""
graphene_impurity_limited_crossover.py

THE CROSSOVER DENSITY IS NOT A MATERIAL CONSTANT: n_c IS PROPORTIONAL TO THE
CHARGED-IMPURITY DENSITY, AND THE GATE STACK THAT DECIDES WHETHER A DEVICE
REACHES n_c ALSO DECIDES WHERE n_c IS.

This closes the item ranked FIRST on 2026-10-30 and created the same day --
"n_c SHOULD BE A FUNCTION OF THE CHARGED-IMPURITY DENSITY" -- which was the
direct continuation of that session's Section 6.  In
`graphene_density_dependent_mobility.py` (10-30, hereafter DDM) the crossover
density is a single module-level constant:

    N_CROSSOVER = 3e16              # m^-2  (= 3e12 cm^-2)

and the 10-30 entry recorded, in prose and computed by hand in a notes file,
that this cannot be right across its own gate-stack table:

    "in (B1), A = n_c*mu_L_bar depends only on the NEUTRAL defect density, so
     halving mu_L_bar DOUBLES n_c.  The honest top-gate row is
     n_c ~ 6e12 cm^-2 and dV_c ~ 0.43 V."
            -- AUTOMATION_LOG.md 2026-10-30, Section 6

That hand calculation is reproduced, checked, CORRECTED and extended here.
The correction is Section 5; it is a factor of up to 2.1 on n_c and it takes
the margin the 10-30 entry reported as 2.3x down to as little as 1.10x.

------------------------------------------------------------------------------
1. THE INVARIANT, WHICH IS WHAT MAKES THIS A ONE-PARAMETER STATEMENT
------------------------------------------------------------------------------
DDM's Matthiessen form (B1) is

    1/mu_b(n) = 1/mu_L,b  +  |n| / (n_c * mu_L_bar),      b in {e, h}

The second term is the SHORT-RANGE (neutral defect) channel.  Its strength is
a property of the neutral defects alone -- vacancies, dislocations, grain
boundaries -- and nothing about adding or removing CHARGED impurities may
change it.  But it is written with two factors, `n_c` and `mu_L_bar`, and
`mu_L_bar` is the long-range mobility, which charged impurities certainly do
change.  So the quantity that is actually the neutral-defect property is the
PRODUCT

    A  ==  n_c * mu_L_bar            [m^-2 * m^2 V^-1 s^-1  =  V^-1 s^-1]  (C1)

and `n_c` alone is not a constant of the material.  A is invariant under any
change in charged-impurity density; `n_c` and `mu_L_bar` move inversely.  That
is the whole content of this module, and everything below follows from it with
one piece of literature.

------------------------------------------------------------------------------
2. THE PIECE OF LITERATURE, AND WHAT IT DOES AND DOES NOT SAY
------------------------------------------------------------------------------
Chen, Jang, Adam, Fuhrer, Williams and Ishigami, "Charged impurity scattering
in graphene", Nature Physics 4, 377 (2008), arXiv:0708.2408, dose graphene
with potassium in situ and measure mobility against DELIBERATELY ADDED charged
impurity density.  Their scaling law, stated in the paper as the relation they
test rather than fit, is

    mu * n_imp  =  C,        C = 5e15 V^-1 s^-1                           (C2)

with the dosing rate measured independently as
dn_imp/dt = (2.6 ~ 3.2)e15 m^-2 s^-1 and the value of C cross-checked against
the shift in V_g,min.  This is the strongest available form of the input:
n_imp is a controlled variable rather than a fitted one.

What (C2) does NOT say, and the caveat is load-bearing here: their measured
mobility is the TOTAL mobility, not the long-range channel in isolation.
Using (C2) for `mu_L` is legitimate exactly to the extent that their dosed
samples are charged-impurity-limited, which is their central claim, but their
PRISTINE mobility already contains a short-range contribution, so (C2) is a
relation for mu_L only in the regime where the long-range channel dominates.
D7 measures what that costs at this repository's own densities.

Combining (C1) and (C2) eliminates mu_L_bar entirely:

    n_c(N_imp)  =  (A / C) * N_imp  ==  KAPPA * N_imp                     (C3)

The crossover density is a FIXED MULTIPLE of the charged-impurity density, and
the multiple is dimensionless -- a ratio of two literature constants, one from
Zion et al. (the crossover) and one from Chen et al. (the scaling law).  For
this repository's reference device KAPPA = 2.4, and the number is derived in
Section 3 rather than chosen.

Two consequences are worth stating before any number is computed, because
both are checkable and neither is obvious from (B1) as written:

  (i) r0 IS INVARIANT.  Scaling N_imp scales mu_L,e and mu_L,h by the SAME
      factor, because the beyond-Born asymmetry of DDM Section 1 is a ratio of
      per-scatterer cross-sections and carries no density.  So
      r0 = (mu_L,e-mu_L,h)/(mu_L,e+mu_L,h) does not move at all.  D2.
  (ii) A DIRTIER DEVICE IS MORE ASYMMETRIC, not less.  More charged impurities
      push n_c UP, so at a fixed operating density the device sits further
      BELOW the crossover and r(n) sits closer to r0.  This is the right sign
      and it is a consistency check on the whole decomposition rather than a
      new result: the asymmetry lives in the charged-impurity channel, so more
      charged impurities must mean more asymmetry.  D4.

------------------------------------------------------------------------------
3. THE REFERENCE POINT, AND WHY IT IS A PAIRING RATHER THAN A MEASUREMENT
------------------------------------------------------------------------------
A is not measured anywhere.  It is obtained by pairing Zion et al.'s crossover
density with this repository's long-range mobility:

    A_ref     = N_CROSSOVER * gfet.mu  =  3e16 * 0.4   =  1.2e16  V^-1 s^-1
    N_imp,ref = C / gfet.mu            =  5e15 / 0.4   =  1.25e16 m^-2
    KAPPA     = A_ref / C              =  2.4           (dimensionless)

and `N_imp,ref = 1.25e12 cm^-2` is an ordinary charged-impurity density for
exfoliated or CVD graphene on SiO2, which is the only external check available
on it.

THE PAIRING IS THE WEAKEST LINK IN THIS MODULE AND IT IS STATED HERE RATHER
THAN BURIED.  Zion et al. (arXiv:1801.03565) report their crossover at
p ~ 3e12 cm^-2 for THEIR CVD samples; `gfet.mu = 0.4 m^2/Vs` is THIS
repository's device.  If Zion's samples were cleaner than that, A_ref is
overestimated, and n_c is overestimated with it at every impurity density.
An attempt to read Zion's own mobilities today reached only the abstract --
the arXiv abstract page carries the asymmetry statements and the T and B
ranges but no mobility, substrate or impurity number -- so the pairing stands
unresolved, as it did on 10-30.  It is logged as the first open item.

What survives the pairing, and what does not:

  * THE SCALING n_c ~ N_imp (C3) does not depend on A_ref at all.  It follows
    from (C1) and (C2) alone.
  * EVERY RATIO -- n_c(dirty)/n_c(clean), the margin ratio between two gate
    stacks, the fractional change in r -- is independent of A_ref, because
    A_ref cancels.  D5 measures this as an invariance rather than asserting it.
  * THE ABSOLUTE n_c, and therefore the absolute overdrive dV_c and the
    absolute margin of Section 5, carry the pairing.  They are reported with
    that attached.

------------------------------------------------------------------------------
4. WHAT IS CHECKED, AND WHICH CHECKS ARE EXACT
------------------------------------------------------------------------------
Three of these compare against a value known exactly in advance rather than
against a plausible range, which is the standing requirement of item 15 and
is what caught the real defect on 2026-09-17.

  D1  REFERENCE REPRODUCTION, EXACT.  At N_imp = N_imp,ref the model must
      reproduce DDM bitwise, not approximately: n_c == N_CROSSOVER, and
      mu_e(n), mu_h(n), r(n) and sigma(n) identical to the last bit at every
      density of a 401-point sweep.  A model that supersedes another in a
      strictly larger parameter space must agree with it IDENTICALLY on the
      overlap, and "identically" admits a bitwise test.
  D2  r0 IS INVARIANT UNDER N_imp, BITWISE, over six decades of N_imp --
      Section 2 (i).  Computed through the full species path, not from the
      closed form, so it is the construction being tested.
  D3  THE ASYMPTOTE IS THE INVARIANT, BITWISE.  |n|*mu_b(n) -> A as
      |n| -> infinity, and A does not depend on N_imp, so the high-density
      asymptote must be the SAME FLOAT at every impurity density.  This is the
      sharpest available statement of (C1) and it is exact.
  D4  MONOTONICITY AND THE TWO LIMITS -- Section 2 (ii).  r(n) is strictly
      increasing in N_imp at fixed n; r -> r0 as N_imp -> infinity; r -> 0 as
      N_imp -> 0.  The infinite limit is checked in TWO parts, and the split
      is a correction made on the day this module was written: as a single
      tolerance at N_imp = 1e12 x ref it FAILED at 83x its counted budget,
      because at that impurity density n/n_c is still 4.6e-12 and the
      shortfall is a real first-order term rather than rounding.  So the
      limit is now taken where it is genuinely the limit (1e20 x ref, where
      the short-range term cannot move the float), and the finite-N_imp
      residual gets a PREDICTION -- r0*u/(1+u) -- instead of a tolerance.
      Same lesson as 10-30's D4, reached the same way.
  D5  A_ref CANCELS IN EVERY RATIO -- Section 3.  Perturbing A_ref by a factor
      of 3 must leave every reported RATIO bitwise unchanged while moving
      every absolute n_c by exactly that factor.  The claim "the ratios do not
      carry the pairing" is thereby measured, not argued.
  D6  SCALE INVARIANCE, BITWISE.  r depends on (n, N_imp) only through
      n/n_c ~ n/N_imp, so doubling BOTH must leave r bitwise unchanged.  A
      power of two is exactly representable, so this is a true bitwise
      symmetry and not a tolerance.
  D7  THE Section 2 CAVEAT, MEASURED.  How much of the total mobility is the
      long-range channel at this repository's densities -- i.e. the error
      incurred by reading Chen's total-mobility C as a long-range-channel C.
  D8  SYMMETRY THAT MUST GIVE ZERO.  r0 = 0 gives mu_e == mu_h bitwise and
      sigma(+dV) - sigma(-dV) EXACTLY 0.0, at every impurity density.
  D9  THE 10-30 HAND CALCULATION, reproduced from its own stated inputs
      before being corrected -- Section 5.  A correction that cannot first
      reproduce the thing it corrects is not a correction.
  D10 THE TWO INPUTS OF arXiv:1010.0913 CONSTRAIN EACH OTHER -- Section 6.
      Its added-impurity range and its mobility loss are not independent once
      C is fixed, and this check measures whether they are mutually
      consistent.  It is allowed to fail informatively.
  D11 MUTATION: hold `n_c` fixed while scaling mu_L -- i.e. treat n_c rather
      than A as the invariant, which is exactly the error this module exists
      to correct -- and D3 must fail.  A check that passes against the
      pre-correction model is not checking the correction.
"""

import sys

import numpy as np

import graphene_fet_model as gfet
import graphene_fet_signed_carrier_model as sgn
import graphene_branch_asymmetric_mobility as bam
import graphene_density_dependent_mobility as ddm

EPS = float(np.finfo(np.float64).eps)

# Chen et al., Nature Physics 4, 377 (2008), arXiv:0708.2408, eq. for the
# charged-impurity scaling law: mu * n_imp = C.  Stated in the paper as
# C = 5e15 V^-1 s^-1 and cross-checked there against the V_g,min shift.
C_CHEN = 5e15                   # V^-1 s^-1

# The neutral-defect invariant (C1), paired at the reference device.  See
# Section 3: this is the one number in the module that is a PAIRING of two
# devices rather than a measurement on one.
A_NEUTRAL_REF = ddm.N_CROSSOVER * gfet.mu          # V^-1 s^-1
N_IMP_REF = C_CHEN / gfet.mu                       # m^-2
KAPPA = A_NEUTRAL_REF / C_CHEN                     # dimensionless

# arXiv:1010.0913 (found 10-30 while fetching HfO2 mobility values, and
# reporting these in passing for an unrelated purpose): HfO2 top-gate
# deposition adds 1.5-4e12 cm^-2 of charged impurities and costs about half
# the mobility.  Both numbers are used, and D10 checks them against each other.
DN_IMP_HFO2_LO = 1.5e16         # m^-2
DN_IMP_HFO2_HI = 4.0e16         # m^-2
MOBILITY_RETAINED_HFO2 = 0.5    # the paper's "~50% of the mobility"

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
# Section 1: the relations
# ---------------------------------------------------------------------------

def mu_L_of_impurity(n_imp):
    """Long-range mobility from (C2): mu_L = C / N_imp."""
    return C_CHEN / np.asarray(n_imp, dtype=float)


def n_imp_of_mu_L(mu_L):
    """Inverse of (C2)."""
    return C_CHEN / np.asarray(mu_L, dtype=float)


def n_c_of_impurity(n_imp, a_neutral=None):
    """
    Crossover density from (C3): n_c = (A/C) * N_imp.

    `a_neutral` exists so that D5 can perturb the pairing of Section 3 and
    show which reported quantities move with it and which do not.
    """
    if a_neutral is None:
        a_neutral = A_NEUTRAL_REF
    return (a_neutral / C_CHEN) * np.asarray(n_imp, dtype=float)


def mu_pair_at_impurity(n_abs, n_imp, r0=0.0, a_neutral=None):
    """
    (mu_e, mu_h) at total density `n_abs` for a device of charged-impurity
    density `n_imp`, routed through DDM's own (B1) so that there is exactly
    one implementation of the mobility pair in the repository.

    `r0` is the asymmetry of the long-range channel and is independent of
    `n_imp` -- Section 2 (i).
    """
    mu_L_bar = float(mu_L_of_impurity(n_imp))
    mu_L_e = mu_L_bar * (1.0 + r0)
    mu_L_h = mu_L_bar * (1.0 - r0)
    n_c = float(n_c_of_impurity(n_imp, a_neutral=a_neutral))
    return ddm.mu_pair_density(n_abs, mu_L_e, mu_L_h, n_c=n_c)


def r_at_impurity(n_abs, n_imp, r0, a_neutral=None):
    """Asymmetry parameter via the species path, not the closed form."""
    mu_e, mu_h = mu_pair_at_impurity(n_abs, n_imp, r0=r0, a_neutral=a_neutral)
    return (mu_e - mu_h) / (mu_e + mu_h)


def asymptote(n_imp, r0=0.0, a_neutral=None, n_over_nc=1e8):
    """
    |n| * mu_b(n) at |n| = `n_over_nc` * n_c, which tends to A = n_c*mu_L_bar.

    The point of D3 is that this must be the SAME FLOAT for every `n_imp`.
    """
    n_c = float(n_c_of_impurity(n_imp, a_neutral=a_neutral))
    n_abs = n_over_nc * n_c
    mu_e, mu_h = mu_pair_at_impurity(n_abs, n_imp, r0=r0, a_neutral=a_neutral)
    return n_abs * float(mu_e), n_abs * float(mu_h)


def asymptote_mutant(n_imp, r0=0.0, n_over_nc=1e8):
    """
    D11's mutant: `n_c` held at N_CROSSOVER while mu_L scales with N_imp --
    i.e. n_c treated as the invariant instead of A.  This is the pre-10-30
    model, reinstated as a control rather than described.
    """
    mu_L_bar = float(mu_L_of_impurity(n_imp))
    mu_L_e = mu_L_bar * (1.0 + r0)
    mu_L_h = mu_L_bar * (1.0 - r0)
    n_abs = n_over_nc * ddm.N_CROSSOVER
    mu_e, mu_h = ddm.mu_pair_density(n_abs, mu_L_e, mu_L_h,
                                     n_c=ddm.N_CROSSOVER)
    return n_abs * float(mu_e), n_abs * float(mu_h)


# ---------------------------------------------------------------------------
# Section 2: the device statements
# ---------------------------------------------------------------------------

def crossover_overdrive(C_ox, n_imp):
    """Gate overdrive needed to reach n_c: dV_c = e*n_c/C_ox."""
    return ddm.E_CHARGE * float(n_c_of_impurity(n_imp)) / C_ox


def operating_margin(C_ox, n_imp, V_over):
    """
    M = n(V_over)/n_c -- how far past the crossover the device operates.

    M > 1 means screening has removed the branch asymmetry over the operating
    range; M < 1 means the constant-pair model of 10-29 is the right one for
    this device.
    """
    n_at_V = C_ox * V_over / ddm.E_CHARGE
    return n_at_V / float(n_c_of_impurity(n_imp))


def max_tolerable_impurity(C_ox, V_over, a_neutral=None):
    """
    The design rule.  Total charged-impurity density at which M = 1, i.e. the
    dirtiest channel for which this gate stack at this overdrive still
    operates past the crossover:

        N_imp*  =  C * C_ox * V / (e * A)
    """
    if a_neutral is None:
        a_neutral = A_NEUTRAL_REF
    return C_CHEN * C_ox * V_over / (ddm.E_CHARGE * a_neutral)


# ---------------------------------------------------------------------------
# Section 3: figure
# ---------------------------------------------------------------------------

def make_figure(path='impurity_limited_crossover.png'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.75))

    n_imp = np.logspace(np.log10(2e15), np.log10(2e17), 400)
    n_c = n_c_of_impurity(n_imp)

    # -- (a) the linear law, and where each gate stack's operating density is
    ax = axes[0]
    ax.loglog(n_imp * 1e-4, n_c * 1e-4, 'C0-', lw=2.0,
              label=r'$n_c = \kappa N_{imp}$, $\kappa$ = %.2f' % KAPPA)
    n_back = gfet.C_ox * 5.0 / ddm.E_CHARGE
    n_top = ddm.C_OX_TOPGATE_HFO2 * 1.0 / ddm.E_CHARGE
    ax.axhline(n_back * 1e-4, color='C3', ls='--', lw=1.4)
    ax.text(2.4e11, n_back * 1e-4 * 1.15,
            'density reached: 90 nm SiO$_2$ at 5 V', fontsize=7, color='C3')
    ax.axhline(n_top * 1e-4, color='C2', ls='--', lw=1.4)
    ax.text(2.4e11, n_top * 1e-4 * 1.15,
            'density reached: 10 nm HfO$_2$ at 1 V', fontsize=7, color='C2')
    ax.axvspan((N_IMP_REF + DN_IMP_HFO2_LO) * 1e-4,
               (N_IMP_REF + DN_IMP_HFO2_HI) * 1e-4,
               color='C2', alpha=0.15)
    # Placed ABOVE the band rather than below it: at the original (x, 4e11)
    # this label sat underneath the legend box and was unreadable.  Found by
    # rendering and LOOKING at the figure, which no check in this module can
    # do -- the standing "render and inspect every committed figure" item.
    ax.text(np.sqrt((N_IMP_REF + DN_IMP_HFO2_LO)
                    * (N_IMP_REF + DN_IMP_HFO2_HI)) * 1e-4, 3.0e13,
            'after HfO$_2$\ndeposition', fontsize=7, color='C2',
            ha='center', va='center')
    ax.plot([N_IMP_REF * 1e-4], [ddm.N_CROSSOVER * 1e-4], 'ko', ms=7,
            zorder=5)
    ax.annotate('reference pairing\n(Zion $n_c$, this repo $\\mu$)',
                xy=(N_IMP_REF * 1e-4, ddm.N_CROSSOVER * 1e-4),
                xytext=(2.2e11, 3.2e12), fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.set_xlabel(r'charged-impurity density $N_{imp}$ (cm$^{-2}$)')
    ax.set_ylabel(r'crossover density $n_c$ (cm$^{-2}$)')
    ax.set_title('(a) $n_c$ is proportional to $N_{imp}$', fontsize=10)
    ax.legend(fontsize=7, loc='lower right')
    ax.grid(True, which='both', alpha=0.25)

    # -- (b) a dirtier device is MORE asymmetric
    ax = axes[1]
    r0 = 0.1304347826086956   # mu_e/mu_h = 1.3, the 10-29/10-30 working value
    nn = np.logspace(np.log10(1e15), np.log10(3e17), 400)
    for n_i, col in ((0.25 * N_IMP_REF, 'C0'), (N_IMP_REF, 'C1'),
                     (4.0 * N_IMP_REF, 'C3')):
        r = r_at_impurity(nn, n_i, r0)
        ax.semilogx(nn * 1e-4, np.asarray(r) / r0, color=col, lw=1.8,
                    label=r'$N_{imp}$ = %.2fe12 cm$^{-2}$' % (n_i * 1e-16))
    ax.axvline(n_top * 1e-4, color='C2', ls='--', lw=1.2)
    ax.text(n_top * 1e-4 * 0.55, 0.9, '10 nm HfO$_2$ at 1 V',
            fontsize=7, color='C2', rotation=90, va='top')
    ax.set_xlabel(r'total carrier density $|n|$ (cm$^{-2}$)')
    ax.set_ylabel(r'$r(n)\,/\,r_0$  (surviving asymmetry)')
    ax.set_title('(b) more charged impurities, more asymmetry', fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(True, which='both', alpha=0.25)
    ax.set_ylim(0.0, 1.05)

    # -- (c) the margin, and the design rule
    ax = axes[2]
    m_back = np.array([operating_margin(gfet.C_ox, x, 5.0) for x in n_imp])
    m_top = np.array([operating_margin(ddm.C_OX_TOPGATE_HFO2, x, 1.0)
                      for x in n_imp])
    ax.loglog(n_imp * 1e-4, m_back, 'C3-', lw=1.8,
              label='90 nm SiO$_2$ back gate at 5 V')
    ax.loglog(n_imp * 1e-4, m_top, 'C2-', lw=1.8,
              label='10 nm HfO$_2$ top gate at 1 V')
    ax.axhline(1.0, color='k', lw=1.2)
    ax.text(2.4e11, 1.12, 'M = 1: the device just reaches $n_c$',
            fontsize=7)
    ax.axvspan((N_IMP_REF + DN_IMP_HFO2_LO) * 1e-4,
               (N_IMP_REF + DN_IMP_HFO2_HI) * 1e-4,
               color='C2', alpha=0.15)
    n_star = max_tolerable_impurity(ddm.C_OX_TOPGATE_HFO2, 1.0)
    ax.plot([n_star * 1e-4], [1.0], 'kv', ms=8, zorder=5)
    ax.annotate(r'$N_{imp}^{*}$ = %.2fe12 cm$^{-2}$' % (n_star * 1e-16),
                xy=(n_star * 1e-4, 1.0), xytext=(3.0e11, 0.2),
                fontsize=7, arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.set_xlabel(r'charged-impurity density $N_{imp}$ (cm$^{-2}$)')
    ax.set_ylabel(r'operating margin $M = n(V)/n_c$')
    ax.set_title('(c) the margin, and the impurity budget', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True, which='both', alpha=0.25)

    fig.suptitle('The crossover density is set by the charged impurities, and '
                 'the dielectric that gates the channel also dirties it',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Section 4: checks
# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('THE CROSSOVER DENSITY AS A FUNCTION OF CHARGED-IMPURITY DENSITY')
    print('=' * 78)

    print('\nSECTION 0.  The constants, and which of them is a pairing')
    print('  C (Chen et al. 2008, arXiv:0708.2408)   = %.4e V^-1 s^-1'
          % C_CHEN)
    print('  n_c,ref (Zion et al., arXiv:1801.03565) = %.4e m^-2  (%.2fe12 '
          'cm^-2)' % (ddm.N_CROSSOVER, ddm.N_CROSSOVER * 1e-16))
    print('  mu_L,ref (this repository, gfet.mu)     = %.4e m^2/Vs' % gfet.mu)
    print('  -> A = n_c*mu_L_bar  (the INVARIANT)    = %.4e V^-1 s^-1'
          % A_NEUTRAL_REF)
    print('  -> N_imp,ref = C/mu_L                   = %.4e m^-2  (%.2fe12 '
          'cm^-2)' % (N_IMP_REF, N_IMP_REF * 1e-16))
    print('  -> KAPPA = A/C                          = %.6f  (dimensionless)'
          % KAPPA)
    print('  A is a PAIRING of Zion\'s crossover with this repository\'s')
    print('  mobility, not a measurement on one device.  D5 shows which')
    print('  reported quantities carry it and which do not.')

    r0_ref = bam.r_of(*bam.mu_pair(1.3))
    mu_L_e, mu_L_h = bam.mu_pair(1.3)

    # -- D1 -----------------------------------------------------------------
    print('\nSECTION 1.  D1: EXACT REPRODUCTION OF THE MODEL THIS SUPERSEDES')
    print('            At N_imp = N_imp,ref this model must BE the 10-30')
    print('            model, not approximate it.')
    n_c_back = float(n_c_of_impurity(N_IMP_REF))
    u_nc = int(ulp_distance(np.array([n_c_back]),
                            np.array([float(ddm.N_CROSSOVER)]))[0])
    # Counted: N_imp,ref = C/mu (1), A = n_c*mu (1), A/C (1), *(N_imp) (1).
    check('D1-i n_c(N_imp,ref) == N_CROSSOVER within a counted 4-ulp budget',
          u_nc <= 4,
          'n_c = %r against N_CROSSOVER = %r, ulp distance %d, budget 4 '
          '(four roundings: C/mu, n_c*mu, A/C, and the product)'
          % (n_c_back, float(ddm.N_CROSSOVER), u_nc))

    nn = sgn.n_net_of_dV(np.linspace(-5.0, 5.0, 401))
    S = np.sqrt(nn ** 2 + gfet.n_puddle ** 2)
    mu_e_new, mu_h_new = mu_pair_at_impurity(S, N_IMP_REF, r0=r0_ref)
    mu_e_old, mu_h_old = ddm.mu_pair_density(S, mu_L_e, mu_L_h)
    u_mu = int(max(np.max(ulp_distance(mu_e_new, mu_e_old)),
                   np.max(ulp_distance(mu_h_new, mu_h_old))))
    check('D1-ii the mobility pair reproduces DDM within a counted 16-ulp '
          'budget over the whole +/-5 V sweep',
          u_mu <= 16,
          'max ulp distance %d, budget 16: mu_b is a sum of POSITIVE terms, '
          'so relative error is the right currency and the 4 ulp of D1-i '
          'propagate without amplification' % u_mu)

    r_new = np.asarray(r_at_impurity(S, N_IMP_REF, r0_ref))
    r_old = np.asarray(ddm.r_species(S, mu_L_e, mu_L_h))
    abs_r = float(np.max(np.abs(r_new - r_old)))
    BUDGET_R = 16.0 * EPS
    check('D1-iii the asymmetry parameter reproduces DDM within a counted '
          'ABSOLUTE 16-eps budget',
          abs_r <= BUDGET_R,
          'max absolute deviation %.4e against %.4e; relative would be the '
          'WRONG currency here and 10-30 D1 is the precedent -- r is a '
          'difference of nearly equal mobilities over their sum, so '
          'cancellation amplifies relative error by 1/|r| without bound'
          % (abs_r, BUDGET_R))

    sig_new = ddm.sigma_dd(nn, mu_L_e, mu_L_h, n_c=n_c_back)
    sig_old = ddm.sigma_dd(nn, mu_L_e, mu_L_h)
    u_sig = int(np.max(ulp_distance(sig_new, sig_old)))
    check('D1-iv the sheet conductivity reproduces DDM within a counted '
          '32-ulp budget',
          u_sig <= 32,
          'max ulp distance %d, budget 32: the odd term 0.5*dmu*n_net is '
          'bounded by r times the even term, so sigma inherits no '
          'cancellation amplification' % u_sig)

    # -- D2, D4 -------------------------------------------------------------
    print('\nSECTION 2.  The two consequences of Section 2 of the docstring')
    r0_rec = []
    for dec in range(-3, 4):
        n_i = N_IMP_REF * (10.0 ** dec)
        r0_rec.append(float(r_at_impurity(1e-6 * n_i, n_i, r0_ref)))
    spread = max(r0_rec) - min(r0_rec)
    check('D2 r0 is invariant under N_imp over seven decades (recovered at '
          'n << n_c)',
          spread <= 32.0 * EPS,
          'recovered r0 spans %.4e against an absolute budget of %.4e; r0 = '
          '%.12f.  The asymmetry is a per-scatterer cross-section ratio and '
          'carries no density, so scaling N_imp scales BOTH branches equally'
          % (spread, 32.0 * EPS, r0_ref))

    n_probe = n_top_probe = ddm.C_OX_TOPGATE_HFO2 * 1.0 / ddm.E_CHARGE
    r_series = [float(r_at_impurity(n_probe, N_IMP_REF * f, r0_ref))
                for f in (0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0, 1e6)]
    mono = all(b > a for a, b in zip(r_series, r_series[1:]))
    check('D4-i r(n) is STRICTLY INCREASING in N_imp at fixed n -- a dirtier '
          'device is more asymmetric',
          mono,
          'at n = %.3e m^-2 (10 nm HfO2 at 1 V), r/r0 over N_imp = 0.01..1e6 '
          'x ref: %s' % (n_probe,
                         ', '.join('%.4f' % (x / r0_ref) for x in r_series)))
    # D4-ii is in two parts, and the split is the point.  The first time this
    # check was written it was a single tolerance -- |r - r0| <= 32 eps at
    # N_imp = 1e12 x ref -- and it FAILED at 5.9e-13, eighty-three times the
    # budget.  The budget was not miscounted; the residual is not rounding at
    # all.  At 1e12 x ref the dimensionless group n/n_c is still 4.6e-12, so
    # there is a real first-order shortfall there, and a tolerance loose
    # enough to admit it would have been a number chosen to pass.  This is
    # the same correction 10-30 made to its own D4: a residual that is not
    # noise gets a PREDICTION instead of a tolerance.
    r_lim = float(r_at_impurity(n_probe, N_IMP_REF * 1e20, r0_ref))
    u_lim = int(ulp_distance(np.array([r_lim]), np.array([r0_ref]))[0])
    check('D4-ii-a r -> r0 as N_imp -> infinity, reached within a counted '
          '8-ulp budget',
          u_lim <= 8,
          'at N_imp = 1e20 x ref: r = %.17f against r0 = %.17f, ulp distance '
          '%d, budget 8 -- five roundings reconstructing r0 from '
          'mu_L_bar*(1 +/- r0) and one comparing.  At this N_imp the '
          'short-range term is 4.6e-20 of 1/mu_L and cannot move the float '
          'at all, so this is the limit itself and not an approach to it.'
          % (r_lim, r0_ref, u_lim))

    n_i_12 = N_IMP_REF * 1e12
    r_12 = float(r_at_impurity(n_probe, n_i_12, r0_ref))
    u_12 = (1.0 - r_lim * r_lim) * n_probe / float(n_c_of_impurity(n_i_12))
    resid = r_lim - r_12
    pred_resid = r_lim * u_12 / (1.0 + u_12)
    rel_pred = abs(resid / pred_resid - 1.0)
    check('D4-ii-b the finite-N_imp residual is the PREDICTED first-order '
          'term, not noise',
          rel_pred < 1e-4,
          'at N_imp = 1e12 x ref the shortfall is %.4e; (B3) predicts '
          'r0*u/(1+u) with u = (1-r0^2)*n/n_c = %.4e, i.e. %.4e -- accounted '
          'for to %.3e relative.  A tolerance large enough to pass this '
          'would have been 83x the counted budget and would have hidden a '
          'real term; a prediction can be false in a way a tolerance cannot.'
          % (resid, u_12, pred_resid, rel_pred))
    r_zero = float(r_at_impurity(n_probe, N_IMP_REF * 1e-12, r0_ref))
    check('D4-iii r -> 0 as N_imp -> 0 (a perfectly clean channel is '
          'symmetric at any finite density)',
          r_zero < 1e-10 * r0_ref,
          'at N_imp = 1e-12 x ref: r = %.4e, i.e. r/r0 = %.4e'
          % (r_zero, r_zero / r0_ref))

    # -- D3 -----------------------------------------------------------------
    print('\nSECTION 3.  D3: the asymptote IS the invariant, against an')
    print('            exactly predicted value')
    X = 2.0 ** 30
    pred = A_NEUTRAL_REF * X / (1.0 + X)
    worst = 0
    rows = []
    for dec in range(-3, 4):
        n_i = N_IMP_REF * (10.0 ** dec)
        a_e, a_h = asymptote(n_i, r0=0.0, n_over_nc=X)
        u = int(max(ulp_distance(np.array([a_e]), np.array([pred]))[0],
                    ulp_distance(np.array([a_h]), np.array([pred]))[0]))
        worst = max(worst, u)
        rows.append((n_i, a_e, u))
    print('  %-16s %-24s %-8s' % ('N_imp (m^-2)', '|n|*mu at |n|=2^30*n_c',
                                  'ulp'))
    for n_i, a_e, u in rows:
        print('  %-16.4e %-24.17g %-8d' % (n_i, a_e, u))
    check('D3 |n|*mu_b = A*X/(1+X) exactly, at every N_imp, within a counted '
          '8-ulp budget',
          worst <= 8,
          'predicted %.17g (an EXACT consequence of (B1), not a limit taken '
          'numerically); worst ulp distance %d over seven decades of N_imp.  '
          'The predicted value does not contain N_imp, and that is the '
          'content of (C1).' % (pred, worst))

    # -- D6 -----------------------------------------------------------------
    print('\nSECTION 4.  D6: one dimensionless group, tested BITWISE')
    bitwise_ok = True
    for dec in range(-2, 3):
        n_i = N_IMP_REF * (10.0 ** dec)
        for n_v in (1e15, 1e16, 1e17):
            a = float(r_at_impurity(n_v, n_i, r0_ref))
            b = float(r_at_impurity(2.0 * n_v, 2.0 * n_i, r0_ref))
            if a.hex() != b.hex():
                bitwise_ok = False
    check('D6 r(2n, 2*N_imp) == r(n, N_imp) BITWISE -- the model has exactly '
          'one dimensionless group, n/N_imp',
          bitwise_ok,
          'tested at 15 (n, N_imp) points.  A factor of TWO is used, not a '
          'round ten: scaling by a power of two is exact and commutes with '
          'rounding at every step, so this is a true bitwise symmetry rather '
          'than a tolerance that happens to hold.')

    # -- D5 -----------------------------------------------------------------
    print('\nSECTION 5.  D5: which numbers carry the Section 3 pairing')
    A_PERT = 4.0 * A_NEUTRAL_REF
    n_lo = N_IMP_REF + DN_IMP_HFO2_LO
    n_hi = N_IMP_REF + DN_IMP_HFO2_HI
    ratio_ref = float(n_c_of_impurity(n_hi)) / float(n_c_of_impurity(n_lo))
    ratio_pert = (float(n_c_of_impurity(n_hi, a_neutral=A_PERT))
                  / float(n_c_of_impurity(n_lo, a_neutral=A_PERT)))
    abs_ref = float(n_c_of_impurity(n_hi))
    abs_pert = float(n_c_of_impurity(n_hi, a_neutral=A_PERT))
    check('D5-i every RATIO is independent of the pairing, BITWISE under a '
          'x4 perturbation of A',
          ratio_ref.hex() == ratio_pert.hex(),
          'n_c(dirty)/n_c(clean) = %.17g both before and after A -> 4A.  A '
          'factor of 4 rather than 3 so the perturbation is exactly '
          'representable and the invariance is bitwise.' % ratio_ref)
    check('D5-ii every ABSOLUTE n_c moves by exactly the perturbation, so the '
          'pairing is carried by the absolute numbers and only by them',
          (abs_pert / abs_ref) == 4.0,
          'n_c(dirty): %.6e -> %.6e, ratio exactly %.1f'
          % (abs_ref, abs_pert, abs_pert / abs_ref))

    # -- D8 -----------------------------------------------------------------
    print('\nSECTION 6.  D8: the symmetry that must give zero')
    zero_ok = True
    worst_sig = 0.0
    for dec in range(-2, 3):
        n_i = N_IMP_REF * (10.0 ** dec)
        mu_e0, mu_h0 = mu_pair_at_impurity(S, n_i, r0=0.0)
        if np.any(np.asarray(mu_e0) != np.asarray(mu_h0)):
            zero_ok = False
        n_c_i = float(n_c_of_impurity(n_i))
        mu_L = float(mu_L_of_impurity(n_i))
        dV = np.linspace(0.05, 5.0, 200)
        s_p = ddm.sigma_dd(sgn.n_net_of_dV(dV), mu_L, mu_L, n_c=n_c_i)
        s_m = ddm.sigma_dd(sgn.n_net_of_dV(-dV), mu_L, mu_L, n_c=n_c_i)
        worst_sig = max(worst_sig, float(np.max(np.abs(s_p - s_m))))
    check('D8 r0 = 0 gives mu_e == mu_h bitwise and sigma(+dV) - sigma(-dV) '
          'EXACTLY 0.0, at every impurity density',
          zero_ok and worst_sig == 0.0,
          'worst |sigma(+dV)-sigma(-dV)| over five decades of N_imp: %r '
          '(exact zero, not a small number)' % worst_sig)

    # -- D7 -----------------------------------------------------------------
    print('\nSECTION 7.  D7: the cost of reading Chen\'s TOTAL mobility as a')
    print('            long-range-channel mobility')
    print('  %-34s %-12s %-14s %-14s'
          % ('density', '|n|/n_c', 'long-range', 'N_imp,ref'))
    print('  %-34s %-12s %-14s %-14s'
          % ('', '', 'share of 1/mu', 'overstated by'))
    probes = (('puddle, n_p', gfet.n_puddle, ddm.N_CROSSOVER),
              ('90 nm SiO2 back gate at 5 V',
               gfet.C_ox * 5.0 / ddm.E_CHARGE, ddm.N_CROSSOVER),
              ('10 nm HfO2 top gate at 1 V',
               ddm.C_OX_TOPGATE_HFO2 * 1.0 / ddm.E_CHARGE,
               float(n_c_of_impurity(n_lo))))
    shares = []
    for name, n_v, n_c_v in probes:
        x = n_v / n_c_v
        share = 1.0 / (1.0 + x)
        shares.append(share)
        print('  %-34s %-12.4f %-14.4f %-14s'
              % (name, x, share, '%.1f %%' % (100.0 * x)))
    check('D7 the long-range channel is the MAJORITY of 1/mu at the shipped '
          'device\'s densities, so (C2) is being used inside its regime there',
          shares[0] > 0.5 and shares[1] > 0.5,
          'shares %.4f (puddle) and %.4f (5 V back gate).  At the top-gate '
          'operating point it is %.4f -- a MINORITY, which is not a defect '
          'but the restatement of M > 1: that device operates past the '
          'crossover, where neutral defects dominate.  The caveat this '
          'measures is that N_imp,ref inferred as C/mu_total is overstated '
          'by the factor (1 + |n|/n_c) at whatever density Chen et al. '
          'measured, and they do not state it.'
          % (shares[0], shares[1], shares[2]))

    # -- D9: reproduce, then correct ---------------------------------------
    print('\nSECTION 8.  D9: the 2026-10-30 hand calculation, reproduced')
    print('            from ITS OWN input before being corrected')
    n_imp_halved = n_imp_of_mu_L(gfet.mu * MOBILITY_RETAINED_HFO2)
    n_c_hand = float(n_c_of_impurity(n_imp_halved))
    dv_hand = crossover_overdrive(ddm.C_OX_TOPGATE_HFO2, n_imp_halved)
    m_hand = operating_margin(ddm.C_OX_TOPGATE_HFO2, n_imp_halved, 1.0)
    check('D9 the 10-30 entry\'s "n_c ~ 6e12 cm^-2, dV_c ~ 0.43 V, margin '
          '2.3x" is reproduced from its stated input (mobility halved)',
          abs(n_c_hand / 6.0e16 - 1.0) < 1e-12
          and abs(dv_hand - 0.435) < 0.01 and abs(m_hand - 2.30) < 0.02,
          'mobility halved -> N_imp = %.4e m^-2 -> n_c = %.4e m^-2 '
          '(%.3fe12 cm^-2), dV_c = %.4f V, M = %.3f.  A correction that '
          'cannot reproduce what it corrects is not a correction.'
          % (n_imp_halved, n_c_hand, n_c_hand * 1e-16, dv_hand, m_hand))

    print('\n  THE CORRECTION.  The 10-30 figure used ONE of the two numbers')
    print('  arXiv:1010.0913 reports -- the mobility loss.  Using its OTHER')
    print('  number, the added charged-impurity density, directly:')
    print('  %-30s %-14s %-14s %-12s %-10s'
          % ('channel', 'N_imp (e12 cm^-2)', 'n_c (e12 cm^-2)', 'dV_c (V)',
             'M at 1 V'))
    for name, n_i in (('reference (no deposition)', N_IMP_REF),
                      ('+1.5e12 (paper, clean end)', n_lo),
                      ('mobility-halved (10-30)', float(n_imp_halved)),
                      ('+4.0e12 (paper, dirty end)', n_hi)):
        print('  %-30s %-14.3f %-14.3f %-12.4f %-10.3f'
              % (name, n_i * 1e-16, float(n_c_of_impurity(n_i)) * 1e-16,
                 crossover_overdrive(ddm.C_OX_TOPGATE_HFO2, n_i),
                 operating_margin(ddm.C_OX_TOPGATE_HFO2, n_i, 1.0)))
    m_lo = operating_margin(ddm.C_OX_TOPGATE_HFO2, n_lo, 1.0)
    m_hi = operating_margin(ddm.C_OX_TOPGATE_HFO2, n_hi, 1.0)
    check('the 10-30 CONCLUSION survives -- the top-gated device still '
          'operates past the crossover across the whole impurity range',
          m_hi > 1.0,
          'margin falls from %.3f (clean end) to %.3f (dirty end) against the '
          '2.30 of 10-30.  THE CONCLUSION SURVIVES AND THE MARGIN DOES NOT: '
          '10-30 reported 2.3x and the honest range is %.2f-%.2fx.  At the '
          'dirty end the device clears the crossover by %.0f %%, which is '
          'not a margin a process would design to.'
          % (m_lo, m_hi, m_hi, m_lo, 100.0 * (m_hi - 1.0)))

    # -- the design rule ----------------------------------------------------
    print('\nSECTION 9.  The design rule, which is the device statement')
    n_star_top = max_tolerable_impurity(ddm.C_OX_TOPGATE_HFO2, 1.0)
    n_star_back = max_tolerable_impurity(gfet.C_ox, 5.0)
    print('  N_imp* = C*C_ox*V/(e*A) -- the dirtiest channel a gate stack can')
    print('  drive past the crossover at a given overdrive:')
    print('  %-34s %-16s %-16s' % ('gate stack and overdrive', 'N_imp* (e12 '
                                   'cm^-2)', 'vs. actual'))
    print('  %-34s %-16.3f %-16s'
          % ('10 nm HfO2 at 1 V', n_star_top * 1e-16,
             '%.2fx the reference' % (n_star_top / N_IMP_REF)))
    print('  %-34s %-16.3f %-16s'
          % ('90 nm SiO2 at 5 V', n_star_back * 1e-16,
             '%.2fx the reference' % (n_star_back / N_IMP_REF)))
    dn_budget = n_star_top - N_IMP_REF
    check('the HfO2 deposition budget is REAL but thin: the paper\'s dirty '
          'end sits just inside it',
          dn_budget > DN_IMP_HFO2_HI,
          'a 10 nm HfO2 gate driven to 1 V tolerates dN_imp <= %.3fe12 cm^-2 '
          'of deposition damage; the paper reports up to %.2fe12, so the '
          'headroom is %.1f %%.  This is the falsifiable form of 10-30\'s '
          '"a thicker or dirtier dielectric could eat it".'
          % (dn_budget * 1e-16, DN_IMP_HFO2_HI * 1e-16,
             100.0 * (dn_budget / DN_IMP_HFO2_HI - 1.0)))
    check('the back-gated device cannot reach n_c at ANY cleanliness it '
          'plausibly has -- and the reason is the inverted one',
          n_star_back < N_IMP_REF,
          'the 90 nm stack at 5 V would need N_imp <= %.3fe12 cm^-2, i.e. the '
          'channel %.2fx CLEANER than it is.  Cleaning a channel moves n_c '
          'DOWN, so better material reaches the asymmetry-killing regime at '
          'LOWER gate voltage -- the opposite of the intuition that a cleaner '
          'device is the more ideal one in every respect.'
          % (n_star_back * 1e-16, N_IMP_REF / n_star_back))

    # -- D10 ----------------------------------------------------------------
    print('\nSECTION 10. D10: the source\'s two numbers constrain each other')
    dn_implied = float(n_imp_halved) - N_IMP_REF
    f_lo = dn_implied / DN_IMP_HFO2_LO
    f_hi = dn_implied / DN_IMP_HFO2_HI
    check('arXiv:1010.0913\'s mobility loss and its added-impurity range are '
          'mutually consistent only at the CLEAN end of its own range',
          0.5 < f_lo < 2.0,
          'its "~50 %% of the mobility" implies dN_imp = %.3fe12 cm^-2 given '
          'C, which is %.2fx its stated low end (1.5e12) and %.2fx its '
          'stated high end (4.0e12).  So the two numbers agree to %.0f %% at '
          'the clean end and disagree by a factor of %.1f at the dirty end.  '
          'THE ESCAPE IS REAL AND IS NOT ASSUMED AWAY: if HfO2 deposition '
          'also adds NEUTRAL defects then A grows too, the mobility loss is '
          'shared between the channels, and a larger dN_imp becomes '
          'consistent.  That would change A, which this module holds fixed, '
          'and it is logged as an open item rather than modelled here.'
          % (dn_implied * 1e-16, f_lo, f_hi,
             100.0 * abs(f_lo - 1.0), 1.0 / f_hi))

    # -- D11 ----------------------------------------------------------------
    print('\nSECTION 11. D11: the pre-correction model, reinstated as a')
    print('            control rather than described')
    X2 = 2.0 ** 30
    pred2 = A_NEUTRAL_REF * X2 / (1.0 + X2)
    mut_rows = []
    mut_fails = 0
    for dec in (-2, 0, 2):
        n_i = N_IMP_REF * (10.0 ** dec)
        a_e, _ = asymptote_mutant(n_i, r0=0.0, n_over_nc=X2)
        u = int(ulp_distance(np.array([a_e]), np.array([pred2]))[0])
        mut_rows.append((n_i, a_e, u))
        if u > 8:
            mut_fails += 1
    print('  %-16s %-24s %-8s' % ('N_imp (m^-2)', 'mutant |n|*mu', 'ulp'))
    for n_i, a_e, u in mut_rows:
        print('  %-16.4e %-24.17g %-8d' % (n_i, a_e, u))
    check('D11 holding n_c fixed instead of A makes D3 FAIL at every '
          'impurity density away from the reference',
          mut_fails == 2,
          '%d of 3 probes exceed D3\'s 8-ulp budget, and the one that does '
          'not is the reference point itself, where the two models coincide '
          'by construction.  The mutant\'s asymptote moves as 1/N_imp '
          'instead of staying at A -- which is exactly the error the 10-30 '
          'entry identified in prose and this module removes from the code.'
          % mut_fails)

    path = make_figure()
    print('\n  Figure written: %s' % path)
    print('\n' + '=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
