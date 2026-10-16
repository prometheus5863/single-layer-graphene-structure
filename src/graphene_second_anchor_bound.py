"""
graphene_second_anchor_bound.py

Chapter 7's TOP PHYSICS ITEM, open since 2026-09-21 and untouched for
twenty-four consecutive sessions:

    "A second anchor for Delta_c, at any separation other than 3.3 A."

The item was written expecting a NUMBER -- a second (d, Delta_c) pair that
would turn the one-parameter form

    Delta_c(d ; ell) = 0.9 eV * exp( -(d - 3.3 A) / ell )              (N1)

of graphene_per_metal_crossover_model into a determined curve with no swept
parameter.  Today's literature pass found that the source does not contain
such a pair, and that what it does contain is better than nothing and worse
than a number:

    KHOMYAKOV ET AL. STATE THE SECOND ANCHOR AS AN INEQUALITY, IN WORDS.

Khomyakov, Giovannetti, Rusu, Brocks, van den Brink, Kelly,
"First-principles study of the interaction and charge transfer between
graphene and metals", Phys. Rev. B 79, 195425 (2009), arXiv:0902.1203 --
the same paper (N1)'s single anchor comes from.  Section III.2, on the
large-separation calculation performed at d = 5.0 A:

    "the crossover point from n- to p-type doping decreases to its expected
     value, W_M ~ W_G, for large separations"

and Section III.3, on the short-range term Delta_c(d) = e^(-gamma d)(a_0 +
a_1 d + a_2 d^2) (their Eq. 4):

    "|Delta_c(d)| becomes negligible for large d anyway"

The four fitted coefficients (gamma, a_0, a_1, a_2) are NOT reported
numerically anywhere in the paper's text -- which is why 2026-09-21 built
(N1) in the first place, and why twenty-four sessions of hoping for a
number were hoping for something that is not there.  What IS there is a
second separation, d = 5.0 A, at which the paper asserts Delta_c is
negligible.  Since W_0(d) = W_G + Delta_c(d), "crossover ~ W_G at 5.0 A" is
exactly the statement Delta_c(5.0 A) <= epsilon for some small epsilon the
paper does not quantify.

So the second anchor converts (N1)'s FREE parameter into a BOUNDED one:

    Delta_c(5.0 A ; ell) <= epsilon
      <=>  0.9 * exp( -(5.0 - 3.3) / ell ) <= epsilon
      <=>  ell <= ell_max(epsilon) = (5.0 A - 3.3 A) / ln(0.9 eV / epsilon)

This module computes that bound, reports it as a FAMILY over epsilon rather
than at one chosen epsilon (the paper quantifies no epsilon, and pretending
otherwise would manufacture precision the source does not carry), and asks
the only question that matters for this repository: DOES THE BOUND MOVE ANY
COMMITTED NUMBER?

RESULT.  At epsilon = 0.08 eV -- the only small energy scale the paper
itself quotes in this context, its "small deviation of ~0.08 eV at the
equilibrium separation" -- ell <= 0.7024 A.  The committed sweep range is
ELL_RANGE = [0.3, 1.5] A, so the source's own large-separation statement
EXCLUDES 66.5 % of the range this repository has been sweeping, and the
committed default ELL_DEFAULT = 0.6 A SURVIVES.  That it survives was not
guaranteed and is not an assertion: it is checked here, and the top of the
swept range is shown to put the 5.0 A crossover 0.29 eV above W_G, which is
not "~ W_G" by any reading.

HONESTY ABOUT WHAT THIS IS.  epsilon = 0.08 eV is borrowed from a sentence
about a DIFFERENT quantity (the non-rigid-shift correction), because it is
the only number the paper attaches to the word "small" nearby.  The bound is
therefore reported over epsilon in [0.02, 0.20] eV throughout, and no
committed number is changed on the strength of one borrowed epsilon.  The
bound is one-sided: nothing here constrains ell from BELOW, and the item
stays open for a two-sided constraint.

Validation policy follows 2026-09-17: at least one check against an exactly
known value.  There are six, including a planted-ell round trip that the
inversion must recover and two degenerate limits that must be exact.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_per_metal_crossover_model as pmc
from graphene_per_metal_crossover_model import (
    D_ANCHOR, DC_ANCHOR, ELL_DEFAULT, ELL_RANGE, delta_c, w_cross_for_metal,
)
from graphene_contact_doping_model import METAL_WORK_FUNCTIONS, W_GRAPHENE
from graphene_contact_doping_nonlinear_model import D_EQ

# ---------------------------------------------------------------------------
# The second anchor, as the source states it
# ---------------------------------------------------------------------------
D_FAR = 5.0e-10          # m  -- the separation of Khomyakov et al.'s
                         #       large-separation calculation (their Sec. III.2
                         #       and Fig. 4's upper point set)
EPS_PAPER = 0.08         # eV -- the only small energy scale the paper quotes
                         #       in this context; see HONESTY note above
EPS_RANGE = (0.02, 0.20) # eV -- the family the bound is reported over

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


def ulps(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    return abs(a - b) / np.spacing(max(abs(a), abs(b)))


# ---------------------------------------------------------------------------
# The bound
# ---------------------------------------------------------------------------

def ell_max(epsilon, d_far=D_FAR):
    """Largest ell consistent with Delta_c(d_far) <= epsilon.

    ell_max(epsilon) = (d_far - D_ANCHOR) / ln(DC_ANCHOR / epsilon).

    epsilon >= DC_ANCHOR constrains nothing and returns +inf EXACTLY (the
    logarithm is <= 0 there); epsilon <= 0 is inadmissible and returns 0.0.
    """
    eps = np.asarray(epsilon, dtype=float)
    span = float(d_far) - D_ANCHOR
    with np.errstate(divide='ignore', invalid='ignore'):
        lg = np.log(DC_ANCHOR / eps)
        out = np.where(lg > 0.0, span / lg, np.inf)
    return out if out.shape else float(out)


def eps_implied(ell, d_far=D_FAR):
    """The forward direction: Delta_c(d_far ; ell).  The inverse of ell_max."""
    return float(delta_c(float(d_far), float(ell)))


def excluded_fraction_of_swept_range(epsilon):
    """Fraction of the committed ELL_RANGE that the bound rules out."""
    lo, hi = ELL_RANGE
    lm = float(ell_max(epsilon))
    if lm >= hi:
        return 0.0
    if lm <= lo:
        return 1.0
    return (hi - lm) / (hi - lo)


# ---------------------------------------------------------------------------
# Validations against exactly known values
# ---------------------------------------------------------------------------

def validations():
    print('VALIDATIONS (exactly known values, not plausible ranges)')

    # (1) A bound at epsilon == the anchor itself constrains nothing: the
    #     logarithm is exactly zero and ell_max must be exactly +inf, not a
    #     large finite number.
    v = float(ell_max(DC_ANCHOR))
    check('epsilon = DC_ANCHOR gives ell_max = +inf EXACTLY',
          v == np.inf, 'ell_max(0.9 eV) = %r' % v)

    # (2) PLANTED-ELL ROUND TRIP.  Pick an ell, compute what Delta_c(5.0 A)
    #     it implies, feed that back as epsilon.  The inversion must return
    #     the planted ell.  This is the check that the algebra is right and
    #     not merely plausible.
    worst = 0.0
    worst_at = None
    for ell_true in np.linspace(0.2e-10, 1.8e-10, 65):
        back = float(ell_max(eps_implied(ell_true)))
        u = ulps(back, float(ell_true))
        if u > worst:
            worst, worst_at = u, float(ell_true)
    check('planted-ell round trip over 65 values, worst error <= 4 ulp',
          worst <= 4.0,
          'worst %.3g ulp at ell = %.4f A' % (worst, (worst_at or 0.0) * 1e10))

    # (3) The bound must not touch the anchor.  Delta_c(3.3 A, ell) == 0.9
    #     BITWISE for every ell in the SURVIVING window, not merely in the
    #     old swept range (which 2026-09-21 already checked).
    lm = float(ell_max(EPS_PAPER))
    ells = np.linspace(ELL_RANGE[0], lm, 61)
    oks = [float(delta_c(D_ANCHOR, e)) == DC_ANCHOR for e in ells]
    check('anchor untouched across the surviving window, bitwise',
          all(oks), 'Delta_c(3.3 A, ell) == 0.9 for %d/%d ell in [%.3f, %.4f] A'
          % (sum(oks), len(oks), ELL_RANGE[0] * 1e10, lm * 1e10))

    # (4) DEGENERATE LIMIT: as the far separation collapses onto the anchor,
    #     the bound becomes vacuous in the other direction -- ell_max -> 0
    #     exactly, because no ell can make Delta_c fall below epsilon over
    #     zero distance.
    v = float(ell_max(EPS_PAPER, d_far=D_ANCHOR))
    check('d_far = D_ANCHOR gives ell_max = 0 EXACTLY',
          v == 0.0, 'ell_max at zero span = %r' % v)

    # (5) STRICT MONOTONICITY in epsilon, as an exact ordering (a looser
    #     bound can never admit fewer ell).
    eps = np.linspace(EPS_RANGE[0], EPS_RANGE[1], 97)
    lms = np.array([float(ell_max(x)) for x in eps])
    check('ell_max strictly increasing in epsilon over 97 values',
          bool(np.all(np.diff(lms) > 0.0)),
          'ell_max spans %.4f -> %.4f A' % (lms[0] * 1e10, lms[-1] * 1e10))

    # (6) POSITIVE CONTROL ON THE EXCLUSION ARITHMETIC.  A bound placed
    #     exactly at the top of the swept range must exclude exactly 0, and
    #     one placed exactly at the bottom exactly 1 -- so a zero exclusion
    #     reported below is a measurement and not a broken counter.
    hi_eps = eps_implied(ELL_RANGE[1])
    lo_eps = eps_implied(ELL_RANGE[0])
    f_hi = excluded_fraction_of_swept_range(hi_eps)
    f_lo = excluded_fraction_of_swept_range(lo_eps)
    check('exclusion counter reads exactly 0 at the top and 1 at the bottom',
          f_hi == 0.0 and f_lo == 1.0,
          'f(eps=%.4f eV) = %r ; f(eps=%.4f eV) = %r'
          % (hi_eps, f_hi, lo_eps, f_lo))


# ---------------------------------------------------------------------------
# What the bound does to this repository's numbers
# ---------------------------------------------------------------------------

def report_bound_family():
    print('\nTHE BOUND, AS A FAMILY OVER epsilon (the paper quantifies none)')
    print('  %-12s %-14s %-12s %s' % ('eps / eV', 'ell_max / A',
                                      'excluded', 'ELL_DEFAULT survives?'))
    for eps in (0.20, 0.15, 0.10, EPS_PAPER, 0.05, 0.02):
        lm = float(ell_max(eps))
        frac = excluded_fraction_of_swept_range(eps)
        mark = 'yes' if ELL_DEFAULT <= lm else 'NO'
        star = '   <- paper' if eps == EPS_PAPER else ''
        print('  %-12.3f %-14.4f %-12s %s%s'
              % (eps, lm * 1e10, '%.1f %%' % (100.0 * frac), mark, star))


def report_swept_range_endpoints():
    print('\nWHAT THE SWEPT RANGE IMPLIES AT 5.0 A (the forward direction)')
    print('  %-16s %-18s %s' % ('ell / A', 'Delta_c(5.0 A) / eV',
                                'crossover W_0(5.0 A) / eV'))
    for ell, tag in ((ELL_RANGE[0], 'bottom of range'),
                     (ELL_DEFAULT, 'COMMITTED DEFAULT'),
                     (float(ell_max(EPS_PAPER)), 'the new bound'),
                     (ELL_RANGE[1], 'top of range')):
        dc = eps_implied(ell)
        print('  %-16s %-18.4f %.4f    (%s)'
              % ('%.4f' % (ell * 1e10), dc, W_GRAPHENE + dc, tag))
    print('  W_G = %.2f eV.  "crossover ~ W_G at 5.0 A" is the source\'s own'
          % W_GRAPHENE)
    print('  description of the large-separation limit; the top of the swept')
    print('  range puts it %.4f eV above W_G, which is not that.'
          % eps_implied(ELL_RANGE[1]))


def report_per_metal_effect():
    """Does restricting ell to the surviving window move any metal's verdict?

    This is the question the item exists to answer.  A bound that excludes
    two thirds of a swept parameter but changes no conclusion is a
    tightening of the error bar, not a correction -- and saying which it is
    requires measuring, not asserting.
    """
    lm = float(ell_max(EPS_PAPER))
    print('\nPER-METAL EFFECT OF RESTRICTING ell TO (0, %.4f A]' % (lm * 1e10))
    print('  %-6s %-8s %-11s %-11s %-11s %s'
          % ('metal', 'd_eq/A', 'W_0 @0.3A', 'W_0 @bound', 'W_0 @1.5A',
             'verdict moves?'))
    moved = []
    for m in sorted(METAL_WORK_FUNCTIONS, key=lambda k: METAL_WORK_FUNCTIONS[k]):
        d = D_EQ.get(m)
        if d is None:
            continue
        w_lo, st = w_cross_for_metal(m, ELL_RANGE[0])
        w_bd, _ = w_cross_for_metal(m, lm)
        w_hi, _ = w_cross_for_metal(m, ELL_RANGE[1])
        wm = METAL_WORK_FUNCTIONS[m]
        # The verdict is the SIGN of (W_M - W_0): positive = p-type doping.
        sgn_bd = np.sign(wm - w_bd)
        sgn_hi = np.sign(wm - w_hi)
        changed = (sgn_bd != sgn_hi)
        if changed:
            moved.append(m)
        print('  %-6s %-8.2f %-11.4f %-11.4f %-11.4f %s%s'
              % (m, d * 1e10, w_lo, w_bd, w_hi,
                 'YES' if changed else 'no',
                 '' if st == 'ok' else '  [%s]' % st))
    if moved:
        print('  Verdict moves for: %s' % ', '.join(moved))
    else:
        print('  No metal\'s n/p verdict changes between the old top of the')
        print('  swept range and the new bound.  The bound tightens the error')
        print('  bar on W_0 and leaves every doping sign where it was.')
    return moved


def plot_bound():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    eps = np.linspace(EPS_RANGE[0], EPS_RANGE[1], 400)
    lms = np.array([float(ell_max(x)) for x in eps]) * 1e10
    ax1.plot(eps, lms, color='navy', lw=2.0,
             label=r'$\ell_{\max}(\epsilon)$ from $\Delta_c(5.0\,$\AA$)\leq\epsilon$')
    ax1.axhspan(ELL_RANGE[0] * 1e10, ELL_RANGE[1] * 1e10, color='grey',
                alpha=0.18, label='swept range [0.3, 1.5] '+ r'\AA')
    ax1.axhline(ELL_DEFAULT * 1e10, color='seagreen', ls='--', lw=1.6,
                label=r'committed default 0.6 \AA')
    ax1.axvline(EPS_PAPER, color='crimson', ls=':', lw=1.6,
                label=r"paper's own 0.08 eV scale")
    ax1.plot([EPS_PAPER], [float(ell_max(EPS_PAPER)) * 1e10], 'o',
             color='crimson', ms=7, zorder=5)
    ax1.set_xlabel(r'$\epsilon$  (eV)  --  how small "negligible" is taken to be')
    ax1.set_ylabel(r'$\ell_{\max}$  (\AA)')
    ax1.set_title('The second anchor is a one-sided bound')
    ax1.legend(fontsize=8, loc='lower right')
    ax1.grid(alpha=0.3)

    d = np.linspace(2.0, 5.6, 400) * 1e-10
    for ell, col, ls, lab in (
            (ELL_RANGE[0], 'steelblue', '-', r'$\ell=0.3$ \AA'),
            (ELL_DEFAULT, 'seagreen', '-', r'$\ell=0.6$ \AA (default)'),
            (float(ell_max(EPS_PAPER)), 'crimson', '--',
             r'$\ell=0.702$ \AA (the bound)'),
            (ELL_RANGE[1], 'darkorange', ':', r'$\ell=1.5$ \AA (EXCLUDED)')):
        ax2.plot(d * 1e10, delta_c(d, ell), color=col, ls=ls, lw=1.9, label=lab)
    ax2.plot([D_ANCHOR * 1e10], [DC_ANCHOR], 'k*', ms=14, zorder=6,
             label=r'anchor 1: $\Delta_c(3.3\,$\AA$)=0.9$ eV')
    ax2.plot([D_FAR * 1e10], [EPS_PAPER], 'kv', ms=10, zorder=6,
             label=r'anchor 2: $\Delta_c(5.0\,$\AA$)\leq\epsilon$')
    ax2.annotate('', xy=(D_FAR * 1e10, 0.0), xytext=(D_FAR * 1e10, EPS_PAPER),
                 arrowprops=dict(arrowstyle='-|>', color='k', lw=1.2))
    ax2.set_xlabel(r'graphene--metal separation $d$  (\AA)')
    ax2.set_ylabel(r'$\Delta_c$  (eV)')
    ax2.set_ylim(-0.03, 1.0)
    ax2.set_title('One number and one inequality')
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.suptitle('The second anchor for $\\Delta_c$: Khomyakov et al. state it '
                 'at $d=5.0$ \\AA\\ as an inequality, not a number',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig('second_anchor_bound.png', dpi=150)
    plt.close(fig)
    print('\nWrote second_anchor_bound.png')


def main():
    plt.rcParams['text.usetex'] = False
    print('=' * 72)
    print('THE SECOND ANCHOR FOR Delta_c -- found, and it is a bound')
    print('=' * 72)
    print('Source: Khomyakov et al., Phys. Rev. B 79, 195425 (2009),')
    print('        arXiv:0902.1203, Sections III.2 and III.3.')
    print('Anchor 1: Delta_c(%.1f A) = %.1f eV      (used since 2026-09-21)'
          % (D_ANCHOR * 1e10, DC_ANCHOR))
    print('Anchor 2: Delta_c(%.1f A) <= epsilon     (NEW -- stated in words)'
          % (D_FAR * 1e10))
    print()
    validations()
    report_bound_family()
    report_swept_range_endpoints()
    moved = report_per_metal_effect()
    plot_bound()
    print('\n' + '=' * 72)
    print('SUMMARY')
    print('  ell <= %.4f A at the paper\'s own 0.08 eV scale.' % (float(ell_max(EPS_PAPER)) * 1e10))
    print('  %.1f %% of the committed swept range is excluded by the source.'
          % (100.0 * excluded_fraction_of_swept_range(EPS_PAPER)))
    print('  The committed default ELL_DEFAULT = %.1f A survives (checked, not assumed).'
          % (ELL_DEFAULT * 1e10))
    print('  Metals whose n/p verdict moves: %s'
          % (', '.join(moved) if moved else 'none'))
    print('  The item stays OPEN for a two-sided constraint: nothing here')
    print('  bounds ell from below.')
    print('  Validations: %d passed, %d failed' % (passed, failed))
    print('=' * 72)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
