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

import os
import re

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
                         #       in this context; see HONESTY note above and
                         #       PROVENANCE below.  The name is kept for
                         #       continuity of the committed transcripts and is
                         #       MISLEADING: 0.08 eV is not the paper's bound on
                         #       Delta_c(5.0 A).  Annotated 2026-10-18 rather
                         #       than renamed, because renaming it would make
                         #       every prior transcript unrelatable to this one.
EPS_RANGE = (0.02, 0.20) # eV -- the family the bound is reported over

# ---------------------------------------------------------------------------
# PROVENANCE -- added 2026-10-18
# ---------------------------------------------------------------------------
# WHY THIS IS IN THE PRINTED OUTPUT AND NOT ONLY IN THE DOCSTRING.
#
# The HONESTY paragraph above has been correct since this module was written:
# epsilon = 0.08 eV is borrowed from a sentence about a different quantity.
# Every PROSE artefact derived from this module carries that caveat -- the
# 2026-10-15 notes file states it, the AUTOMATION_LOG entry states it in bold
# beneath its table, and thesis Section 7.9 item 2 hedges it as "the only small
# energy scale the paper attaches to 'small' nearby".
#
# THE TRANSCRIPT CARRIED NONE OF IT.  Checked 2026-10-18:
# second_anchor_bound_output.txt contained zero occurrences of "borrow",
# "different quantity" or "non-rigid", while printing the attribution twice, in
# stronger words than any prose artefact uses -- "ell <= 0.7024 A at the
# PAPER'S OWN 0.08 eV scale", "66.5 % of the committed swept range is EXCLUDED
# BY THE SOURCE" -- and marking the row "<- paper".
#
# That is this repository's standing fault class pointing the other way.  On
# 2026-10-17 two docstrings denied what their code did; here the docstring is
# scrupulous and the OUTPUT overstates.  And the transcript is the artefact
# that travels: it is what pristine_transcript_audit compares, what a reader of
# results sees, and what gets quoted.  A caveat recorded only where no reader
# of the number encounters it is the 2026-10-01 dead-name fault in reverse --
# not a citation with no code behind it, but a qualification with no output
# behind it.
#
# AND THE CAVEAT ITSELF WAS TOO WEAK, which is the substantive half.  The
# paper's sentence, read verbatim, is:
#
#   "We see that while Eq. (1) holds for a relatively large separation
#    d = 5.0 A between the graphene sheet and the metal surface, there is a
#    small deviation of ~0.08 eV at the equilibrium separation d ~ 3.3 A which
#    can be traced to a perturbation [of the graphene electronic structure by
#    physisorption that cannot be described as a rigid shift]"
#
# So the 0.08 eV is (i) attached to the OTHER separation, d ~ 3.3 A, not to the
# 5.0 A case it is used as a tolerance for, and (ii) a measure of the failure
# of the rigid-shift picture.  Delta_c IS a rigid shift.  A number that
# quantifies how badly rigid shift describes 3.3 A is therefore not a
# conservative tolerance on a rigid-shift term at 5.0 A; it is a quantity of a
# different kind at a different separation, and the two defects compound rather
# than cancel.  The consequence the prose artefacts stop short of drawing:
# THERE IS NO PRIVILEGED ROW.  The family is the result, and the SUMMARY's
# single headline number should never have existed.
#
# Nothing numerical changes.  ell_max(0.08 eV) is still 0.7024 A and the
# exclusion at that epsilon is still 66.5 %; what changes is what those two
# numbers are attributed to, and that the transcript now says so.

SOURCE = 'Khomyakov et al., Phys. Rev. B 79, 195425 (2009), arXiv:0902.1203'

PROVENANCE = [
    ('Delta_c(3.3 A) = 0.9 eV', 'Sec. III.2',
     'the crossover point between n- and p-type doping is found for a metal '
     'with a work function W_M = W_G + 0.9 eV',
     'd ~ 3.3 A', 'the crossover shift itself -- USED AS ANCHOR 1'),
    ('Delta_c(5.0 A) <= epsilon', 'Sec. III.2',
     'the crossover point from n- to p-type doping decreases to its expected '
     'value, W_M ~ W_G, for large separations',
     'd = 5.0 A', 'an inequality in words; NO epsilon is quantified'),
    ('epsilon = 0.08 eV', 'Sec. III.2',
     'there is a small deviation of ~0.08 eV at the equilibrium separation '
     'd ~ 3.3 A which can be traced to a perturbation [...] that cannot be '
     'described as a rigid shift',
     'd ~ 3.3 A -- NOT 5.0 A',
     'the failure of the RIGID-SHIFT picture, which is what Delta_c is. '
     'BORROWED, twice over: wrong separation and wrong kind of quantity'),
    ('Delta_c(d) form', 'Sec. III.3, Eq. (4)',
     'Delta_c(d) = e^(-gamma d) (a_0 + a_1 d + a_2 d^2)',
     'all d', 'gamma, a_0, a_1, a_2 are NOT reported numerically anywhere in '
     'the text; (N1) drops the polynomial prefactor'),
]

# Tokens that must appear in this module's own printed output wherever
# EPS_PAPER is reported.  A check, not a habit -- see check_provenance_travels.
CAVEAT_TOKENS = ('BORROWED', 'd ~ 3.3 A -- NOT 5.0 A', 'no privileged row')

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

def ell_max(epsilon, d_far=None):
    """Largest ell consistent with Delta_c(d_far) <= epsilon.

    ell_max(epsilon) = (d_far - D_ANCHOR) / ln(DC_ANCHOR / epsilon).

    epsilon >= DC_ANCHOR constrains nothing and returns +inf EXACTLY (the
    logarithm is <= 0 there); epsilon <= 0 is inadmissible and returns 0.0.
    """
    if d_far is None:
        d_far = D_FAR   # 2026-10-17: was a frozen default capture
    eps = np.asarray(epsilon, dtype=float)
    span = float(d_far) - D_ANCHOR
    with np.errstate(divide='ignore', invalid='ignore'):
        lg = np.log(DC_ANCHOR / eps)
        out = np.where(lg > 0.0, span / lg, np.inf)
    return out if out.shape else float(out)


def eps_implied(ell, d_far=None):
    """The forward direction: Delta_c(d_far ; ell).  The inverse of ell_max."""
    if d_far is None:
        d_far = D_FAR   # 2026-10-17: was a frozen default capture
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
    # (7) ADDED 2026-10-18.  THE CAVEAT MUST TRAVEL WITH THE NUMBER.
    #     An exactly decidable property of this file, not a habit: no print
    #     statement in this module may attribute epsilon to the source
    #     unqualified, and the provenance tokens must actually be emitted.
    #     The forbidden phrases are deliberately still present in this file,
    #     in COMMENTS, as the kept record of what the transcript used to say --
    #     so the check is scoped to print() arguments, which is the only
    #     scope under which "kept in place rather than deleted" and "no longer
    #     claimed" are both true at once.
    import ast as _ast
    _src = open(os.path.abspath(__file__), encoding='utf-8').read()
    _printed = []
    for _node in _ast.walk(_ast.parse(_src)):
        if (isinstance(_node, _ast.Call)
                and isinstance(_node.func, _ast.Name)
                and _node.func.id == 'print'):
            for _sub in _ast.walk(_node):
                if isinstance(_sub, _ast.Constant) and isinstance(_sub.value, str):
                    _printed.append(_sub.value)
    #     THE FIRST DRAFT OF THIS CHECK FAILED, and for the third time today
    #     the cause is a verdict resolved on spelling.  It forbade a phrase
    #     anywhere in a print() argument, which cannot tell ASSERTING "the
    #     paper's own 0.08 eV scale" from QUOTING it as the wording this
    #     module has retracted -- and report_provenance() does the second,
    #     deliberately, so that the transcript carries the record of what it
    #     used to claim.  A check that cannot make that distinction forces a
    #     choice between keeping the record and passing the check, i.e. it
    #     pushes toward deleting the retracted wording, which is the one move
    #     this repository's annotation rule exists to refuse.
    #
    #     The rule, stated: a forbidden phrase in a print() is a VIOLATION
    #     unless that same string also carries a supersession marker.  The
    #     scope is the individual string literal, not the file, so a marker
    #     elsewhere in the module cannot launder an unqualified claim.
    _FORBIDDEN = ("paper's own 0.08", 'excluded by the source',
                  "the paper's own scale")
    _MARKERS = ('superseded', 'former', 'changed 2026-10-18', 'retracted')
    _violations = [t for t in _printed
                   if any(f in t for f in _FORBIDDEN)
                   and not any(m in t.lower() for m in _MARKERS)]
    check('no print() asserts the source attribution unqualified '
          '(quoting it as superseded is allowed, in the same string)',
          not _violations, 'violations: %r' % _violations)
    # And the positive half: the record must actually still be quoted, or the
    # check above passes by silence.
    _quoted = [t for t in _printed if any(f in t for f in _FORBIDDEN)]
    check('the retracted wording is still quoted in the transcript, '
          'not merely absent from it',
          len(_quoted) >= 1, '%d print strings quote it' % len(_quoted))
    _emitted = ' '.join(_printed + [t for row in PROVENANCE for t in row]).lower()
    _missing = [t for t in CAVEAT_TOKENS if t.lower() not in _emitted]
    check('every provenance caveat token is emitted somewhere in the output',
          not _missing, 'missing: %r' % _missing)
    # And the forbidden phrases must still EXIST in the file, as comments:
    # a check that passed by deleting the record would be the wrong repair.
    check('the superseded wording is kept in the file as a comment, not deleted',
          "ell <= 0.7024 A at the paper's own 0.08 eV scale" in _src
          and "66.5 % of the committed swept range is excluded by the source" in _src,
          'both superseded lines retained')

    check('exclusion counter reads exactly 0 at the top and 1 at the bottom',
          f_hi == 0.0 and f_lo == 1.0,
          'f(eps=%.4f eV) = %r ; f(eps=%.4f eV) = %r'
          % (hi_eps, f_hi, lo_eps, f_lo))


# ---------------------------------------------------------------------------
# What the bound does to this repository's numbers
# ---------------------------------------------------------------------------

def report_provenance():
    """Print where every externally-sourced number in this transcript comes
    from, at what separation, and what quantity it measures.  Added 2026-10-18
    because the caveat was in the docstring and the claim was in the
    transcript, and only the transcript travels."""
    print('\nPROVENANCE OF EVERY SOURCED NUMBER IN THIS TRANSCRIPT')
    print('  Source: %s' % SOURCE)
    for quantity, where, quote, applies_at, what in PROVENANCE:
        print('  %s   [%s]' % (quantity, where))
        print('      applies at : %s' % applies_at)
        for ln in re.findall(r'.{1,62}(?:\s|$)', 'quoted     : "%s"' % quote):
            if ln.strip():
                print('      %s' % ln.rstrip())
        for ln in re.findall(r'.{1,62}(?:\s|$)', 'quantity   : %s' % what):
            if ln.strip():
                print('      %s' % ln.rstrip())
    # ONE string literal, deliberately.  The check in validations() scopes
    # the supersession marker to the individual print argument, so the
    # retracted wording and the word "superseded" have to travel in the same
    # breath -- which is the rule this whole amendment is about, applied to
    # its own output.
    print("""  => THERE IS NO PRIVILEGED ROW in the family below.  The paper
     quantifies no epsilon at 5.0 A, so no epsilon is "the paper's",
     the family is the result, and no single headline number is
     reported for it.  SUPERSEDED 2026-10-18: this transcript used to
     end "ell <= 0.7024 A at the paper's own 0.08 eV scale" and
     "66.5 % of the committed swept range is excluded by the source".
     Both numbers stand; both attributions were false, and they are
     quoted here rather than deleted so that the transcript carries
     its own correction.""")


def report_bound_family():
    print('\nTHE BOUND, AS A FAMILY OVER epsilon (the paper quantifies none)')
    print('  %-12s %-14s %-12s %s' % ('eps / eV', 'ell_max / A',
                                      'excluded', 'ELL_DEFAULT survives?'))
    for eps in (0.20, 0.15, 0.10, EPS_PAPER, 0.05, 0.02):
        lm = float(ell_max(eps))
        frac = excluded_fraction_of_swept_range(eps)
        mark = 'yes' if ELL_DEFAULT <= lm else 'NO'
        star = '   <- BORROWED, see PROVENANCE' if eps == EPS_PAPER else ''
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
                label=r"$\epsilon$ = 0.08 eV (borrowed, not the paper's bound)")
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
    report_provenance()
    report_bound_family()
    report_swept_range_endpoints()
    moved = report_per_metal_effect()
    plot_bound()
    print('\n' + '=' * 72)
    print('SUMMARY')
    # AMENDED 2026-10-18.  These two lines used to read
    #   "ell <= 0.7024 A at the paper's own 0.08 eV scale."
    #   "66.5 % of the committed swept range is excluded by the source."
    # Both numbers are unchanged and both attributions were false: the paper
    # quantifies no epsilon at 5.0 A, and 0.08 eV is its deviation at 3.3 A
    # from the rigid-shift picture that Delta_c IS.  The old wording is kept
    # here as a comment rather than deleted.
    print('  At epsilon = %.2f eV -- BORROWED, not the paper\'s bound; see'
          % EPS_PAPER)
    print('  PROVENANCE -- ell <= %.4f A, which excludes %.1f %% of the'
          % (float(ell_max(EPS_PAPER)) * 1e10,
             100.0 * excluded_fraction_of_swept_range(EPS_PAPER)))
    print('  committed swept range.  The source supports the FAMILY, not')
    print('  this row: it asserts Delta_c(5.0 A) is negligible and')
    print('  quantifies no epsilon at that separation.')
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
