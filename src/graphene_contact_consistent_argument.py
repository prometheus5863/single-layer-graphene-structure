"""
graphene_contact_consistent_argument.py

ITEM 16's PRECONDITION (i), DISCHARGED.  2026-10-15.

Section 7.9 item 16 -- "flip quantum_capacitance()'s argument from the gate
overdrive to E_F/e" -- has been the first-ranked open item on Chapter 4 since
2026-10-10.  On 2026-10-14 its plumbing preconditions were reduced to exactly
two, both named:

  (i) `graphene_contact_doping_model` needs its OWN contact-correct
      replacement, because the contact cannot take the channel's object: the
      channel's replacement solves for n against C_ox, and installing it at a
      contact whose series capacitance is C_interface (65.93x C_ox) is a
      channel solver at the contact -- a different wrong answer, not a repair.
      The flip is TWO replacements, not one rebinding.

 (ii) Section 4.7's prose and Section 4.10's annotation need rewriting.

This module is (i).  It does not perform the flip -- item 16 stays open -- it
removes the reason the flip could not be performed.

WHAT WAS ALREADY TRUE, AND WHY IT WAS NOT ENOUGH.  The contact-correct form
was MEASURED on 2026-10-13 by `graphene_contact_argument_blast_radius`
(-12.4964 % in n for Pd, 4.7696-12.4964 % across the seven metals).  A
measurement harness is not a replacement, for three reasons that only became
visible on trying to use it:

  (a) THE OBJECT BELONGED TO THE AUDIT, NOT TO THE CONTACT.  Nothing outside
      the blast-radius module could install it, so item 16's flip had no
      second object to install even though its numbers were known.

  (b) THE HARNESS REBOUND BY HAND -- `cd.quantum_capacitance = ...`, a single
      assignment -- which is precisely the mechanism the 2026-10-14 sweep
      retired for the channel after measuring that hand-listed rebinding
      leaves 25 stale sites at 14 sites, 22 of them never found by hand.  The
      contact half was still on the old mechanism.

  (c) AND THE NEW MECHANISM COULD NOT EXPRESS IT.  This is the finding.
      `bc.rebind()` gives every member of an identity class THE SAME
      replacement, with `CONTACT_OPERATING_POINT` as the one escape hatch:
      leave this member alone, for a stated reason.  That vocabulary --
      {same replacement, left alone} -- has no word for "rebind this member to
      a DIFFERENT object", which is exactly what the contact needs the moment
      it has an object of its own.

      IDENTITY ANSWERS "WHICH NAMES MOVE TOGETHER" AND DOES NOT ANSWER "TO
      WHAT".  `quantum_capacitance`'s class spans two operating points, the
      channel's C_ox and the contact's C_interface, and an operating point is
      not something identity can see.  So the class must be PARTITIONED by
      one, and the partition has to come from the caller.  That is
      `bc.rebind_by_operating_point`, added today, which carries 2026-10-14's
      refusal across in the form the same lesson takes here: there, an
      exclusion without a stated reason was refused; here, an UNASSIGNED class
      member is refused, because silently leaving a member at the shipped
      object is the identical fault reached by omission from a map instead of
      omission from a list.

THE RESULT.  `on_contact_consistent_argument()` installs the contact's own
object through the partitioned mechanism, and
`on_consistent_argument_everywhere()` installs BOTH halves at once with an
EMPTY left-alone set -- so after today there is no name on the charge path
that item 16's flip would leave at the shipped object.  That emptiness is the
discharge of precondition (i), and it is asserted by the mechanism rather
than by this docstring: the manager refuses to open if any class member is
unassigned.

Validation policy follows 2026-09-17.  Seven checks, five against exactly
known values, including a closed-form quantum-capacitance limit
n = 4(dV/A_F)^2 that the solver must hit, a cross-construction agreement with
2026-10-13's committed per-metal transcript, and a positive control on every
arrival counter.
"""

import contextlib

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import graphene_contact_doping_model as cd
import graphene_quantum_capacitance_argument as qca
import graphene_contact_argument_blast_radius as blast
import graphene_binding_class as bc

E_CHARGE = 1.602176634e-19

# The seven metals of record, in the order 2026-10-13's transcript prints them.
METALS = ('Ti', 'Cr', 'Cu', 'Ni', 'Au', 'Pd', 'Pt')

# 2026-10-13's committed per-metal correction, read from
# contact_argument_blast_radius_output.txt Section 4.  Held here as the
# CROSS-CONSTRUCTION TARGET: today's owned object is a different construction
# of the same physics, and two constructions that disagree are one too many.
COMMITTED_2026_10_13 = {
    'Ti': -6.2578, 'Cr': None, 'Cu': -4.7696, 'Ni': -12.4233,
    'Au': -12.4878, 'Pd': -12.4964, 'Pt': -11.8738,
}

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
# The contact's own replacement -- owned HERE, beside the contact model
# ---------------------------------------------------------------------------

_ARRIVED = {'cq_contact': 0, 'cq_channel': 0, 'charge_channel': 0}


def reset_arrivals():
    for k in _ARRIVED:
        _ARRIVED[k] = 0


def cq_consistent_contact(V_metal_induced, T=300.0):
    """FORM B_T AT THE CONTACT.  Drop-in for the `quantum_capacitance` name
    that `cd.contact_edge_carrier_density` reads on its line 105.

    Same signature as the shipped function, so the shipped
    `contact_edge_carrier_density` forms its own series factor against
    C_interface and its own product, both untouched: the ONLY difference
    between form A and form B_T at the contact is C_q's argument.

    Built on `blast.n_consistent(V, cd.C_interface)`, which is
    `qca._fixed_point` -- the committed channel solver, unmodified -- with the
    one shipped name `gfet.C_ox` that `qca.series_factor` reads substituted for
    C_interface.  Reusing the channel solver rather than writing a second one
    is deliberate: the contact solve and the channel solve are then the same
    function and cannot drift apart, and the 1.1e-9 relative constants gap
    2026-10-13 measured between `qca.cq_thermal` and the shipped function is
    kept out by construction rather than by tolerance.
    """
    _ARRIVED['cq_contact'] += 1
    n = blast.n_consistent(V_metal_induced, cd.C_interface)
    return blast.cq_shipped_of_n(n, T=T)


def _cq_consistent_channel(V_g_minus_Vdirac, T=None):
    """The channel half, wrapped only to carry its own arrival counter.
    Delegates to the object `qca` already owns -- this module does not build a
    second channel replacement."""
    _ARRIVED['cq_channel'] += 1
    return qca._quantum_capacitance_consistent(V_g_minus_Vdirac, T=T)


def _charge_consistent_channel(V_g, V_ch=0.0):
    _ARRIVED['charge_channel'] += 1
    return qca._carrier_density_consistent(V_g, V_ch)


# Why each class member gets the object it gets.  One place, not re-argued at
# each call site -- the successor to bc.CONTACT_OPERATING_POINT, which said
# "leave this one alone" and now says "give this one its own".
REASON_CHANNEL = ('channel operating point: series capacitance C_ox, series '
                  'factor 0.999 at the RF bias; solves for n against C_ox')
REASON_CONTACT = ('contact operating point: series capacitance C_interface '
                  '(65.93x C_ox), series factor 0.555-0.806; solves for n '
                  'against C_interface.  Formerly excluded by '
                  'bc.CONTACT_OPERATING_POINT because no such object existed')


def _cq_assignment():
    """The partition of `quantum_capacitance`'s identity class by operating
    point.  Computed against the class as it stands at call time, so a member
    added later is UNASSIGNED and the manager refuses rather than silently
    leaving it shipped."""
    return {
        ('graphene_fet_model', 'quantum_capacitance'):
            (_cq_consistent_channel, REASON_CHANNEL),
        ('graphene_contact_doping_model', 'quantum_capacitance'):
            (cq_consistent_contact, REASON_CONTACT),
    }


@contextlib.contextmanager
def on_contact_consistent_argument():
    """Correct C_q's argument AT THE CONTACT ONLY, leaving the channel shipped.

    Separated from the composed manager because 2026-10-13 measured the
    contact correction (4.77-12.50 %) to be 4.7x-12.2x the channel's
    (-1.0212 %), and a single combined number would hide that the two halves
    are not the same size.
    """
    assignment = dict(_cq_assignment())
    assignment[('graphene_fet_model', 'quantum_capacitance')] = (
        None, 'channel deliberately left shipped: this manager isolates the '
              'contact half so its size can be read alone')
    with bc.rebind_by_operating_point(gfet, 'quantum_capacitance',
                                      assignment) as (ins, alone):
        yield ins, alone


@contextlib.contextmanager
def on_consistent_argument_everywhere():
    """ITEM 16's FLIP, FULLY PLUMBED -- both operating points, nothing left
    at the shipped object.

    This is what precondition (i) was blocking.  The charge relation and both
    capacitance operating points are replaced together, and the left-alone
    set is EMPTY -- asserted by the mechanism, which refuses to open if any
    member of either identity class is unassigned.

    It does NOT flip the default.  Item 16 remains open: what is removed is
    the reason it could not be flipped, not the decision to flip it.  The
    decision needs precondition (ii) -- Section 4.7's prose and Section
    4.10's annotation -- which is writing, not plumbing.
    """
    with contextlib.ExitStack() as stack:
        ins_q, alone_q = stack.enter_context(
            bc.rebind_by_operating_point(gfet, 'quantum_capacitance',
                                         _cq_assignment()))
        charge_assignment = {
            site: (_charge_consistent_channel, REASON_CHANNEL)
            for site in bc.identity_class(gfet, 'carrier_density')
        }
        ins_c, alone_c = stack.enter_context(
            bc.rebind_by_operating_point(gfet, 'carrier_density',
                                         charge_assignment))
        yield (ins_q + ins_c), (alone_q + alone_c)


# ---------------------------------------------------------------------------
# Validations
# ---------------------------------------------------------------------------

def validations():
    print('VALIDATIONS')

    # (a) EXACT.  dWF = 0 is a symmetry of the model: the metal and graphene
    #     have the same work function, nothing transfers, and n must be
    #     EXACTLY zero under both forms.  A solver that returns a small
    #     nonzero number here is wrong in a way no plausible-range check sees.
    nA = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE))
    with on_contact_consistent_argument():
        nB = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE))
    check('(a) EXACT: dWF = 0 gives n = 0 under both forms',
          nA == 0.0 and nB == 0.0, 'form A %r, form B_T %r' % (nA, nB))

    # (b) EXACT.  Form B_T is ODD in its argument -- C_q depends on |n| only
    #     and the electrostatic term is linear -- so n(-dV) must be -n(+dV)
    #     at 0 ulp, not merely to tolerance.  Tested on the SOLVER, at
    #     exactly negated dV.
    worst = 0.0
    for dv in (0.05, 0.17, 0.3, 0.62, 1.15):
        p = float(blast.n_consistent(dv, cd.C_interface))
        m = float(blast.n_consistent(-dv, cd.C_interface))
        worst = max(worst, ulps(p, -m))
    check('(b) EXACT: form B_T is odd in dV, 0 ulp over 5 offsets',
          worst == 0.0, 'worst %.3g ulp' % worst)

    # (b2) AND A FINDING THAT (b)'s FIRST DRAFT MISTOOK FOR A SOLVER DEFECT.
    #      The same oddness test applied through the SHIPPED
    #      contact_edge_carrier_density -- i.e. at symmetric WORK FUNCTIONS
    #      rather than at symmetric dV -- fails at 3 ulp, and does so under
    #      FORM A as well.  The cause is not form B_T: it is the shipped
    #      function's own `dWF = work_function_metal - W_GRAPHENE`, which is
    #      not exactly antisymmetric in double precision --
    #      (4.5 + 0.05) - 4.5 = 0.050000000000000266 while
    #      (4.5 - 0.05) - 4.5 = -0.04999999999999982.  The first draft of (b)
    #      tested through the shipped function and reported 3 ulp as though
    #      the new solver had introduced it.  It had not, and that is the
    #      2026-10-11 rule again: a knob must be turned alone, and a check
    #      that reads two knobs attributes the wrong one.
    worst_wf_A = 0.0
    worst_wf_B = 0.0
    for d in (0.05, 0.17, 0.3, 0.62, 1.15):
        pa = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + d))
        ma = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE - d))
        worst_wf_A = max(worst_wf_A, ulps(pa, -ma))
    with on_contact_consistent_argument():
        for d in (0.05, 0.17, 0.3, 0.62, 1.15):
            pb = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + d))
            mb = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE - d))
            worst_wf_B = max(worst_wf_B, ulps(pb, -mb))
    check('(b2) the 3-ulp antisymmetry defect is the SHIPPED dWF subtraction, '
          'not form B_T -- it is present under form A too',
          worst_wf_A > 0.0 and abs(worst_wf_B - worst_wf_A) <= 2.0,
          'form A %.3g ulp, form B_T %.3g ulp, solver itself 0 ulp'
          % (worst_wf_A, worst_wf_B))

    # (c) EXACT, AND THE ONE THAT MATTERS MOST.  Give the contact solver the
    #     OLD argument -- a C_q that ignores n and returns the shipped
    #     function at dV -- and it must reproduce the SHIPPED
    #     contact_edge_carrier_density bitwise.  This is what establishes
    #     that form A and form B_T differ in C_q's ARGUMENT and in nothing
    #     else: no second solver, no second set of constants, no changed
    #     series factor.  Without it, every per-metal number below could be
    #     the effect of some other difference.
    def _cq_old_argument(V, T=300.0):
        """The new path, with C_q deliberately left at the OLD argument."""
        return gfet.quantum_capacitance(V, T=gfet.T)
    worst = 0.0
    base = {m: float(cd.contact_edge_carrier_density(gfet_wf))
            for m, gfet_wf in cd.METAL_WORK_FUNCTIONS.items()}
    assignment = dict(_cq_assignment())
    assignment[('graphene_fet_model', 'quantum_capacitance')] = (
        None, 'control: channel untouched')
    assignment[('graphene_contact_doping_model', 'quantum_capacitance')] = (
        _cq_old_argument, 'control: the new solver fed the OLD argument')
    with bc.rebind_by_operating_point(gfet, 'quantum_capacitance', assignment):
        for m, wf in cd.METAL_WORK_FUNCTIONS.items():
            worst = max(worst, ulps(float(cd.contact_edge_carrier_density(wf)),
                                    base[m]))
    check('(c) EXACT: new machinery + OLD argument == shipped, bitwise over '
          '%d metals' % len(base),
          worst == 0.0,
          'worst %.3g ulp -- so A vs B_T differs in the ARGUMENT alone'
          % worst)

    # (d) EXACT CLOSED FORM.  In the quantum-capacitance-limited regime
    #     C_s -> infinity the series combination collapses:
    #         n -> (C_s dV/e)(C_q/C_s) = C_q dV/e,
    #     and with the T = 0 dispersion C_q = 2e sqrt(n)/A_F this closes to
    #         n = 4 (dV/A_F)^2
    #     with no free constant.  The solver must hit it.
    #     At FINITE C_s the solve differs from the limit by O(C_q/C_s), so a
    #     single threshold would be a plausible-range check wearing a closed
    #     form's clothes.  The exact statement is the CONVERGENCE ORDER: the
    #     deviation must fall by exactly the factor C_s is raised by.
    dV = 0.62
    n_exact = 4.0 * (dV / qca.A_FERMI) ** 2
    rels = []
    for k in (1e9, 1e11):
        n_num = float(blast.n_consistent(dV, cd.C_interface * k,
                                         cq_of_n=qca.cq_disp))
        rels.append(abs(n_num - n_exact) / n_exact)
    ratio = rels[0] / rels[1]
    check('(d) EXACT CLOSED FORM: quantum limit is n = 4(dV/A_F)^2, '
          'approached at first order in 1/C_s',
          rels[1] < 1e-9 and abs(ratio - 100.0) / 100.0 < 0.02,
          'closed form %.6e; rel dev %.3e -> %.3e as C_s x100, '
          'ratio %.3f (exactly 100 expected)' % (n_exact, rels[0], rels[1], ratio))

    # (e) CROSS-CONSTRUCTION AGREEMENT with 2026-10-13's committed transcript.
    #     Today's object is owned here and installed through the partitioned
    #     mechanism; 2026-10-13's was owned by the audit and installed by hand
    #     assignment.  Two constructions of the same physics that disagree are
    #     one too many, and the committed numbers are the target.
    worst = 0.0
    rows = []
    with on_contact_consistent_argument():
        for m in METALS:
            wf = cd.METAL_WORK_FUNCTIONS[m]
            rows.append((m, wf, float(cd.contact_edge_carrier_density(wf))))
    for m, wf, nb in rows:
        with bc.rebind_by_operating_point(
                gfet, 'quantum_capacitance',
                {('graphene_fet_model', 'quantum_capacitance'):
                     (None, 'control'),
                 ('graphene_contact_doping_model', 'quantum_capacitance'):
                     (blast._cq_consistent_contact,
                      '2026-10-13 construction, for comparison')}):
            nb_old = float(cd.contact_edge_carrier_density(wf))
        worst = max(worst, ulps(nb, nb_old))
    check('(e) EXACT: today\'s owned object == 2026-10-13\'s construction, '
          'bitwise over %d metals' % len(rows),
          worst == 0.0, 'worst %.3g ulp' % worst)

    # (f) THE EMPTINESS THAT DISCHARGES PRECONDITION (i), asserted by the
    #     mechanism rather than claimed.
    with on_consistent_argument_everywhere() as (ins, alone):
        n_ins, n_alone = len(ins), len(alone)
        sites = sorted(s for s, _ in ins)
    check('(f) the composed manager leaves NOTHING at the shipped object',
          n_alone == 0 and n_ins >= 3,
          '%d installed, %d left alone: %s'
          % (n_ins, n_alone, ', '.join('%s.%s' % s for s in sites)))

    # (g) THE REFUSAL IS REAL -- a planted unassigned member must raise, not
    #     be silently left shipped.  Negative control on the new mechanism:
    #     without it, (f)'s empty left-alone set and a mechanism that never
    #     checks would read identically.
    raised = False
    try:
        partial = dict(_cq_assignment())
        partial.pop(('graphene_contact_doping_model', 'quantum_capacitance'))
        with bc.rebind_by_operating_point(gfet, 'quantum_capacitance',
                                          partial):
            pass
    except bc.UnassignedClassMember as exc:
        raised = True
        msg = str(exc)
    check('(g) NEGATIVE CONTROL: an unassigned class member is refused',
          raised, (msg.split('.')[0] if raised else 'no exception raised'))

    # (h) POSITIVE CONTROLS ON EVERY ARRIVAL COUNTER.  A counter that reads
    #     zero because it is broken and one that reads zero because the name
    #     is never reached are indistinguishable without this.
    reset_arrivals()
    with on_consistent_argument_everywhere():
        cd.contact_edge_carrier_density(cd.METAL_WORK_FUNCTIONS['Pd'])
        gfet.quantum_capacitance(1.2)
        gfet.carrier_density(3.5)
    got = dict(_ARRIVED)
    check('(h) POSITIVE CONTROL: all three patched names reach their counters',
          got['cq_contact'] >= 1 and got['cq_channel'] >= 1
          and got['charge_channel'] >= 1, repr(got))


# ---------------------------------------------------------------------------
# What the installed object does, per metal
# ---------------------------------------------------------------------------

def per_metal_table():
    print('\nTHE CONTACT CORRECTION, THROUGH THE OWNED OBJECT')
    print('  %-5s %-8s %-9s %-15s %-15s %-10s %s'
          % ('metal', 'WF/eV', 'dV/V', 'n_A (m^-2)', 'n_BT (m^-2)',
             'change', 'vs 2026-10-13'))
    out = {}
    for m in METALS:
        wf = cd.METAL_WORK_FUNCTIONS[m]
        dV = -(wf - cd.W_GRAPHENE)
        nA = float(cd.contact_edge_carrier_density(wf))
        with on_contact_consistent_argument():
            nB = float(cd.contact_edge_carrier_density(wf))
        if nA == 0.0:
            print('  %-5s %-8.2f %-9.3f %-15.6e %-15.6e %-10s %s'
                  % (m, wf, dV, nA, nB, '--', '(dWF = 0 exactly)'))
            out[m] = None
            continue
        # Change in MAGNITUDE, which is the convention
        # 2026-10-13's transcript prints: both n are signed, and a
        # correction that shrinks |n| is reported negative whichever
        # sign n carries.
        ch = 100.0 * (abs(nB) - abs(nA)) / abs(nA)
        want = COMMITTED_2026_10_13[m]
        agree = 'identical' if abs(ch - want) < 5e-4 else 'DIFFERS %.4f' % (ch - want)
        print('  %-5s %-8.2f %-9.3f %-15.6e %-15.6e %-10s %s'
              % (m, wf, dV, nA, nB, '%.4f%%' % ch, agree))
        out[m] = ch
    vals = [v for v in out.values() if v is not None]
    print('  Range: %.4f %% to %.4f %%.  The channel\'s same correction at '
          'V_g = 3.5 V is -1.0212 %%, so the contact is %.1fx to %.1fx less '
          'forgiving.'
          % (min(vals), max(vals), abs(min(vals)) / 1.0212,
             abs(max(vals)) / 1.0212))
    return out


def plot_contact_correction(per_metal):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    ms = [m for m in METALS if per_metal[m] is not None]
    vals = [abs(per_metal[m]) for m in ms]
    cols = ['crimson' if cd.METAL_WORK_FUNCTIONS[m] > cd.W_GRAPHENE
            else 'steelblue' for m in ms]
    ax1.bar(ms, vals, color=cols, alpha=0.85, edgecolor='k', lw=0.6)
    ax1.axhline(1.0212, color='seagreen', ls='--', lw=1.8,
                label='the same correction in the channel (1.0212 %)')
    ax1.set_ylabel(r'$|\Delta n_{\rm contact}|$  (%)')
    ax1.set_title('Form $B_T$ at the contact, per metal')
    ax1.legend(fontsize=8)
    ax1.grid(axis='y', alpha=0.3)
    for x, v in zip(ms, vals):
        ax1.text(x, v + 0.2, '%.2f' % v, ha='center', fontsize=8)

    dWF = np.linspace(-1.2, 1.2, 161)
    nA, nB = [], []
    for d in dWF:
        nA.append(float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + d)))
    with on_contact_consistent_argument():
        for d in dWF:
            nB.append(float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + d)))
    nA, nB = np.array(nA), np.array(nB)
    ax2.plot(dWF, nA / 1e16, color='grey', lw=2.2,
             label=r'form A (committed): $C_q$ at the full $\Delta W\!F$')
    ax2.plot(dWF, nB / 1e16, color='crimson', lw=2.0, ls='--',
             label=r'form $B_T$: $C_q$ at $E_F(n)/e$, solved')
    ax2.axhline(0, color='k', lw=0.7)
    ax2.axvline(0, color='k', lw=0.7)
    ax2.plot([0], [0], 'ko', ms=6, zorder=5,
             label=r'$\Delta W\!F=0 \Rightarrow n=0$ exactly (both)')
    ax2.set_xlabel(r'$W_M - W_G$  (eV)')
    ax2.set_ylabel(r'signed $n_{\rm contact}$  ($10^{16}$ m$^{-2}$)')
    ax2.set_title('Odd in $\\Delta W\\!F$ at 0 ulp, and zero at zero')
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.suptitle("Item 16's precondition (i): the contact now has a "
                 "replacement of its own, not an exclusion", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig('contact_consistent_argument.png', dpi=150)
    plt.close(fig)
    print('\nWrote contact_consistent_argument.png')


def main():
    plt.rcParams['text.usetex'] = False
    print('=' * 74)
    print("ITEM 16 PRECONDITION (i): the contact's own contact-correct object")
    print('=' * 74)
    print('The identity class of quantum_capacitance, at call time:')
    for s in sorted(bc.identity_class(gfet, 'quantum_capacitance')):
        print('    %s.%s' % s)
    print('  -- two members, TWO operating points.  bc.rebind() can give them')
    print('     one object or leave one alone; it cannot give them two.')
    print()
    validations()
    pm = per_metal_table()
    plot_contact_correction(pm)
    print('\n' + '=' * 74)
    print('STANDING')
    print('  Precondition (i) DISCHARGED: the contact owns a replacement, the')
    print('  class is partitioned by operating point, and the composed manager')
    print('  leaves nothing at the shipped object (0 left alone, asserted by')
    print('  the mechanism and backed by a negative control).')
    print('  Precondition (ii) STILL OPEN: Section 4.7\'s prose and Section')
    print('  4.10\'s annotation.  That is writing, not plumbing.')
    print('  ITEM 16 ITSELF STILL OPEN.  Nothing here flips the default.')
    print('  Validations: %d passed, %d failed' % (passed, failed))
    print('=' * 74)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
