"""graphene_installed_surface_audit.py -- every leverage number in this
repository was taken with the replacement's own flip NOT entered, and until
today nothing had checked that the two configurations agree.

WHY THIS EXISTS.

2026-10-25's ranked-first item asked for the obvious continuation of the
value-read surface work:

    "COMPLETE THE SURFACE FOR EVERY REPLACEMENT IN `rra.FLIPS`, NOT FOUR
     HOLDERS ... walk every rebind-set member and adjudicate."

Enumerating the population exposes a question the item did not ask.  A
replacement is a function a flip manager INSTALLS at a site; it only ever
runs with that manager entered.  Every leverage number committed here --
2026-10-22's, 2026-10-24's and 2026-10-25's alike -- was taken by calling the
replacement DIRECTLY, in isolation, with no flip entered at all.  Both
instruments involved are configuration-relative: `csg.capture_surface` and
`vrs.value_reads` are static walks that resolve each name to THE OBJECT BOUND
AT THE MOMENT OF THE WALK, and `csg.leverage` calls the holder in whatever
namespace is current.  So a verdict measured in isolation is a verdict about
a configuration in which the replacement never runs, and nothing said whether
that mattered.

THE POPULATION.  `rra.FLIPS` x `rra.rebind_set` gives **8 (flip, site) pairs
and 6 distinct installed objects**.  2026-10-25 audited four holders; two of
the six -- the arrival-counting wrappers `flip._counting` builds -- had never
had a read surface computed at all.

THE ANSWER, AND IT IS MOSTLY GOOD NEWS, MEASURED RATHER THAN ASSUMED.

  * **All six objects are BITWISE flip-invariant under their own flip.**
    Every name on the surface, every counterfactual, and the holder's own
    value: identical to the isolated measurement, bit for bit.  THAT is the
    proposition which licenses every leverage figure this repository has
    committed, and it had never been stated, let alone checked.  It is now a
    check (A1-A3) rather than a silent assumption.
  * **Five of the six are invariant under all three flips.**  One is not.
    `blast._cq_consistent_contact` -- the live, unrepaired contact
    replacement -- is bitwise invariant under its own flip (15/15 names) and
    under NEITHER of the other two: under `qca.on_consistent_argument()` and
    under `flip.on_default_flip(live_contact=False)` its value moves from
    0x1.0c835a6063492p-4 to 0x1.4ef081d146707p-7 (a factor of 6.4) and **11
    of its 15 names cross the exact-zero boundary** -- 5 go bitwise dead
    (`gfet.e` +278.4911 %, `gfet.hbar` and `gfet.v_F` -63.5273 %, `gfet.kB`
    +0.0242 %, `gfet._thermal_rounding` +58.3421 %) and 6 come alive from
    bitwise zero (`qca.cq_thermal`, `qca.cq_thermal_of_n` and
    `qca.CQ_PREFACTOR` at +52.4103 %, `qca.KB` +35.8863 %,
    `qca.n_arg_EF_thermal` +9.0738 %).
  * **The value change is NOT a new result** -- it is 2026-10-21's "the two
    halves do not compose", which measured exactly this composition and
    costed it at -51 % to -62 % in n_contact.  What is new is that the same
    non-composition shows up as a change in the READ SURFACE, which is the
    object every guard in this repository demands adjudication over.  The
    standing item "ADJUDICATE THE RELATION FOR THE FLIPS 10-23 DID NOT
    COMPOSE" is partly discharged here: the relation is measured over all
    6 x 3 cells, and exactly one cell row is non-composing.
  * **A surface change is not the same thing as a verdict change**, and the
    table separates them.  Four of the six objects gain names on their
    surface in some configuration (the counting wrapper for
    `carrier_density` goes 15 -> 17 under `qca.on_consistent_argument()`)
    while every one of their counterfactuals stays bitwise identical.  Only
    `_cq_consistent_contact` has a cell where the VERDICTS would change.

THE REASON THAT OVER-CLAIMS.  `gfet.C_ox` is adjudicated PROTECTED on
`_cq_consistent_contact` with the stated reason "the contact solve closes
against C_interface, so the oxide capacitance reaches this value not at all
... the zero is a property of the contact model".  The verdict is correct in
the configuration it was measured in and stays correct here.  The REASON is
not: under the default flip the same holder gives `gfet.C_ox` **+7.8110 %**,
so the exact zero is a property of the contact model COMPOSED WITH NOTHING
ELSE, not of the contact model.  This repository's own rule is that the
stated reason is part of the verdict, so the reason is corrected in place and
the old wording kept as the superseded text.

A SECOND, SMALLER FINDING: `flip._counting` IS A CAPTURE.  Each counting
wrapper's surface contains exactly one name the underlying replacement's does
not -- the underlying replacement's OWN name (`qca._carrier_density_
consistent`, `qca._quantum_capacitance_consistent`) -- and a rebind of that
name leaves the wrapper's value BITWISE UNCHANGED, because `_counting` holds
the function object in a cell.  So the objects the default flip installs are
mediated by a third, unadjudicated capture, and the channel replacement NAME
is immune under the default flip.  The verdict on it is UNREACHED and not
PROTECTED, by csg's own 2026-10-24 argument for `CDQC`: the wrapper never
looks the name up (closure cell, not a global read), so the bitwise zero is
not a capture neutralising a read and PROTECTED would credit the capture with
work it did not do.  This module's first draft said PROTECTED; the 2026-10-24
guard refused it, which is the fourth time a planted rule written on an
earlier day caught a claim written on a later one.  The wrapper's value is
bitwise identical to the underlying replacement's (W3) and every shared name
has bitwise identical leverage on both (W4), which is what lets 2026-10-25's
verdicts transfer onto the wrappers WITH PROOF rather than by assumption.

WHAT THIS MODULE ADDS.

  1. `installed_population()` -- the population, DERIVED from `rra.FLIPS` via
     `rra.rebind_set`, with the object actually installed at each site.  A
     flip that grew a fourth rebinding would be audited without editing this
     file.
  2. A four-valued dependence axis per (holder, configuration, name):
     `BITWISE_SAME`, `BOTH_LIVE`, `DEAD_AS_INSTALLED`,
     `LIVE_ONLY_AS_INSTALLED`.  Four-valued and not a bool because only the
     last two reverse a verdict, and a two-valued axis puts "moved" and
     "moved in a way that makes the committed verdict false" in one bucket.
  3. `require_installed_adjudication` -- a STRICT strengthening of
     `vrs.require_complete_capture_adjudication`, which it enters unchanged
     at the end, so every refusal of the 2026-10-24 and 2026-10-25 guards
     still applies.  Added: a verdict is required for every name on the
     INSTALLED surface; `PROTECTED` must be bitwise exact AS INSTALLED and
     `EXPOSED` must move AS INSTALLED; `UNREACHED` is refused for a name on
     the installed surface; and the caller must DECLARE the holder's
     flip-invariance, which is measured.  The declaration has no default on
     purpose -- 2026-10-25 passed with the question unasked, and a default
     would let the next session do the same.

WHAT THIS MODULE DOES NOT CLAIM.  `value_reads` is still over-approximate and
its negatives are still unsound (2026-10-25); nothing here repairs that.
Three configurations are measured, not all of them: `live_contact=True`, the
simultaneous composition of two managers, and the `n`-side holders under a
contact flip are not covered, and S7.9 item 23 records that.  "Invariant under
its own flip" is measured at one probe point (`csg.PROBE_V`) with one
perturbation (`csg.PERTURBATION`), so it is a bitwise statement about that
cell and not a theorem.

Usage:  PYTHONPATH=... python3 graphene_installed_surface_audit.py
"""


import contextlib
import os
import sys

import graphene_fet_model as gfet
import graphene_quantum_capacitance_argument as qca
import graphene_contact_argument_blast_radius as blast
import graphene_default_flip as flip
import graphene_replacement_read_audit as rra
import graphene_capture_surface_guard as csg
import graphene_value_read_surface as vrs

from graphene_capture_surface_guard import (
    PROTECTED, UNREACHED, EXPOSED, GOVERNED, VERDICTS, _dotted)


class UnadjudicatedInstalledRead(RuntimeError):
    """A name on the holder's INSTALLED read surface carries no verdict, or
    the caller did not declare the holder's flip-invariance.  Raised before
    anything is measured, like `vrs.UnadjudicatedValueRead`."""


class FalseInstalledVerdict(RuntimeError):
    """A verdict, or a declared flip-invariance, disagrees with the
    counterfactual taken IN THE CONFIGURATION THE HOLDER RUNS IN."""


# The four-valued dependence axis.  Not a bool, for the reason the sibling
# repository's census rows keep re-earning: a two-valued axis forces
# "different" and "different in a way that reverses the verdict" into one
# bucket, and only the second is a defect.
BITWISE_SAME = 'BITWISE_SAME'
BOTH_LIVE = 'BOTH_LIVE'
DEAD_AS_INSTALLED = 'DEAD_AS_INSTALLED'
LIVE_ONLY_AS_INSTALLED = 'LIVE_ONLY_AS_INSTALLED'
DEPENDENCE = (BITWISE_SAME, BOTH_LIVE, DEAD_AS_INSTALLED,
              LIVE_ONLY_AS_INSTALLED)

# Dependence classes that reverse an exact-zero verdict, i.e. the ones a
# PROTECTED/EXPOSED adjudication taken in the wrong configuration gets wrong.
REVERSING = (DEAD_AS_INSTALLED, LIVE_ONLY_AS_INSTALLED)


# ---------------------------------------------------------------------------
# The population, DERIVED from rra.FLIPS
# ---------------------------------------------------------------------------

def installed_population():
    """((flip_label, site, installed_object), ...) over every rebind.

    `rra.rebind_set` is measured OUTSIDE any flip -- entering a manager while
    already inside it is a different measurement (it records the manager's
    writes over an already-substituted namespace, and reports two extra
    sites), and the population wanted here is what a flip installs into the
    PRISTINE tree.
    """
    mods = rra.repo_modules()
    pop = []
    for label, factory in rra.FLIPS:
        sites = sorted(rra.rebind_set(factory))
        with factory():
            for site in sites:
                pop.append((label, site, getattr(mods[site[0]], site[1])))
    return tuple(pop)


def factory_for(label):
    return dict(rra.FLIPS)[label]


# ---------------------------------------------------------------------------
# The same instruments, taken inside the flip
# ---------------------------------------------------------------------------

def installed_surface(holder, factory):
    """`vrs.complete_read_surface(holder)` computed with the flip ENTERED."""
    with factory():
        return vrs.complete_read_surface(holder)


def installed_read_surface(holder, factory):
    """The names the holder READS live, as installed.

    `csg.capture_surface`'s second element plus `vrs.value_reads`: a bare
    scalar read is a live read, and the 2026-10-25 module is where that half
    of the read surface lives.  Kept separate from the full surface because
    the verdict vocabulary distinguishes them -- `UNREACHED` is a statement
    about a name the holder never LOOKS UP, which a name held in a cell does
    not contradict (csg's own 2026-10-24 argument for `CDQC`).
    """
    with factory():
        return csg.capture_surface(holder)[1] | vrs.value_reads(holder)


def installed_leverage(holder, site, factory):
    with factory():
        return csg.leverage(holder, site)


def installed_identity_leverage(holder, site, factory):
    with factory():
        return csg.identity_leverage(holder, site)


def base_invariance(holder, factory, V=csg.PROBE_V):
    """(isolated, installed, bitwise_same) for the holder's own value."""
    iso = float(holder(V))
    with factory():
        inst = float(holder(V))
    return iso, inst, iso.hex() == inst.hex()


def classify(iso_exact, inst_exact, iso, iso_m, inst, inst_m):
    if iso_exact and inst_exact:
        return BITWISE_SAME
    if iso_exact and not inst_exact:
        return LIVE_ONLY_AS_INSTALLED
    if inst_exact and not iso_exact:
        return DEAD_AS_INSTALLED
    a = (iso_m - iso) / iso
    b = (inst_m - inst) / inst
    return BITWISE_SAME if float(a).hex() == float(b).hex() else BOTH_LIVE


def dependence_rows(holder, factory):
    """One row per name on the UNION of the two surfaces.

    The union, not the installed surface: a name that LEAVES the surface when
    the flip is entered is exactly as interesting as one that enters, and it
    is the class the 2026-10-25 contact verdicts live in.
    """
    iso_s = vrs.complete_read_surface(holder)[0]
    inst_s = installed_surface(holder, factory)[0]
    rows = []
    for site in sorted(iso_s | inst_s):
        a, am, ae = csg.leverage(holder, site)
        b, bm, be = installed_leverage(holder, site, factory)
        rows.append({
            'site': site,
            'on_isolated': site in iso_s,
            'on_installed': site in inst_s,
            'iso': (a, am, ae),
            'inst': (b, bm, be),
            'dependence': classify(ae, be, a, am, b, bm),
        })
    return rows


# ---------------------------------------------------------------------------
# THE STRENGTHENED GUARD
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def require_installed_adjudication(label, holder, factory, adjudication,
                                   flip_invariant, report=None):
    """The 2026-10-25 guard, plus the configuration the holder runs in.

    `flip_invariant` is not optional and has no default: a caller that has
    not decided whether its holder's value survives the flip has not done the
    measurement this guard exists for, and a default would let it skip the
    decision silently.  2026-10-25's own transcript is the argument -- it
    passed with the question unasked.
    """
    if flip_invariant not in (True, False):
        raise UnadjudicatedInstalledRead(
            '%s: flip_invariant must be declared True or False, not %r'
            % (label, flip_invariant))

    inst_complete, _, _ = installed_surface(holder, factory)
    inst_read = installed_read_surface(holder, factory)
    flip_rebinds = set()
    for _, fac in rra.FLIPS:
        flip_rebinds |= rra.rebind_set(fac)

    missing = sorted(inst_complete - set(adjudication))
    if missing:
        raise UnadjudicatedInstalledRead(
            '%s: %d name(s) reach this value IN THE CONFIGURATION THE '
            'REPLACEMENT IS INSTALLED IN and carry no verdict: %s.  A surface '
            'computed with the flip not entered is a surface of a '
            'configuration this holder never runs in (2026-10-26).'
            % (label, len(missing), ', '.join(_dotted(m) for m in missing)))

    iso, inst, same = base_invariance(holder, factory)
    if same != bool(flip_invariant):
        raise FalseInstalledVerdict(
            '%s: declared flip_invariant=%s, but the holder returns %s in '
            'isolation and %s as installed'
            % (label, bool(flip_invariant), iso.hex(), inst.hex()))

    rows, wrong = installed_verdict_errors(label, holder, factory,
                                           adjudication, sorted(inst_complete),
                                           inst_read, flip_rebinds)
    if wrong:
        raise FalseInstalledVerdict(
            '%s: %d verdict(s) disagree with the counterfactual taken IN THE '
            'INSTALLED CONFIGURATION:\n%s'
            % (label, len(wrong), '\n'.join(
                '    %-46s %-10s %s' % (_dotted(s), v, why)
                for s, v, why in wrong)))

    with vrs.require_complete_capture_adjudication(
            label, holder, adjudication, report) as (surface, old_rows):
        yield inst_complete, list(old_rows) + rows


def installed_verdict_errors(label, holder, factory, adjudication, sites,
                             inst_read=None, flip_rebinds=None):
    """(rows, wrong) for `sites` only, measured AS INSTALLED.

    Split out of the guard so that the refusal branches added on 2026-10-26
    have POSITIVE CONTROLS THAT REACH THEM.  The first draft's G6 planted a
    false PROTECTED on a holder whose adjudication was incomplete in that
    configuration, so the missing-name refusal fired first and the PROTECTED
    branch was never entered: the control scored green and could not have
    scored anything else, and mutation N4 (delete the branch) went
    undetected.  Sixth instance in this repository of *a verdict needs a
    planted control that reaches it*, and the first in which the unreachable
    control was in the same file as the branch.
    """
    if inst_read is None:
        inst_read = installed_read_surface(holder, factory)
    if flip_rebinds is None:
        flip_rebinds = set()
        for _, fac in rra.FLIPS:
            flip_rebinds |= rra.rebind_set(fac)
    rows, wrong = [], []
    for site in sites:
        verdict, reason = adjudication[site]
        if verdict not in VERDICTS:
            raise UnadjudicatedInstalledRead(
                '%s: adjudication for %s has verdict %r, not one of %s'
                % (label, _dotted(site), verdict, '/'.join(VERDICTS)))
        if not reason:
            raise UnadjudicatedInstalledRead(
                '%s: adjudication for %s has no stated reason'
                % (label, _dotted(site)))

        b0, b1, ident_exact = installed_identity_leverage(holder, site,
                                                          factory)
        if not ident_exact:
            wrong.append((site, verdict,
                          'the IDENTITY rebind of %s moved the value AS '
                          'INSTALLED (%s -> %s): the probe perturbs what it '
                          'measures' % (_dotted(site), b0.hex(), b1.hex())))
            continue
        base, moved, exact = installed_leverage(holder, site, factory)
        rows.append((site, verdict, base, moved, exact))
        if verdict == UNREACHED and site in inst_read:
            wrong.append((site, verdict,
                          'claimed UNREACHED, but %s is on the INSTALLED READ '
                          'surface -- the holder looks this name up in the '
                          'configuration it runs in' % _dotted(site)))
        elif verdict == UNREACHED and not exact:
            wrong.append((site, verdict,
                          'claimed UNREACHED but a rebind of %s moves the '
                          'value by %+.4f %% AS INSTALLED'
                          % (_dotted(site), 100.0 * (moved - base) / base)))
        elif verdict == PROTECTED and not exact:
            wrong.append((site, verdict,
                          'claimed PROTECTED but a rebind of %s moves the '
                          'value BY %+.4f %% AS INSTALLED (it is bitwise '
                          'exact only in isolation)'
                          % (_dotted(site), 100.0 * (moved - base) / base)))
        elif verdict == EXPOSED and exact:
            wrong.append((site, verdict,
                          'claimed EXPOSED but a rebind of %s leaves the '
                          'value BITWISE UNCHANGED AS INSTALLED'
                          % _dotted(site)))
        elif verdict == GOVERNED and site not in flip_rebinds:
            wrong.append((site, verdict,
                          'claimed GOVERNED but %s is in no rebind set'
                          % _dotted(site)))
    return rows, wrong


# ---------------------------------------------------------------------------
# Planted controls -- the mechanism with no physics in it
# ---------------------------------------------------------------------------
# `_CTRL_A` is read bare inside `_ctrl_route_a`, `_CTRL_B` bare inside
# `_ctrl_route_b`, and `_ctrl_route` is the name the control flip rebinds.
# A holder that calls `_ctrl_route` therefore has `_CTRL_A` on its surface in
# isolation and `_CTRL_B` on it as installed -- the EXACT shape of the real
# `blast._cq_consistent_contact` finding, in nine lines with no solver.

_CTRL_A = 3.0
_CTRL_B = 5.0
_CTRL_C = 7.0


def _ctrl_route_a(V):
    return _CTRL_A * float(V)


def _ctrl_route_b(V):
    return _CTRL_B * float(V)


_ctrl_route = _ctrl_route_a

_HERE_MOD = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]
CTRL_ROUTE = (_HERE_MOD, '_ctrl_route')
CTRL_A = (_HERE_MOD, '_CTRL_A')
CTRL_B = (_HERE_MOD, '_CTRL_B')
CTRL_C = (_HERE_MOD, '_CTRL_C')


def _this_module():
    return vrs._modules()[_HERE_MOD]


@contextlib.contextmanager
def control_flip(identity=False):
    """A flip with no physics in it: rebinds `_ctrl_route` to the other route,
    or (identity=True) to itself -- the instrument-is-not-the-perturbation
    control, which must leave every number bitwise unchanged."""
    mod = _this_module()
    old = mod._ctrl_route
    mod._ctrl_route = old if identity else _ctrl_route_b
    try:
        yield
    finally:
        mod._ctrl_route = old


def ctrl_holder_routed(V):
    """Reads the rebound name.  Must be flip-DEPENDENT."""
    return float(_ctrl_route(V))


def ctrl_holder_fixed(V):
    """Reads a name no flip touches.  Must be flip-INVARIANT."""
    return _CTRL_C * float(V)


def _ctrl_surface(holder, installed):
    if installed:
        with control_flip():
            return vrs.complete_read_surface(holder)[0]
    return vrs.complete_read_surface(holder)[0]


def _ctrl_leverage(holder, site, installed, identity=False):
    fn = vrs.identity_leverage_here if identity else vrs.leverage_here
    if installed:
        with control_flip():
            return fn(holder, site)
    return fn(holder, site)


def _ctrl_base(holder, installed):
    if installed:
        with control_flip():
            return float(holder(csg.PROBE_V))
    return float(holder(csg.PROBE_V))


def _pct(base, moved, exact):
    return 'exact 0' if exact else '%+10.4f %%' % (100.0 * (moved - base)
                                                   / base)


# ---------------------------------------------------------------------------
# Adjudications: the 2026-10-25 dicts, extended where the population is wider
# ---------------------------------------------------------------------------

CDC = ('graphene_quantum_capacitance_argument', '_carrier_density_consistent')
CQC = ('graphene_quantum_capacitance_argument',
       '_quantum_capacitance_consistent')
COX = ('graphene_fet_model', 'C_ox')
QCT = ('graphene_quantum_capacitance_argument', 'cq_thermal')
QCTN = ('graphene_quantum_capacitance_argument', 'cq_thermal_of_n')
QNAE = ('graphene_quantum_capacitance_argument', 'n_arg_EF_thermal')
CQP = ('graphene_quantum_capacitance_argument', 'CQ_PREFACTOR')
QKB = ('graphene_quantum_capacitance_argument', 'KB')

_COUNTING_CAPTURE = (
    'the name of the replacement this counting wrapper wraps.  THE SECOND '
    'FINDING OF 2026-10-26: `flip._counting` holds the function object in a '
    'CELL, so this is a THIRD capture in this repository and nothing had '
    'noticed it -- a rebind of the name leaves the installed object BITWISE '
    'unchanged, which means the channel replacement name is immune under the '
    'default flip.  UNREACHED and not PROTECTED, by csg\'s own 2026-10-24 '
    'argument for CDQC: the wrapper never LOOKS THIS NAME UP (it is a closure '
    'cell, not a global read), so the bitwise zero is not the capture '
    'neutralising a read, and calling it PROTECTED would credit the capture '
    'with work it did not do.  The wrapper returns a value bitwise identical '
    'to the replacement it wraps (W3) and every shared name has bitwise '
    'identical leverage on both (W4), which is what transfers the 2026-10-25 '
    'verdicts onto this holder with proof rather than by assumption.')

# The two counting wrappers: their 2026-10-25 holder's adjudication, plus the
# one name their surface adds.  Built by `dict(...)` and not by mutating the
# imported dicts, because `vrs.ADJ_CHARGE` is also the live argument of the
# 2026-10-25 transcript and a module that edited it in place would change
# that transcript's meaning by being imported.
ADJ_WRAP_CHARGE = dict(vrs.ADJ_CHARGE)
ADJ_WRAP_CHARGE[CDC] = (UNREACHED, _COUNTING_CAPTURE)

ADJ_WRAP_CQ = dict(vrs.ADJ_CQ)
ADJ_WRAP_CQ[CQC] = (UNREACHED, _COUNTING_CAPTURE)

# The capture holder, under the DEFAULT flip.  Five qca names enter its
# surface there (because `gfet.quantum_capacitance` resolves to the channel
# replacement and the walk descends into qca's thermal path) and every one is
# BITWISE ZERO on it -- the capture is what makes them unreachable.  The same
# five names carry +52.4103 %/+35.8863 %/+9.0738 % on the UNREPAIRED
# replacement in the same configuration, so this dict and ADJ_CONTACT_LIVE
# below are the two halves of one control pair.
_CAPTURED_AWAY = (
    'on the INSTALLED surface under the default flip and BITWISE ZERO: the '
    '2026-10-21 repair captured the shipped C_q, so the channel path this '
    'name belongs to cannot reach this value however the flip rebinds it.  '
    'PROTECTED and not UNREACHED because the name IS on the installed '
    'surface -- the walk finds it, the counterfactual says it does not '
    'arrive, and those are different statements.')

# NOTE the base dict: `csg.ADJUDICATION_REPAIRED` is the capture holder's own
# adjudication (2026-10-24), not `vrs.ADJ_CONTACT`, which belongs to the
# UNREPAIRED `blast._cq_consistent_contact`.  The two holders are bitwise
# indistinguishable in isolation, so using the wrong one here is a mistake
# that measures clean -- it was made in this module's first draft and the
# inner guard refused it on `cd.quantum_capacitance`.
ADJ_CAPTURE = dict(csg.ADJUDICATION_REPAIRED)
# `csg.ADJUDICATION_REPAIRED` predates the bare-scalar half of the read
# surface, so it carries none of the five names 2026-10-25 added for the
# contact holders; that dict is explicitly stated there to cover BOTH contact
# holders, so it is reused rather than re-reasoned.
ADJ_CAPTURE.update(vrs.ADDED_CONTACT)
for _s in (CQP, QKB, QCT, QCTN, QNAE):
    ADJ_CAPTURE[_s] = (PROTECTED, _CAPTURED_AWAY)

# The live contact replacement under its OWN flip, where every 2026-10-25
# verdict holds bitwise.  Only one entry changes, and it is a REASON and not
# a verdict.
_COX_SUPERSEDED = (
    'EXACT ZERO under this replacement\'s own flip and in isolation, and the '
    'verdict PROTECTED is correct in both.  THE REASON IS CORRECTED '
    '2026-10-26.  The superseded wording was: "the contact solve closes '
    'against C_interface, so the oxide capacitance reaches this value not at '
    'all ... the zero is a property of the contact model and not of the '
    '2026-10-21 repair."  The first clause is right and the generalisation is '
    'wrong: compose this replacement with `qca.on_consistent_argument()` or '
    'with `flip.on_default_flip(live_contact=False)` and the SAME holder '
    'gives gfet.C_ox +7.8110 %.  The exact zero is a property of the contact '
    'model COMPOSED WITH NOTHING ELSE.  The 2026-10-25 number is not '
    'withdrawn -- it is correct, and its scope is now stated.')

ADJ_CONTACT_LIVE = dict(vrs.ADJ_CONTACT)
ADJ_CONTACT_LIVE[COX] = (PROTECTED, _COX_SUPERSEDED)


def _holder_table():
    """((label, holder, factory, adjudication, flip_invariant), ...).

    One row per DISTINCT installed object, each paired with the flip that
    installs it -- which is the configuration its committed verdicts are
    about.  Built from `installed_population()` so a new rebinding appears
    here as a KeyError on a missing adjudication rather than as silence.
    """
    pop = installed_population()
    by_obj = {}
    for label, site, obj in pop:
        by_obj.setdefault(id(obj), (obj, []))[1].append((label, site))
    adj_for = {
        '_carrier_density_consistent': vrs.ADJ_CHARGE,
        '_quantum_capacitance_consistent': vrs.ADJ_CQ,
        '_cq_consistent_contact': ADJ_CONTACT_LIVE,
        'replacement': ADJ_CAPTURE,
    }
    rows = []
    for obj, where in by_obj.values():
        name = obj.__code__.co_name
        if name == 'wrapper':
            held = csg.held_objects(obj)
            inner = [h for h in held if getattr(h, '__name__', '')
                     in ('_carrier_density_consistent',
                         '_quantum_capacitance_consistent')]
            adj = (ADJ_WRAP_CHARGE
                   if inner and inner[0].__name__ ==
                   '_carrier_density_consistent' else ADJ_WRAP_CQ)
            label = '_counting(%s)' % (inner[0].__name__ if inner else '?')
        else:
            adj = adj_for[name]
            label = name
        own = sorted(set(l for l, _ in where))[0]
        rows.append((label, obj, factory_for(own), adj, own,
                     sorted(set('%s.%s' % s for _, s in where))))
    rows.sort(key=lambda r: r[0])
    return rows


# ---------------------------------------------------------------------------

def main():
    checks = []

    def check(name, ok, detail=''):
        detail = detail if isinstance(detail, str) else repr(detail)
        checks.append((name, bool(ok), detail))
        print('    %-4s %s%s' % ('PASS' if ok else 'FAIL', name,
                                 ('   -- ' + detail) if detail else ''))

    print('=' * 78)
    print('THE INSTALLED SURFACE: a verdict measured in isolation is a')
    print('  verdict about a configuration the replacement never runs in')
    print('=' * 78)
    print()

    table = _holder_table()
    pop = installed_population()

    # -- 0.  The population ------------------------------------------------
    print('-' * 78)
    print('0.  THE POPULATION, DERIVED FROM rra.FLIPS (the ranked-first item)')
    print('-' * 78)
    print()
    for label, site, obj in pop:
        print('    %-44s %-46s -> %s'
              % (label[:44], '%s.%s' % site, obj.__code__.co_name))
    print()
    print('    %d (flip, site) pairs, %d distinct installed objects'
          % (len(pop), len(table)))
    check('P1 the population is DERIVED: every pair comes from rra.FLIPS via '
          'rra.rebind_set, so a flip that grew a rebinding is audited '
          'without editing this file', len(pop) == 8 and len(table) == 6,
          '%d pairs, %d objects' % (len(pop), len(table)))
    audited_2026_10_25 = {'_carrier_density_consistent',
                          '_quantum_capacitance_consistent',
                          '_cq_consistent_contact', 'replacement'}
    never = sorted(r[0] for r in table
                   if r[1].__code__.co_name not in audited_2026_10_25)
    check('P2 two of the six installed objects had NEVER had a read surface '
          'computed: the arrival-counting wrappers the default flip installs',
          len(never) == 2 and all('_counting' in n for n in never), never)
    print()

    # -- 1.  Controls ------------------------------------------------------
    print('-' * 78)
    print('1.  THE MECHANISM ON CASES WHOSE ANSWER IS KNOWN IN ADVANCE')
    print('-' * 78)
    print()
    routed, fixed = ctrl_holder_routed, ctrl_holder_fixed
    s_iso = _ctrl_surface(routed, False)
    s_ins = _ctrl_surface(routed, True)
    check('C1 a holder reading a name the flip rebinds: _CTRL_A is on the '
          'ISOLATED surface and off the INSTALLED one, and _CTRL_B the other '
          'way round -- the real contact finding with no solver in it',
          CTRL_A in s_iso and CTRL_A not in s_ins
          and CTRL_B in s_ins and CTRL_B not in s_iso,
          'iso=%s inst=%s' % (sorted(n for _, n in s_iso),
                              sorted(n for _, n in s_ins)))
    a_iso = _ctrl_leverage(routed, CTRL_A, False)
    a_ins = _ctrl_leverage(routed, CTRL_A, True)
    b_iso = _ctrl_leverage(routed, CTRL_B, False)
    b_ins = _ctrl_leverage(routed, CTRL_B, True)
    check('C2 and the counterfactuals reverse with the surface: _CTRL_A '
          'moves the value in isolation and is BITWISE DEAD as installed; '
          '_CTRL_B is bitwise dead in isolation and moves as installed',
          (not a_iso[2]) and a_ins[2] and b_iso[2] and (not b_ins[2]),
          'A %s/%s  B %s/%s' % (_pct(*a_iso), _pct(*a_ins),
                                _pct(*b_iso), _pct(*b_ins)))
    cls_a = classify(a_iso[2], a_ins[2], a_iso[0], a_iso[1],
                     a_ins[0], a_ins[1])
    cls_b = classify(b_iso[2], b_ins[2], b_iso[0], b_iso[1],
                     b_ins[0], b_ins[1])
    check('C3 the four-valued axis classes them DEAD_AS_INSTALLED and '
          'LIVE_ONLY_AS_INSTALLED -- the two classes that reverse a verdict',
          cls_a == DEAD_AS_INSTALLED and cls_b == LIVE_ONLY_AS_INSTALLED,
          '%s / %s' % (cls_a, cls_b))
    c_iso = _ctrl_leverage(fixed, CTRL_C, False)
    c_ins = _ctrl_leverage(fixed, CTRL_C, True)
    cls_c = classify(c_iso[2], c_ins[2], c_iso[0], c_iso[1],
                     c_ins[0], c_ins[1])
    check('C4 NEGATIVE control: a holder reading a name no flip touches is '
          'BITWISE_SAME, and its value is bitwise equal in both '
          'configurations',
          cls_c == BITWISE_SAME
          and _ctrl_base(fixed, False).hex() == _ctrl_base(fixed, True).hex(),
          cls_c)
    check('C3b every class the axis returns is a declared member of '
          'DEPENDENCE, and exactly two of the four are in REVERSING -- the '
          'ones that make a verdict taken in the wrong configuration false',
          cls_a in DEPENDENCE and cls_b in DEPENDENCE and cls_c in DEPENDENCE
          and set(REVERSING) == {DEAD_AS_INSTALLED, LIVE_ONLY_AS_INSTALLED}
          and len(DEPENDENCE) == 4
          and sorted(c in REVERSING for c in (cls_a, cls_b, cls_c)).count(True)
          == 2,
          '%s/%s/%s' % (cls_a, cls_b, cls_c))
    check('C3c the control reaches the value through the NAME the flip '
          'rebinds -- rebinding CTRL_ROUTE is what swaps _CTRL_A for _CTRL_B '
          'on the surface -- yet CTRL_ROUTE is itself INVISIBLE to '
          'csg.capture_surface, because the script module is outside '
          'rra.repo_modules().  The same gap vrs._modules() was written for '
          'on 2026-10-25, and the reason the controls here must measure '
          'through vrs.leverage_here rather than csg.leverage',
          CTRL_ROUTE not in csg.capture_surface(routed)[0]
          and CTRL_A in vrs.value_reads(routed)
          and _this_module()._ctrl_route is _ctrl_route_a,
          _dotted(CTRL_ROUTE))
    with control_flip(identity=True):
        id_base = float(routed(csg.PROBE_V))
        id_surf = vrs.complete_read_surface(routed)[0]
    check('C5 IDENTITY control flip: rebinding the route to ITSELF leaves the '
          'value and the surface bitwise unchanged -- the instrument is not '
          'the perturbation',
          id_base.hex() == float(routed(csg.PROBE_V)).hex()
          and id_surf == s_iso, id_base.hex())
    ia = _ctrl_leverage(routed, CTRL_A, False, identity=True)
    check('C6 the per-name identity rebind is bitwise neutral in the control '
          'too', ia[2], ia[0].hex())
    print()

    # -- 2.  The cross-configuration table ---------------------------------
    print('-' * 78)
    print('2.  EVERY INSTALLED OBJECT x EVERY FLIP')
    print('-' * 78)
    print()
    own_invariant, all_rows, nonrow = True, {}, []
    for label, obj, own_fac, adj, own_label, sites in table:
        print('    %s' % label)
        print('      installed at: %s' % ', '.join(sites))
        print('      by          : %s' % own_label)
        iso = vrs.complete_read_surface(obj)[0]
        print('      %-46s %-7s %-9s %s'
              % ('configuration', '|surf|', 'base', 'dependence classes'))
        print('      %-46s %-7d %-9s' % ('(isolation)', len(iso), '--'))
        for flab, ffac in rra.FLIPS:
            rows = dependence_rows(obj, ffac)
            inst = installed_surface(obj, ffac)[0]
            b_iso, b_ins, b_same = base_invariance(obj, ffac)
            counts = {}
            for r in rows:
                counts[r['dependence']] = counts.get(r['dependence'], 0) + 1
            same_all = all(r['dependence'] == BITWISE_SAME for r in rows)
            tag = '  <== its own flip' if flab == own_label else ''
            print('      %-46s %-7d %-9s %s%s'
                  % (flab, len(inst), 'bitwise' if b_same else 'MOVES',
                     ', '.join('%s=%d' % kv for kv in sorted(counts.items())),
                     tag))
            all_rows[(label, flab)] = rows
            if flab == own_label:
                own_invariant = own_invariant and same_all and b_same
            elif not same_all:
                nonrow.append((label, flab, counts))
        print()
    check('A1 EVERY one of the six installed objects is BITWISE flip-'
          'invariant under ITS OWN flip -- every name, every counterfactual, '
          'and the value itself.  This is the proposition every committed '
          'leverage number rests on, measured for the first time',
          own_invariant)
    offenders = sorted(set(l for l, _, _ in nonrow))
    check('A2 exactly ONE object is not invariant under the other flips, and '
          'it is the live unrepaired contact replacement',
          offenders == ['_cq_consistent_contact'], offenders)
    print()

    # -- 3.  The one non-composing row -------------------------------------
    print('-' * 78)
    print('3.  THE ONE ROW THAT MOVES, NAME BY NAME')
    print('-' * 78)
    print()
    live = blast._cq_consistent_contact
    dfac = factory_for('flip.on_default_flip(live_contact=False)')
    rows = dependence_rows(live, dfac)
    b_iso, b_ins, _ = base_invariance(live, dfac)
    print('    blast._cq_consistent_contact under '
          'flip.on_default_flip(live_contact=False)')
    print('      value  isolation %s   as installed %s   ratio %.4f'
          % (b_iso.hex(), b_ins.hex(), b_iso / b_ins))
    print()
    print('      %-52s %-5s %-12s %-12s %s'
          % ('name', 'surf', 'ISOLATED', 'INSTALLED', 'class'))
    for r in rows:
        w = ('both' if r['on_isolated'] and r['on_installed']
             else ('iso' if r['on_isolated'] else 'inst'))
        print('      %-52s %-5s %-12s %-12s %s'
              % (_dotted(r['site']), w, _pct(*r['iso']).strip(),
                 _pct(*r['inst']).strip(), r['dependence']))
    dead = [r for r in rows if r['dependence'] == DEAD_AS_INSTALLED]
    alive = [r for r in rows if r['dependence'] == LIVE_ONLY_AS_INSTALLED]
    check('A3 11 of 15 names cross the exact-zero boundary: 5 go bitwise '
          'dead and 6 come alive from bitwise zero',
          len(dead) == 5 and len(alive) == 6,
          'dead=%d alive=%d' % (len(dead), len(alive)))
    cap_obj0 = [r for r in table if r[0] == 'replacement'][0][1]
    check('A4 and this reproduces 2026-10-21 rather than contradicting it: '
          'the composition it costed at -51..-62 %% in n_contact is the same '
          'one, here measured as a change in the READ SURFACE.  The '
          'invariance FLAG itself is asserted, bitwise and in both '
          'directions -- False for the live replacement, True for the '
          'repaired one in the same configuration -- because a tolerant '
          'comparison here would report a factor of 6.4 as invariant (N8)',
          not b_iso.hex() == b_ins.hex() and b_iso / b_ins > 1.0
          and base_invariance(live, dfac)[2] is False
          and base_invariance(cap_obj0, dfac)[2] is True,
          'ratio %.4f' % (b_iso / b_ins))
    cox = [r for r in rows if r['site'] == COX][0]
    check('A5 gfet.C_ox: bitwise ZERO in isolation and under this '
          "replacement's own flip, +7.8110 %% under the default flip.  The "
          'PROTECTED verdict stands; its stated REASON generalised from one '
          'configuration and is corrected in place',
          cox['iso'][2] and not cox['inst'][2]
          and abs(100.0 * (cox['inst'][1] - cox['inst'][0]) / cox['inst'][0]
                  - 7.8110) < 1e-3,
          _pct(*cox['inst']).strip())
    print()

    # -- 4.  The counting wrappers are captures ----------------------------
    print('-' * 78)
    print('4.  flip._counting IS A CAPTURE (the second finding)')
    print('-' * 78)
    print()
    wraps = [r for r in table if r[0].startswith('_counting')]
    inner_of = {'_counting(_carrier_density_consistent)':
                (qca._carrier_density_consistent, CDC),
                '_counting(_quantum_capacitance_consistent)':
                (qca._quantum_capacitance_consistent, CQC)}
    ok_extra, ok_zero, ok_same, ok_rest = True, True, True, True
    for label, obj, own_fac, adj, own_label, sites in wraps:
        inner, site = inner_of[label]
        w_s = vrs.complete_read_surface(obj)[0]
        i_s = vrs.complete_read_surface(inner)[0]
        extra = w_s - i_s
        b, m, e = csg.leverage(obj, site)
        bw, bi = float(obj(csg.PROBE_V)), float(inner(csg.PROBE_V))
        print('    %s' % label)
        print('      surface %d vs wrapped %d; adds %s'
              % (len(w_s), len(i_s), sorted(_dotted(s) for s in extra)))
        print('      rebinding the wrapped name: %s' % _pct(b, m, e).strip())
        print('      wrapper value %s   wrapped value %s' % (bw.hex(),
                                                             bi.hex()))
        ok_extra = ok_extra and extra == {site}
        ok_zero = ok_zero and e
        ok_same = ok_same and bw.hex() == bi.hex()
        shared = sorted(w_s & i_s)
        for s in shared:
            a1 = csg.leverage(obj, s)
            a2 = csg.leverage(inner, s)
            ok_rest = ok_rest and a1[2] == a2[2] \
                and a1[1].hex() == a2[1].hex()
    print()
    check('W1 each counting wrapper adds EXACTLY ONE name to the surface of '
          'the replacement it wraps: that replacement\'s own name', ok_extra)
    check('W2 and a rebind of that name is BITWISE ZERO on the wrapper -- '
          '`_counting` holds the function in a cell, so this is a third '
          'capture in this repository and the channel replacement NAME is '
          'immune under the default flip', ok_zero)
    check('W3 EXACT ANCHOR: the wrapper returns a value bitwise identical to '
          'the replacement it wraps, so the wrapper is not a perturbation',
          ok_same)
    check('W4 and every name the two surfaces SHARE has bitwise identical '
          'leverage on both -- which is what transfers the 2026-10-25 '
          'verdicts onto the wrappers with proof', ok_rest)
    print()

    # -- 5.  Exact anchors on the (hbar v_F)^2 coincidence -----------------
    print('-' * 78)
    print('5.  EXACT ANCHORS: TWO NAMES THAT MUST MOVE THE VALUE IDENTICALLY')
    print('-' * 78)
    print()
    HB = ('graphene_fet_model', 'hbar')
    VF = ('graphene_fet_model', 'v_F')
    cap = [r for r in table if r[0] == 'replacement'][0]
    cap_obj, cap_fac = cap[1], cap[2]
    pairs = []
    for cfg, fac in (('isolation', None),) + tuple(
            (l, f) for l, f in rra.FLIPS):
        if fac is None:
            h1, h2 = csg.leverage(cap_obj, HB), csg.leverage(cap_obj, VF)
        else:
            h1 = installed_leverage(cap_obj, HB, fac)
            h2 = installed_leverage(cap_obj, VF, fac)
        pairs.append((cfg, h1, h2))
        print('    %-46s hbar %s   v_F %s'
              % (cfg, h1[1].hex(), h2[1].hex()))
    check('X1 hbar and v_F move the capture holder to the BITWISE IDENTICAL '
          'value in every configuration, because the C_q prefactor depends '
          'on them only through (hbar v_F)^2.  A probe reaching the value by '
          'any other route would not reproduce this',
          all(h1[1].hex() == h2[1].hex() for _, h1, h2 in pairs))
    check('X2 and both are BITWISE DEAD on the live replacement under the '
          'default flip, where the channel path supplies C_q instead',
          installed_leverage(live, HB, dfac)[2]
          and installed_leverage(live, VF, dfac)[2])
    print()

    # -- 6.  The guard -----------------------------------------------------
    print('-' * 78)
    print('6.  THE STRENGTHENED GUARD, ON EVERY HOLDER')
    print('-' * 78)
    print()
    passed = []
    for label, obj, own_fac, adj, own_label, sites in table:
        _, _, inv = base_invariance(obj, own_fac)
        try:
            with require_installed_adjudication(label, obj, own_fac, adj,
                                                flip_invariant=inv) as (s, r):
                passed.append((label, len(s), len(r)))
        except (UnadjudicatedInstalledRead, FalseInstalledVerdict,
                vrs.UnadjudicatedValueRead, vrs.FalseValueReadVerdict,
                csg.UnadjudicatedCaptureSurface, csg.FalseCaptureVerdict) as e:
            print('    REFUSED %s: %s' % (label, e))
            passed.append((label, -1, -1))
    for label, ns, nr in passed:
        print('    %-46s installed surface %2d   adjudicated rows %2d'
              % (label, ns, nr))
    check('G1 the guard accepts all six holders once every name on the '
          'INSTALLED surface carries a verdict that agrees with the '
          'installed counterfactual', all(ns >= 0 for _, ns, _ in passed))

    # the refusals, each with a planted positive control
    def refuses(fn, exc):
        try:
            fn()
        except exc:
            return True
        except Exception:
            return False
        return False

    wl, wo, wf, wadj, _, _ = [r for r in table
                              if r[0].startswith('_counting')][0]
    bad = dict(wadj)
    del bad[CDC if CDC in bad else CQC]

    def _missing():
        with require_installed_adjudication('ctl', wo, wf, bad,
                                            flip_invariant=True):
            pass
    check('G2 a name on the installed surface with NO verdict is refused',
          refuses(_missing, UnadjudicatedInstalledRead))

    badv = dict(wadj)
    badv[COX] = (UNREACHED, 'claims the oxide capacitance is not read')

    def _unreached():
        with require_installed_adjudication('ctl', wo, wf, badv,
                                            flip_invariant=True):
            pass
    check('G3 UNREACHED on a name that IS on the installed READ surface is '
          'refused -- and the rule is the READ surface, not the full one, so '
          'a name held in a closure cell can still be UNREACHED',
          refuses(_unreached, FalseInstalledVerdict))

    def _badinv():
        with require_installed_adjudication('ctl', wo, wf, wadj,
                                            flip_invariant=False):
            pass
    check('G4 a FALSE flip-invariance declaration is refused -- the '
          'declaration is measured, not taken on trust',
          refuses(_badinv, FalseInstalledVerdict))

    def _nodecl():
        with require_installed_adjudication('ctl', wo, wf, wadj,
                                            flip_invariant=None):
            pass
    check('G5 an UNDECLARED flip-invariance is refused: there is no default, '
          'because 2026-10-25 passed with the question unasked',
          refuses(_nodecl, UnadjudicatedInstalledRead))

    # G6/G7: the two branches added today, exercised through the extracted
    # entry so that the missing-name refusal cannot fire first and hide them.
    # The first draft tested them through the whole guard on an incompletely
    # adjudicated holder; the control was green and unreachable, and mutation
    # N4 proved it.
    GE_ = ('graphene_fet_model', 'e')
    _, w6 = installed_verdict_errors(
        'ctl', live, dfac, {COX: (PROTECTED, 'bitwise zero in isolation, and '
                                             'that is the whole claim')},
        [COX])
    check('G6 THE REAL SHAPE, on the real name: gfet.C_ox is bitwise exact '
          'in isolation, so PROTECTED passes every 2026-10-24 and 2026-10-25 '
          'check; AS INSTALLED under the default flip it moves +7.8110 %, '
          'and the branch added today refuses it',
          len(w6) == 1 and w6[0][0] == COX and 'AS INSTALLED' in w6[0][2],
          w6[0][2] if w6 else 'NOT REFUSED')
    _, w7 = installed_verdict_errors(
        'ctl', live, dfac, {GE_: (EXPOSED, 'claims the elementary charge '
                                           'moves, which it does -- in '
                                           'isolation')},
        [GE_])
    check('G7 and the mirror, also on a real name: gfet.e carries +278.4911 '
          '% in isolation and is BITWISE DEAD as installed, so EXPOSED is '
          'refused there',
          len(w7) == 1 and w7[0][0] == GE_
          and 'BITWISE UNCHANGED' in w7[0][2],
          w7[0][2] if w7 else 'NOT REFUSED')
    _, w7b = installed_verdict_errors(
        'ctl', live, dfac, {COX: (EXPOSED, 'the corrected verdict for this '
                                           'configuration')},
        [COX])
    check('G8 NEGATIVE control for both: the verdict that IS true as '
          'installed passes the same branch -- so G6 and G7 are not a '
          'function that refuses everything', not w7b)

    # G9: no assignment of the four verdicts to gfet.C_ox on this holder is
    # true in BOTH configurations.  Exhaustive, four cases.
    flip_reb = set()
    for _, _fac in rra.FLIPS:
        flip_reb |= rra.rebind_set(_fac)
    verdict_ok = {}
    for v in VERDICTS:
        _, wi = installed_verdict_errors('ctl', live, dfac,
                                         {COX: (v, 'exhaustive probe')}, [COX])
        b, m, e = csg.leverage(live, COX)
        iso_read = csg.capture_surface(live)[1] | vrs.value_reads(live)
        iso_ok = not ((v == PROTECTED and not e)
                      or (v == PROTECTED and COX not in iso_read)
                      or (v == EXPOSED and e)
                      or (v == UNREACHED and not e)
                      or (v == UNREACHED and COX in iso_read)
                      or (v == GOVERNED and COX not in flip_reb))
        verdict_ok[v] = (not wi, iso_ok)
    print()
    print('    gfet.C_ox on blast._cq_consistent_contact, all four verdicts:')
    for v in VERDICTS:
        ins, iso = verdict_ok[v]
        print('      %-10s isolated: %-8s as installed: %s'
              % (v, 'ok' if iso else 'refused', 'ok' if ins else 'refused'))
    check('G9 EXHAUSTIVE: none of the four verdicts is true in BOTH '
          'configurations for gfet.C_ox on this holder, so the vocabulary '
          'cannot express the state of this name and "the flips do not '
          'compose" is not dischargeable by adjudication alone',
          not any(ins and iso for ins, iso in verdict_ok.values()),
          repr({k: v for k, v in verdict_ok.items()}))
    print()

    # -- 7.  Summary -------------------------------------------------------
    print('=' * 78)
    ok = sum(1 for _, g, _ in checks if g)
    print('VALIDATION: %d/%d' % (ok, len(checks)))
    print('%d/%d checks passed' % (ok, len(checks)))
    print('=' * 78)
    return 0 if ok == len(checks) else 1


if __name__ == '__main__':
    sys.exit(main())
