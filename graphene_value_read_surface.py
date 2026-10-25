"""graphene_value_read_surface.py -- the read surface's NEGATIVES are NOT
sound, and that is the one property every `UNREACHED` verdict rests on.

WHY THIS EXISTS.

2026-10-24 built `graphene_capture_surface_guard.py`, which adjudicates a
capture by COUNTERFACTUAL and separates `PROTECTED` from `UNREACHED` using
the static read surface.  Its Section 0 states the ground for that split:

    "an over-approximate read set has reliable NEGATIVES, so 'not in the
     static read surface' is trustworthy even though 'in it' is not."

THAT SENTENCE IS FALSE, and this module is the measurement that shows it.
`rra.static_reads` records a name in exactly two situations: a module object
loaded and then an attribute taken from it (`qca.A_FERMI`), or a bare global
that resolves to a REPO FUNCTION (`cq_thermal`).  A bare global that resolves
to a module-level SCALAR -- read from inside the very module that defines it,
which is how every constant in this repository is read by its own functions --
matches neither branch and is recorded NOWHERE.  The surface's negatives
therefore do not mean "not read"; they mean "not read in a spelling the
walker can see".

The standing item asked whether `gfet._thermal_rounding` and `qca.A_FERMI`
carry 2026-10-24's +58 % onto the other three replacements, "which would make
them repository-wide levers rather than a contact fact".  Both halves of that
question are answered here, and they answer OPPOSITE WAYS:

    gfet._thermal_rounding   a CONTACT FACT.  Bitwise EXACTLY ZERO on both
                             qca replacements -- and that zero is measured by
                             counterfactual, not inferred from the name's
                             absence from a surface, which after this module
                             would prove nothing (V3).

    qca.A_FERMI              REPOSITORY-WIDE, and INVISIBLE where it matters.
                             +29.8900 % on qca._quantum_capacitance_consistent
                             and +0.0121 % on qca._carrier_density_consistent,
                             while being OFF the surface of both -- because
                             inside `qca` it is spelled bare.  The same name
                             is ON the surface of the two contact holders,
                             where `blast` spells it `qca.A_FERMI`.  One name,
                             one value, one mechanism, and its visibility to
                             every audit in this repository is decided by the
                             SPELLING at the read site (V4).

AND THE LARGER ONE IT TURNED UP.  `qca.CQ_PREFACTOR` carries +50.4145 % on
the capacitance replacement on the same off-surface footing, and `gfet.e`
carries +278.4911 % on both contact holders.  Across the four holders audited
here the completed surface adds 18 (holder, name) pairs, and EVERY ONE OF
THEM MOVES THE VALUE.  Not one was ever adjudicated, because the guard's
`missing` check can only require a verdict for a name its surface names.

SO 2026-10-24's "the capture surface, fully adjudicated" IS SUPERSEDED.  It
adjudicated 9 names on `flip.make_contact_cq(live=False)`; 5 more reach that
holder's value, the largest at +278.4911 %.  Yesterday's nine verdicts are
not wrong -- every one is re-measured here and every one still holds -- but
the claim of COMPLETENESS around them was false, and §4.6.20 is annotated in
place rather than rewritten.

REPOSITORY-WIDE EXTENT.  75 of the 90 module-level scalars in the imported
repository are read bare inside their own module and are therefore invisible
to every surface this repository computes -- `rra.static_reads`, the
2026-10-23 write/read guard, `bc.identity_class`, and yesterday's capture
surface alike.  83 %, not an edge case.  The counterfactual is what settles
whether any one of them is a lever on a given holder; the 18 measured here
are the ones reachable from the four replacements actually under audit.

A DEFECT IN THIS MODULE'S OWN CHECKS, CAUGHT BY ITS OWN MUTATION TEST.  The
first draft of V10c read `(blind - b_ex) == set(_PLANTED_APPARATUS)`, and
since `b_ex` is `blind - _PLANTED_APPARATUS` by construction that is TRUE FOR
ANY VALUE of the exclusion -- a tautology.  M7 widens the exclusion to swallow
`qca.A_FERMI` and `qca.CQ_PREFACTOR`, i.e. hides the two findings of the day
behind the apparatus carve-out the 2026-10-17 standing item warned about, and
the first draft scored 31/31 against it.  V10c now compares the exclusion with
the census's own result for this module and V10f requires every excluded name
to belong to this module, whose name is taken from `__file__`; neither can be
satisfied by editing the exclusion.  This is the fourth consecutive session in
which a check that passed its plausible reading was false under mutation.

WHAT THIS MODULE DOES NOT CLAIM.  `value_reads` is still a STATIC reader and
still over-approximate: a name on a branch never taken is reported.  Its
negatives are better than `static_reads`' but they are not proven sound
either -- a scalar reached through a dict, a getattr or a comprehension cell
could still escape it.  The honest statement is narrower and is the one made
here: the completed surface is a STRICT SUPERSET of the old one (V1, exact),
it adds 18 pairs with measured leverage, and no verdict anywhere should be
read as "this name cannot reach the value" when its only evidence is absence
from a static surface.

Usage:  PYTHONPATH=... python3 graphene_value_read_surface.py
"""

import contextlib
import dis
import os
import sys
import types

import graphene_fet_model as gfet
import graphene_quantum_capacitance_argument as qca
import graphene_contact_argument_blast_radius as blast
import graphene_default_flip as flip
import graphene_replacement_read_audit as rra
import graphene_capture_surface_guard as csg

from graphene_capture_surface_guard import (
    PROTECTED, UNREACHED, EXPOSED, GOVERNED, VERDICTS, _dotted)


class UnadjudicatedValueRead(RuntimeError):
    """A name on the COMPLETED read surface -- reached by a bare scalar read
    the static surface cannot see -- carries no verdict.  Raised before
    anything is measured, like `csg.UnadjudicatedCaptureSurface`."""


class FalseValueReadVerdict(RuntimeError):
    """A verdict on such a name disagrees with the counterfactual."""


# The planted constants of Section 5's controls.  Declared here, at module
# level, because the 2026-10-17 standing item ("a module must not plant the
# pattern it censuses in the tree other audits read") applies to this module
# with full force: `_PLANTED_BARE` read bare below IS a member of the blind
# class this module censuses.  Rather than hide that, the census reports the
# unfiltered number AND the number with these names removed, and V9c checks
# that the removed set is EXACTLY this tuple -- so the exclusion cannot be
# widened to drop anything else.
_PLANTED_BARE = 3.0
_PLANTED_APPARATUS = (('graphene_value_read_surface', '_PLANTED_BARE'),)



_HERE = os.path.dirname(os.path.abspath(__file__))
_THIS_MODULE = os.path.splitext(
    os.path.basename(os.path.abspath(__file__)))[0]


def _modules():
    """`rra.repo_modules()` plus the module being RUN AS A SCRIPT.

    `rra.repo_modules()` keys on `sys.modules` and skips any name starting
    with `__`, so the module executing as `__main__` is absent from its own
    table.  That is harmless for an audit that only ever looks at OTHER
    modules, and it is not harmless here: Section 1's control has to read a
    scalar from inside the module that defines it, and that module is this
    one.  Without this, the control could not be written at all -- and a
    mechanism whose isolated control cannot be written is a mechanism
    demonstrated only on cases that also contain physics.
    """
    mods = dict(rra.repo_modules())
    main = sys.modules.get('__main__')
    f = getattr(main, '__file__', None)
    if f and os.path.dirname(os.path.abspath(f)) == _HERE:
        stem = os.path.splitext(os.path.basename(os.path.abspath(f)))[0]
        mods.setdefault(stem, main)
    return mods


@contextlib.contextmanager
def _rebound_here(site, factory):
    """`csg.rebound` over `_modules()` instead of `rra.repo_modules()`.

    Used ONLY for the planted controls, whose site lives in the script
    module.  Every measurement on a real holder below goes through
    `csg.leverage` unchanged, so what is reported about this repository is
    reported by the instrument that was audited yesterday and not by a new
    one written today.
    """
    mod = _modules()[site[0]]
    old = getattr(mod, site[1])
    setattr(mod, site[1], factory(old))
    try:
        yield old
    finally:
        setattr(mod, site[1], old)


def leverage_here(holder, site, k=csg.PERTURBATION, V=csg.PROBE_V):
    """(base, moved, exact) for a site in the script module."""
    base = float(holder(V))
    with _rebound_here(site, lambda old: csg.perturb(old, k)):
        moved = float(holder(V))
    return base, moved, base.hex() == moved.hex()


def identity_leverage_here(holder, site, V=csg.PROBE_V):
    base = float(holder(V))
    with _rebound_here(site, lambda old: old):
        again = float(holder(V))
    return base, again, base.hex() == again.hex()


# ---------------------------------------------------------------------------
# The reader the static surface is missing
# ---------------------------------------------------------------------------

def value_reads(fn, maxdepth=8):
    """{(module, attr)} for module-level SCALARS reached by a BARE global read.

    Walks the same call graph as `rra.static_reads` -- same load opcodes, same
    closure resolution, same depth bound, same visited set on code identity --
    and differs in exactly one place: where `static_reads` discards a bare
    global that is not a repo function, this records it if it is an `int` or
    `float` attribute of the module whose globals the read happened in.

    `bool` is excluded for the reason `csg.perturb` excludes it: there is no
    perturbation of a bool that is not also a different kind of value, so a
    bool could carry no counterfactual and a verdict on it would be untestable.
    """
    mods = _modules()
    # Keyed on the TABLE's key, not on `m.__name__`: the script module is in
    # the table under its file stem while its `__name__` is `__main__`, and
    # keying on the latter would hand back a name the table cannot resolve.
    mod_by_id = {id(m): k for k, m in mods.items()}
    glob_owner = {id(m.__dict__): k for k, m in mods.items()}
    vals, seen = set(), set()

    def walk(code, g, freemap, depth):
        if depth > maxdepth or id(code) in seen:
            return
        seen.add(id(code))
        owner = glob_owner.get(id(g))
        ins = list(dis.get_instructions(code))
        for i, x in enumerate(ins):
            if x.opname not in rra._LOAD_OPS:
                continue
            name = x.argval
            if x.opname in ('LOAD_DEREF', 'LOAD_CLASSDEREF'):
                base = freemap.get(name, rra._MISSING)
            else:
                base = g.get(name, rra._MISSING)
            if base is rra._MISSING:
                continue
            nxt = ins[i + 1] if i + 1 < len(ins) else None
            if id(base) in mod_by_id and nxt is not None \
                    and nxt.opname in rra._ATTR_OPS:
                mname, attr = mod_by_id[id(base)], nxt.argval
                tgt = getattr(mods[mname], attr, rra._MISSING)
                if isinstance(tgt, types.FunctionType) \
                        and id(tgt.__globals__) in glob_owner:
                    walk(tgt.__code__, tgt.__globals__,
                         rra._closure_map(tgt), depth + 1)
            elif isinstance(base, types.FunctionType) \
                    and id(base.__globals__) in glob_owner:
                walk(base.__code__, base.__globals__,
                     rra._closure_map(base), depth + 1)
            elif owner is not None and isinstance(base, (int, float)) \
                    and not isinstance(base, bool) \
                    and getattr(mods[owner], name, rra._MISSING) is base:
                vals.add((owner, name))
        for const in code.co_consts:
            if isinstance(const, types.CodeType):
                walk(const, g, freemap, depth)

    walk(fn.__code__, fn.__globals__, rra._closure_map(fn), 0)
    return vals


def complete_read_surface(fn):
    """(complete, old_surface, added) for a holder `fn`.

    `complete` is `csg.capture_surface(fn)[0] | value_reads(fn)`, so it is a
    superset of the old surface BY CONSTRUCTION -- which is why V1 checks it
    anyway: the construction is a union of two independent walks and a bug in
    either would show up there and nowhere else.
    """
    old = csg.capture_surface(fn)[0]
    vals = value_reads(fn)
    return old | vals, old, vals - old


def blind_class(exclude=()):
    """Repository-wide: module-level scalars read BARE inside their own module.

    Returns (total, blind) as sets of (module, attr).  `total` is every
    module-level int/float attribute of every imported repo module; `blind` is
    the subset some function of that same module loads as a bare global.  This
    is computed from the module table and not from the four holders, so its
    agreement with `value_reads` (V9a) is a cross-check of two separate walks.
    """
    mods = _modules()
    exclude = set(exclude)
    total = set()
    for mname, mod in mods.items():
        for attr, val in vars(mod).items():
            if attr.startswith('__') or isinstance(val, bool) \
                    or not isinstance(val, (int, float)):
                continue
            total.add((mname, attr))
    blind = set()
    for mname, mod in mods.items():
        for _, fn in vars(mod).items():
            if not isinstance(fn, types.FunctionType) \
                    or fn.__globals__ is not mod.__dict__:
                continue
            stack, seen = [fn.__code__], set()
            while stack:
                co = stack.pop()
                if id(co) in seen:
                    continue
                seen.add(id(co))
                ins = list(dis.get_instructions(co))
                for i, x in enumerate(ins):
                    if x.opname not in ('LOAD_GLOBAL', 'LOAD_NAME'):
                        continue
                    nxt = ins[i + 1] if i + 1 < len(ins) else None
                    if nxt is not None and nxt.opname in rra._ATTR_OPS:
                        continue
                    if (mname, x.argval) in total:
                        blind.add((mname, x.argval))
                for c in co.co_consts:
                    if isinstance(c, types.CodeType):
                        stack.append(c)
    return total - exclude, blind - exclude


# ---------------------------------------------------------------------------
# THE STRENGTHENED GUARD
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def require_complete_capture_adjudication(label, holder, adjudication,
                                          report=None):
    """`csg.require_capture_adjudication` plus the names it cannot see.

    Strictly stronger: every refusal of the 2026-10-24 guard still applies,
    because that guard is entered unchanged at the end.  What is added is a
    refusal for the names on the COMPLETED surface that the old one omitted,
    with the same four verdicts and the same counterfactual rules.

    `UNREACHED` is not available on an added name and is refused outright: an
    added name is one some reachable code loads, so the one thing `UNREACHED`
    asserts -- that the holder never looks it up -- is exactly what the walk
    that produced the name contradicts.  A name that genuinely is not reached
    does not appear in `value_reads` and so is never asked about.
    """
    complete, old, added = complete_read_surface(holder)
    flip_rebinds = set()
    for _, fac in rra.FLIPS:
        flip_rebinds |= rra.rebind_set(fac)

    missing = sorted(added - set(adjudication))
    if missing:
        raise UnadjudicatedValueRead(
            '%s: %d name(s) reach this value through a BARE SCALAR read that '
            'the static surface cannot see, and carry no verdict: %s.  The '
            "read surface's negatives are not sound (2026-10-25); absence "
            'from it is not evidence that a name is unreached.'
            % (label, len(missing), ', '.join(_dotted(m) for m in missing)))

    rows, wrong = [], []
    for site in sorted(added):
        spec = adjudication[site]
        if not (isinstance(spec, tuple) and len(spec) == 2):
            raise UnadjudicatedValueRead(
                '%s: adjudication for %s must be a (verdict, reason) pair'
                % (label, _dotted(site)))
        verdict, reason = spec
        if verdict not in VERDICTS:
            raise UnadjudicatedValueRead(
                '%s: adjudication for %s has verdict %r, not one of %s'
                % (label, _dotted(site), verdict, '/'.join(VERDICTS)))
        if not reason:
            raise UnadjudicatedValueRead(
                '%s: adjudication for %s has no stated reason'
                % (label, _dotted(site)))

    for site in sorted(added):
        verdict, _ = adjudication[site]
        base, ident, ident_exact = csg.identity_leverage(holder, site)
        if not ident_exact:
            wrong.append((site, verdict,
                          'the IDENTITY rebind of %s moved the value (%s -> '
                          '%s): the probe perturbs what it measures'
                          % (_dotted(site), base.hex(), ident.hex())))
            continue
        base, moved, exact = csg.leverage(holder, site)
        rows.append((site, verdict, base, moved, exact))
        if verdict == UNREACHED:
            wrong.append((site, verdict,
                          'claimed UNREACHED, but %s was found by a walk of '
                          'the code that produces this value, so it IS on '
                          'the (completed) read surface' % _dotted(site)))
        elif verdict == PROTECTED and not exact:
            wrong.append((site, verdict,
                          'claimed PROTECTED but a rebind of %s moved the '
                          'value by %+.4f %%'
                          % (_dotted(site), 100.0 * (moved - base) / base)))
        elif verdict == EXPOSED and exact:
            wrong.append((site, verdict,
                          'claimed EXPOSED but a rebind of %s left the value '
                          'bitwise unchanged -- the stated reason describes '
                          'nothing' % _dotted(site)))
        elif verdict == GOVERNED and site not in flip_rebinds:
            wrong.append((site, verdict,
                          'claimed GOVERNED but %s is in no rebind set'
                          % _dotted(site)))
    if wrong:
        raise FalseValueReadVerdict(
            '%s: %d verdict(s) on bare-scalar reads disagree with the '
            'counterfactual:\n%s'
            % (label, len(wrong), '\n'.join(
                '    %-46s %-10s %s' % (_dotted(s), v, why)
                for s, v, why in wrong)))

    with csg.require_capture_adjudication(label, holder, adjudication,
                                          report) as (surface, old_rows):
        yield complete, list(old_rows) + rows


# ---------------------------------------------------------------------------
# Planted controls
# ---------------------------------------------------------------------------

def control_bare_scalar():
    """Reads `_PLANTED_BARE` BARE, from inside the module that defines it.

    The known answer, decided before anything is measured: the name must be
    ABSENT from `csg.capture_surface`, PRESENT in `value_reads`, and a rebind
    of it must move the value.  This is the mechanism in isolation, with no
    physics in it.
    """
    def holder(V):
        return _PLANTED_BARE * float(V)
    return holder


def control_qualified_scalar():
    """Reads `gfet.T` QUALIFIED.  The other half: the old surface DOES see
    this one, so `value_reads` must add nothing the old surface lacks here.
    A reader that simply reported every scalar would pass the control above
    and fail this one."""
    def holder(V):
        return float(gfet.T) * float(V)
    return holder


def control_no_scalar():
    """Reads nothing.  `value_reads` must be EMPTY -- the negative control
    that makes the two above mean something."""
    def holder(V):
        return 2.0 * float(V)
    return holder


# ---------------------------------------------------------------------------
# The holders under audit, and their adjudication
# ---------------------------------------------------------------------------

AF = ('graphene_quantum_capacitance_argument', 'A_FERMI')
CQP = ('graphene_quantum_capacitance_argument', 'CQ_PREFACTOR')
QEC = ('graphene_quantum_capacitance_argument', 'E_CHARGE')
QKB = ('graphene_quantum_capacitance_argument', 'KB')
GE = ('graphene_fet_model', 'e')
GHBAR = ('graphene_fet_model', 'hbar')
GKB = ('graphene_fet_model', 'kB')
GVF = ('graphene_fet_model', 'v_F')
TR = ('graphene_fet_model', '_thermal_rounding')

_CONSTANT = ('a physical constant imported from scipy.constants.  No flip '
             'would rebind it and no replacement writes it, so EXPOSED is '
             'not a warning -- it is the record that this name is a live '
             'lever on this value, which is what the counterfactual '
             'measures and what the static surface could not say at all.')

_PREFACTOR = ('a DERIVED model prefactor, not a constant of nature: it is '
              'the kind of name item 16\'s three charge forms differ in, it '
              'is in no rebind set, it is written by no replacement, and it '
              'was off every surface this repository computes.  The same '
              'shape as 2026-10-24\'s gfet._thermal_rounding finding, one '
              'level further in.')

# The five names that reach the CONTACT holders -- `blast._cq_consistent_
# contact` and `flip.make_contact_cq(live=False)` share them, because the
# second's captured object is the shipped C_q the first also calls.
ADDED_CONTACT = {
    GE: (EXPOSED, 'the elementary charge, read bare by gfet; +278.4911 % -- '
                  'the largest single lever found anywhere in this '
                  'repository, and it was on no surface.  ' + _CONSTANT),
    GHBAR: (EXPOSED, 'hbar, read bare by gfet; -63.5273 %.  ' + _CONSTANT),
    GKB: (EXPOSED, 'Boltzmann, read bare by gfet; +0.0242 %.  ' + _CONSTANT),
    GVF: (EXPOSED, 'the Fermi velocity: a MEASURED material parameter of '
                   'graphene with a literature spread, read bare by gfet, '
                   'carrying -63.5273 %.  Of the five names here this is '
                   'the one a thesis chapter could legitimately vary, and '
                   'it had no verdict.  ' + _PREFACTOR),
    QEC: (EXPOSED, 'qca\'s own copy of the elementary charge; -21.2518 %.  '
                   'A SECOND name for the same physical quantity as '
                   'graphene_fet_model.e above, with a different leverage on '
                   'the same value because the two enter at different '
                   'points.  ' + _CONSTANT),
}

ADDED_QCA_CHARGE = {
    AF: (EXPOSED, 'E_F/e = A_FERMI*sqrt(n); +0.0121 % here and +29.8900 % '
                  'on the capacitance replacement.  THE ANSWER to the '
                  '2026-10-24 ranked-first item: repository-wide, and '
                  'off-surface in both places it matters, because inside qca '
                  'it is spelled bare.  ' + _PREFACTOR),
    CQP: (EXPOSED, 'the C_q prefactor; +0.0227 % here and +50.4145 % on the '
                   'capacitance replacement.  ' + _PREFACTOR),
    QEC: (EXPOSED, '-0.5361 %.  ' + _CONSTANT),
    QKB: (EXPOSED, '+0.0173 %; the thermal scale kT/e the smoothing uses.  '
                   + _CONSTANT),
}

ADDED_QCA_CQ = {
    AF: (EXPOSED, '+29.8900 % -- the headline.  A derived prefactor with '
                  'nearly a third of the value riding on it, in no rebind '
                  'set, written by nothing, and invisible to every surface '
                  'this repository computed before today.  ' + _PREFACTOR),
    CQP: (EXPOSED, '+50.4145 % -- LARGER than A_FERMI, and on the same '
                   'off-surface footing.  ' + _PREFACTOR),
    QEC: (EXPOSED, '-25.3101 %.  ' + _CONSTANT),
    QKB: (EXPOSED, '+27.6751 %.  ' + _CONSTANT),
}

_UNREACHED_TR = (
    'the 2026-10-24 item asked whether this name carries its +58.3421 % '
    'onto the qca replacements.  It does NOT: a rebind of it leaves this '
    'value BITWISE unchanged, so it is a contact fact and not a '
    'repository-wide lever.  The verdict is UNREACHED and not PROTECTED '
    'because nothing is captured here, and it is asserted on the strength '
    'of the bitwise zero -- NOT on the name\'s absence from the static '
    'surface, which this module shows proves nothing.')


def _fmt(site, verdict, base, moved, exact):
    pct = 0.0 if exact else 100.0 * (moved - base) / base
    return ('    %-52s %-9s %s %s %s'
            % (_dotted(site), verdict, base.hex(), moved.hex(),
               'exact 0' if exact else '%+10.4f %%' % pct))


XFT = ('graphene_exact_charge_rederivation', 'floor_total')
COX = ('graphene_fet_model', 'C_ox')
GT = ('graphene_fet_model', 'T')
VD = ('graphene_fet_model', 'V_dirac')
NP = ('graphene_fet_model', 'n_puddle')
QFP = ('graphene_quantum_capacitance_argument', '_fixed_point')
QCT = ('graphene_quantum_capacitance_argument', 'cq_thermal')
QCTN = ('graphene_quantum_capacitance_argument', 'cq_thermal_of_n')
QNAE = ('graphene_quantum_capacitance_argument', 'n_arg_EF_thermal')
QSF = ('graphene_quantum_capacitance_argument', 'series_factor')
BCQS = ('graphene_contact_argument_blast_radius', 'cq_shipped_of_n')
BNC = ('graphene_contact_argument_blast_radius', 'n_consistent')
CI = ('graphene_contact_doping_model', 'C_interface')
GQC = ('graphene_fet_model', 'quantum_capacitance')
CDQC = ('graphene_contact_doping_model', 'quantum_capacitance')

_SOLVE = ('an ingredient of the self-consistency this replacement IS, so the '
          'dependence is meant; recorded as a measured number so that a '
          'future flip over this name is a refusal rather than a silent '
          'change of a committed figure.')

# -- the ranked-first item: the OLD surface of each of the three ------------
ADJ_CHARGE = {
    XFT: (EXPOSED, 'the puddle floor applied to the solved n; +50.0000 % -- '
                   'the largest lever on this replacement, and the standing '
                   '"whether the floor belongs inside C_q\'s argument" item '
                   'now has a number on the outside of it too.'),
    COX: (EXPOSED, 'the gate oxide capacitance sets n_es; ' + _SOLVE),
    GT: (EXPOSED, 'the lattice temperature, +0.0173 %; ' + _SOLVE),
    VD: (EXPOSED, 'the Dirac-point offset enters the overdrive; ' + _SOLVE),
    NP: (EXPOSED, '+49.2125 %, and THE FIRST NUMBER the standing "what '
                  'n_puddle actually means" item has ever had: at V=0.5 the '
                  'puddle floor is comparable to the solved density itself, '
                  'so what this name means is not a detail of the low-bias '
                  'limit.  In no rebind set.'),
    QFP: (EXPOSED, 'the form-B_T fixed-point solver; ' + _SOLVE),
    QCT: (EXPOSED, 'the thermal C_q spelling; ' + _SOLVE),
    QCTN: (EXPOSED, 'C_q at the corrected argument -- the whole of item 16; '
                    + _SOLVE),
    QNAE: (EXPOSED, 'form B_T itself; ' + _SOLVE),
    QSF: (EXPOSED, 'the series factor C_q/(C_q+C_ox); ' + _SOLVE),
    TR: (UNREACHED, _UNREACHED_TR),
}
ADJ_CHARGE.update(ADDED_QCA_CHARGE)

ADJ_CQ = {
    COX: (EXPOSED, 'reaches through the fixed point\'s n_es; ' + _SOLVE),
    GT: (EXPOSED, '+27.6751 % -- this replacement is far more '
                  'temperature-sensitive than the charge one (+0.0173 %), '
                  'because here T enters C_q directly and there it only '
                  'moves the solved n; ' + _SOLVE),
    QFP: (EXPOSED, 'the fixed-point solver; ' + _SOLVE),
    QCT: (EXPOSED, '+50.4145 %; ' + _SOLVE),
    QCTN: (EXPOSED, '+50.4145 %, identical to cq_thermal because the one '
                    'calls the other with no other route; ' + _SOLVE),
    QNAE: (EXPOSED, 'form B_T; ' + _SOLVE),
    QSF: (EXPOSED, 'the series factor; ' + _SOLVE),
    TR: (UNREACHED, _UNREACHED_TR),
}
ADJ_CQ.update(ADDED_QCA_CQ)

ADJ_CONTACT = {
    BCQS: (EXPOSED, 'the shipped C_q at the contact\'s own n; +58.3421 %, '
                    'bitwise the same move as gfet._thermal_rounding and '
                    'gfet.quantum_capacitance, because all three are the '
                    'same route reached by three names.'),
    BNC: (EXPOSED, 'the contact self-consistency; +22.4740 %; ' + _SOLVE),
    CI: (EXPOSED, 'C_interface is what the contact solve is referenced TO; '
                  + _SOLVE),
    COX: (PROTECTED, 'EXACT ZERO, and on the read surface: the contact solve '
                     'closes against C_interface, so the oxide capacitance '
                     'reaches this value not at all.  The same bitwise zero '
                     '2026-10-24 measured on the capture, re-measured here '
                     'on the replacement that has no capture -- so the zero '
                     'is a property of the contact model and not of the '
                     '2026-10-21 repair.'),
    GT: (EXPOSED, 'the lattice temperature, +0.0034 %; ' + _SOLVE),
    TR: (EXPOSED, 'the 2026-10-24 finding, re-measured on the UNCAPTURED '
                  'replacement: +58.3421 %, bitwise identical to the move it '
                  'makes on the captured one.  The capture was never what '
                  'exposed this name.'),
    GQC: (GOVERNED, 'the shipped C_q name, in flip\'s rebind set, so an '
                    'existing mechanism moves it deliberately; its +58.3421 % '
                    'here is that mechanism\'s own leverage.'),
    AF: (EXPOSED, '+58.3417 % -- and here, uniquely among the four holders, '
                  'this name IS on the old surface, because blast spells it '
                  '`qca.A_FERMI`.  Section 4 is that coincidence measured.'),
    QFP: (EXPOSED, 'the shared fixed-point solver; ' + _SOLVE),
    QSF: (EXPOSED, 'the series factor; ' + _SOLVE),
}
ADJ_CONTACT.update(ADDED_CONTACT)


def main():
    checks = []

    def check(name, ok, detail=''):
        detail = detail if isinstance(detail, str) else repr(detail)
        checks.append((name, bool(ok), detail))
        print('    %-4s %s%s' % ('PASS' if ok else 'FAIL', name,
                                 ('   -- ' + detail) if detail else ''))

    print('=' * 78)
    print("THE READ SURFACE'S NEGATIVES ARE NOT SOUND")
    print('  answers the 2026-10-24 ranked-first item, and supersedes the')
    print('  completeness claim that item was written under')
    print('=' * 78)
    print()

    shipped = gfet.quantum_capacitance
    cap = flip.make_contact_cq(shipped, live=False)
    HOLDERS = (
        ('qca._carrier_density_consistent', qca._carrier_density_consistent,
         ADJ_CHARGE),
        ('qca._quantum_capacitance_consistent',
         qca._quantum_capacitance_consistent, ADJ_CQ),
        ('blast._cq_consistent_contact', blast._cq_consistent_contact,
         ADJ_CONTACT),
    )

    # -- 0.  The false sentence -------------------------------------------
    print('-' * 78)
    print('0.  THE SENTENCE 2026-10-24 RESTED `UNREACHED` ON')
    print('-' * 78)
    print()
    print('    "an over-approximate read set has reliable NEGATIVES, so')
    print('     \'not in the static read surface\' is trustworthy."')
    print()
    h = qca._quantum_capacitance_consistent
    live, _ = rra.static_reads(h)
    b, m, ex = csg.leverage(h, AF)
    print('    qca._quantum_capacitance_consistent:')
    print('      static read surface              : %d names' % len(live))
    print('      qca.A_FERMI in it                : %s' % (AF in live))
    print('      rebinding qca.A_FERMI moves it by: %+.4f %%'
          % (100.0 * (m - b) / b))
    check('V0 a name ABSENT from the static read surface moves the value by '
          'nearly a third: the negatives are not sound',
          AF not in live and not ex,
          '%s -> %s' % (b.hex(), m.hex()))
    print()

    # -- 1.  Controls ------------------------------------------------------
    print('-' * 78)
    print('1.  THE READER, ON CASES WHOSE ANSWER IS KNOWN IN ADVANCE')
    print('-' * 78)
    print()
    cb, cq_, cn = (control_bare_scalar(), control_qualified_scalar(),
                   control_no_scalar())
    PB = _PLANTED_APPARATUS[0]
    vb, oldb = value_reads(cb), csg.capture_surface(cb)[0]
    _, _, pb_ident = identity_leverage_here(cb, PB)
    b, m, ex = leverage_here(cb, PB)
    check('V5a BARE scalar read from inside its own module: absent from the '
          'old surface, present in value_reads, and it MOVES the value',
          PB not in oldb and PB in vb and not ex and pb_ident,
          '%+.4f %%' % (100.0 * (m - b) / b))
    vq, oldq = value_reads(cq_), csg.capture_surface(cq_)[0]
    check('V5b QUALIFIED scalar read: the old surface already has it, and '
          'value_reads adds nothing -- a reader that reported every scalar '
          'would pass V5a and fail this',
          GT in oldq and not (vq - oldq), sorted(vq - oldq))
    check('V5c a holder that reads no module scalar at all: value_reads is '
          'EMPTY (the negative control)', not value_reads(cn),
          sorted(value_reads(cn)))
    print()

    # -- 2.  The gap -------------------------------------------------------
    print('-' * 78)
    print('2.  THE COMPLETED SURFACE, AND WHAT IT ADDS')
    print('-' * 78)
    print()
    all_added, super_ok, ident_ok, all_move = {}, True, True, True
    for label, fn, _ in HOLDERS + (('flip.make_contact_cq(live=False)',
                                    cap, None),):
        complete, old, added = complete_read_surface(fn)
        all_added[label] = added
        super_ok = super_ok and old <= complete
        print('    %s' % label)
        print('      old surface %2d   completed %2d   ADDED %d'
              % (len(old), len(complete), len(added)))
        for site in sorted(added):
            bb, ii, ie = csg.identity_leverage(fn, site)
            ident_ok = ident_ok and ie
            bb, mm, ee = csg.leverage(fn, site)
            all_move = all_move and not ee
            print(_fmt(site, '(added)', bb, mm, ee))
        print()
    n_added = sum(len(v) for v in all_added.values())
    check('V1 the completed surface is a SUPERSET of the old one for every '
          'holder (a union of two independent walks; a bug in either shows '
          'up only here)', super_ok)
    check('V2 the IDENTITY rebind of every added name moves the value by '
          'EXACTLY ZERO -- the probe does not perturb what it measures',
          ident_ok)
    check('V2b every one of the %d added (holder, name) pairs MOVES the '
          'value: not one is a false positive of the reader' % n_added,
          all_move)
    print()

    # -- 3.  The standing question, answered ------------------------------
    print('-' * 78)
    print('3.  THE 2026-10-24 QUESTION: ARE THESE TWO REPOSITORY-WIDE?')
    print('-' * 78)
    print()
    tr_zero, af_moves = True, True
    for label, fn, _ in HOLDERS:
        for site in (TR, AF):
            bb, mm, ee = csg.leverage(fn, site)
            print('    %-38s %-50s %s'
                  % (label, _dotted(site),
                     'exact 0' if ee else '%+10.4f %%'
                     % (100.0 * (mm - bb) / bb)))
        print()
    for label, fn, _ in HOLDERS[:2]:
        tr_zero = tr_zero and csg.leverage(fn, TR)[2]
        af_moves = af_moves and not csg.leverage(fn, AF)[2]
    check('V3 gfet._thermal_rounding is BITWISE EXACTLY ZERO on both qca '
          'replacements: a CONTACT FACT, not a repository-wide lever -- and '
          'the zero is measured, not inferred from absence', tr_zero)
    check('V3b qca.A_FERMI moves BOTH qca replacements: repository-wide, '
          'and the opposite answer to the same question', af_moves)
    print()

    # -- 4.  The same name, decided by spelling ---------------------------
    print('-' * 78)
    print('4.  ONE NAME, ONE VALUE, VISIBILITY DECIDED BY SPELLING')
    print('-' * 78)
    print()
    spell = []
    for label, fn in (('qca._carrier_density_consistent',
                       qca._carrier_density_consistent),
                      ('qca._quantum_capacitance_consistent',
                       qca._quantum_capacitance_consistent),
                      ('blast._cq_consistent_contact',
                       blast._cq_consistent_contact),
                      ('flip.make_contact_cq(live=False)', cap)):
        old = csg.capture_surface(fn)[0]
        comp = complete_read_surface(fn)[0]
        bb, mm, ee = csg.leverage(fn, AF)
        spell.append((label, AF in old, AF in comp))
        print('    %-38s old=%-5s completed=%-5s %+10.4f %%'
              % (label, AF in old, AF in comp, 100.0 * (mm - bb) / bb))
    print()
    print('    `qca` reads A_FERMI bare; `blast` spells it `qca.A_FERMI`.')
    check('V4 qca.A_FERMI is OFF the old surface for both qca holders and ON '
          'it for both contact holders -- an exact binary prediction of the '
          'mechanism, on one name',
          [s[1] for s in spell] == [False, False, True, True])
    check('V4b and it is on the COMPLETED surface for all four',
          all(s[2] for s in spell))
    print()

    # -- 5.  What the old guard accepts -----------------------------------
    print('-' * 78)
    print('5.  THE OLD GUARD ACCEPTS AN ADJUDICATION MISSING EVERY LEVER')
    print('-' * 78)
    print()
    for label, fn, adj in HOLDERS:
        added = all_added[label]
        thin = {k: v for k, v in adj.items() if k not in added}
        try:
            with csg.require_capture_adjudication(label, fn, thin):
                pass
            ok_old = True
            why = 'accepted, with %d unadjudicated lever(s)' % len(added)
        except Exception as exc:
            ok_old, why = False, '%s: %s' % (type(exc).__name__, exc)
        check('V6 %-36s the 2026-10-24 guard ACCEPTS it' % label,
              ok_old, why)
    print()
    for label, fn, adj in HOLDERS:
        added = all_added[label]
        thin = {k: v for k, v in adj.items() if k not in added}
        exc = None
        try:
            with require_complete_capture_adjudication(label, fn, thin):
                pass
        except UnadjudicatedValueRead as e:
            exc = e
        check('V7 %-36s the completed guard REFUSES it' % label,
              exc is not None and all(_dotted(a) in str(exc) for a in added),
              str(exc).splitlines()[0][:66] if exc else 'no exception')
    print()
    for label, fn, adj in HOLDERS:
        exc = None
        try:
            with require_complete_capture_adjudication(label, fn, adj) as (
                    surf, rows):
                n = len(surf)
        except Exception as e:
            exc = e
        check('V8 %-36s the FULL adjudication is accepted (%s)'
              % (label, 'surface %d' % n if exc is None else 'raised'),
              exc is None,
              '' if exc is None else '%s: %s' % (type(exc).__name__,
                                                 str(exc)[:200]))
    print()

    # -- 6.  Refusal branches of the new code -----------------------------
    print('-' * 78)
    print('6.  A PLANTED CONTROL FOR EVERY NEW REFUSAL BRANCH')
    print('-' * 78)
    print()
    fn, adj = blast._cq_consistent_contact, ADJ_CONTACT

    def expect(name, mutate, exctype):
        bad = dict(adj)
        mutate(bad)
        exc = None
        try:
            with require_complete_capture_adjudication('control', fn, bad):
                pass
        except Exception as e:
            exc = e
        check(name, isinstance(exc, exctype),
              type(exc).__name__ if exc else 'no exception')

    expect('V9a a malformed (verdict, reason) spec on an added name is '
           'refused', lambda d: d.__setitem__(GE, EXPOSED),
           UnadjudicatedValueRead)
    expect('V9b a verdict outside the vocabulary is refused',
           lambda d: d.__setitem__(GE, ('SAFE', 'a word that is not one')),
           UnadjudicatedValueRead)
    expect('V9c an empty reason is refused',
           lambda d: d.__setitem__(GE, (EXPOSED, '')),
           UnadjudicatedValueRead)
    expect('V9d UNREACHED on a name the walk FOUND is refused -- the one '
           'verdict an added name cannot carry',
           lambda d: d.__setitem__(GE, (UNREACHED, 'claims it is not read')),
           FalseValueReadVerdict)
    expect('V9e PROTECTED on an added name that moves is refused',
           lambda d: d.__setitem__(GE, (PROTECTED, 'claims a zero')),
           FalseValueReadVerdict)
    expect('V9f GOVERNED on an added name in no rebind set is refused',
           lambda d: d.__setitem__(GE, (GOVERNED, 'claims a mechanism')),
           FalseValueReadVerdict)
    print()
    print('    ONE BRANCH HERE HAS NO CONTROL, stated rather than omitted.')
    print('    `EXPOSED` on an added name whose rebind leaves the value')
    print('    BITWISE unchanged is refused by this module and that refusal')
    print('    is never reached above, because all 18 added pairs move.  To')
    print('    reach it needs an added name with exact-zero leverage, and the')
    print('    only ones available are the planted controls, whose sites')
    print('    `csg.leverage` cannot resolve (Section 1).  This is a new')
    print('    instance of the standing "give every verdict in every')
    print('    vocabulary a planted control" item, found in this session\'s')
    print('    own code, and it is the sixth.')
    print()

    # -- 7.  Repository-wide -----------------------------------------------
    print('-' * 78)
    print('7.  REPOSITORY-WIDE EXTENT OF THE BLIND CLASS')
    print('-' * 78)
    print()
    total, blind = blind_class()
    t_ex, b_ex = blind_class(exclude=_PLANTED_APPARATUS)
    print('    module-level scalars in imported repo modules : %d' % len(total))
    print('    of those read BARE inside their own module    : %d  (%.0f %%)'
          % (len(blind), 100.0 * len(blind) / len(total)))
    print('    with this module\'s own planted constants out  : %d of %d'
          % (len(b_ex), len(t_ex)))
    print()
    by_mod = {}
    for mname, attr in sorted(b_ex):
        by_mod.setdefault(mname, []).append(attr)
    for mname in sorted(by_mod, key=lambda k: (-len(by_mod[k]), k))[:8]:
        print('      %-50s %2d' % (mname, len(by_mod[mname])))
    print()
    measured = set()
    for v in all_added.values():
        measured |= v
    check('V10a all %d distinct names behind the 18 added (holder, name) '
          'pairs are in the repository-wide blind class -- the per-holder '
          'walk and the module-table census agree exactly' % len(measured),
          measured <= blind, sorted(measured - blind))
    TOX = ('graphene_fet_model', 't_ox')
    EPS = ('graphene_fet_model', 'epsilon_0')
    check('V10d the census is not just reporting every scalar: gfet.t_ox and '
          'gfet.epsilon_0 are module-level scalars and are NOT in the blind '
          'class, because they are read only at IMPORT time to build C_ox '
          'and never inside a function body',
          {TOX, EPS} <= total and not ({TOX, EPS} & blind))
    check('V10e so the blind class is a STRICT subset: %d of %d scalars are '
          'outside it' % (len(total - blind), len(total)),
          blind < total)
    check('V10b the census FINDS this module\'s own planted constant, so it '
          'is not blind to the pattern it censuses',
          set(_PLANTED_APPARATUS) <= blind)
    own = {site for site in blind if site[0] == _THIS_MODULE}
    check('V10c the apparatus exclusion is exactly the blind names the census '
          'found IN THIS MODULE -- checked against the census, not against '
          'itself',
          set(_PLANTED_APPARATUS) == own,
          'exclusion %s vs own %s' % (sorted(_PLANTED_APPARATUS), sorted(own)))
    check('V10f and every excluded name BELONGS to this module, whose name '
          'comes from __file__ and not from the exclusion: the exclusion '
          'cannot reach a name in a model module, so it cannot be used to '
          'hide a finding',
          all(m == _THIS_MODULE for m, _ in _PLANTED_APPARATUS),
          sorted(m for m, _ in _PLANTED_APPARATUS if m != _THIS_MODULE))
    print()

    print('=' * 78)
    npass = sum(1 for _, ok, _ in checks if ok)
    print('VALIDATION: %d/%d' % (npass, len(checks)))
    for name, ok, detail in checks:
        if not ok:
            print('   FAILED: %s   %s' % (name, detail))
    print('=' * 78)
    return 0 if npass == len(checks) else 1


if __name__ == '__main__':
    sys.exit(main())
