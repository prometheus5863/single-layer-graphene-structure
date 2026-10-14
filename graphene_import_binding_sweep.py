"""
graphene_import_binding_sweep.py

Chapter 7 Section 7.9 -- "SWEEP EVERY IMPORT-TIME BINDING OF EVERY PATCHED
NAME", created 2026-10-13 and ranked SECOND there, immediately below item 16,
on the stated ground that it is "what makes item 16's flip verifiable rather
than hopeful".

The defect it generalises has now appeared three times, each time found by
hand and repaired as a one-off:

  2026-10-03  the frozen-default fault, in its monkey-patch form.
  2026-10-10  `crx.carrier_density` -- contact_resistance_crossover.py does
              `from graphene_fet_model import carrier_density`, so rebinding
              only `gfet.carrier_density` left Section 4.4 running the
              committed form while the harness reported it as re-derived.
              Repaired in on_exact_charge() and on_consistent_argument().
  2026-10-13  `cd.quantum_capacitance` -- the same shape, one module over,
              found by hand while scoping item 16.  Not repaired; the
              contact model was simply moved outside the harness.

The 2026-10-13 entry states the rule the third occurrence forced:

    A REPAIR TO ONE BINDING IS NOT A REPAIR TO THE CLASS.

and the obvious consequence, which is this module: nobody has ever enumerated
the bindings of `sheet_conductivity`, `n_puddle`, `C_ox`, `mu`, `lambda_decay`
or `W_GRAPHENE`, any one of which a future audit will want to substitute.
Three hand-found instances in ten days is a rate, not a coincidence.

WHAT THIS MODULE DOES

  Section 1  An AST pass over every module in the repository, emitting the
             complete table of module-level bindings that alias a name out of
             another repository module.  Static, no import required, and it
             sees modules this session never runs.

  Section 2  The same table built a second way, at RUNTIME, by object
             identity: import the modules and find every attribute of every
             module that IS the owner's object.  The two methods have
             different blind spots -- the AST pass cannot see a binding made
             by exec or by a conditional import, and identity cannot
             distinguish a shared binding from an independently-created equal
             object -- so each is the other's control, and every disagreement
             is reported rather than reconciled.

  Section 3  Every rebinding site in the repository (`mod.name = ...`), its
             enclosing function, and for each one the REQUIRED rebind set
             from Sections 1-2.  A rebinding that does not cover its required
             set is a MISS.

  Section 4  Whether each miss is LIVE or LATENT, measured and not argued:
             enter the manager, run it, and count arrivals at the stale
             binding with a counter installed INSIDE the stale holder.  The
             2026-10-13 rule applies to this module's own controls -- a
             counter that lives in the patch cannot see whether the shipped
             object was reached, so each counter is installed in the object
             that would be reached if the patch had missed.

  Section 5  Validations against exactly known values, not plausible ranges
             (2026-09-17).  Six, listed at the top of that section, of which
             three are facts the repository's own history fixes exactly:
             the 2026-10-10 repair must show ZERO misses for its pair, the
             2026-10-13 hand-found miss must be REPRODUCED, and a planted
             binding in a throwaway module must be found exactly once.
"""

import ast
import contextlib
import importlib
import io
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.abspath(__file__))

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


# ---------------------------------------------------------------------------
# Section 1: the static binding table
# ---------------------------------------------------------------------------

def repo_modules():
    """Every importable module name in the repository root.

    The first draft of this function skipped every file whose name begins
    with an underscore, on the unexamined ground that such a file is
    private.  Section 5's planted-binding control FAILED on it: the plant
    was written to `_binding_sweep_planted_probe.py` and the sweep reported
    zero bindings in a file it had itself just created with one.  A module
    whose name begins with an underscore is as importable as any other and
    can hold exactly the bindings this module exists to enumerate, so the
    filter was a blind spot and not a convention.  Only dunder names are
    excluded now, and the planted control is what found it -- the fourth
    time in this repository that a control has failed on the instrument
    rather than on the thing measured.
    """
    out = []
    for fn in sorted(os.listdir(REPO)):
        if fn.endswith('.py') and not fn.startswith('__'):
            out.append(fn[:-3])
    return out


def _is_module_level(node, tree):
    return node in tree.body


def static_bindings(modnames):
    """Module-level `from M import n [as a]` where M is a repo module.

    Returns (bindings, deferred, module_aliases, parse_errors) where

      bindings       [(holder, local_name, owner, owner_name)]  -- the
                     dangerous kind: a reference COPIED at import time.
      deferred       the same, but inside a function body: re-evaluated on
                     every call, so a patch of the owner's attribute DOES
                     reach it.  Reported separately and not counted as a
                     binding to rebind.
      module_aliases [(holder, alias, owner)] from `import M [as a]`: the
                     SAFE form, since the attribute is looked up at call time.
    """
    bindings = []
    deferred = []
    aliases = []
    errors = []
    known = set(modnames)
    for m in modnames:
        path = os.path.join(REPO, m + '.py')
        try:
            with open(path, 'r', encoding='utf-8') as f:
                tree = ast.parse(f.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError) as exc:
            errors.append((m, repr(exc)))
            continue
        # Mark every node that is lexically at module level.
        at_top = set()

        def mark(body):
            for n in body:
                at_top.add(id(n))
                for sub in ('body', 'orelse', 'finalbody'):
                    if isinstance(getattr(n, sub, None), list) and not isinstance(
                            n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                        mark(getattr(n, sub))

        mark(tree.body)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module in known:
                for al in node.names:
                    rec = (m, al.asname or al.name, node.module, al.name)
                    if id(node) in at_top:
                        bindings.append(rec)
                    else:
                        deferred.append(rec)
            elif isinstance(node, ast.Import):
                for al in node.names:
                    if al.name in known:
                        aliases.append((m, al.asname or al.name, al.name))
    return bindings, deferred, aliases, errors


# ---------------------------------------------------------------------------
# Section 3 (static half): the rebinding sites
# ---------------------------------------------------------------------------

def rebind_sites(modnames):
    """Every `<alias>.<name> = ...` assignment, with its enclosing function.

    Returns [(holder_module, enclosing_function, alias, name, lineno)].
    """
    out = []
    for m in modnames:
        path = os.path.join(REPO, m + '.py')
        try:
            with open(path, 'r', encoding='utf-8') as f:
                tree = ast.parse(f.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError):
            continue
        parents = {}
        for node in ast.walk(tree):
            for ch in ast.iter_child_nodes(node):
                parents[id(ch)] = node

        def enclosing(node):
            cur = node
            while id(cur) in parents:
                cur = parents[id(cur)]
                if isinstance(cur, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    return cur.name
            return '<module>'

        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            for tgt in node.targets:
                if (isinstance(tgt, ast.Attribute)
                        and isinstance(tgt.value, ast.Name)):
                    out.append((m, enclosing(node), tgt.value.id, tgt.attr,
                                node.lineno))
    return out



def class_rebind_sites(modnames, aliases):
    """Every `<bc>.rebind(<module alias>, '<name>', ...)` call site.

    ADDED 2026-10-14, in the same session that created the thing it looks for,
    and for a reason worth stating: the first version of this sweep found the
    repository's rebinding sites by looking for `alias.attr = ...`.  Replacing
    three of those with graphene_binding_class.rebind() made them INVISIBLE to
    the sweep, and validation (2) -- which exists to keep 2026-10-13's
    hand-found miss on the record -- failed because the manager it names had
    stopped being a rebinding site in the only sense the sweep understood.

    A repair that blinds the audit is not a repair.  So the sweep now sees
    both forms, and reports a class-complete site as covered BY CONSTRUCTION
    with its exclusions listed, rather than as absent.
    """
    out = []
    alias_of = {}
    for holder, alias, owner in aliases:
        alias_of.setdefault(holder, {})[alias] = owner
    for m in modnames:
        path = os.path.join(REPO, m + '.py')
        try:
            with open(path, 'r', encoding='utf-8') as f:
                tree = ast.parse(f.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError):
            continue
        parents = {}
        for node in ast.walk(tree):
            for ch in ast.iter_child_nodes(node):
                parents[id(ch)] = node

        def enclosing(node):
            cur = node
            while id(cur) in parents:
                cur = parents[id(cur)]
                if isinstance(cur, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    return cur.name
            return '<module>'

        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            if not (isinstance(fn, ast.Attribute) and fn.attr == 'rebind'
                    and isinstance(fn.value, ast.Name)):
                continue
            if alias_of.get(m, {}).get(fn.value.id) != 'graphene_binding_class':
                continue
            if len(node.args) < 2:
                continue
            a0, a1 = node.args[0], node.args[1]
            owner = (alias_of.get(m, {}).get(a0.id)
                     if isinstance(a0, ast.Name) else None)
            name = a1.value if isinstance(a1, ast.Constant) else None
            if owner is None or not isinstance(name, str):
                continue
            excl = [k.arg for k in node.keywords if k.arg == 'exclude_reasons']
            out.append((m, enclosing(node), owner, name, bool(excl),
                        node.lineno))
    return out

# ---------------------------------------------------------------------------
# Section 2: the runtime binding table, by object identity
# ---------------------------------------------------------------------------

# Modules imported for the runtime half.  Deliberately NOT every module in the
# repository: this list is exactly the set that the static table of Section 1
# names as a holder or an owner of a binding, plus the three modules that own
# a patch manager.  Importing the rest would add import-time cost without
# adding a binding, and an import failure in an unrelated module would be
# reported as a blind spot in this one.
def runtime_universe(bindings, aliases):
    names = set()
    for holder, _local, owner, _oname in bindings:
        names.add(holder)
        names.add(owner)
    for holder, _alias, owner in aliases:
        names.add(holder)
        names.add(owner)
    return sorted(names)


def import_universe(names):
    mods = {}
    failures = []
    for n in names:
        try:
            mods[n] = importlib.import_module(n)
        except BaseException as exc:          # noqa: BLE001 -- reported, not hidden
            failures.append((n, type(exc).__name__ + ': ' + str(exc)[:120]))
    return mods, failures


def runtime_bindings(mods, owner_names):
    """Find every module attribute that IS one of the owners' objects.

    owner_names: {owner_module: [attribute names to trace]}
    Returns [(holder, local_name, owner, owner_name)].
    """
    targets = []
    for owner, attrs in owner_names.items():
        if owner not in mods:
            continue
        for a in attrs:
            if hasattr(mods[owner], a):
                targets.append((owner, a, getattr(mods[owner], a)))
    out = []
    for holder, mod in sorted(mods.items()):
        d = vars(mod)
        for local, val in list(d.items()):
            if local.startswith('__'):
                continue
            for owner, oname, obj in targets:
                if holder == owner and local == oname:
                    continue              # the owner's own definition
                if val is obj:
                    out.append((holder, local, owner, oname))
    return out


# ---------------------------------------------------------------------------
# Section 4: arrival counting, installed in the STALE holder
# ---------------------------------------------------------------------------

_ARRIVALS = {}


def _counting(key, fn):
    """Wrapper that records an arrival under `key` and then calls `fn`.

    2026-10-13's rule, applied to this module: a counter that lives in the
    REPLACEMENT cannot see whether the shipped object was reached.  Every
    counter here is installed in the object that would be reached if the
    patch had missed a binding, so a nonzero count IS the miss being used.
    """
    _ARRIVALS.setdefault(key, 0)

    def wrapped(*a, **k):
        _ARRIVALS[key] += 1
        return fn(*a, **k)

    wrapped.__name__ = getattr(fn, '__name__', 'wrapped')
    wrapped.__wrapped_by_sweep__ = True
    return wrapped


@contextlib.contextmanager
def count_at(mods, holder, local, fn):
    """Install a counting wrapper at mods[holder].<local> for the duration."""
    key = '%s.%s' % (holder, local)
    saved = getattr(mods[holder], local)
    _ARRIVALS[key] = 0
    try:
        setattr(mods[holder], local, _counting(key, fn))
        yield key
    finally:
        setattr(mods[holder], local, saved)


@contextlib.contextmanager
def identity_patch_all(mods, bindings, owner_names):
    """Wrap EVERY traced binding (and every owner definition) in a wrapper
    that returns the shipped value unchanged.

    The neutrality control of 2026-10-13 Section 6: if any number moves under
    this, the instrument is the perturbation and nothing it measures is
    trustworthy.  Only callables are wrapped; a float binding has nothing to
    wrap and is left alone.
    """
    saved = []
    try:
        sites = set()
        for owner, attrs in owner_names.items():
            for a in attrs:
                sites.add((owner, a))
        for holder, local, _owner, _oname in bindings:
            sites.add((holder, local))
        for holder, local in sorted(sites):
            if holder not in mods or not hasattr(mods[holder], local):
                continue
            val = getattr(mods[holder], local)
            if not callable(val):
                continue
            saved.append((holder, local, val))

            def make(f):
                def ident(*a, **k):
                    return f(*a, **k)
                return ident

            setattr(mods[holder], local, make(val))
        yield len(saved)
    finally:
        for holder, local, val in saved:
            setattr(mods[holder], local, val)


# ---------------------------------------------------------------------------
# Probes: the numbers a miss can move, taken from the STALE HOLDER itself
# ---------------------------------------------------------------------------

def _flat(obj):
    """Flatten whatever a probe returns into a tuple of floats."""
    out = []

    def walk(o):
        if isinstance(o, dict):
            for k in sorted(o.keys(), key=str):
                walk(o[k])
        elif isinstance(o, (list, tuple)):
            for x in o:
                walk(x)
        elif isinstance(o, np.ndarray):
            for x in np.ravel(o):
                out.append(float(x))
        elif isinstance(o, (int, float, np.floating, np.integer)):
            out.append(float(o))
    walk(obj)
    return tuple(out)


def probe_registry(mods):
    """One cheap, deterministic numeric probe per holder module.

    A holder with no entry here is reported UNMEASURED rather than guessed
    at: the point of Section 4 is a measured verdict, and a probe invented to
    fill a row would be a verdict about the probe.
    """
    P = {}

    cdm = mods.get('graphene_contact_doping_model')
    if cdm is not None:
        def p_cd():
            ws = sorted(cdm.METAL_WORK_FUNCTIONS.values())
            vals = [cdm.contact_edge_carrier_density(w) for w in ws]
            vals.append(cdm._sheet_conductivity_puddle_regularized(1.0e15))
            vals.append(cdm.junction_extra_resistance(5.1, 2.0e16, n_points=60))
            return _flat(vals)
        P['graphene_contact_doping_model'] = p_cd

    crx = mods.get('contact_resistance_crossover')
    if crx is not None:
        P['contact_resistance_crossover'] = lambda: _flat(
            [crx.channel_resistance_at_length(200e-9, 2.0),
             crx.channel_resistance_at_length(1.0e-6, 3.5)])

    pcm = mods.get('graphene_photodetector_collection_model')
    if pcm is not None:
        P['graphene_photodetector_collection_model'] = lambda: _flat(
            pcm.modeled_eqe_by_metal(n_points=201))

    psc = mods.get('graphene_photodetector_signed_carrier_model')
    if psc is not None:
        P['graphene_photodetector_signed_carrier_model'] = lambda: _flat(
            psc.net_response(4.5, 5.4, n_points=401))

    ptc = mods.get('graphene_photodetector_two_contact_model')
    if ptc is not None:
        P['graphene_photodetector_two_contact_model'] = lambda: _flat(
            ptc.net_response(4.5, 5.4, n_points=401))

    return P


@contextlib.contextmanager
def close_misses(mods, misses):
    """Rebind each stale holder to whatever the OWNER currently holds.

    Must be entered INSIDE the manager under test: the replacement is read
    out of the owner's attribute at entry, so this closes the miss with
    exactly the object the manager installed, without this module needing to
    know what that object is.
    """
    saved = []
    try:
        for owner, oname, holder, local in misses:
            if owner not in mods or holder not in mods:
                continue
            repl = getattr(mods[owner], oname)
            saved.append((holder, local, getattr(mods[holder], local)))
            setattr(mods[holder], local, repl)
        yield len(saved)
    finally:
        for holder, local, val in saved:
            setattr(mods[holder], local, val)


# ---------------------------------------------------------------------------
# The managers whose rebinds can be ENTERED and therefore measured
# ---------------------------------------------------------------------------

MANAGERS = [
    ('graphene_quantum_capacitance_argument.on_consistent_argument',
     'graphene_quantum_capacitance_argument', 'on_consistent_argument', {}),
    ('graphene_exact_charge_rederivation.on_exact',
     'graphene_exact_charge_rederivation', 'on_exact', {}),
    ('graphene_double_regularisation.regularisers(floor=False)',
     'graphene_double_regularisation', 'regularisers', {'floor': False}),
    ('graphene_contact_argument_blast_radius.contact_on_consistent_argument',
     'graphene_contact_argument_blast_radius', 'contact_on_consistent_argument',
     {}),
]


def sections_1_to_3():
    print('=' * 78)
    print('IMPORT-BINDING SWEEP -- every module-level alias of every patchable name')
    print('Section 7.9, ranked second on 2026-10-13: the precondition for item 16')
    print('=' * 78)
    print()

    mods_all = repo_modules()
    print('Repository modules parsed: %d' % len(mods_all))

    # -- Section 1 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('SECTION 1 -- THE STATIC BINDING TABLE (AST, no import)')
    print('-' * 78)
    bindings, deferred, aliases, perr = static_bindings(mods_all)
    for m, e in perr:
        print('  [PARSE ERROR, reported not hidden] %s: %s' % (m, e))

    by_owner = {}
    for holder, local, owner, oname in bindings:
        by_owner.setdefault((owner, oname), []).append((holder, local))

    print('  module-level from-imports of repo names : %d' % len(bindings))
    print('  function-level (deferred, SAFE)         : %d' % len(deferred))
    print('  `import M as a` module aliases (SAFE)   : %d' % len(aliases))
    print('  distinct (owner, name) pairs aliased    : %d' % len(by_owner))
    print()
    print('  Every owner.name with MORE THAN ONE live binding -- i.e. every name')
    print('  for which a single `setattr(owner, name, ...)` is NOT sufficient:')
    print()
    print('  %-34s %-26s %s' % ('owner.name', 'stale holders', 'n'))
    multi = []
    for (owner, oname), hs in sorted(by_owner.items()):
        multi.append(((owner, oname), hs))
    for (owner, oname), hs in multi:
        print('  %-34s %-26s %d' % ('%s.%s' % (owner.replace('graphene_', 'g_'), oname),
                                    ', '.join(h.replace('graphene_', 'g_') for h, _l in hs)[:26],
                                    len(hs)))
    print()
    print('  Total names requiring a multi-site rebind: %d' % len(multi))
    print('  Total stale sites across those names     : %d' % len(bindings))

    # -- Section 2 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('SECTION 2 -- THE SAME TABLE BY RUNTIME OBJECT IDENTITY (the control)')
    print('-' * 78)
    universe = runtime_universe(bindings, aliases)
    mods, ifail = import_universe(universe)
    print('  modules imported: %d of %d' % (len(mods), len(universe)))
    for n, e in ifail:
        print('  [IMPORT FAILED, reported not hidden] %s: %s' % (n, e))

    owner_names = {}
    for (owner, oname) in by_owner:
        owner_names.setdefault(owner, []).append(oname)
    rt = runtime_bindings(mods, owner_names)
    s_static = set((h, l, o, n) for h, l, o, n in bindings if o in mods and h in mods)
    s_rt = set(rt)
    only_static = sorted(s_static - s_rt)
    only_rt = sorted(s_rt - s_static)
    print('  static (importable subset): %d   runtime: %d   agreed: %d'
          % (len(s_static), len(s_rt), len(s_static & s_rt)))
    if only_static:
        print('  seen STATICALLY but not at runtime (rebound already, or shadowed):')
        for h, l, o, n in only_static:
            print('      %s.%s  <- %s.%s' % (h, l, o, n))
    if only_rt:
        print('  seen at RUNTIME but not statically (identity coincidence, or a')
        print('  binding this AST pass cannot see -- each one is a blind spot):')
        for h, l, o, n in only_rt:
            print('      %s.%s  is  %s.%s' % (h, l, o, n))
    if not only_static and not only_rt:
        print('  the two methods agree exactly.')

    # -- Section 3 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('SECTION 3 -- EVERY REBINDING SITE, AND WHAT IT FAILS TO COVER')
    print('-' * 78)
    sites = rebind_sites(mods_all)
    alias_map = {}
    for holder, alias, owner in aliases:
        alias_map[(holder, alias)] = owner
    groups = {}
    for holder, func, alias, name, lineno in sites:
        owner = alias_map.get((holder, alias))
        if owner is None:
            continue                 # not a repo-module attribute assignment
        groups.setdefault((holder, func), set()).add((owner, name))

    class_sites = class_rebind_sites(mods_all, aliases)
    print('  CLASS-COMPLETE rebinding sites (graphene_binding_class.rebind):')
    print('  these cover their required set by construction -- the class is')
    print('  computed at entry -- so they contribute no miss, and what matters')
    print('  about them is their EXCLUSIONS, listed with the reason in')
    print('  graphene_binding_class.CONTACT_OPERATING_POINT.')
    print()
    print('  %-52s %-26s %s' % ('site', 'name', 'exclusions?'))
    for holder, func, owner, name, excl, lineno in sorted(class_sites):
        print('  %-52s %-26s %s'
              % (('%s.%s' % (holder.replace('graphene_', 'g_'), func))[:52],
                 '%s.%s' % (owner.replace('graphene_', 'g_'), name),
                 'yes' if excl else 'none'))
    print('  class-complete sites: %d' % len(class_sites))
    print()
    miss_table = []
    print('  HAND-LISTED rebinding sites (`alias.name = ...`), and their misses:')
    print('  %-52s %5s %5s' % ('rebinding site (module.function)', 'rebnd', 'MISS'))
    for (holder, func), pairs in sorted(groups.items()):
        covered = set(pairs)
        req = set()
        for (owner, name) in pairs:
            req.add((owner, name))
            for h, l in by_owner.get((owner, name), []):
                req.add((h, l))
        missing = sorted(req - covered)
        label = '%s.%s' % (holder.replace('graphene_', 'g_'), func)
        print('  %-52s %5d %5d' % (label[:52], len(covered), len(missing)))
        for (h, l) in missing:
            owners = [(o, n) for (o, n) in pairs if n == l]
            owner, oname = owners[0] if owners else ('?', l)
            miss_table.append((holder, func, owner, oname, h, l))
    print()
    print('  TOTAL STATIC MISSES: %d, across %d rebinding sites'
          % (len(miss_table), len(set((a, b) for a, b, _c, _d, _e, _f in miss_table))))
    print('  ANNOTATION, same day: this total is an UPPER BOUND and the number')
    print('  is retained rather than replaced.  A static import graph asks')
    print('  whether the repository contains a second binding; a miss needs the')
    print('  holder to be IMPORTED IN THE PATCHING PROCESS.  Section 4c measures')
    print('  that per site, and it removes rows -- including the one this section')
    print('  calls the third unrepaired instance of the 2026-10-10 fault.')
    print()
    print('  %-44s %s' % ('site', 'stale binding it leaves behind'))
    for holder, func, owner, oname, h, l in miss_table:
        print('  %-44s %s.%s  (shadows %s.%s)'
              % (('%s.%s' % (holder.replace('graphene_', 'g_'), func))[:44],
                 h.replace('graphene_', 'g_'), l, owner.replace('graphene_', 'g_'), oname))
    return (mods, bindings, by_owner, miss_table, groups, alias_map,
            class_sites)


# ---------------------------------------------------------------------------
# Section 1b: transitive aliases -- `LOCAL = imported_name` at module level
# ---------------------------------------------------------------------------

def transitive_aliases(modnames, bindings, aliases):
    """Module-level `NEW = OLD` where OLD is itself a binding or an owner name.

    Added after Section 2's identity control found fifteen bindings this
    file's first AST rule could not see.  Every one of them was a RENAMING
    alias of exactly this shape (`W_CROSS_VACUUM = W_GRAPHENE`,
    `C_OX = gfet.C_ox`), and every one of them is a real second-order binding:
    a patch of the owner's attribute does not reach it, and neither does a
    patch of the first-order holder's.  The control is kept and this rule is
    the repair -- 2026-10-11's "let the structure set the criterion", applied
    to the instrument.
    """
    first = {}
    for holder, local, owner, oname in bindings:
        first.setdefault(holder, {})[local] = (owner, oname)
    alias_of = {}
    for holder, alias, owner in aliases:
        alias_of.setdefault(holder, {})[alias] = owner
    out = []
    for m in modnames:
        path = os.path.join(REPO, m + '.py')
        try:
            with open(path, 'r', encoding='utf-8') as f:
                tree = ast.parse(f.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in tree.body:
            if not isinstance(node, ast.Assign) or len(node.targets) != 1:
                continue
            tgt = node.targets[0]
            if not isinstance(tgt, ast.Name):
                continue
            v = node.value
            if isinstance(v, ast.Name):
                src = first.get(m, {}).get(v.id)
                if src is not None:
                    out.append((m, tgt.id, src[0], src[1]))
            elif (isinstance(v, ast.Attribute) and isinstance(v.value, ast.Name)):
                owner = alias_of.get(m, {}).get(v.value.id)
                if owner is not None:
                    out.append((m, tgt.id, owner, v.attr))
    return out


# ---------------------------------------------------------------------------
# Section 4a: is the stale binding USED by its holder?
# ---------------------------------------------------------------------------

def perturb(val, factor=1.37):
    if callable(val):
        def pert(*a, **k):
            r = val(*a, **k)
            return r * factor
        return pert
    try:
        return type(val)(val * factor) if val != 0 else factor
    except (TypeError, ValueError):
        return None


def usedness(mods, probes, holder, local, factor=1.37):
    """Change ONLY the holder's own copy and see whether its probe moves.

    This is what makes a miss a defect rather than a cosmetic finding, and it
    needs no manager to be entered: if the holder's numbers do not depend on
    its local copy, no patch that misses that copy can change them.
    """
    if holder not in probes or holder not in mods:
        return None
    base = probes[holder]()
    saved = getattr(mods[holder], local)
    pert = perturb(saved, factor)
    if pert is None:
        return None
    try:
        setattr(mods[holder], local, pert)
        now = probes[holder]()
    finally:
        setattr(mods[holder], local, saved)
    if len(base) != len(now):
        return float('inf')
    worst = 0.0
    for a, b in zip(base, now):
        if a == b:
            continue
        den = abs(a) if a != 0 else 1.0
        worst = max(worst, abs(b - a) / den)
    return worst


# ---------------------------------------------------------------------------
# Section 4b: does the patched run REACH the stale holder?
# ---------------------------------------------------------------------------

def reachability(mods, mod_name, mgr_name, kwargs, stale, workload):
    """Enter the manager, run its own workload, and count arrivals at the
    STALE binding.

    The counter is installed in the stale holder, i.e. in the object that is
    reached only if the patch missed -- 2026-10-13's rule, which this module
    is required to obey about itself.  Returns (arrivals, ran_ok).
    """
    m = mods.get(mod_name)
    if m is None or not hasattr(m, mgr_name):
        return None, False
    holder, local = stale
    if holder not in mods or not hasattr(mods[holder], local):
        return None, False
    shipped = getattr(mods[holder], local)
    if not callable(shipped):
        return None, False           # a float cannot be counted; see 4a
    key = '%s.%s' % (holder, local)
    _ARRIVALS[key] = 0
    try:
        with getattr(m, mgr_name)(**kwargs):
            setattr(mods[holder], local, _counting(key, shipped))
            try:
                workload()
            finally:
                setattr(mods[holder], local, shipped)
    except BaseException as exc:     # noqa: BLE001 -- reported, not hidden
        print('        [workload raised, reported not hidden] %s: %s'
              % (type(exc).__name__, str(exc)[:90]))
        return _ARRIVALS.get(key, 0), False
    return _ARRIVALS.get(key, 0), True


def sections_1_to_4():
    (mods, bindings, by_owner, miss_table, groups, alias_map,
     class_sites) = sections_1_to_3()

    # -- Section 1b / 2 reconciliation -------------------------------------
    print()
    print('-' * 78)
    print('SECTION 1b -- THE BLIND SPOT SECTION 2 FOUND, AND ITS REPAIR')
    print('-' * 78)
    mods_all = repo_modules()
    _b, _d, aliases, _e = static_bindings(mods_all)
    trans = transitive_aliases(mods_all, bindings, aliases)
    print('  transitive (renaming) aliases the first AST rule could not see: %d'
          % len(trans))
    for h, l, o, n in sorted(trans):
        print('      %-58s <- %s.%s' % ('%s.%s' % (h, l), o, n))
    print()
    print('  These are REAL bindings, not identity coincidences: each is a')
    print('  module-level `NEW = OLD` whose value a patch of either the owner')
    print('  or the first-order holder leaves untouched.')
    trans_set = set((h, l, o, n) for h, l, o, n in trans)
    owner_names = {}
    for (owner, oname) in by_owner:
        owner_names.setdefault(owner, []).append(oname)
    for h, l, o, n in trans_set:
        owner_names.setdefault(o, [])
        if n not in owner_names[o]:
            owner_names[o].append(n)
    rt = set(runtime_bindings(mods, owner_names))
    static_union = set((h, l, o, n) for h, l, o, n in bindings
                       if o in mods and h in mods) | set(
        (h, l, o, n) for h, l, o, n in trans_set if o in mods and h in mods)
    print()
    print('  reconciliation, with the traced-name set made symmetric:')
    print('    static directed edges %d, runtime identity pairs %d, agreed %d'
          % (len(static_union), len(rt), len(static_union & rt)))
    resid_st = sorted(static_union - rt)
    resid_rt = sorted(rt - static_union)
    for h, l, o, n in resid_st:
        print('      static-only (no longer shares the object -- already rebound,')
        print('        shadowed, or copied by value): %s.%s <- %s.%s' % (h, l, o, n))
    print()
    print('  The %d remaining runtime-only pairs are not a defect in either' % len(resid_rt))
    print('  method.  Runtime identity does not yield a DIRECTED edge; it yields')
    print('  the EQUIVALENCE CLASS of every name that currently shares one')
    print('  object.  Where a class has more members than the import graph has')
    print('  edges into it, the extra members are still sites a patch must')
    print('  rebind -- so for the question this module exists to answer, the')
    print('  class is the safe criterion and the import graph is the unsafe one.')
    print()
    classes = {}
    for h, l, o, n in sorted(rt):
        classes.setdefault((o, n), set()).add((h, l))
    for h, l, o, n in sorted(static_union):
        classes.setdefault((o, n), set()).add((h, l))
    wide = [(k, v) for k, v in classes.items() if len(v) >= 2]
    print('  %-44s %s' % ('patch target', 'complete rebind set (identity class)'))
    for (o, n), members in sorted(wide, key=lambda kv: (-len(kv[1]), kv[0])):
        print('  %-44s %d sites: %s'
              % ('%s.%s' % (o.replace('graphene_', 'g_'), n), len(members),
                 ', '.join(sorted('%s.%s' % (a.replace('graphene_', 'g_'), b)
                                  for a, b in members))[:200]))

    # -- Section 4a --------------------------------------------------------
    print()
    print('-' * 78)
    print('SECTION 4a -- IS THE STALE BINDING USED?  (measured, per holder)')
    print('-' * 78)
    probes = probe_registry(mods)
    print('  holders with a probe: %d (%s)'
          % (len(probes), ', '.join(sorted(k.replace('graphene_', 'g_')
                                           for k in probes))))
    print()
    print('  %-46s %-10s %s' % ('stale binding', 'max rel dn', 'verdict'))
    used = {}
    for holder, func, owner, oname, h, l in miss_table:
        if (h, l) in used:
            continue
        u = usedness(mods, probes, h, l)
        used[(h, l)] = u
        if u is None:
            verdict = 'UNMEASURED (no probe for this holder)'
            shown = '--'
        elif u == 0.0:
            verdict = 'UNUSED -- the miss is cosmetic'
            shown = '0'
        else:
            verdict = 'USED -- the miss changes this holder''s numbers'
            shown = '%.4e' % u
        print('  %-46s %-10s %s' % ('%s.%s' % (h.replace('graphene_', 'g_'), l),
                                    shown, verdict))

    # -- Section 4b --------------------------------------------------------
    print()
    print('-' * 78)
    print('SECTION 4b -- DOES THE PATCHED RUN REACH IT?  (arrival counters)')
    print('-' * 78)
    print('  Only the four rebinding sites that are enterable context managers')
    print('  can be measured this way; the other ten are inline try/finally')
    print('  blocks inside a function and are reported UNENTERABLE rather than')
    print('  assumed either way.')
    print()
    workloads = {}
    qca = mods.get('graphene_quantum_capacitance_argument')
    if qca is not None and hasattr(qca, 'crossover_lengths'):
        workloads['graphene_quantum_capacitance_argument'] = (
            lambda: qca.crossover_lengths())
    X = mods.get('graphene_exact_charge_rederivation')
    if X is not None and hasattr(X, 'peak_rf'):
        workloads['graphene_exact_charge_rederivation'] = (
            lambda: X.peak_rf(n_vg=40))
    dr = mods.get('graphene_double_regularisation')
    if dr is not None and hasattr(dr, 'cq'):
        workloads['graphene_double_regularisation'] = (
            lambda: dr.cq(np.array([0.0, 1.2])))
    blast = mods.get('graphene_contact_argument_blast_radius')
    if blast is not None and hasattr(blast, 'residual'):
        workloads['graphene_contact_argument_blast_radius'] = (
            lambda: blast.residual(1.2, 2.529768e-02))

    reach_rows = []
    for label, mod_name, mgr_name, kwargs in MANAGERS:
        misses_here = [(h, l) for holder, func, owner, oname, h, l in miss_table
                       if holder == mod_name and func == mgr_name]
        if not misses_here:
            print('  %-56s no static miss' % label[:56])
            continue
        wl = workloads.get(mod_name)
        for (h, l) in misses_here:
            if wl is None:
                print('  %-56s %s.%s  NO WORKLOAD' % (label[:56], h, l))
                continue
            n, ok = reachability(mods, mod_name, mgr_name, kwargs, (h, l), wl)
            if n is None:
                print('  %-50s %-22s not countable (float)'
                      % (label[:50], '%s.%s' % (h.replace('graphene_', 'g_'), l)))
                reach_rows.append((label, '%s.%s' % (h, l), None))
            else:
                verdict = ('LIVE: the shipped object was reached %d time(s) '
                           'inside the patch' % n) if n > 0 else (
                    'LATENT: never reached by this workload')
                print('  %-50s %-22s %s'
                      % (label[:50], '%s.%s' % (h.replace('graphene_', 'g_'), l),
                         verdict))
                reach_rows.append((label, '%s.%s' % (h, l), n))
    return (mods, bindings, by_owner, miss_table, probes, used, trans,
            reach_rows, class_sites)


# ---------------------------------------------------------------------------
# Section 4c: is the stale holder even IMPORTED where the patch happens?
# ---------------------------------------------------------------------------

def imported_universe(mod_name):
    """The repository modules present after `import mod_name`, in a fresh
    interpreter.

    A subprocess rather than this process, because this one has already
    imported half the repository for Sections 2 and 4a-b; asking it which
    modules a patch site can see would give the sweep's own import list.
    That is the 2026-09-30 arrival rule applied to a static question: the
    instrument must not be what makes the answer true.
    """
    import subprocess
    code = ('import sys, importlib\n'
            'importlib.import_module(%r)\n'
            'import os\n'
            'REPO = %r\n'
            'out = []\n'
            'for n, m in list(sys.modules.items()):\n'
            '    f = getattr(m, "__file__", None)\n'
            '    if f and os.path.dirname(os.path.abspath(f)) == REPO:\n'
            '        out.append(n)\n'
            'print(" ".join(sorted(out)))\n') % (mod_name, REPO)
    try:
        r = subprocess.run([sys.executable, '-c', code], cwd=REPO,
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return None
    if r.returncode != 0:
        return None
    return set(r.stdout.split())


def section_4c(miss_table):
    print()
    print('-' * 78)
    print('SECTION 4c -- IS THE STALE HOLDER IMPORTED WHERE THE PATCH HAPPENS?')
    print('-' * 78)
    print('  Section 3 counts a miss whenever the REPOSITORY holds a second')
    print('  binding.  A miss can only fire if the holder is in the patching')
    print('  process\'s sys.modules.  Measured here in a fresh interpreter per')
    print('  rebinding-site module, so the sweep\'s own imports cannot supply')
    print('  the answer.')
    print()
    sites = sorted(set(r[0] for r in miss_table))
    present = {}
    for m in sites:
        u = imported_universe(m)
        present[m] = u
        if u is None:
            print('  %-46s import failed or timed out (reported, not hidden)' % m)
        else:
            print('  %-46s %d repo modules in process' % (m, len(u)))
    print()
    print('  %-44s %-30s %s' % ('rebinding site', 'stale holder', 'verdict'))
    real, absent, unknown = [], [], []
    for holder, func, owner, oname, h, l in miss_table:
        u = present.get(holder)
        if u is None:
            unknown.append((holder, func, h, l))
            v = 'UNKNOWN'
        elif h in u:
            real.append((holder, func, h, l))
            v = 'PRESENT -- a real candidate'
        else:
            absent.append((holder, func, h, l))
            v = 'ABSENT -- cannot fire in this process'
        print('  %-44s %-30s %s'
              % (('%s.%s' % (holder.replace('graphene_', 'g_'), func))[:44],
                 ('%s.%s' % (h.replace('graphene_', 'g_'), l))[:30], v))
    print()
    print('  Of Section 3\'s %d static misses: %d have the holder PRESENT, %d'
          % (len(miss_table), len(real), len(absent)))
    print('  have it ABSENT, %d unknown.  Section 3\'s total stands as written'
          % len(unknown))
    print('  and is an upper bound; this is the subset that can fire.')
    return real, absent, unknown


# ---------------------------------------------------------------------------
# Section 5: validations against exactly known values
# ---------------------------------------------------------------------------

PLANT_NAME = '_binding_sweep_planted_probe'


def validations(mods, bindings, by_owner, miss_table, probes, trans,
                class_sites):
    print()
    print('-' * 78)
    print('SECTION 5 -- VALIDATIONS AGAINST EXACTLY KNOWN VALUES')
    print('-' * 78)
    print('  Six, and three of them are fixed exactly by this repository\'s own')
    print('  history rather than by a tolerance chosen here.')
    print()

    # (1) The 2026-10-10 repair must still be covered -- now by construction.
    MGRS = (('graphene_quantum_capacitance_argument', 'on_consistent_argument'),
            ('graphene_exact_charge_rederivation', 'on_exact'))
    cd_misses = [r for r in miss_table
                 if r[5] == 'carrier_density' and (r[0], r[1]) in MGRS]
    cd_class = [c for c in class_sites
                if (c[0], c[1]) in MGRS and c[3] == 'carrier_density']
    check('(1) the 2026-10-10 repair of crx.carrier_density: 0 misses in both '
          'managers, and both now cover it class-completely',
          len(cd_misses) == 0 and len(cd_class) == 2,
          'misses %d (expected 0), class-complete sites %d (expected 2)'
          % (len(cd_misses), len(cd_class)))

    # (2) The 2026-10-13 hand-found miss, kept on the record after its
    #     manager stopped being a hand-listed site.  REWRITTEN 2026-10-14:
    #     this check used to assert that on_consistent_argument left exactly
    #     cd.quantum_capacitance stale.  It now asserts the same FACT in the
    #     form the repair left it in -- the manager rebinds the class, and
    #     that one member is excluded by name with a stated reason.  The
    #     original form is not deleted: it is the first clause below, which
    #     still requires zero SILENT misses.
    import graphene_binding_class as bc_mod
    qca_misses = [(r[4], r[5]) for r in miss_table
                  if (r[0], r[1]) == MGRS[0]]
    qca_class_qc = [c for c in class_sites
                    if (c[0], c[1]) == MGRS[0] and c[3] == 'quantum_capacitance']
    key = ('graphene_contact_doping_model', 'quantum_capacitance')
    reason = bc_mod.CONTACT_OPERATING_POINT.get(key, '')
    check('(2) 2026-10-13\'s hand-found miss is on the record as a REASONED '
          'exclusion, not as a silent miss',
          (qca_misses == []
           and len(qca_class_qc) == 1 and qca_class_qc[0][4] is True
           and list(bc_mod.CONTACT_OPERATING_POINT) == [key]
           and len(reason) > 80),
          'silent misses %r; class-complete C_q sites %d with exclusions=%r; '
          'CONTACT_OPERATING_POINT keys %r; reason %d chars'
          % (qca_misses, len(qca_class_qc),
             qca_class_qc[0][4] if qca_class_qc else None,
             list(bc_mod.CONTACT_OPERATING_POINT), len(reason)))

    # (3) A planted binding must be found exactly once, by both methods.
    path = os.path.join(REPO, PLANT_NAME + '.py')
    before_s = len([b for b in bindings if b[0] == PLANT_NAME])
    with open(path, 'w', encoding='utf-8') as f:
        f.write('from graphene_fet_model import quantum_capacitance\n'
                'ALIASED = quantum_capacitance\n')
    try:
        mods_all2 = repo_modules()
        b2, _d2, a2, _e2 = static_bindings(mods_all2)
        t2 = transitive_aliases(mods_all2, b2, a2)
        found_s = [b for b in b2 if b[0] == PLANT_NAME]
        found_t = [t for t in t2 if t[0] == PLANT_NAME]
        planted = importlib.import_module(PLANT_NAME)
        mods2 = dict(mods)
        mods2[PLANT_NAME] = planted
        rt2 = runtime_bindings(
            mods2, {'graphene_fet_model': ['quantum_capacitance']})
        found_r = [r for r in rt2 if r[0] == PLANT_NAME]
        ok = (before_s == 0 and len(found_s) == 1
              and len(found_t) == 1 and len(found_r) == 2)
        check('(3) a planted binding is found exactly once by the AST rule, '
              'once by the transitive rule, and twice by identity',
              ok,
              'static %d (expected 1), transitive %d (expected 1), '
              'identity %d (expected 2: the import and the rename)'
              % (len(found_s), len(found_t), len(found_r)))
    finally:
        os.unlink(path)
        sys.modules.pop(PLANT_NAME, None)

    # (4) Identity-patch neutrality: the instrument must not be the signal.
    owner_names = {}
    for (owner, oname) in by_owner:
        owner_names.setdefault(owner, []).append(oname)
    base = {k: p() for k, p in probes.items()}
    worst_ulp = 0.0
    nwrapped = 0
    with identity_patch_all(mods, bindings, owner_names) as nwrapped:
        for k, p in probes.items():
            now = p()
            for a, b in zip(base[k], now):
                if a != b:
                    worst_ulp = max(worst_ulp,
                                    abs(a - b) / np.spacing(max(abs(a), abs(b))))
    check('(4) identity-patch neutrality over %d wrapped sites: every probe '
          'bit-identical' % nwrapped,
          worst_ulp == 0.0,
          'worst departure %.3g ulp; required exactly 0' % worst_ulp)

    # (5) A null that must be exactly zero: perturbing by unity.
    worst_null = 0.0
    for holder in sorted(probes):
        for local in ('mu', 'n_puddle', 'W', 'quantum_capacitance',
                      'carrier_density', 'C_ox'):
            if not hasattr(mods[holder], local):
                continue
            u = usedness(mods, probes, holder, local, factor=1.0)
            if u is not None:
                worst_null = max(worst_null, u)
    check('(5) the usedness instrument\'s own null: a perturbation by a factor '
          'of exactly 1.0 moves nothing',
          worst_null == 0.0,
          'worst %.3g; required exactly 0' % worst_null)

    # (6) A POSITIVE control for the arrival counter, so that every LATENT
    #     verdict in Section 4b is a measurement and not a broken counter.
    cdm = mods.get('graphene_contact_doping_model')
    n_metals = len(cdm.METAL_WORK_FUNCTIONS)
    key = 'graphene_contact_doping_model.quantum_capacitance'
    shipped = cdm.quantum_capacitance
    _ARRIVALS[key] = 0
    try:
        cdm.quantum_capacitance = _counting(key, shipped)
        for w in cdm.METAL_WORK_FUNCTIONS.values():
            cdm.contact_edge_carrier_density(w)
    finally:
        cdm.quantum_capacitance = shipped
    got = _ARRIVALS[key]
    check('(6) positive control: the SAME counter that read 0 in Section 4b '
          'reads exactly %d when the holder is exercised' % n_metals,
          got == n_metals,
          'counted %d of an exactly known %d' % (got, n_metals))
    return worst_ulp


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(by_owner, trans, miss_table, used, reach_rows,
                fname='import_binding_sweep.png'):
    classes = {}
    for (owner, oname), hs in by_owner.items():
        classes.setdefault((owner, oname), set()).update(hs)
    for h, l, o, n in trans:
        classes.setdefault((o, n), set()).add((h, l))
    wide = sorted(((k, v) for k, v in classes.items() if len(v) >= 2),
                  key=lambda kv: -len(kv[1]))[:10]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))

    ax = axes[0]
    labs = ['%s.%s' % (o.replace('graphene_', '').replace('_model', ''), n)
            for (o, n), _v in wide]
    vals = [len(v) + 1 for _k, v in wide]          # +1 for the owner itself
    ax.barh(range(len(vals))[::-1], vals, color='#4C72B0')
    ax.set_yticks(range(len(vals))[::-1])
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_xlabel('sites a single setattr() FAILS to reach, +1 for the owner')
    ax.set_title('(a) the rebind set is not one site\n%d names need more than one'
                 % len([1 for _k, v in classes.items() if len(v) >= 1]))
    ax.grid(alpha=0.3, axis='x')

    ax = axes[1]
    rows = [(k, v) for k, v in used.items() if v is not None]
    rows.sort(key=lambda kv: -(kv[1] if kv[1] is not None else -1))
    labs2 = ['%s.%s' % (h.replace('graphene_', '').replace('_model', ''), l)
             for (h, l), _v in rows]
    vals2 = [v for _k, v in rows]
    cols = ['#C44E52' if v > 0 else '#999999' for v in vals2]
    ax.barh(range(len(vals2))[::-1], [max(v, 1e-4) for v in vals2], color=cols)
    ax.set_xscale('log')
    ax.set_yticks(range(len(vals2))[::-1])
    ax.set_yticklabels(labs2, fontsize=8)
    ax.set_xlabel('max relative change in the holder\'s own numbers (log)')
    ax.set_title('(b) which misses are USED\nred = a defect waiting; grey = cosmetic')
    ax.grid(alpha=0.3, axis='x')

    ax = axes[2]
    per_site = {}
    for holder, func, owner, oname, h, l in miss_table:
        per_site['%s.%s' % (holder.replace('graphene_', ''), func)] = \
            per_site.get('%s.%s' % (holder.replace('graphene_', ''), func), 0) + 1
    items = sorted(per_site.items(), key=lambda kv: -kv[1])
    ax.barh(range(len(items))[::-1], [v for _k, v in items], color='#DD8452')
    ax.set_yticks(range(len(items))[::-1])
    ax.set_yticklabels([k for k, _v in items], fontsize=7)
    ax.set_xlabel('stale bindings left behind')
    n_live = len([r for r in reach_rows if r[2]])
    ax.set_title('(c) %d misses at %d sites\n%d measured LIVE, %d LATENT, rest '
                 'unenterable'
                 % (len(miss_table), len(items), n_live,
                    len([r for r in reach_rows if r[2] == 0])))
    ax.grid(alpha=0.3, axis='x')

    fig.suptitle('Import-binding sweep: a repair to one binding is not a repair '
                 'to the class', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fname, dpi=150)
    print('  figure written: %s' % fname)


def main():
    (mods, bindings, by_owner, miss_table, probes, used, trans,
     reach_rows, class_sites) = sections_1_to_4()
    real, absent, unknown = section_4c(miss_table)
    validations(mods, bindings, by_owner, miss_table, probes, trans,
                class_sites)
    print()
    print('-' * 78)
    print('SECTION 6 -- WHAT THIS MEANS FOR ITEM 16')
    print('-' * 78)
    live_used = [(h, l) for (h, l), v in used.items() if v]
    print('  Three axes, and a miss has to clear all three to be a defect:')
    print('    PRESENT (4c)  the holder is imported where the patch runs: %d of %d'
          % (len(real), len(miss_table)))
    print('    USED    (4a)  the holder\'s own numbers depend on its copy')
    print('    REACHED (4b)  the patched run actually calls it')
    print('  Section 3\'s %d is the static count and is kept as the upper'
          ' bound.' % len(miss_table))
    print()
    print('  Item 16 flips the default to form B_T.  Whatever installs it must')
    print('  reach every member of the identity class of BOTH charge-path names,')
    print('  not the two that history happens to have found by hand.  Measured')
    print('  here: %d of the %d distinct stale bindings are USED by their own'
          % (len(live_used), len(used)))
    print('  holder, so each is a number that would silently not move.')
    print()
    print('  And one of them cannot be repaired by rebinding at all.')
    print('  cd.quantum_capacitance is reached at a contact whose capacitance is')
    print('  C_interface, while the replacement the owner holds solves for n with')
    print('  C_ox.  Closing that miss with the owner\'s object would install a')
    print('  channel solver at the contact -- a DIFFERENT wrong answer, not a')
    print('  repair.  2026-10-13 measured the contact-correct form separately')
    print('  (-12.4964 % in n for Pd); it is not the object this sweep would')
    print('  copy.  So the sweep\'s output is two rules, not one:')
    print()
    print('    (i)  every member of the identity class must be rebound; and')
    print('    (ii) a holder that uses the name at a DIFFERENT operating point')
    print('         needs its OWN replacement, not the owner\'s.')
    print()
    make_figure(by_owner, trans, miss_table, used, reach_rows)
    print()
    print('=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    buf = io.StringIO()

    class Tee:
        def __init__(self, *s):
            self.s = s

        def write(self, d):
            for x in self.s:
                x.write(d)

        def flush(self):
            for x in self.s:
                x.flush()

    old = sys.stdout
    sys.stdout = Tee(old, buf)
    try:
        rc = main()
    finally:
        sys.stdout = old
    with open('import_binding_sweep_output.txt', 'w') as f:
        f.write(buf.getvalue())
    sys.exit(rc)
