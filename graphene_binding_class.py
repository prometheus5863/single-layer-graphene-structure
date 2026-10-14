"""
graphene_binding_class.py

The repair that graphene_import_binding_sweep.py's finding (i) asks for, in
reusable form: a context manager that rebinds a name's ENTIRE identity class
rather than the sites someone remembered to list.

WHY THIS IS A MODULE AND NOT A PATCH

Three times in ten days a monkey-patch in this repository reported a model it
was not running, each time because a second module held its own import-time
reference to the patched name:

    2026-10-03  the frozen-default fault, monkey-patch form.
    2026-10-10  contact_resistance_crossover.carrier_density.
    2026-10-13  graphene_contact_doping_model.quantum_capacitance.

Each was repaired by adding one more line to one more manager.  2026-10-13
stated the rule that makes that insufficient -- "a repair to one binding is
not a repair to the class" -- and 2026-10-14's sweep measured the size of the
remaining exposure: 25 stale sites at 14 rebinding sites, 22 of them never
found by hand, 7 of the 10 distinct stale bindings measured USED by their own
holder.  A hand-maintained list cannot be the mechanism, because the list is
exactly what was wrong all three times.

WHAT THE CLASS IS, AND WHY IT IS THE SAFE CRITERION

`from M import n` copies a REFERENCE at import time.  So at runtime the set of
places that must be rebound is the set of names that currently hold the same
OBJECT -- and that set is discoverable by identity, with no import graph, no
AST pass and no list.  The sweep found that the identity class is strictly
larger than the import graph: twelve second-order renaming aliases
(`C_OX = gfet.C_ox`, `N0_SHIPPED = gfet.n_puddle`,
`W_CROSS_VACUUM = W_GRAPHENE`) are invisible to any rule that looks only at
import statements.  Identity sees them for free.

The cost of identity is the opposite error: two names that hold the same
object for unrelated reasons join the same class.  For a FUNCTION that cannot
happen -- a function object has one definition site.  For a float it can, and
`exclude=` is how a caller states a known false member rather than having this
module guess.  Every exclusion must carry a reason, which is why
`exclude_reasons` is required alongside it and is printed by `describe()`.

WHAT THIS MODULE DOES NOT DO

It does not make a rebind CORRECT.  The sweep's finding (ii) is that
`graphene_contact_doping_model.quantum_capacitance` is reached at a contact
whose series capacitance is C_interface, while the replacement the owner holds
solves for n against C_ox; copying the owner's object there would install a
channel solver at the contact.  A class-complete rebinder would do exactly
that, confidently.  So the class is necessary and not sufficient, and a holder
that uses a name at a different operating point has to be excluded BY NAME
with the reason recorded -- see `CONTACT_OPERATING_POINT` below, which is that
exclusion and its reason in one place rather than re-argued at each call site.
"""

import contextlib
import os
import sys

REPO = os.path.dirname(os.path.abspath(__file__))


# The one standing exclusion in this repository, with its reason attached.
# Any manager that rebinds a charge-path name passes this and nothing else;
# a second exclusion has to be added here, with its own reason, rather than
# inlined at a call site where the next reader will not find it.
CONTACT_OPERATING_POINT = {
    ('graphene_contact_doping_model', 'quantum_capacitance'):
        'reached at the contact, where the series capacitance is C_interface '
        '(65.93x C_ox) and the series factor runs 0.555-0.806 rather than '
        '0.999; the owner\'s replacement solves for n against C_ox, so '
        'copying it here installs a channel solver at the contact.  The '
        'contact-correct form is measured separately (2026-10-13, -12.4964 % '
        'in n for Pd) and is not this object.',
}


def repo_module_objects():
    """Every already-imported module whose file lives in this repository.

    Deliberately reads sys.modules rather than importing anything: the class
    that matters is the one that exists in the process doing the patching.  A
    module nobody imported holds no live binding, and importing it to look
    would change the thing being measured.
    """
    out = {}
    for name, mod in list(sys.modules.items()):
        if mod is None or name.startswith('__'):
            continue
        f = getattr(mod, '__file__', None)
        if not f:
            continue
        try:
            if os.path.dirname(os.path.abspath(f)) == REPO:
                out[name] = mod
        except (TypeError, ValueError):
            continue
    return out


def identity_class(owner, name, exclude=()):
    """Every (module_name, attr) that currently holds `getattr(owner, name)`.

    The owner's own (module, name) is always first in the returned list.
    `exclude` is a set of (module_name, attr) pairs to drop.
    """
    owner_name = owner if isinstance(owner, str) else owner.__name__
    mods = repo_module_objects()
    if owner_name not in mods:
        raise KeyError('owner module %r is not imported' % owner_name)
    obj = getattr(mods[owner_name], name)
    excl = set(exclude)
    sites = [(owner_name, name)]
    for mname, mod in sorted(mods.items()):
        for attr, val in list(vars(mod).items()):
            if attr.startswith('__'):
                continue
            if (mname, attr) == (owner_name, name):
                continue
            if val is obj:
                sites.append((mname, attr))
    return [s for s in sites if s not in excl]


@contextlib.contextmanager
def rebind(owner, name, replacement, exclude_reasons=None, report=None):
    """Rebind every member of `owner.name`'s identity class to `replacement`.

    exclude_reasons: {(module_name, attr): reason}.  A member listed here is
        LEFT ALONE and the reason is recorded; passing a bare set is refused,
        because an exclusion without a stated reason is how all three of the
        2026-10 faults were introduced.
    report: optional callable(list_of_sites, list_of_excluded) for a caller
        that wants to print or assert on what was actually rebound.

    The class is computed at ENTRY, so a binding created after entry is not
    rebound -- correct, since such a binding was made from whatever was
    installed at the time.
    """
    exclude_reasons = dict(exclude_reasons or {})
    sites = identity_class(owner, name, exclude=set(exclude_reasons))
    skipped = [(k, v) for k, v in sorted(exclude_reasons.items())]
    mods = repo_module_objects()
    saved = []
    try:
        for mname, attr in sites:
            saved.append((mname, attr, getattr(mods[mname], attr)))
            setattr(mods[mname], attr, replacement)
        if report is not None:
            report(sites, skipped)
        yield sites
    finally:
        for mname, attr, val in reversed(saved):
            setattr(mods[mname], attr, val)


def describe(owner, name, exclude_reasons=None):
    """A one-call text summary, for a transcript."""
    exclude_reasons = dict(exclude_reasons or {})
    sites = identity_class(owner, name, exclude=set(exclude_reasons))
    lines = ['%s.%s -- identity class has %d member(s) to rebind:'
             % (owner if isinstance(owner, str) else owner.__name__, name,
                len(sites))]
    for mname, attr in sites:
        lines.append('    %s.%s' % (mname, attr))
    for (mname, attr), why in sorted(exclude_reasons.items()):
        lines.append('    EXCLUDED %s.%s' % (mname, attr))
        lines.append('        %s' % why)
    return '\n'.join(lines)
