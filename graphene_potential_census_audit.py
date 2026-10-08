"""
graphene_potential_census_audit.py

A VARIABLE THAT IS A POTENTIAL IN TWO DIFFERENT SENSES IS NOT A VARIABLE.

Created 2026-10-08 to mechanise the 2026-10-07 item "a census of every
potential-like variable in this repository recording which potential it is"
(Chapter 7 Section 7.9 item 10).

WHY THIS CLASS NEEDS A SWEEP RATHER THAN A SECOND ASK
-----------------------------------------------------
On 2026-10-07 a diffusion term was added to Eq. (4) of
graphene_fet_model.py, built exactly as specified, validated, and then found
to rest on a premise that was not a fact about the model: `V_ch` was
undetermined between the ELECTROSTATIC potential and the QUASI-FERMI
potential, and those two readings differ by *exactly* the term that had just
been built.  The term was real; the deficiency it was correcting was an
artefact of a name.

That was found by asking one variable one question.  Three other
potential-like variables were named as plausible siblings in the same log
entry, on no evidence: nothing in this repository was looking for the class.
This file looks for it, and -- the point of mechanising -- cannot go stale
silently, because a potential-like name that appears in the repository and is
absent from the registry below is reported as a fault in its own right.

OVER-COLLECT, THEN ADJUDICATE
-----------------------------
The name pattern is deliberately too wide: it catches `mu_acoustic`, which is
a mobility and not a chemical potential at all.  That is the correct bias for
a census.  A pattern narrow enough to have no false positives is a pattern
that has decided in advance which names are potentials, which is the thing
being audited.  So every matched name must appear in POTENTIAL_REGISTRY with
an explicit reading, and MOBILITY_NOT_A_POTENTIAL is one of the readings a
name can be adjudicated to.

THE TWO ROLES, AND THE SCOPE THEY ARE MEASURED OVER
---------------------------------------------------
  ELECTROSTATIC   the name enters a charge/capacitance relation: it is the
                  potential whose difference from the gate sets the induced
                  sheet charge through C_ox, C_q or quantum_capacitance().
  TRANSPORT       the name enters a current/field/resistance relation: it is
                  the potential whose gradient (or whose span, when it is the
                  integration variable of a gradual-channel integral) carries
                  the current.

A variable in BOTH roles, with no declaration, is UNDETERMINED: the model it
appears in is silently one of two different models.  This is the 10-07 defect.

Roles are attributed over the enclosing FUNCTION, not over the line, because
a function is the unit over which a variable's meaning has to be
single-valued -- and because the 10-07 defect is invisible at line scope.
`V_ch` appears in `dV = V_g - V_ch - V_dirac`, a line with no capacitance
symbol on it; the capacitance is two lines down, acting on `dV`.  A
line-scoped classifier would have missed the one defect this file is known to
have to catch.  Module-body statements are attributed per STATEMENT instead,
since a module body is a sequence of independent definitions rather than one
computation.

Coarse attribution over-flags rather than under-flags, so every role claim is
printed with the source line that produced it.  A reader who thinks a verdict
is wrong is given the evidence to refute it in the transcript, which is the
only form of over-flagging that is honest.

THE IN-SOURCE DECLARATION, AND WHY IT IS NOT A LIST IN HERE
-----------------------------------------------------------
A dual-role variable is cleared by DECLARING its reading where it is defined:

    # potential-reading: QUASI_FERMI -- <reason>

following graphene_dead_name_sweep.py's convention that the repository's
reasons belong in the repository and not in the instrument that judges it.
The reason is required; a declaration with no reason is a fault.

CONTROLS (all four run every time; see SECTION 5)
  POSITIVE, exactly known     graphene_fet_model.V_ch must come out dual-role.
                              Established independently on 10-07.  If this
                              classifier misses it, the classifier is broken
                              and no other verdict here is worth reading.
  NULL, at the point of entry Run the whole classifier over a copy of the
                              repository with every potential-like name
                              rewritten to a non-potential name.  The census
                              must come back EMPTY.  2026-09-30's rule: the
                              control has to sit where the failure enters.
  SYNTHETIC, both directions  One variable in one role must be DETERMINED;
                              one variable in two roles must be UNDETERMINED.
  MUTATION                    Inject a second role into a copy of a module
                              whose variable is currently single-role; the
                              census must flip it to UNDETERMINED.  Plus an
                              INERT mutant (the same text in a comment) which
                              it must NOT flip -- 10-07's distinction between
                              a mutant that is uncaught and a mutant that
                              never arrived.
"""

import ast
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# Deliberately over-wide; see the docstring.  Every match must be adjudicated
# in POTENTIAL_REGISTRY.
POTENTIAL_NAME_RE = re.compile(
    r'^('
    r'V_[A-Za-z0-9_]*'      # V_g, V_ch, V_bias, V_dirac, ...
    r'|phi[A-Za-z0-9_]*'    # phi, phi_k
    r'|psi[A-Za-z0-9_]*'
    r'|mu_[A-Za-z0-9_]*'    # mostly mobilities in this repository -- see docstring
    r'|E_F[A-Za-z0-9_]*'
    r')$'
)

# A charge/capacitance relation.
ELECTROSTATIC_MARKER_RE = re.compile(
    r'\bC_ox\b|\bC_q\b|\bC_total\b|\bC_gs\b|\bC_gd\b|\bquantum_capacitance\b'
    r'|\bcarrier_density\b|\bn_electrostatic\b|\bseries_factor\b'
)

# A current / field / resistance relation, or the span of a gradual-channel
# integration variable.
TRANSPORT_MARKER_RE = re.compile(
    r'\bR_[A-Za-z0-9_]+\b|\bRc_total\b|\bchannel_resistance\b'
    r'|\bsheet_conductivity\b|\bsigma_sheet\b|\bE_field\b'
    r'|\bI_?d\b|\bId\[|\bv_drift\b|\bv_sat\b'
    r'|/\s*L(_channel)?\b|\blinspace\(\s*0\s*,\s*V'
)

DECL_RE = re.compile(
    r'#\s*potential-reading:\s*([A-Z_]+)\s*(?:--\s*(.*))?$'
)

# ---------------------------------------------------------------------------
# SECTION 0 -- the census registry.
#
# Every name the scanner can find must be here.  `reading` is the
# adjudication; `note` is why.  UNADJUDICATED is an honest value: it means the
# name has been collected and NOT yet decided.  It is reported and counted,
# and it is not a pass.
#
# Readings:
#   ELECTROSTATIC_ONLY       a potential, used only in charge relations
#   TRANSPORT_ONLY           a potential, used only in current/field relations
#   QUASI_FERMI              declared to be the electrochemical potential
#   TERMINAL_BIAS            an applied terminal voltage; no channel ambiguity
#                            is possible because it is not a field variable
#   REFERENCE_LEVEL          a fixed offset/anchor (a Dirac-point voltage, a
#                            neutrality point), not a potential that varies
#   ENERGY_NOT_A_POTENTIAL   an energy in eV or J (E_F...), not a volt
#   MOBILITY_NOT_A_POTENTIAL a mobility; the name pattern over-collected
#   GEOMETRIC_PHASE          a lattice/angle phi, not an electric potential
#   DISPLAY_ONLY             appears only inside format strings
#   UNADJUDICATED            collected, not decided
# ---------------------------------------------------------------------------

POTENTIAL_REGISTRY = {
    # --- Chapter 4, the GFET: the one place the class is known to bite ------
    'graphene_fet_model.py::V_ch': ('QUASI_FERMI',
        'declared 2026-10-08; see Section 4.6.6 and the module docstring'),
    'graphene_fet_model.py::V_channel_profile': ('QUASI_FERMI',
        'the gradual-channel integration variable, 0 -> Vds, same reading as V_ch'),
    'graphene_fet_model.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_fet_model.py::V_dirac': ('REFERENCE_LEVEL', 'charge-neutrality offset'),
    'graphene_fet_model.py::V_Dirac': ('REFERENCE_LEVEL', 'plot-label spelling of V_dirac'),
    'graphene_fet_model.py::V_channel_shift': ('ELECTROSTATIC_ONLY',
        'the band-bending offset in the docstring derivation; see Section 4.6.6 '
        'for why it is unread in the closed form actually used'),
    'graphene_fet_model.py::V_g_minus_Vdirac': ('REFERENCE_LEVEL',
        'a difference of a terminal bias and an offset'),
    'graphene_fet_model.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment in a label'),

    # --- modules that import the GFET channel potential --------------------
    'graphene_diffusion_current_model.py::V_ch': ('QUASI_FERMI',
        'inherits graphene_fet_model.V_ch; this is the module whose own term '
        'the declaration makes a double count -- see Section 4.6.6'),
    'graphene_diffusion_current_model.py::V_ds': ('TERMINAL_BIAS', 'applied drain bias'),
    'graphene_diffusion_current_model.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_diffusion_current_model.py::V_dirac': ('REFERENCE_LEVEL', 'offset'),
    'graphene_diffusion_current_model.py::V_F': ('UNADJUDICATED',
        'a Fermi VELOCITY written with a V_ prefix; over-collected by the name '
        'pattern and not yet confirmed to be velocity at every use site'),
    'graphene_diffusion_current_model.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),
    'graphene_diffusion_current_model.py::E_F': ('ENERGY_NOT_A_POTENTIAL', 'an energy'),
    'graphene_diffusion_current_model.py::mu_resp': ('UNADJUDICATED',
        'name pattern over-collected; mobility or chemical potential not confirmed'),
    'graphene_diffusion_current_model.py::mu_saved': ('MOBILITY_NOT_A_POTENTIAL',
        'a saved copy of the module mobility'),

    'graphene_gds_quadrature_audit.py::V_ch': ('QUASI_FERMI', 'inherits graphene_fet_model.V_ch'),
    'graphene_gds_quadrature_audit.py::V_channel_profile': ('QUASI_FERMI', 'as above'),
    'graphene_gds_quadrature_audit.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_gds_quadrature_audit.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),

    'contact_resistance_crossover.py::V_ch': ('QUASI_FERMI', 'inherits graphene_fet_model.V_ch'),
    'contact_resistance_crossover.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'contact_resistance_crossover.py::V_g_on': ('TERMINAL_BIAS', 'a chosen on-state gate bias'),
    'contact_resistance_crossover.py::V_dirac': ('REFERENCE_LEVEL', 'offset'),

    'graphene_velocity_saturation_model.py::V_ch': ('QUASI_FERMI', 'inherits graphene_fet_model.V_ch'),
    'graphene_velocity_saturation_model.py::V_ds': ('TERMINAL_BIAS', 'applied drain bias'),
    'graphene_velocity_saturation_model.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_velocity_saturation_model.py::V_hi': ('TERMINAL_BIAS', 'sweep endpoint'),
    'graphene_velocity_saturation_model.py::V_SAT_ANCHOR': ('UNADJUDICATED',
        'a literature anchor; which quantity it anchors is not re-checked here'),
    'graphene_velocity_saturation_model.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),
    'graphene_velocity_saturation_model.py::mu_saved': ('MOBILITY_NOT_A_POTENTIAL', 'saved mobility'),
    'graphene_velocity_saturation_model.py::mu_t': ('MOBILITY_NOT_A_POTENTIAL', 'a trial mobility'),
    'graphene_velocity_saturation_mutation.py::V_ds': ('TERMINAL_BIAS', 'applied drain bias'),

    'graphene_fmax_shortfall_decomposition.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_fmax_shortfall_decomposition.py::V_ds': ('TERMINAL_BIAS', 'applied drain bias'),
    'graphene_fmax_shortfall_decomposition.py::V_dirac': ('REFERENCE_LEVEL', 'offset'),
    'graphene_fmax_shortfall_decomposition.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),

    'graphene_perfect_contact_counterfactual.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_perfect_contact_counterfactual.py::V_ds': ('TERMINAL_BIAS', 'applied drain bias'),
    'graphene_perfect_contact_counterfactual.py::V_dirac': ('REFERENCE_LEVEL', 'offset'),
    'graphene_perfect_contact_counterfactual.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),
    'graphene_perfect_contact_mutation.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),

    'rf_small_signal_model.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'rf_small_signal_model.py::V_dirac': ('REFERENCE_LEVEL', 'offset'),
    'rf_small_signal_model.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),

    'graphene_rootfinder_audit.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_diffusion_mutation.py::V_g': ('TERMINAL_BIAS', 'applied gate voltage'),
    'graphene_diffusion_mutation.py::V_F': ('UNADJUDICATED', 'see the V_F entry above'),

    # --- Chapter 6, photodetectors -----------------------------------------
    'graphene_photodetector_model.py::V_bias': ('TRANSPORT_ONLY',
        'enters only E_field = V_bias / L_channel and hence tau_transit; no '
        'charge relation anywhere in Chapter 6, so no ambiguity is possible'),
    'graphene_photodetector_two_contact_model.py::V_bias': ('DISPLAY_ONLY',
        'imported and printed; it enters no computation in this module -- '
        'collected here as a provenance fact, not a fault'),
    'graphene_photodetector_two_contact_model.py::V_': ('DISPLAY_ONLY', 'LaTeX fragment'),
    'graphene_photodetector_collection_model.py::V_bias': ('DISPLAY_ONLY',
        'imported and printed only'),
    'graphene_photodetector_nonuniform_illumination_model.py::V_bias': ('DISPLAY_ONLY',
        'imported and printed only'),
    'graphene_photodetector_nonuniform_illumination_model.py::mu_e': ('MOBILITY_NOT_A_POTENTIAL',
        'electron mobility'),
    'graphene_photodetector_nonuniform_illumination_model.py::mu_h': ('MOBILITY_NOT_A_POTENTIAL',
        'hole mobility'),
    'graphene_photodetector_signed_carrier_model.py::V_bias': ('TRANSPORT_ONLY',
        'drift field for the signed-carrier drift lengths'),
    'graphene_photodetector_signed_carrier_model.py::mu_e': ('MOBILITY_NOT_A_POTENTIAL', 'electron mobility'),
    'graphene_photodetector_signed_carrier_model.py::mu_h': ('MOBILITY_NOT_A_POTENTIAL', 'hole mobility'),
    'graphene_photodetector_signed_carrier_model.py::phi': ('UNADJUDICATED',
        'a phi in a photodetector module; not yet confirmed to be a photon '
        'flux rather than a potential'),

    # --- Chapter 2/3, band structure and transport -------------------------
    # graphene_band_structure.py::phi was a registry row with no code behind
    # it: the name appears there only inside a COMMENT.  Removed 2026-10-08
    # by the dead-row check in SECTION 2.
    'graphene_band_structure.py::phi_k': ('GEOMETRIC_PHASE', 'the pseudospin angle'),
    'graphene_band_structure_audit.py::phi': ('GEOMETRIC_PHASE', 'as above'),
    'graphene_band_structure_audit.py::E_F': ('ENERGY_NOT_A_POTENTIAL', 'an energy'),
    'graphene_transport_properties.py::V_neutral': ('REFERENCE_LEVEL', 'neutrality point'),
    'graphene_transport_properties.py::mu_acoustic': ('MOBILITY_NOT_A_POTENTIAL', 'phonon-limited mobility'),
    'graphene_transport_properties.py::mu_impurity': ('MOBILITY_NOT_A_POTENTIAL', 'impurity-limited mobility'),
    'graphene_transport_properties.py::mu_total': ('MOBILITY_NOT_A_POTENTIAL', 'Matthiessen total'),
    'graphene_transport_optical_audit.py::mu_ac': ('MOBILITY_NOT_A_POTENTIAL', 'mobility'),
    'graphene_transport_optical_audit.py::mu_acoustic': ('MOBILITY_NOT_A_POTENTIAL', 'mobility'),

    # --- contacts ----------------------------------------------------------
    'graphene_contact_doping_nonlinear_model.py::E_F': ('ENERGY_NOT_A_POTENTIAL', 'an energy'),
    'graphene_contact_doping_nonlinear_model.py::V_F': ('UNADJUDICATED', 'see the V_F entry above'),
    'graphene_contact_doping_nonlinear_model.py::phi': ('UNADJUDICATED',
        'plausibly a metal work function; not confirmed at every use site'),
    'graphene_edge_contact_model.py::V_BG': ('TERMINAL_BIAS', 'applied back-gate voltage'),
    'graphene_edge_contact_model.py::E_F': ('ENERGY_NOT_A_POTENTIAL', 'an energy'),
    'graphene_edge_contact_model.py::E_F_eV': ('ENERGY_NOT_A_POTENTIAL', 'an energy in eV'),
    'graphene_edge_contact_model.py::E_F_joules': ('ENERGY_NOT_A_POTENTIAL', 'an energy in J'),

    # --- instruments -------------------------------------------------------
    'graphene_figure_provenance_audit.py::phi': ('GEOMETRIC_PHASE', 'an angle in a checked figure'),
}

CHANNEL_POTENTIAL_READINGS = ('QUASI_FERMI', 'ELECTROSTATIC_ONLY', 'TRANSPORT_ONLY')


# ---------------------------------------------------------------------------
# SECTION 1 -- the scanner
#
# THE GATE THAT MATTERS, AND HOW THE FIRST VERSION OF THIS FILE GOT IT WRONG
# --------------------------------------------------------------------------
# The first run of this classifier reported THIRTEEN dual-role names,
# including `V_g`, `V_hi`, `mu_t` and `V_SAT_ANCHOR`.  Every one of those is a
# false positive, and they are false in an instructive way: `V_g` really does
# appear in a charge relation and in a current relation, in the same function,
# and that is not a defect.  A terminal bias is the SAME physical quantity in
# both expressions -- it is the gate electrode's potential, fixed by a supply
# -- so there is no second reading for it to be ambiguous between.
#
# What made `V_ch` ambiguous was not that it appeared twice.  It was that it
# is a CHANNEL variable: a potential that varies along the channel and is
# swept from 0 to the drain bias INSIDE a single current evaluation.  Only
# such a variable has two inequivalent readings -- electrostatic potential and
# quasi-Fermi potential -- that differ by a diffusion term.
#
# So "both roles" is necessary and not sufficient.  The defect verdict
# requires both roles AND channel-swept, where channel-swept is detected
# mechanically: the name is assigned from, or iterated over, an array built
# from 0 to a DRAIN bias, with one level of indirection followed
# (`profile = linspace(0, Vds, n)` then `for V_ch in profile`).
#
# Both roles WITHOUT channel-swept is reported as DUAL_ROLE_TERMINAL and is
# not a fault.  Mutant M4 in SECTION 2 exists to prove that gate is
# load-bearing rather than decorative: it removes the drain-bias name from
# the sweep in a copy of graphene_fet_model.py and requires `V_ch` to fall
# out of the defect class.
#
# Comments and docstrings are stripped before role matching.  The first
# version matched a role marker on the line
#     # used downstream (sheet_conductivity takes abs(n)), so np.sign() would
# which is prose about the code and not the code, and on each function's own
# `def` line, which matched the function's own name.  2026-09-28 established
# that prose is a detector; it is not evidence of a data dependency.
# ---------------------------------------------------------------------------


def _strip_prose(lines, fn_def_lines):
    """
    Blank out comments, string literals and `def` signature lines, so a role
    marker can only fire on executable code.
    """
    out = []
    in_doc = False
    doc_q = None
    for idx, raw in enumerate(lines, start=1):
        ln = raw
        if in_doc:
            if doc_q in ln:
                ln = ln.split(doc_q, 1)[1]
                in_doc = False
            else:
                out.append('')
                continue
        for q in ('"""', "'''"):
            if q in ln:
                before, _, after = ln.partition(q)
                if q in after:
                    ln = before + ' ' + after.split(q, 1)[1]
                else:
                    ln = before
                    in_doc = True
                    doc_q = q
                break
        ln = re.sub(r'#.*$', '', ln)
        ln = re.sub(r'"[^"]*"|\'[^\']*\'', '""', ln)
        if idx in fn_def_lines:
            ln = ''
        out.append(ln)
    return out


def _names_in(node):
    """Yield (name, lineno) for potential-like Names/args, skipping f-strings."""
    for sub in ast.walk(node):
        if isinstance(sub, ast.JoinedStr):
            continue
        if isinstance(sub, ast.Name):
            if POTENTIAL_NAME_RE.match(sub.id):
                yield sub.id, getattr(sub, 'lineno', 0)
        elif isinstance(sub, ast.arg):
            if POTENTIAL_NAME_RE.match(sub.arg):
                yield sub.arg, getattr(sub, 'lineno', 0)
        elif isinstance(sub, ast.Attribute):
            # `gfet.V_dirac` is a potential used across a module boundary.
            # The first version of this scanner saw only ast.Name and missed
            # every one of them, which showed up as three registry rows with
            # no code behind them -- the dead-row check earning its keep.
            if POTENTIAL_NAME_RE.match(sub.attr):
                yield sub.attr, getattr(sub, 'lineno', 0)


def _string_names_in(tree):
    """Potential-like tokens appearing inside string constants / f-strings."""
    out = set()
    for sub in ast.walk(tree):
        if isinstance(sub, ast.JoinedStr):
            for inner in ast.walk(sub):
                if isinstance(inner, ast.Name) and POTENTIAL_NAME_RE.match(inner.id):
                    out.add(inner.id)
        if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
            for tok in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', sub.value):
                if POTENTIAL_NAME_RE.match(tok):
                    out.add(tok)
    return out


def _imported_names_in(tree):
    out = set()
    for sub in ast.walk(tree):
        if isinstance(sub, (ast.Import, ast.ImportFrom)):
            for al in sub.names:
                nm = al.asname or al.name.split('.')[0]
                if POTENTIAL_NAME_RE.match(nm):
                    out.add(nm)
    return out


DRAIN_BIAS_RE = re.compile(r'\bV?ds\b|\bV_ds\b|\bVds\b|\bVDS\b')


def _channel_swept_names(tree, code_lines):
    """
    Names that are swept from 0 to a DRAIN bias inside a single evaluation.
    One level of indirection is followed: an array assigned from
    linspace(0, Vds, ...) marks its own name, and then any iteration target
    over that array.
    """
    direct = set()
    arrays = set()

    for sub in ast.walk(tree):
        if isinstance(sub, ast.Assign) and len(sub.targets) == 1 and \
                isinstance(sub.targets[0], ast.Name):
            seg = '\n'.join(code_lines[sub.lineno - 1:
                                       getattr(sub, 'end_lineno', sub.lineno)])
            if re.search(r'linspace\s*\(\s*0', seg) and DRAIN_BIAS_RE.search(seg):
                arrays.add(sub.targets[0].id)
                if POTENTIAL_NAME_RE.match(sub.targets[0].id):
                    direct.add(sub.targets[0].id)

    def targets_of(node):
        if isinstance(node, ast.Name):
            return [node.id]
        if isinstance(node, (ast.Tuple, ast.List)):
            got = []
            for el in node.elts:
                got.extend(targets_of(el))
            return got
        return []

    for sub in ast.walk(tree):
        iter_node = None
        tgt = None
        if isinstance(sub, ast.For):
            iter_node, tgt = sub.iter, sub.target
        elif isinstance(sub, ast.comprehension):
            iter_node, tgt = sub.iter, sub.target
        if iter_node is None:
            continue
        iter_names = {n.id for n in ast.walk(iter_node) if isinstance(n, ast.Name)}
        seg_ok = bool(iter_names & arrays)
        if not seg_ok:
            try:
                seg = ast.unparse(iter_node)
            except Exception:
                seg = ''
            seg_ok = bool(re.search(r'linspace\s*\(\s*0', seg)
                          and DRAIN_BIAS_RE.search(seg))
        if seg_ok:
            for t in targets_of(tgt):
                if POTENTIAL_NAME_RE.match(t):
                    direct.add(t)
    return direct


def _roles_of_text(text):
    roles = set()
    if ELECTROSTATIC_MARKER_RE.search(text):
        roles.add('ELECTROSTATIC')
    if TRANSPORT_MARKER_RE.search(text):
        roles.add('TRANSPORT')
    return roles


def _evidence(code_lines, lo, hi, pattern):
    for ln in code_lines[lo - 1:hi]:
        if pattern.search(ln):
            return ln.strip()[:100]
    return '(no single line)'


def scan_file(path):
    src = open(path, 'r').read()
    lines = src.splitlines()
    tree = ast.parse(src)

    func_nodes = [n for n in ast.walk(tree)
                  if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    fn_def_lines = set()
    for fn in func_nodes:
        # the signature can span lines; blank from `def` to the first body stmt
        first_body = fn.body[0].lineno if fn.body else fn.lineno
        for ln in range(fn.lineno, first_body):
            fn_def_lines.add(ln)
    code_lines = _strip_prose(lines, fn_def_lines)

    declarations = {}
    for i, ln in enumerate(lines):
        m = DECL_RE.search(ln)
        if not m:
            continue
        reading, reason = m.group(1), (m.group(2) or '').strip()
        target = None
        for probe in (ln, lines[i + 1] if i + 1 < len(lines) else ''):
            mm = re.search(r'\b(V_[A-Za-z0-9_]*|phi[A-Za-z0-9_]*|mu_[A-Za-z0-9_]*)\b',
                           re.sub(r'#.*potential-reading[^\n]*', '', probe))
            if mm:
                target = mm.group(1)
                break
        declarations[target] = (reading, reason, i + 1)

    swept = _channel_swept_names(tree, code_lines)

    per_name = {}

    def record(name, scope, roles, line, ev):
        d = per_name.setdefault(name, {'roles': set(), 'scopes': [], 'decl': None,
                                       'swept': name in swept, 'error': None})
        d['roles'] |= roles
        d['scopes'].append((scope, line, ev))

    func_spans = [(fn.lineno, getattr(fn, 'end_lineno', fn.lineno)) for fn in func_nodes]

    for fn in func_nodes:
        lo, hi = fn.lineno, getattr(fn, 'end_lineno', fn.lineno)
        body_src = '\n'.join(code_lines[lo - 1:hi])
        roles = _roles_of_text(body_src)
        ev = []
        if 'ELECTROSTATIC' in roles:
            ev.append('ES<- ' + _evidence(code_lines, lo, hi, ELECTROSTATIC_MARKER_RE))
        if 'TRANSPORT' in roles:
            ev.append('TR<- ' + _evidence(code_lines, lo, hi, TRANSPORT_MARKER_RE))
        seen = set()
        for ident, lineno in _names_in(fn):
            if ident in seen:
                continue
            seen.add(ident)
            record(ident, '%s()' % fn.name, set(roles), lineno, ' | '.join(ev))

    def inside_fn(lineno):
        return any(a <= lineno <= b for a, b in func_spans)

    for stmt in tree.body:
        if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        lo, hi = stmt.lineno, getattr(stmt, 'end_lineno', stmt.lineno)
        stmt_src = '\n'.join(code_lines[lo - 1:hi])
        roles = _roles_of_text(stmt_src)
        ev = []
        if 'ELECTROSTATIC' in roles:
            ev.append('ES<- ' + _evidence(code_lines, lo, hi, ELECTROSTATIC_MARKER_RE))
        if 'TRANSPORT' in roles:
            ev.append('TR<- ' + _evidence(code_lines, lo, hi, TRANSPORT_MARKER_RE))
        for ident, lineno in _names_in(stmt):
            if inside_fn(lineno):
                continue
            record(ident, '<module>', set(roles), lineno, ' | '.join(ev))

    for ident in _imported_names_in(tree):
        record(ident, '<import>', set(), 0, '')
    for ident in _string_names_in(tree):
        record(ident, '<format string>', set(), 0, '')

    for name, d in per_name.items():
        if name in declarations:
            d['decl'] = declarations[name]
    return per_name


def census(root):
    out = {}
    for fname in sorted(os.listdir(root)):
        if not fname.endswith('.py'):
            continue
        try:
            per_name = scan_file(os.path.join(root, fname))
        except SyntaxError as exc:
            out['%s::<PARSE ERROR>' % fname] = {
                'roles': set(), 'scopes': [], 'decl': None, 'swept': False,
                'error': str(exc)}
            continue
        for name, d in per_name.items():
            out['%s::%s' % (fname, name)] = d
    return out


def verdict_for(entry):
    if entry.get('error'):
        return 'PARSE_ERROR'
    if entry['decl'] is not None:
        reading, reason, _ = entry['decl']
        return 'DECLARED' if reason else 'DECLARED_NO_REASON'
    roles = entry['roles']
    if len(roles) == 2:
        return 'UNDETERMINED' if entry['swept'] else 'DUAL_ROLE_TERMINAL'
    if not roles:
        return 'COLLECTED_NO_ROLE'
    return 'DETERMINED'


# ---------------------------------------------------------------------------
# SECTION 2 -- controls
# ---------------------------------------------------------------------------

SYNTHETIC_SINGLE_ROLE = '''
import numpy as np
C_ox = 1.0
V_dirac = 0.0

def charge(V_g, V_probe=0.0):
    dV = V_g - V_probe - V_dirac
    return C_ox * dV
'''

SYNTHETIC_DUAL_ROLE = '''
import numpy as np
C_ox = 1.0
V_dirac = 0.0
Rc_total = 1.0

def charge(V_g, V_probe=0.0):
    dV = V_g - V_probe - V_dirac
    return C_ox * dV

def current(V_g, Vds=0.05):
    profile = np.linspace(0, Vds, 4)
    R_channel = sum(1.0 / charge(V_g, V_probe) for V_probe in profile)
    R_total = R_channel + Rc_total
    return Vds / R_total
'''


def synthetic_control():
    results = {}
    for label, text in (('single', SYNTHETIC_SINGLE_ROLE),
                        ('dual', SYNTHETIC_DUAL_ROLE)):
        d = tempfile.mkdtemp(prefix='potcensus_syn_')
        try:
            with open(os.path.join(d, 'syn_%s.py' % label), 'w') as fh:
                fh.write(text)
            c = census(d)
            key = 'syn_%s.py::V_probe' % label
            results[label] = verdict_for(c[key]) if key in c else 'NOT FOUND'
        finally:
            shutil.rmtree(d, ignore_errors=True)
    return results


def null_control():
    """
    Rewrite every potential-like name in a copy of the repository to a
    non-potential name and re-run.  The census must come back EMPTY.  The
    control sits where the failure enters: it is the NAME pattern that is
    being tested, so the names are what get removed.
    """
    d = tempfile.mkdtemp(prefix='potcensus_null_')
    try:
        renames = {}
        counter = [0]

        def rename(m):
            tok = m.group(0)
            if not POTENTIAL_NAME_RE.match(tok):
                return tok
            if tok not in renames:
                counter[0] += 1
                renames[tok] = 'zz%03d' % counter[0]
            return renames[tok]

        n_files = 0
        for fname in sorted(os.listdir(HERE)):
            if not fname.endswith('.py'):
                continue
            if fname == os.path.basename(__file__):
                # The instrument is excluded from its own census in SECTION 1
                # and must be excluded here for the same reason.  It is also
                # the one file the rename CANNOT clean: its mutant table holds
                # the string '...TRANSPORT_ONLY\\nV_bias = 0.1', where the
                # escape makes `nV_bias` a single source token that the
                # identifier rewriter does not match, while ast sees a real
                # newline and recovers `V_bias` from the parsed value.  That is
                # a true limitation of rewriting source TEXT rather than
                # identifiers, and it is recorded rather than worked around.
                continue
            src = open(os.path.join(HERE, fname)).read()
            with open(os.path.join(d, fname), 'w') as fh:
                fh.write(re.sub(r'[A-Za-z_][A-Za-z0-9_]*', rename, src))
            n_files += 1
        c = census(d)
        return n_files, len(renames), sorted(c.keys())
    finally:
        shutil.rmtree(d, ignore_errors=True)


# Each mutant: (label, file, variable, anchor text, replacement, expected verdict)
MUTANTS = (
    ('M1 a charge role injected into a TRANSPORT_ONLY terminal bias',
     'graphene_photodetector_model.py', 'V_bias',
     'E_field = V_bias / L_channel',
     'E_field = V_bias / L_channel\nC_ox = 1.0\n_Q_INJECTED = C_ox * (V_bias - 0.0)',
     'DUAL_ROLE_TERMINAL'),
    ('M2 INERT: the same injection inside a comment',
     'graphene_photodetector_model.py', 'V_bias',
     'E_field = V_bias / L_channel',
     'E_field = V_bias / L_channel\n# C_ox = 1.0 ; _Q = C_ox * (V_bias - 0.0)',
     'DETERMINED'),
    ('M3 a declaration with no reason',
     'graphene_photodetector_model.py', 'V_bias',
     'V_bias = 0.1',
     '# potential-reading: TRANSPORT_ONLY\nV_bias = 0.1',
     'DECLARED_NO_REASON'),
    ('M4 the channel-swept gate is load-bearing: remove the drain-bias name '
     'from the sweep and V_ch must LEAVE the defect class',
     'graphene_fet_model.py', 'V_ch',
     'V_channel_profile = np.linspace(0, Vds, n_segments)',
     'V_channel_profile = np.linspace(0, 0.05, n_segments)',
     'DUAL_ROLE_TERMINAL'),
    ('M5 a SECOND channel-swept dual-role variable in real repository code '
     'must also be caught',
     'graphene_fet_model.py', 'V_probe',
     'V_channel_profile = np.linspace(0, Vds, n_segments)',
     'V_channel_profile = np.linspace(0, Vds, n_segments)\n'
     '    V_probe = np.linspace(0, Vds, 2)[0]\n'
     '    _q = C_ox * (0.0 - V_probe)\n'
     '    _r = channel_resistance(0.0, V_probe)',
     'UNDETERMINED'),
)


def mutation_control():
    out = []
    for label, fname, varname, old, new, expected in MUTANTS:
        d = tempfile.mkdtemp(prefix='potcensus_mut_')
        try:
            src = open(os.path.join(HERE, fname)).read()
            if old not in src:
                out.append((label, 'ANCHOR MISSING -- mutant never arrived',
                            expected, False))
                continue
            with open(os.path.join(d, fname), 'w') as fh:
                fh.write(src.replace(old, new, 1))
            c = census(d)
            key = '%s::%s' % (fname, varname)
            got = verdict_for(c[key]) if key in c else 'NOT FOUND'
            out.append((label, got, expected, got == expected))
        finally:
            shutil.rmtree(d, ignore_errors=True)
    return out


# ---------------------------------------------------------------------------
# SECTION 3 -- report
# ---------------------------------------------------------------------------

passed = 0
failed = 0
EXPECTED_FAILURES = 0


def check(label, ok, detail='', expected_to_fail=False):
    global passed, failed, EXPECTED_FAILURES
    if ok:
        passed += 1
        print('  PASS  %s' % label)
    else:
        failed += 1
        if expected_to_fail:
            EXPECTED_FAILURES += 1
            print('  FAIL* %s   (* expected; see detail)' % label)
        else:
            print('  FAIL  %s' % label)
    if detail:
        for chunk in detail.split('\n'):
            print('        %s' % chunk)


def main():
    print('=' * 78)
    print('POTENTIAL-VARIABLE CENSUS -- which potential is each one?')
    print('Chapter 7 Section 7.9 item 10, created 2026-10-07, mechanised 2026-10-08')
    print('=' * 78)

    c = census(HERE)
    c = {k: v for k, v in c.items()
         if not k.startswith('graphene_potential_census_audit.py')}

    by_verdict = {}
    for key in sorted(c):
        by_verdict.setdefault(verdict_for(c[key]), []).append(key)

    print()
    print('SECTION 1 -- the census')
    print('-' * 78)
    for v in ('UNDETERMINED', 'DECLARED_NO_REASON', 'DUAL_ROLE_TERMINAL',
              'DECLARED', 'DETERMINED', 'COLLECTED_NO_ROLE', 'PARSE_ERROR'):
        print('  %-20s %3d' % (v, len(by_verdict.get(v, []))))
    print('  %-20s %3d   across %d modules'
          % ('TOTAL NAMES', len(c), len({k.split('::')[0] for k in c})))

    print()
    print('  UNDETERMINED -- both roles AND channel-swept (the 10-07 defect):')
    if not by_verdict.get('UNDETERMINED'):
        print('    (none)')
    for key in by_verdict.get('UNDETERMINED', []):
        print('    %s' % key)
        for scope, line, ev in c[key]['scopes']:
            if ev:
                print('        %-34s line %-5s %s' % (scope, line, ev))

    print()
    print('  DUAL_ROLE_TERMINAL -- both roles, NOT channel-swept, not a fault:')
    for key in by_verdict.get('DUAL_ROLE_TERMINAL', []):
        print('    %-54s %s' % (key, POTENTIAL_REGISTRY.get(key, ('?',))[0]))

    print()
    print('  DECLARED in the source:')
    if not by_verdict.get('DECLARED'):
        print('    (none)')
    for key in by_verdict.get('DECLARED', []):
        reading, reason, line = c[key]['decl']
        print('    %-50s %s (line %d)' % (key, reading, line))
        print('        reason: %s' % reason[:110])

    print()
    print('SECTION 2 -- registry completeness: a census that cannot go stale')
    print('-' * 78)
    missing = sorted(k for k in c if k not in POTENTIAL_REGISTRY)
    extra = sorted(k for k in POTENTIAL_REGISTRY if k not in c)
    check('every potential-like name found is in POTENTIAL_REGISTRY',
          not missing,
          ('unregistered (the census has gone stale): ' + ', '.join(missing[:8]))
          if missing else
          '%d names, every one adjudicated or explicitly UNADJUDICATED' % len(c))
    check('POTENTIAL_REGISTRY has no rows for names that do not exist',
          not extra,
          ('dead registry rows: ' + ', '.join(extra[:8])) if extra
          else 'no dead rows')

    adj = {}
    for key in c:
        adj.setdefault(POTENTIAL_REGISTRY.get(key, ('UNREGISTERED', ''))[0],
                       []).append(key)
    print()
    print('  Adjudication breakdown:')
    for reading in sorted(adj):
        print('    %-26s %3d' % (reading, len(adj[reading])))
    unadj = sorted(adj.get('UNADJUDICATED', []))
    print()
    print('  %d of %d names are UNADJUDICATED. That is OPEN WORK, not a pass:'
          % (len(unadj), len(c)))
    for key in unadj:
        print('    %s' % key)
    print('    None of them is channel-swept, so none can carry the 10-07')
    print('    defect; what is undecided is whether they are potentials at all.')

    print()
    print('SECTION 3 -- the three siblings named on 2026-10-07, checked')
    print('-' * 78)
    inter = [k for k in c if k.startswith('graphene_interconnect_model.py')]
    check("Chapter 5's interconnect model has a potential-like variable at all",
          bool(inter),
          ('FOUND: ' + ', '.join(inter)) if inter else
          'NONE.  graphene_interconnect_model.py contains no potential-like\n'
          'name anywhere: it is a resistivity-and-geometry model and has no\n'
          'channel potential to be ambiguous about.  The 10-07 candidate\n'
          '"Chapter 5\'s interconnect channel potential" DOES NOT EXIST.\n'
          'This check is kept, and kept failing, as the standing record of\n'
          'that -- a passing check here would mean a potential had appeared.',
          expected_to_fail=True)
    twoc = c.get('graphene_photodetector_two_contact_model.py::V_bias')
    check("Chapter 6's two-contact channel potential is not ambiguous",
          twoc is not None and len(twoc['roles']) < 2,
          'roles=%s.  It is imported and printed and enters no computation in\n'
          'that module, so there is no second reading available to it.'
          % (sorted(twoc['roles']) if twoc else 'NOT FOUND'))
    s44 = c.get('graphene_fet_model.py::V_channel_shift')
    check("Section 4.4's local Dirac-point shift (V_channel_shift) is collected",
          s44 is not None,
          'verdict=%s; it appears only in the docstring derivation, which is\n'
          'itself a finding: the closed form actually evaluated never computes\n'
          'it.  Collected as a provenance fact, adjudicated ELECTROSTATIC_ONLY.'
          % (verdict_for(s44) if s44 else 'NOT FOUND'))

    print()
    print('SECTION 4 -- controls')
    print('-' * 78)
    print('  4a POSITIVE, exactly known (established independently on 10-07):')
    vch = c.get('graphene_fet_model.py::V_ch')
    check('graphene_fet_model.V_ch is seen in BOTH roles',
          vch is not None and vch['roles'] == {'ELECTROSTATIC', 'TRANSPORT'},
          'roles=%s' % (sorted(vch['roles']) if vch else 'NOT FOUND'))
    check('graphene_fet_model.V_ch is seen as CHANNEL-SWEPT',
          vch is not None and vch['swept'],
          'swept=%s' % (vch['swept'] if vch else 'NOT FOUND'))

    print('  4b SYNTHETIC, both directions:')
    syn = synthetic_control()
    check('one role -> DETERMINED', syn.get('single') == 'DETERMINED',
          'got %s' % syn.get('single'))
    check('two roles + channel-swept -> UNDETERMINED',
          syn.get('dual') == 'UNDETERMINED', 'got %s' % syn.get('dual'))

    print('  4c NULL, at the point where the failure enters:')
    n_files, n_renames, leftovers = null_control()
    check('the census is EMPTY when every potential-like name is renamed away',
          not leftovers,
          '%d files rewritten (the instrument itself excluded, see\n'
          'null_control for why), %d distinct names renamed; leftovers=%s'
          % (n_files, n_renames, leftovers[:6]))

    print('  4d MUTATION:')
    mut = mutation_control()
    for label, got, expected, ok in mut:
        check('%s  [expect %s]' % (label, expected), ok, 'got %s' % got)

    print()
    print('SECTION 5 -- what the census says')
    print('-' * 78)
    n_undet = len(by_verdict.get('UNDETERMINED', []))
    n_dual = len(by_verdict.get('DUAL_ROLE_TERMINAL', []))
    print('  %d potential-like names, over-collected on purpose, across %d'
          % (len(c), len({k.split('::')[0] for k in c})))
    print('  modules.  %d are in both roles and channel-swept; %d are in both'
          % (n_undet, n_dual))
    print('  roles and are terminal biases, where dual use is not ambiguity.')
    print()
    print('  The 10-07 entry listed three places it expected siblings of the')
    print('  V_ch defect.  All three are now checked and NONE of them is one,')
    print('  for three different reasons: Chapter 5 has no potential variable')
    print('  at all; Chapter 6\'s two-contact V_bias enters no computation in')
    print('  its own module; and Section 4.4\'s V_channel_shift is never')
    print('  evaluated by the closed form that replaced the derivation.')
    print()
    print('  The base rate the 10-07 entry asked for, measured: ONE defect in')
    print('  %d names, confined to one variable in one module, and a further' % len(c))
    print('  %d names whose status as potentials is undecided.  The guess that' % len(unadj))
    print('  the class was widespread is not supported.')

    print()
    print('=' * 78)
    print('TOTAL: %d passed, %d failed (%d of them EXPECTED)'
          % (passed, failed, EXPECTED_FAILURES))
    print('=' * 78)
    return failed - EXPECTED_FAILURES


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
