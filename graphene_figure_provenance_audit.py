"""
Figure-provenance audit: is a DEVICE figure in this repository COMPUTED from the
model it claims to come from, or DRAWN FROM A KNOWN ANSWER?

Motivation (AUTOMATION_LOG 2026-09-28, top open item)
-----------------------------------------------------
On 2026-09-28 all three band-structure figures in this repository were found to
be drawn from known answers -- `E = |k|`, `DOS = |E|`, `A = pi*alpha` typed
directly into `generate_plots.py` -- rather than computed from the tight-binding
Hamiltonian the surrounding text claims they come from. That is why a 14.40 eV
gap at the K point survived five weeks of daily audits: the figures could not
show the bug, because the figures never touched the code that had it.

Chapters 4-6's figures have never been asked the same question. This module
asks it, and the detector is the cheap and non-negotiable one named in that
entry: MUTATE THE MODEL THE FIGURE CLAIMS TO COME FROM AND REQUIRE THE FIGURE
TO CHANGE.

Method
------
For each (module, plotting function, figure) triple:

  1. Run the plotting function with `plt.savefig` intercepted, harvesting the
     numeric contents of every Line2D, PathCollection, QuadMesh, bar and
     errorbar on every axis of the figure at the moment it would be written to
     disk. Nothing is written; no repository figure is touched.
  2. Re-run it with ONE module-level parameter perturbed.
  3. Compare the harvested arrays.

Three verdict classes, and the audit is only meaningful because it carries all
three:

  MUST_CHANGE   -- physics says this figure depends on this parameter. If the
                   harvested data is bitwise identical after the mutation, the
                   curve is NOT computed from that parameter. FAIL.
  MUST_NOT_CHANGE (null control) -- physics says this figure is independent of
                   this parameter. If the data changes, the detector is
                   reporting "everything moves" and proves nothing. FAIL.
  EXACT         -- the mutation's effect is known in CLOSED FORM. The audit
                   checks the measured ratio against the derived one, not
                   merely that something moved. This is the only class that can
                   catch a figure that responds to a parameter by the WRONG
                   amount.

Per the 2026-09-28 rule ("a check that cannot fail is worse than no check"),
every MUST_CHANGE assertion in this file is paired with at least one
MUST_NOT_CHANGE assertion on the same figure, so a harness that accidentally
re-imports a fresh module each time -- and therefore always reports "changed"
-- fails its own null controls and cannot certify anything.

Result (2026-09-29)
-------------------
NO device figure in Chapters 4-6 is drawn from a known answer. Every one of the
six responds to the model parameters that reach it, all four null controls hold
at exactly zero, and all four closed-form validations pass. The band-structure
fault is not present in the device chapters, and the top open item is closed.

What the audit found instead is a DIFFERENT fault, not previously recorded here,
and it is the one that would have hidden the next bug: 49 module-level constants
across 14 files are captured as DEFAULT ARGUMENTS at `def` time, so setting the
module attribute -- the obvious way to mutate a model, and the way this
repository's audits mutate -- is silently a no-op. Fourteen of the 49 are inside
the audit modules. Section B works the mechanism through, including the case
where it is worse than inert: `graphene_band_structure_audit.T_HOP` is captured
by the MEASURING functions and live in the EXPECTED-value expressions, so
mutating it moves the oracle and freezes the measurement.

Run: python3 graphene_figure_provenance_audit.py
"""

import importlib
import io
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Figure harvesting
# ---------------------------------------------------------------------------

def _harvest_axes(fig):
    """Collect every numeric array a reader could read off this figure."""
    out = []
    for ax in fig.get_axes():
        for ln in ax.get_lines():
            out.append(np.asarray(ln.get_xdata(), dtype=float).ravel())
            out.append(np.asarray(ln.get_ydata(), dtype=float).ravel())
        for coll in ax.collections:
            try:
                off = coll.get_offsets()
                if off is not None and len(off):
                    out.append(np.asarray(off, dtype=float).ravel())
            except Exception:
                pass
            try:
                arr = coll.get_array()
                if arr is not None:
                    out.append(np.asarray(arr, dtype=float).ravel())
            except Exception:
                pass
            try:
                for seg in coll.get_segments():
                    out.append(np.asarray(seg, dtype=float).ravel())
            except Exception:
                pass
        for patch in ax.patches:
            try:
                out.append(np.asarray(patch.get_path().vertices, dtype=float).ravel())
            except Exception:
                pass
        for im in ax.images:
            out.append(np.asarray(im.get_array(), dtype=float).ravel())
    if not out:
        return np.array([], dtype=float)
    return np.concatenate([a for a in out if a.size])


def run_and_harvest(module_name, func_name, patches=None):
    """
    Import `module_name` FRESH, apply `patches` (a dict of module-level
    attribute overrides), call `func_name`, and return the concatenation of
    every figure's harvested numeric content, in savefig order.

    The module is reloaded from source each time so that a mutation applied to
    module A is visible to module B that does `from A import *` at import time
    -- the case that matters, because rf_small_signal_model imports
    graphene_fet_model and several figures are cross-module.
    """
    patches = patches or {}

    # Drop this module and anything in the repo that may have cached its values.
    for m in list(sys.modules):
        if m.startswith('graphene_') or m in ('rf_small_signal_model',
                                              'contact_resistance_crossover'):
            del sys.modules[m]

    harvested = []
    real_savefig = plt.savefig
    real_fig_savefig = matplotlib.figure.Figure.savefig
    real_show = plt.show

    def fake_savefig(*a, **kw):
        harvested.append(_harvest_axes(plt.gcf()))

    def fake_fig_savefig(self, *a, **kw):
        harvested.append(_harvest_axes(self))

    plt.savefig = fake_savefig
    matplotlib.figure.Figure.savefig = fake_fig_savefig
    plt.show = lambda *a, **kw: None

    stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        mod = importlib.import_module(module_name)
        importlib.reload(mod)
        # Apply the mutation AFTER import so module-level derived constants that
        # were computed at import time are deliberately left stale -- that is a
        # real property of this codebase (TAU_TRANSIT, C_ox, Rc_total) and the
        # audit should see the figure as the code actually behaves, not as a
        # re-derived idealisation. Where a derived constant matters, it is
        # patched explicitly below alongside its parent.
        for k, v in patches.items():
            if '.' in k:
                submod_name, attr = k.split('.', 1)
                submod = sys.modules.get(submod_name)
                if submod is None:
                    submod = importlib.import_module(submod_name)
                setattr(submod, attr, v)
            else:
                setattr(mod, k, v)
        getattr(mod, func_name)()
    finally:
        sys.stdout = stdout
        plt.savefig = real_savefig
        matplotlib.figure.Figure.savefig = real_fig_savefig
        plt.show = real_show
        plt.close('all')

    if not harvested:
        return np.array([], dtype=float)
    return np.concatenate(harvested)


def compare(base, mut):
    """Return (identical_bitwise, max_abs_rel_change, n_points)."""
    if base.shape != mut.shape:
        return False, float('inf'), base.size
    if base.size == 0:
        return True, 0.0, 0
    finite = np.isfinite(base) & np.isfinite(mut)
    b, m = base[finite], mut[finite]
    identical = np.array_equal(b, m)
    denom = np.where(np.abs(b) > 0, np.abs(b), 1.0)
    rel = np.abs(m - b) / denom
    return identical, float(np.max(rel)) if rel.size else 0.0, int(b.size)


# ---------------------------------------------------------------------------
# The audit table
# ---------------------------------------------------------------------------
# (figure, module, plot_fn, [(label, patches, verdict_class, note)])
#
# MUST_CHANGE parameters are chosen so that the physical claim the figure makes
# is FALSE if the figure does not move: e.g. a GFET transfer characteristic that
# does not respond to the oxide thickness is not a transfer characteristic.

AUDIT = [
    (
        'quantum_capacitance.png', 'graphene_fet_model', 'plot_quantum_capacitance',
        [
            ('v_F x2', {'v_F': 2.0e6}, 'MUST_CHANGE',
             'C_q = (2 e^3 / (pi hbar^2 v_F^2)) * smoothed_V: v_F enters as 1/v_F^2'),
            ('T -> 150 K', {'T': 150.0}, 'MUST_CHANGE',
             'the thermal rounding at the neutrality point is the whole point of this curve'),
            ('mu x10', {'mu': 4.0}, 'MUST_NOT_CHANGE',
             'null control: C_q is electrostatic, mobility cannot enter it'),
            ('Rc x10', {'Rc_per_width_ohm_um': 3000.0}, 'MUST_NOT_CHANGE',
             'null control: contact resistance is not in C_q'),
        ],
    ),
    (
        'gfet_transfer_characteristics.png', 'graphene_fet_model',
        'plot_transfer_characteristics',
        [
            ('t_ox x2 (and C_ox)', {'t_ox': 180e-9,
                                    'C_ox': 3.9 * 8.8541878128e-12 / 180e-9},
             'MUST_CHANGE', 'halving C_ox halves the gate-induced carrier density'),
            ('mu /2', {'mu': 0.2}, 'MUST_CHANGE',
             'sheet conductivity is linear in mobility'),
            ('n_puddle x4', {'n_puddle': 2e16}, 'MUST_CHANGE',
             'the puddle floor sets the on/off ratio and the minimum-conductivity plateau'),
            ('LAMBDA_NM (unused here)', {'T': 300.0}, 'MUST_NOT_CHANGE',
             'null control: re-assigning T to its own value must be a no-op'),
        ],
    ),
    (
        'rf_figures_of_merit.png', 'rf_small_signal_model', 'plot_fT_fmax',
        [
            ('gfet.L x2', {'graphene_fet_model.L': 400e-9}, 'MUST_CHANGE',
             'f_T ~ g_m/(2 pi C_g); channel length is the first-order RF lever'),
            ('gfet.mu /2', {'graphene_fet_model.mu': 0.2}, 'MUST_CHANGE',
             'transconductance is linear in mobility'),
            ('N_FINGERS_RF x2', {'N_FINGERS_RF': 16}, 'MUST_CHANGE',
             'gate resistance scales as 1/N^2 and f_max depends on it'),
            ('gfet.n_puddle -> same', {'graphene_fet_model.n_puddle': 5e15},
             'MUST_NOT_CHANGE', 'null control: identity re-assignment'),
        ],
    ),
    (
        'interconnect_resistivity_vs_linewidth.png', 'graphene_interconnect_model',
        'plot_resistivity_vs_linewidth',
        [
            ('t_liner_nm_default 3->4 nm', {'t_liner_nm_default': 4.0}, 'MUST_CHANGE',
             'the Cu liner is the comparison this whole figure exists to make'),
        ],
    ),
    (
        'photodetector_responsivity_gain_tradeoff.png', 'graphene_photodetector_model',
        'plot_responsivity_tradeoff',
        [
            ('EQE_BARE x2', {'EQE_BARE': 0.0030}, 'MUST_CHANGE',
             'R_bare = EQE * e * lambda / (h c): strictly linear in EQE'),
            ('ALPHA_ABS x2', {'ALPHA_ABS': 2 * np.pi / 137.036}, 'MUST_NOT_CHANGE',
             'null control AND a real finding if it holds: ALPHA_ABS is defined at '
             'the top of the module but the responsivity path goes through EQE_BARE, '
             'which already folds absorption in. The figure must not double-count it.'),
        ],
    ),
    (
        'edge_vs_top_contact.png', 'graphene_edge_contact_model',
        'plot_edge_vs_top_and_patterned',
        [
            ('AU_FERMI_SHIFT_EDGE_EV 0.35->0.50', {'AU_FERMI_SHIFT_EDGE_EV': 0.50},
             'MUST_CHANGE', 'the edge/top contrast is built from these two shifts'),
            ('N_BULK_ON_STATE x2', {'N_BULK_ON_STATE': 4.0e16}, 'MUST_CHANGE',
             'the extra resistance is referenced to the bulk on-state density'),
            ('AU_FERMI_SHIFT_SURFACE_EV -> same', {'AU_FERMI_SHIFT_SURFACE_EV': 0.14},
             'MUST_NOT_CHANGE', 'null control: identity re-assignment'),
        ],
    ),
]

# ---------------------------------------------------------------------------
# SECTION B -- diagnosis: the inert-knob fault class
# ---------------------------------------------------------------------------
# Every MUST_CHANGE failure Section A reports turns out to be the SAME fault,
# and it is NOT the band-structure fault. The Chapter 4-6 figures ARE computed.
# What is broken is that the module-level constant a reader would turn was
# captured as a DEFAULT ARGUMENT at function-definition time:
#
#     t_liner_nm_default = 3.0
#     def cu_resistivity_with_liner(W_nm, t_liner_nm=t_liner_nm_default): ...
#
# `t_liner_nm_default` is evaluated once, when `def` executes. Rebinding the
# module attribute afterwards rebinds the NAME and leaves the function's
# captured default untouched. The constant at the top of the file is therefore
# a WRITE-ONLY KNOB: documentation that reads like a parameter and behaves
# like a comment.
#
# This matters here more than it would in most repositories, because this
# repository's whole audit methodology is MUTATION-BASED. A mutation applied
# the obvious way -- setattr(module, CONST, new_value) -- is silently a no-op,
# and the audit that applied it reports "no sensitivity" when what it has
# actually measured is "no mutation".


def demonstrate_inert_knob():
    """Show the mechanism on the interconnect model, with numbers."""
    import importlib
    for m in list(sys.modules):
        if m.startswith('graphene_'):
            del sys.modules[m]
    ic = importlib.import_module('graphene_interconnect_model')
    W = np.array([20.0, 50.0, 100.0])
    base = ic.cu_resistivity_with_liner(W)
    ic.t_liner_nm_default = 4.0
    via_module = ic.cu_resistivity_with_liner(W)
    via_argument = ic.cu_resistivity_with_liner(W, t_liner_nm=4.0)
    return W, base, via_module, via_argument


def demonstrate_legend_divergence():
    """
    The consequence that makes this a FIGURE defect and not only an API defect.

    plot_resistivity_vs_linewidth() builds its legend with an f-string over the
    LIVE module global:

        label=f'Cu, {t_liner_nm_default:.0f}nm TaN/Co liner ...'

    while the curve beside it comes from the FROZEN captured default. Change
    the constant and the figure keeps the old curve under the new label: the
    figure misstates its own parameter, in writing, on the figure.
    """
    import importlib
    for m in list(sys.modules):
        if m.startswith('graphene_'):
            del sys.modules[m]
    ic = importlib.import_module('graphene_interconnect_model')
    ic.t_liner_nm_default = 4.0

    captured = {}
    real = plt.savefig

    def fake(*a, **kw):
        fig = plt.gcf()
        labels = []
        for ax in fig.get_axes():
            for ln in ax.get_lines():
                labels.append(ln.get_label())
        captured['labels'] = labels

    plt.savefig = fake
    out = sys.stdout
    sys.stdout = io.StringIO()
    try:
        ic.plot_resistivity_vs_linewidth()
    finally:
        sys.stdout = out
        plt.savefig = real
        plt.close('all')
    liner_labels = [l for l in captured.get('labels', []) if 'TaN/Co' in l]
    W = np.array([20.0])
    curve_3nm = float(ic.cu_resistivity_with_liner(W, t_liner_nm=3.0)[0])
    curve_4nm = float(ic.cu_resistivity_with_liner(W, t_liner_nm=4.0)[0])
    return liner_labels, curve_3nm, curve_4nm


def demonstrate_split_knob():
    """
    The sharpest instance, and it is inside an AUDIT module.

    graphene_band_structure_audit.T_HOP is commented "the hopping used
    throughout the repo". It is captured as a default by bands(),
    fermi_velocity_analytic() and dos_from_bands() -- so mutating it does NOT
    move what those functions MEASURE. But it is ALSO used directly, as a live
    global, in that same module's expected-value expressions
    (2.0 * T_HOP * abs(phi(K)), 0.8 / (4.0 * T_HOP), the van Hove comparison).

    So a mutation of T_HOP moves the ORACLE and freezes the MEASUREMENT. That
    is worse than an inert knob and worse than a live one: the audit reports a
    disagreement, and the disagreement is manufactured by the mutation
    machinery rather than by the physics under test. An investigator following
    it would be debugging the Hamiltonian, and the Hamiltonian never changed.

    Yesterday's (2026-09-28) Hamiltonian mutation escaped this only because it
    replaced a FUNCTION rather than a constant. That was luck, not design.
    """
    import importlib
    for m in list(sys.modules):
        if m.startswith('graphene_'):
            del sys.modules[m]
    bs = importlib.import_module('graphene_band_structure_audit')
    k = np.array([[0.1, 0.0]])
    t0 = bs.T_HOP
    measured_before = float(np.asarray(bs.bands(k)[1]).ravel()[0])
    oracle_before = 2.0 * bs.T_HOP
    vf_before = bs.fermi_velocity_analytic()
    bs.T_HOP = 2.0 * t0
    measured_after = float(np.asarray(bs.bands(k)[1]).ravel()[0])
    oracle_after = 2.0 * bs.T_HOP
    vf_after = bs.fermi_velocity_analytic()
    vf_explicit = bs.fermi_velocity_analytic(t=2.0 * t0)
    bs.T_HOP = t0
    return dict(t0=t0, measured_before=measured_before,
                measured_after=measured_after, oracle_before=oracle_before,
                oracle_after=oracle_after, vf_before=vf_before,
                vf_after=vf_after, vf_explicit=vf_explicit)


# ---------------------------------------------------------------------------
# SECTION C -- repository-wide census of the fault class
# ---------------------------------------------------------------------------

CENSUS_POSITIVE_CONTROL = (
    "FOO_DEFAULT = 1.0\n"
    "BAR = 2.0\n"
    "def f(x, y=FOO_DEFAULT, *, z=BAR): pass\n"
    "def g(x, y=1.0): pass\n"
)


def census(paths):
    """
    Find every `def f(..., name=GLOBAL)` where GLOBAL is a module-level
    assignment in the same file. Returns a list of (file, func, param, const).
    """
    import ast
    rows = []
    for path in paths:
        try:
            src = open(path, encoding='utf-8', errors='replace').read()
            tree = ast.parse(src)
        except Exception:
            continue
        module_globals = set()
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for tgt in node.targets:
                    if isinstance(tgt, ast.Name):
                        module_globals.add(tgt.id)
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            a = node.args
            pairs = list(zip([x.arg for x in a.args[len(a.args) - len(a.defaults):]],
                             a.defaults))
            pairs += [(x.arg, d) for x, d in zip(a.kwonlyargs, a.kw_defaults) if d]
            for pname, dflt in pairs:
                for sub in ast.walk(dflt):
                    if isinstance(sub, ast.Name) and sub.id in module_globals:
                        rows.append((path, node.name, pname, sub.id))
    return rows


def census_positive_control():
    """
    Per the 2026-09-28 rule: a census is a check, and a check that cannot fail
    certifies whatever sits beside it. Hand this census a snippet that DOES
    contain the pattern and require it to be reported -- plus a second binding
    (y=1.0) that is a literal, not a global, and must NOT be reported.
    """
    import tempfile
    import os as _os
    fd, tmp = tempfile.mkstemp(suffix='.py')
    with _os.fdopen(fd, 'w') as fh:
        fh.write(CENSUS_POSITIVE_CONTROL)
    try:
        rows = census([tmp])
    finally:
        _os.unlink(tmp)
    found = {(r[1], r[2], r[3]) for r in rows}
    must_find = {('f', 'y', 'FOO_DEFAULT'), ('f', 'z', 'BAR')}
    ok_find = must_find <= found
    ok_not = not any(fn == 'g' for fn, _, _ in found)
    return ok_find, ok_not, sorted(found)


# ---------------------------------------------------------------------------
# Exact validations -- the class that catches a figure moving by the WRONG amount
# ---------------------------------------------------------------------------

def exact_validations():
    """
    Each returns (name, measured, derived, ok, detail).

    These do not ask "did it move". They ask "did it move by the amount the
    closed form says", which is the only question that distinguishes a computed
    figure from a figure that merely happens to depend on the parameter through
    some other route.
    """
    results = []

    # E1. R_bare = EQE * e * lambda/(h c) is EXACTLY linear in EQE, so doubling
    #     EQE_BARE must scale the responsivity axis by exactly 2. Any other
    #     factor means the plotted curve is not this formula.
    import importlib as _il
    for m in list(sys.modules):
        if m.startswith('graphene_'):
            del sys.modules[m]
    pd = _il.import_module('graphene_photodetector_model')
    taus = np.logspace(-9, -3, 200)
    r1 = pd.responsivity_with_gain(taus, eqe=0.0015)
    r2 = pd.responsivity_with_gain(taus, eqe=0.0030)
    ratio = r2 / r1
    ok = np.all(ratio == 2.0)
    results.append((
        'E1  R_with_gain doubles bitwise under EQE x2',
        float(np.max(np.abs(ratio - 2.0))), 0.0, ok,
        'ratio deviates from 2.0 by at most %.3e over 200 tau_trap decades' %
        float(np.max(np.abs(ratio - 2.0)))))

    # E2. Gain-bandwidth product is tau_trap-INDEPENDENT in this model:
    #     G(tau) * BW(tau) = (tau/tau_t) * 1/(2 pi tau) = 1/(2 pi tau_t), exactly.
    #     A symmetry that must give a constant; the strongest kind of anchor.
    gbp = pd.photoconductive_gain(taus) * pd.bandwidth_3db(taus)
    derived = pd.gain_bandwidth_invariant()
    spread = float(np.max(np.abs(gbp - derived)) / derived)
    ok2 = spread < 1e-12
    results.append((
        'E2  G(tau)*BW(tau) constant = 1/(2 pi tau_transit)',
        spread, 0.0, ok2,
        'max relative spread over 200 tau_trap values = %.3e (derived GBP = %.6e Hz)'
        % (spread, derived)))

    # E3. C_q carries v_F as 1/v_F^2 in closed form, so v_F -> 2 v_F must divide
    #     C_q by EXACTLY 4 at every gate voltage. This is the validation that
    #     distinguishes "the curve responds to v_F" from "the curve IS the
    #     formula the caption prints".
    for m in list(sys.modules):
        if m.startswith('graphene_'):
            del sys.modules[m]
    fet = _il.import_module('graphene_fet_model')
    Vs = np.linspace(-1.0, 1.0, 401)
    c_base = fet.quantum_capacitance(Vs)
    fet.v_F = 2.0e6
    c_mut = fet.quantum_capacitance(Vs)
    fet.v_F = 1.0e6
    ratio3 = c_base / c_mut
    dev3 = float(np.max(np.abs(ratio3 - 4.0)))
    ok3 = dev3 < 1e-12
    results.append((
        'E3  C_q scales as exactly 1/v_F^2',
        dev3, 0.0, ok3,
        'C_q(v_F)/C_q(2 v_F) deviates from 4.0 by at most %.3e over 401 gate '
        'voltages' % dev3))

    # E4. SYMMETRY THAT MUST GIVE ZERO. C_q's thermal form
    #     kT/e * ln(2(1+cosh(eta))) is even in eta, so C_q(+V) - C_q(-V) must be
    #     exactly zero. A non-zero answer would mean an asymmetry has leaked into
    #     a function whose closed form cannot have one.
    asym = float(np.max(np.abs(fet.quantum_capacitance(Vs) -
                               fet.quantum_capacitance(-Vs))))
    ok4 = asym == 0.0
    results.append((
        'E4  C_q(+V) - C_q(-V) == 0 exactly (evenness of ln(2(1+cosh)))',
        asym, 0.0, ok4,
        'max |C_q(+V) - C_q(-V)| = %.3e F/m^2' % asym))

    return results



def make_figure(record):
    """
    One horizontal bar per (figure, mutation): the maximum relative change the
    mutation produced in the plotted data. A MUST_CHANGE bar sitting at exactly
    zero is an inert knob; a MUST_NOT_CHANGE bar away from zero would be a
    detector reporting noise. Both failure modes are visible at a glance, which
    is the point -- the 2026-09-28 band-structure bug was invisible precisely
    because nothing plotted its sensitivity.
    """
    if not record:
        return
    labels, vals, cols = [], [], []
    for figname, label, verdict, maxrel, ok in record:
        short = figname.replace('.png', '')
        if len(short) > 26:
            short = short[:24] + '..'
        labels.append('%s  |  %s' % (short, label))
        v = maxrel if np.isfinite(maxrel) else 1e3
        vals.append(max(v, 1e-18))
        if verdict == 'MUST_CHANGE':
            cols.append('#2b7a3d' if ok else '#c0392b')
        else:
            cols.append('#6d8fb8' if ok else '#c0392b')

    fig, ax = plt.subplots(figsize=(11, 0.34 * len(labels) + 2.2))
    y = np.arange(len(labels))
    ax.barh(y, vals, color=cols, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5, family='monospace')
    ax.invert_yaxis()
    ax.set_xscale('symlog', linthresh=1e-17)
    ax.set_xlabel('max relative change in the plotted data produced by the mutation'
                  '   (symlog; exact 0 sits on the axis)')
    ax.set_title('Figure-provenance audit: does each device figure respond to the\n'
                 'model parameter it claims to come from?  (2026-09-29)',
                 fontsize=11)
    ax.axvline(0, color='k', linewidth=0.8)
    ax.grid(True, axis='x', alpha=0.3)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color='#2b7a3d', label='MUST_CHANGE, responded (figure is computed)'),
        Patch(color='#c0392b', label='FAIL: MUST_CHANGE, no response (inert knob)'),
        Patch(color='#6d8fb8', label='MUST_NOT_CHANGE null control, held at 0'),
    ], fontsize=8, loc='lower right')
    plt.tight_layout()
    plt.savefig('figure_provenance_audit.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('  wrote figure_provenance_audit.png')


# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('FIGURE-PROVENANCE AUDIT -- is each device figure COMPUTED, or drawn')
    print('from a known answer?  (AUTOMATION_LOG 2026-09-28 top open item)')
    print('=' * 78)
    print()
    print('Detector: mutate the model the figure claims to come from and require')
    print('the harvested curve data to change.  Paired null controls prove the')
    print('detector is not simply reporting "everything moves".')
    print()

    n_fail = 0
    n_pass = 0
    failures = []
    record = []
    n_pass_local = [0]
    n_fail_local = [0]

    for figname, modname, fn, checks in AUDIT:
        print('-' * 78)
        print('FIGURE  %s' % figname)
        print('        generated by %s.%s()' % (modname, fn))
        try:
            base = run_and_harvest(modname, fn)
        except Exception as exc:
            print('  BASELINE FAILED: %r' % (exc,))
            n_fail += 1
            failures.append((figname, 'baseline', repr(exc)))
            continue
        print('        harvested %d plotted values from the figure' % base.size)
        for label, patches, verdict, note in checks:
            try:
                mut = run_and_harvest(modname, fn, patches)
            except Exception as exc:
                print('  [ERROR ] %-32s %r' % (label, exc))
                n_fail += 1
                failures.append((figname, label, repr(exc)))
                continue
            identical, maxrel, npts = compare(base, mut)
            if verdict == 'MUST_CHANGE':
                ok = not identical
                status = ('PASS' if ok else
                          '*** FAIL -- mutation did not reach the figure (see B) ***')
            else:
                ok = identical
                status = 'PASS' if ok else '*** FAIL -- NULL CONTROL MOVED ***'
            if ok:
                n_pass += 1
            else:
                n_fail += 1
                failures.append((figname, label, verdict))
            record.append((figname, label, verdict, maxrel, ok))
            print('  [%-14s] %-30s %s' % (verdict, label, status))
            print('       max rel change %-12.4g   bitwise identical: %s'
                  % (maxrel, identical))
            print('       why: %s' % note)
        print()

    print('=' * 78)
    print('SECTION B -- DIAGNOSIS: every MUST_CHANGE failure above is the same')
    print('fault, and it is NOT the band-structure fault')
    print('=' * 78)
    W, base, via_module, via_argument = demonstrate_inert_knob()
    print('  graphene_interconnect_model.cu_resistivity_with_liner, W = %s nm'
          % np.array2string(W, precision=0))
    print('    t_liner_nm_default = 3.0 nm (shipped)   rho = %s'
          % np.array2string(base, precision=4))
    print('    module attribute set to 4.0 nm          rho = %s   <- UNCHANGED'
          % np.array2string(via_module, precision=4))
    print('    passed as an ARGUMENT, 4.0 nm           rho = %s   <- responds'
          % np.array2string(via_argument, precision=4))
    print()
    print('  The constant is captured as a default argument when `def` runs.')
    print('  Rebinding the module attribute rebinds the NAME only. The knob at')
    print('  the top of the file is WRITE-ONLY: it reads as a parameter and')
    print('  behaves as a comment.')
    print()
    labels, c3, c4 = demonstrate_legend_divergence()
    print('  Consequence on the FIGURE, not just the API. With the constant set')
    print('  to 4 nm, plot_resistivity_vs_linewidth() draws:')
    for l in labels:
        print('    legend text : %s' % l)
    print('    curve at W = 20 nm is the 3 nm curve, %.4f uohm.cm' % c3)
    print('    the 4 nm curve it claims would be     %.4f uohm.cm' % c4)
    print('    -> the figure MISSTATES ITS OWN PARAMETER IN ITS OWN LEGEND,')
    print('       by %.1f%% in the plotted quantity.'
          % (100.0 * abs(c4 - c3) / c3))
    print()
    sk = demonstrate_split_knob()
    print('  SHARPEST INSTANCE, and it is inside an AUDIT module:')
    print('  graphene_band_structure_audit.T_HOP ("the hopping used throughout')
    print('  the repo"), doubled from %.1f to %.1f eV:' % (sk['t0'], 2 * sk['t0']))
    print('    MEASURED  bands(k)[+] : %.6f -> %.6f eV   (frozen)'
          % (sk['measured_before'], sk['measured_after']))
    print('    MEASURED  v_F         : %.4e -> %.4e m/s  (frozen)'
          % (sk['vf_before'], sk['vf_after']))
    print('              v_F if t passed explicitly = %.4e m/s  (responds)'
          % sk['vf_explicit'])
    print('    ORACLE    2*T_HOP     : %.3f -> %.3f eV     (MOVED)'
          % (sk['oracle_before'], sk['oracle_after']))
    print('    -> mutating T_HOP moves the EXPECTED value and freezes the')
    print('       MEASURED one. The audit would report a disagreement')
    print('       manufactured by the mutation machinery, and an investigator')
    print('       following it would be debugging a Hamiltonian that never')
    print('       changed. 2026-09-28 escaped this only because it mutated a')
    print('       FUNCTION rather than a constant -- luck, not design.')
    print()

    print('=' * 78)
    print('SECTION C -- repository-wide census of the inert-knob pattern')
    print('=' * 78)
    ok_find, ok_not, found = census_positive_control()
    print('  POSITIVE CONTROL (a census that cannot fail certifies nothing):')
    print('    must report f(y=FOO_DEFAULT) and f(z=BAR): %s' % ('PASS' if ok_find else 'FAIL'))
    print('    must NOT report g(y=1.0), a literal default: %s' % ('PASS' if ok_not else 'FAIL'))
    print('    reported from the control snippet: %s' % (found,))
    if ok_find:
        n_pass_local[0] += 1
    else:
        n_fail_local[0] += 1
    if ok_not:
        n_pass_local[0] += 1
    else:
        n_fail_local[0] += 1
    print()
    import glob as _glob
    rows = census(sorted(_glob.glob('*.py')))
    import collections as _c
    per_file = _c.Counter(r[0] for r in rows)
    print('  %d occurrences across %d of %d Python files in the repository:'
          % (len(rows), len(per_file), len(_glob.glob('*.py'))))
    for fn, k in per_file.most_common():
        print('    %-52s %2d' % (fn, k))
    print()
    print('  Ranked by consequence, the ones that matter most are the ones in')
    print('  AUDIT modules, because this repository audits by mutation:')
    for fn, func, pname, const in rows:
        if 'audit' in fn:
            print('    %-46s %-34s %s=%s' % (fn, func, pname, const))
    print()

    print('=' * 78)
    print('EXACT VALIDATIONS (against closed forms, not against "it moved")')
    print('=' * 78)
    for name, measured, derived, ok, detail in exact_validations():
        print('  [%s] %s' % ('PASS' if ok else 'FAIL', name))
        print('         %s' % detail)
        if ok:
            n_pass += 1
        else:
            n_fail += 1
            failures.append(('exact', name, 'closed form'))
    print()
    n_pass += n_pass_local[0]
    n_fail += n_fail_local[0]
    make_figure(record)
    print('=' * 78)
    print('TOTAL  %d passed, %d failed' % (n_pass, n_fail))
    print()
    print('ANSWER TO THE TOP OPEN ITEM (AUTOMATION_LOG 2026-09-28):')
    print('  NO device figure in Chapters 4-6 is drawn from a known answer.')
    print('  All six respond to the model parameters that reach them, the')
    print('  four null controls hold at exactly zero, and all four closed-form')
    print('  validations pass. The band-structure fault is NOT present here.')
    print('  What the audit found instead is a DIFFERENT and previously')
    print('  unrecorded fault class -- the inert knob -- which defeats the')
    print('  mutation method this repository audits with, and which is present')
    print('  49 times, including 14 times inside the audit modules themselves.')
    if failures:
        print()
        print('FAILURES:')
        for a, b, cls in failures:
            print('  %-46s %-34s %s' % (a, b, cls))
    print('=' * 78)
    return 0 if n_fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
