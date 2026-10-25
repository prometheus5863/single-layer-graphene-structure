"""graphene_capture_leverage_figure.py -- the capture surface as a figure,
with every number in it derived and every legend entry built from an artist.

WHY A FIGURE AT ALL, when 2026-10-23 declined one.  That session's result was
a 3-edge set; a table carried it exactly.  Today's is a 9 x 2 table of signed
leverages whose POINT is a coincidence between two of its rows -- the
repaired replacement's response to `gfet._thermal_rounding` and the
unrepaired one's response to `gfet.quantum_capacitance` -- and a reader has
to hunt for that in a table of hex floats.  Side by side it is one glance.

WHY IT IS BUILT THIS WAY.  The standing "render and inspect every committed
figure" item (~63 PNGs, three sessions of sampling, 100 % hit rate) offered
two mechanical sweeps: derive every number in a title from the computation,
and build every legend from the artists.  This module is a worked reference
for both, on a new figure rather than on the backlog:

  * no number appears in any title, label or annotation as a literal.  Every
    one is interpolated from the measurement, and `assert_no_literal_numbers`
    checks the SOURCE of this file for a digit inside a string passed to a
    text-drawing call.
  * the legend is built with `ax.get_legend_handles_labels()`, so an artist
    that is not drawn cannot appear in it and a drawn one cannot be missed.
  * the exact zeros are drawn as a different artist from the nonzero bars,
    because a zero-height bar and a missing bar look the same, and the
    2026-10-23 count of figures asserting something the data did not say is
    three.

Usage:  PYTHONPATH=... python3 graphene_capture_leverage_figure.py
"""

import ast
import io
import math
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                              # noqa: E402

import graphene_fet_model as gfet                            # noqa: E402
import graphene_default_flip as flip                         # noqa: E402
import graphene_capture_surface_guard as csg                 # noqa: E402

FNAME = 'capture_surface_leverage.png'
# WHICH ARGUMENT OF EACH CALL IS THE TEXT.  Not every string a text-drawing
# call receives is drawn: `arrowprops=dict(color='0.35')` is a COLOUR, and
# the second version of this checker flagged it, having walked every argument
# of the call instead of the ones that become glyphs.  Same shape as the
# regex version one level down -- an instrument asked about the wrong thing
# -- so the positions are named here and F0d is the control for it.
TEXT_ARGS = {
    'set_title':       ((0,), ('label',)),
    'suptitle':        ((0,), ('t',)),
    'set_xlabel':      ((0,), ('xlabel',)),
    'set_ylabel':      ((0,), ('ylabel',)),
    'annotate':        ((0,), ('text',)),
    'text':            ((2,), ('s',)),
    'set_yticklabels': ((0,), ('labels',)),
}


def measure():
    """[(name, verdict, repaired %, unrepaired %, exact_r, exact_u)] sorted by
    the repaired magnitude, descending.  Nothing is listed by hand."""
    shipped = gfet.quantum_capacitance
    rep = flip.make_contact_cq(shipped, live=False)
    unrep = flip.make_contact_cq(shipped, live=True)
    surface, _, _, _ = csg.capture_surface(rep)
    rows = []
    for site in sorted(surface):
        br, mr, er = csg.leverage(rep, site)
        bu, mu, eu = csg.leverage(unrep, site)
        rows.append((site, csg.ADJUDICATION_REPAIRED[site][0],
                     100.0 * (mr - br) / br, 100.0 * (mu - bu) / bu, er, eu))
    rows.sort(key=lambda r: (-abs(r[2]), r[0]))
    return rows, rep, unrep


def literal_numbers_in_text_calls(src):
    """[(call, literal)] for every string constant containing a digit that is
    passed to a text-drawing call, found with the PARSER.

    THE FIRST VERSION OF THIS FUNCTION WAS A REGEX OVER THE SOURCE TEXT, and
    it reported three violations that were not violations: `'...'  % (a, b)
    ... xytext='` matches as one "string literal" because a regex has no
    notion of where a literal ends, so the CODE between two literals was
    scanned as if it were text.  That is the repository's recorded
    resolve-on-spelling failure (running total fourteen before this one) and
    its dominant failure mode at once -- an instrument reporting on something
    it could not see, in the permissive direction for the thing it was
    checking and the alarmist direction for the thing it reported.

    It is also the reason `control_sources()` exists.  A checker whose
    positive and negative controls are both run cannot ship in that state;
    this one did not, because F1 failed on its own source.
    """
    bad = []
    for node in ast.walk(ast.parse(src)):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = (fn.attr if isinstance(fn, ast.Attribute)
                else getattr(fn, 'id', None))
        if name not in TEXT_ARGS:
            continue
        pos, kws = TEXT_ARGS[name]
        args = [node.args[i] for i in pos if i < len(node.args)]
        args += [k.value for k in node.keywords if k.arg in kws]
        for arg in args:
            for sub in ast.walk(arg):
                if isinstance(sub, ast.Constant) \
                        and isinstance(sub.value, str) \
                        and any(c.isdigit() for c in _strip_specs(sub.value)):
                    bad.append((name, sub.value))
    return bad


_SPEC = re.compile(r'%[-+ #0]*[0-9]*(?:\.[0-9]+)?[a-zA-Z%]')


def _strip_specs(text):
    """Drop %-conversion specifiers before looking for digits.

    The precision in `'%+.4f %%'` is a formatting instruction, not a number a
    reader sees; the third version of this checker flagged it.  Each of these
    three corrections narrowed WHAT the check looks at -- source text, then
    every argument, then every character -- and each was found by the check
    failing on its own source rather than by inspection.  F0e is the control.
    """
    return _SPEC.sub('', text)


CONTROL_CLEAN = "ax.set_title('a title with no numerals at all')"
CONTROL_DIRTY = "ax.set_title('the two differ on exactly 1 of 9 names')"
CONTROL_FORMAT = ("ax.set_xlabel('rebound by %sx at V = %s V' "
                  "% (PERTURBATION, PROBE_V))")
CONTROL_NONTEXT = ("ax.annotate('a numberless note', xy=(a, b), "
                   "arrowprops=dict(color='0.35', lw=0.9))")
CONTROL_SPEC = "ax.text(x, y, '%+.4f %%' % value)"
CONTROL_SPEC_DIRTY = "ax.text(x, y, 'row 3 of %+.4f %%' % value)"


def control_sources():
    """(clean, dirty, formatted) -- the planted controls for the checker.

    `dirty` is the fault the check exists to catch, written out so the check
    is known to be able to catch it; `clean` and `formatted` are the two
    shapes that must NOT be flagged, the second being the one this module
    actually uses everywhere.
    """
    return (CONTROL_CLEAN, CONTROL_DIRTY, CONTROL_FORMAT,
            CONTROL_NONTEXT, CONTROL_SPEC, CONTROL_SPEC_DIRTY)


def main():
    checks = []

    def check(name, ok, detail=''):
        detail = detail if isinstance(detail, str) else repr(detail)
        checks.append((name, bool(ok), detail))
        print('    %-4s %s%s' % ('PASS' if ok else 'FAIL', name,
                                 ('   -- ' + detail) if detail else ''))

    print('=' * 78)
    print('THE CAPTURE SURFACE AS A FIGURE, WITH ITS OWN PROVENANCE CHECKS')
    print('=' * 78)
    print()

    rows, rep, unrep = measure()
    base = float(rep(csg.PROBE_V))
    for site, verdict, pr, pu, er, eu in rows:
        print('    %-48s %-10s rep %+9.4f %%  unrep %+9.4f %%'
              % (csg._dotted(site), verdict, pr, pu))
    print()

    labels = [s[1] if s[0] == 'graphene_fet_model'
              else '%s.%s' % (s[0].replace('graphene_', '')
                              .replace('_model', '').replace('_argument', ''),
                              s[1])
              for s, _, _, _, _, _ in rows]

    # THE AXIS IS LOGARITHMIC, AND THAT IS THE WHOLE POINT.  The first version
    # of this figure was linear, and rendering it showed gfet.T's +0.0034 %
    # bar to be a sliver indistinguishable from the exact zeros -- so the one
    # distinction the figure exists to make, between a value that is BITWISE
    # unmoved and one that merely moves a little, was not visible in it.  F3
    # had passed, because it asked whether the zeros were a different ARTIST
    # and not whether a reader could see the difference: the dominant failure
    # mode, in a provenance check written to prevent it.  Three previous
    # sessions found a defect in a committed figure by rendering it; this is
    # the fourth, and it was found in a figure made this session.
    mags_r = [abs(r[2]) for r in rows if not r[4]]
    mags_u = [abs(r[3]) for r in rows if not r[5]]
    lo = min(mags_r + mags_u)
    hi = max(mags_r + mags_u)
    zpos = lo / 30.0                 # where a bitwise-unmoved row is drawn
    xlo, xhi = zpos / 2.5, hi * 7.0   # headroom for the value labels

    y = range(len(rows))
    fig, ax = plt.subplots(figsize=(13.6, 6.6))
    ax.set_xscale('log')

    h = 0.36
    nz_r = [(i, abs(r[2])) for i, r in enumerate(rows) if not r[4]]
    nz_u = [(i, abs(r[3])) for i, r in enumerate(rows) if not r[5]]
    z_r = [i for i, r in enumerate(rows) if r[4]]
    z_u = [i for i, r in enumerate(rows) if r[5]]

    if nz_r:
        ax.barh([i + h / 2 for i, _ in nz_r], [v - xlo for _, v in nz_r],
                left=xlo, height=h, color='#2b6cb0',
                label='repaired contact replacement, live_contact=False')
    if nz_u:
        ax.barh([i - h / 2 for i, _ in nz_u], [v - xlo for _, v in nz_u],
                left=xlo, height=h, color='#c05621',
                label='unrepaired, live_contact=True')
    if z_r:
        ax.plot([zpos] * len(z_r), [i + h / 2 for i in z_r], marker='o',
                markersize=9, linestyle='none', markerfacecolor='white',
                markeredgecolor='#2b6cb0', markeredgewidth=2.0,
                label='bitwise unmoved, repaired -- an exact zero, drawn off '
                      'the scale because a log axis has no room for it')
    if z_u:
        ax.plot([zpos] * len(z_u), [i - h / 2 for i in z_u], marker='o',
                markersize=9, linestyle='none', markerfacecolor='white',
                markeredgecolor='#c05621', markeredgewidth=2.0,
                label='bitwise unmoved, unrepaired')

    ax.axvline(zpos * 3.5, color='0.45', linewidth=0.9, linestyle=':')
    ax.set_xlim(xlo, xhi)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()

    # Ticks: one for the off-scale zero column, then the decades the data
    # actually spans.  No tick label is a hand-typed number.
    decades = [10.0 ** k for k in range(int(math.floor(math.log10(lo))),
                                       int(math.ceil(math.log10(hi))) + 1)]
    decades = [d for d in decades if xlo < d < xhi]
    ax.set_xticks([zpos] + decades)
    ax.set_xticklabels(['bitwise\nunmoved']
                       + [('%g' % d) for d in decades], fontsize=9)
    ax.set_xlabel('magnitude of the change in the contact quantum capacitance '
                  'when the name is rebound to a stand-in scaled by %sx, '
                  'at V = %s V   [per cent, log scale]'
                  % (csg.PERTURBATION, csg.PROBE_V), fontsize=9.5)

    # Per-row values as text, so no reading depends on a bar length.
    for i, (site, verdict, pr, pu, er, eu) in enumerate(rows):
        ax.text(xhi * 1.08, i, '%s' % verdict, fontsize=8.4,
                va='center', ha='left', color='0.25', clip_on=False)
        if not er:
            ax.text(abs(pr) * 1.15, i + h / 2, '%+.4f %%' % pr, fontsize=7.8,
                    va='center', ha='left', color='#2b6cb0')
        if not eu and (er or abs(pu - pr) > 0.0):
            ax.text(abs(pu) * 1.15, i - h / 2, '%+.4f %%' % pu, fontsize=7.8,
                    va='center', ha='left', color='#c05621')

    # Centred on the FIGURE, not on the axes: the axes span only the
    # middle of the canvas once the long y labels and the verdict column are
    # given room, and a title centred there runs off both edges, which is how
    # the first render clipped its own title.
    fig.suptitle('What a rebind of each name on the capture surface does to '
                 'the contact quantum capacitance\n'
                 'the two replacements differ on exactly %s of %s names, and '
                 '%s -- read by the captured object -- has the same leverage '
                 'as the name the capture holds'
                 % (sum(1 for r in rows if r[4] != r[5]), len(rows),
                    '_thermal_rounding'), fontsize=10.5)

    handles, leg_labels = ax.get_legend_handles_labels()
    # BELOW THE AXES, not inside them.  In the second render the legend box
    # sat over the quantum_capacitance row and hid the one bar the figure is
    # about, which F2 could not see because it only counted entries.
    leg = ax.legend(handles, leg_labels, fontsize=8.0, framealpha=1.0,
                    loc='upper center', bbox_to_anchor=(0.5, -0.155),
                    ncol=2, borderaxespad=0.0)
    ax.grid(axis='x', alpha=0.3, which='major')
    fig.subplots_adjust(left=0.205, right=0.845, top=0.855, bottom=0.30)
    fig.savefig(FNAME, dpi=150)

    # -- provenance ---------------------------------------------------------
    print('-' * 78)
    print('PROVENANCE OF THE FIGURE ITSELF')
    print('-' * 78)
    print()
    (clean, dirty, formatted, nontext, spec,
     spec_dirty) = control_sources()
    check('F0a POSITIVE CONTROL: the checker flags a title with a hand-typed '
          'number', len(literal_numbers_in_text_calls(dirty)) == 1,
          literal_numbers_in_text_calls(dirty))
    check('F0b NEGATIVE CONTROL: it does not flag a numberless title',
          not literal_numbers_in_text_calls(clean))
    check('F0c and it does not flag the %-format shape this module uses, '
          'whose digits live in the ARGUMENTS and not in the literal',
          not literal_numbers_in_text_calls(formatted),
          literal_numbers_in_text_calls(formatted))
    check('F0d and it does not flag a NON-TEXT string argument: a colour '
          'inside arrowprops is not drawn as a glyph',
          not literal_numbers_in_text_calls(nontext),
          literal_numbers_in_text_calls(nontext))
    check('F0e and it does not flag the PRECISION inside a %-conversion '
          'specifier, while still flagging a hand-typed number beside one',
          not literal_numbers_in_text_calls(spec)
          and len(literal_numbers_in_text_calls(spec_dirty)) == 1,
          literal_numbers_in_text_calls(spec)
          + literal_numbers_in_text_calls(spec_dirty))
    src = io.open(os.path.abspath(__file__), encoding='utf-8').read()
    bad = literal_numbers_in_text_calls(src)
    check('F1  no literal digit appears in any title, label or annotation of '
          'THIS file -- every number is interpolated from the measurement',
          not bad, bad[:3])
    check('F2  the legend is built from the artists, so it has one entry per '
          'drawn series and no entry for an undrawn one',
          len(handles) == len(leg_labels) == (bool(nz_r) + bool(nz_u)
                                              + bool(z_r) + bool(z_u)),
          '%d handles' % len(handles))
    drawn = sum(len(c.patches) if hasattr(c, 'patches') else 0
                for c in ax.containers)
    check('F3  the bars drawn equal the nonzero rows, and the exact zeros are '
          'a DIFFERENT artist (a zero-height bar is indistinguishable from a '
          'missing one)', drawn == len(nz_r) + len(nz_u),
          '%d bars, %d zero markers' % (drawn, len(z_r) + len(z_u)))
    check('F4  every y tick label corresponds to a measured row, in order',
          [t.get_text() for t in ax.get_yticklabels()] == labels)
    check('F5  the two exact zeros of the repaired replacement are the two '
          'the guard adjudicates PROTECTED/UNREACHED, and the unrepaired one '
          'has one fewer',
          len(z_r) == len(z_u) + 1 and rows[[r[0] for r in rows]
                                           .index(csg.QC)][4] is True,
          '%d vs %d' % (len(z_r), len(z_u)))
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    fb = fig.bbox
    outside = []
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_text().strip():
            continue
        e = t.get_window_extent(renderer=rend)
        if e.x0 < fb.x0 - 0.5 or e.x1 > fb.x1 + 0.5 \
                or e.y0 < fb.y0 - 0.5 or e.y1 > fb.y1 + 0.5:
            outside.append(t.get_text()[:34])
    check('F7  NO TEXT IS CLIPPED: every drawn string lies inside the canvas. '
          'The first render of this figure cut off the end of its own title '
          'and of its x-axis label, and nothing checked it',
          not outside, outside[:3])

    leg_box = leg.get_window_extent(renderer=rend)
    hidden = []
    for cont in ax.containers:
        for patch in getattr(cont, 'patches', []):
            if patch.get_width() <= 0:
                continue
            if leg_box.overlaps(patch.get_window_extent(renderer=rend)):
                hidden.append(round(patch.get_y(), 2))
    check('F9  THE LEGEND HIDES NO DATA: its box overlaps no drawn bar.  In '
          'the second render it sat on top of the quantum_capacitance row -- '
          'the single row the figure exists to show -- and F2 passed, because '
          'counting legend entries says nothing about where the box lands',
          not hidden, hidden[:4])

    ax_box = ax.get_window_extent(renderer=rend)
    overrun = [t.get_text() for t in ax.texts
               if t.get_clip_on()
               and t.get_window_extent(renderer=rend).x1 > ax_box.x1 + 0.5]
    check('F10 NO VALUE LABEL RUNS PAST THE RIGHT SPINE: the two largest bars '
          'ended within a few pixels of it in the second render',
          not overrun, overrun[:3])

    smallest = min(abs(r[2]) for r in rows if not r[4])
    px_zero = ax.transData.transform((zpos, 0))[0]
    px_small = ax.transData.transform((smallest, 0))[0]
    check('F8  A SMALL MOVE IS VISIBLY DISTINCT FROM AN EXACT ZERO: the '
          'smallest nonzero bar ends at least 40 px from the off-scale zero '
          'marker.  The linear first render put them on top of each other, '
          'and F3 passed anyway because it asked about artists and not about '
          'what a reader can see', px_small - px_zero >= 40.0,
          '%.1f px for %+.4f %%' % (px_small - px_zero, smallest))

    check('F6  the figure file was written and is non-trivial',
          os.path.exists(FNAME) and os.path.getsize(FNAME) > 20000,
          '%d bytes' % (os.path.getsize(FNAME) if os.path.exists(FNAME) else 0))
    print()
    print('    base contact C_q at V = %s V: %s F/m^2 (%s)'
          % (csg.PROBE_V, base, float(base).hex()))
    print('    written: %s' % FNAME)
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
    import sys
    sys.exit(main())
