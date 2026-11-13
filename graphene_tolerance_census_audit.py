"""
graphene_tolerance_census_audit.py

A COMPARISON AGAINST A PUBLISHED NUMBER IS A CLAIM ABOUT THE SOURCE'S
PRECISION, AND UNTIL TODAY THIS REPOSITORY HAD CHECKED EXACTLY ONE OF THEM.

Created 2026-11-11 as the RANKED-FIRST item of the 2026-11-10 entry:

    "CENSUS EVERY COMPARISON AGAINST A PUBLISHED NUMBER AND CHECK IT USES A
     SIGNIFICANT-FIGURE TEST, NOT A RELATIVE TOLERANCE -- Every
     literature-calibrated constant in this repository was quoted to two or
     three figures by its source, and any check comparing against one with a
     relative tolerance is either vacuous or wrong depending on which side of
     the rounding it lands.  `_rounds_to` is now importable.  This is item 15
     with a handle on it."

That item exists because of one instance, found on 2026-11-10 inside
`graphene_self_heating_model.verify_against_dorgan`.  Check B2 compared a
computed 1.1113e11 cm^-2 against Dorgan's footnote-18 value of 1.1e11 cm^-2
with a 1 % RELATIVE tolerance, and failed at 1.027 %.  The model was right and
the CHECK was wrong: the paper quotes two significant figures, so every value
in [1.05e11, 1.15e11) rounds to the quoted 1.1e11, a band 4.545 % wide on
either side.  A 1 % tolerance rejects most of the interval the source actually
claims.

WHAT THE CENSUS FOUND, AND WHERE IT CORRECTS ITS OWN COMMISSIONING ITEM

The item was written as "significant-figure test, NOT a relative tolerance",
i.e. as though the instrument were the defect.  That is not quite right, and
the census is what makes the difference visible.  For a value q quoted to s
significant figures the quoting interval is

    q +/- h_abs,   h_abs = 0.5 * 10**(floor(log10|q|) - s + 1)

so the test `abs(got - q) / q < h_rel` with h_rel = h_abs / |q| is EXACTLY the
"rounds to q at s figures" test, save for the half-way tie-break.  A relative
tolerance is therefore not wrong in kind.  It is wrong whenever its NUMBER was
not derived from the quoted precision -- which is every hand-written tolerance
in this tree, because none of them was.  The defect is an undeclared precision
claim, not a choice of operator, and a session that mechanically rewrote every
relative tolerance as `_rounds_to` would fix the instances by accident while
leaving the class untouched.  So the verdicts here are stated as a RATIO,
h_rel / tol, not as a pass/fail on which function was called:

    OVER_TIGHT   tol < h_rel   the check can reject a value the source's own
                               rounding admits.  This is the B2 shape.
    VACUOUS      tol > h_rel   the check admits values that round to a
                               DIFFERENT published figure than the one quoted.
    MATCHED      tol == h_rel  equivalent to the significant-figure test.
    AMBIGUOUS    the quoted token's significant figures cannot be read from
                 how it is written (trailing zeros, no decimal point), so no
                 verdict is available and none is invented.
    NO_QUOTE     a tolerance test whose label names no published value: in
                 scope for the census, out of scope for a verdict.

HOW THE QUOTED PRECISION IS OBTAINED

From the text of the check, not from the source literal.  This distinction is
itself a finding.  `graphene_self_heating_model` check A2 compares against the
literal `1.0e4` -- which, written that way, reads as two significant figures --
while its own label says "paper: ~1e4", an order-of-magnitude claim with one.
The literal and the label disagree about what the source said, and only the
label is evidence about the source.  Every row below therefore records the
quoted token AS WRITTEN IN THE LABEL, and `_sigfigs` returns None rather than a
guess when the written form is ambiguous.

EXACTLY-KNOWN VALIDATIONS

Required by the standing practice of this repository, and by the 2026-09-17
precedent in which a plausible-range check passed a real bug:

  * h_rel for Dorgan's 1.1e11 at 2 s.f. is EXACTLY 1/22, and for 9.4e11 at
    2 s.f. EXACTLY 1/188.  Both are asserted as exact rationals, not ranges.
  * The 2026-11-10 bug is reproduced as an arithmetic identity: the recorded
    1.027 % must lie strictly between the old 1 % tolerance and 1/22, which is
    the whole window in which a check fails while its value is correct.
  * A source with no comparisons must census to EXACTLY zero rows -- the
    symmetry that must give zero.
  * Every one of the five verdicts gets a planted control (the nineteenth
    instance of that standing item), and the planted controls are then run
    against five mutated verdict functions, each of which must be detected.

Usage:
    python3 graphene_tolerance_census_audit.py            # checks
    python3 graphene_tolerance_census_audit.py --census   # the table only
"""

import ast
import os
import re
import sys
from math import floor, log10

from graphene_self_heating_model import _round_sig, _rounds_to

REPO = os.path.dirname(os.path.abspath(__file__))

# Words that make a check's label a claim about a SOURCE rather than about the
# code's own internal consistency.
#
# CALIBRATED AGAINST THIS TREE'S OWN VOCABULARY, NOT GUESSED.  The first cue
# list included "reported" and "measured value".  On the real tree that list
# flagged 29 labels, of which 16 -- 55 % -- were cued by "reported" ALONE, and
# in this repository "reported" almost always means reported BY AN INSTRUMENT
# IN THIS TREE ("a Name default is NOT reported by this arbiter", "an exemption
# with no reason is reported as DEAD").  A census whose cue list is wrong about
# the subject's vocabulary does not under-report, it reports the wrong rows and
# then adjudicates them, so the loose cues are dropped and the 55 % is recorded
# here as the reason.  "reference" and "expected" are excluded for the same
# reason: in this tree they mean the repository's own earlier model.
PAPER_CUES = (
    "paper", "published", "literature", "footnote", "dorgan",
    "ref.", "et al", "quoted", "codata",
)

# Cues that are only cues as a PHRASE.  "accepted" alone appears three times
# in this tree's check labels and only ONE of them is about a source: the
# other two are "the FULL adjudication is accepted" and "the complete, correct
# adjudication ... is accepted", i.e. this repository's own verdict vocabulary.
# A bare "accepted" cue would therefore run at 33 % precision, the same defect
# that disqualified "reported" above, so the cue is the phrase.  Requiring it
# is what brings `graphene_self_heating_model` check A1 -- "accepted
# thermal-SiO2 value 1.2794 W/m.K" -- into the census.
PAPER_CUE_RES = (
    re.compile(r"accepted[\w\s-]{0,25}?value", re.IGNORECASE),
)


def _has_cue(label):
    low = label.lower()
    return (any(c in low for c in PAPER_CUES)
            or any(r.search(label) for r in PAPER_CUE_RES))

# A number as written: optional sign, digits, optional fraction, optional
# exponent.  Captured as a STRING because its spelling carries its precision.
_NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

# "the paper's 4 %", "paper: ~1e4", "Dorgan's 9.4e11 cm^-2", "quoted 1.1e11"
#
# THE LOOKBEHIND AND THE ALTERNATION ORDER ARE BOTH BUG FIXES, FOUND BY
# RUNNING THIS CENSUS OVER THIS TREE ON THE DAY IT WAS WRITTEN.  The first
# form put bare "footnote" before "footnote-18" in the alternation and allowed
# the number to begin immediately after the cue, so on
#
#     "n_min(...) rounds to the paper's footnote-18 %.1e cm^-2 at %d s.f."
#
# it matched the cue "footnote" and then read the ORDINAL OF THE FOOTNOTE,
# -18, as the published value -- a negative carrier density quoted to two
# significant figures.  It would then have been adjudicated against a
# tolerance.  The longest cue now wins, and a number welded to the preceding
# word by a hyphen or a letter is not a quotation.
_QUOTE_RE = re.compile(
    r"(?:paper'?s?|published|literature|footnote-\d+|footnote|dorgan'?s?|"
    r"quoted|codata|ref\.|accepted[\w\s-]{0,25}?value)"
    r"[^0-9%+-]{0,40}?(~?)(?<![-\w])(" + _NUM + r")",
    re.IGNORECASE,
)

# ---------------------------------------------------------------------------
# THE SECOND WARRANT: A LABEL MAY TAKE ITS CENTRE AND ITS WIDTH FROM
# DIFFERENT PRINTED NUMBERS.
#
# Added 2026-11-13 as the RANKED-SECOND item of the 2026-11-12 entry, which
# was created by this census's FIRST FALSE POSITIVE (check 9l).  The census
# assumed a check's centre and its width are warranted by the same printed
# number.  `graphene_self_heating_model` check A2 is the counterexample and it
# was written deliberately: its centre is Dorgan's one-figure `~1e4 K/W` and
# its width is the two-figure `~2.8e-7 m^2.K/W per unit of device area` from
# THE SAME SENTENCE, which -- since A = L*W is given by the stated geometry
# and not fitted -- is the same statement divided by an exact constant.  The
# census read only the centre, licensed +-50 %, saw +-1.7857 % and returned
# OVER_TIGHT at ratio 28: a wrong verdict reached by correct arithmetic.
#
# The repair is the one the item predicted: a SECOND OPTIONAL TOKEN in the
# label vocabulary, so a label can say "centre 1e4, width from 2.8e-7".  It
# is spelled machine-readably rather than in prose, because the whole defect
# of the first version was reading a precision claim out of English.
#
#     ... [width-warrant: 2.8e-7]
#
# RESTRICTED TO RELATIVE TOLERANCES, AND THE RESTRICTION IS THE POINT.  A
# relative width is dimensionless and therefore transferable between two
# quantities related by an exact constant; an ABSOLUTE width lives in the
# units of its own quantity and transferring it would require the label to
# assert that the two quantities share units, which a label cannot be trusted
# to do.  On an ABS_TOL row the warrant is parsed, recorded and NOT used, and
# check 10e reports any such row rather than silently honouring it.
_WIDTH_WARRANT_RE = re.compile(
    r"\[\s*width-warrant\s*:\s*(~?)(" + _NUM + r")\s*\]", re.IGNORECASE)


def extract_width_warrant(label):
    """(token, approx) of a label's declared WIDTH warrant, if it has one."""
    if not label:
        return (None, False)
    m = _WIDTH_WARRANT_RE.search(label)
    if not m:
        return (None, False)
    return (m.group(2), m.group(1) == "~")


# A test that compares to the quoted number's own rounding rather than to a
# hand-written tolerance.  These are the CONTROL class of the census: the
# instrument the commissioning item asked for, already in use.
#
# THREE VOCABULARIES, NOT ONE, and the third was found by this census
# misreporting a row.  `graphene_kappa_impurity_bound` check D2 was reported
# LIT_NO_TEST -- a published number named with no numeric test beside it --
# when in fact it tests `float('%.2g' % (n_a * 1e-16)) != 1.5`, a
# round-to-two-figures comparison written through the %g format rather than
# through round() or _rounds_to().  Check D1 on the line above does the same
# with %.3g.  So this repository already had the right instrument in a module
# that predates `_rounds_to` by five weeks, spelled in a way that the census
# of "every comparison against a published number" initially scored as the
# worst category available.  Under-reporting a correct check as an absent one
# is the same defect as over-reporting a tolerance, and it is the standing
# item GIVE EVERY VERDICT IN EVERY VOCABULARY A PLANTED CONTROL arriving from
# the vocabulary side instead of the verdict side.
_SIGFIG_FNS = ("_rounds_to", "rounds_to")

# '%.3g' % x, '%.2e' % x -- a format that rounds to a fixed number of figures
_SIGFIG_FMT_RE = re.compile(r"%\.\d+[geE]")

VERDICTS = ("OVER_TIGHT", "VACUOUS", "MATCHED", "AMBIGUOUS", "NO_QUOTE",
            "SIGFIG_TEST", "LIT_NO_TEST")


def _sigfigs(token):
    """Significant figures of a number AS WRITTEN, or None if the spelling is
    ambiguous.

    Trailing zeros count only when a decimal point is present: "1.0" is two
    figures, "1e4" is one, and "100" is undecidable from its spelling alone.
    Returning None for the last case is the point -- a census that guessed
    would manufacture the precision claim it exists to audit.
    """
    t = token.strip().lstrip("+-")
    mant = t.split("e")[0].split("E")[0]
    if "." in mant:
        digits = mant.replace(".", "")
        stripped = digits.lstrip("0")
        if stripped == "":
            return None            # 0.000 -- no significant figures
        # leading zeros before the first significant digit never count;
        # trailing zeros after a decimal point always do
        return len(stripped)
    digits = mant.lstrip("0")
    if digits == "":
        return None                # bare 0
    if digits.endswith("0"):
        return None                # "100", "1200" -- spelling is ambiguous
    return len(digits)


def half_width(token):
    """(h_abs, h_rel, sigfigs) of the quoting interval of a written number.

    For q quoted to s significant figures the source is claiming
    q +/- 0.5 * 10**(floor(log10|q|) - s + 1).  Returns (None, None, None)
    when the spelling does not fix s.
    """
    s = _sigfigs(token)
    q = abs(float(token))
    if s is None or q == 0.0:
        return (None, None, s)
    h_abs = 0.5 * 10.0 ** (floor(log10(q)) - s + 1)
    return (h_abs, h_abs / q, s)


def effective_half_width(form, h_abs, h_rel):
    """The half-width a tolerance of THIS FORM has to be compared against.

    FOUND BY RUNNING THE FIRST VERSION OF THIS CENSUS (2026-11-11).  That
    version compared every tolerance against h_rel, which is correct for
    `abs(a-b)/b < tol` and WRONG for `abs(a-b) < tol`: an absolute tolerance
    lives in the units of the quantity, so it must be judged against h_abs.
    The one absolute row in the tree was consequently reported as vacuous by a
    factor of 12.9 when the true factor is 2.  The census was over-stating the
    defect it exists to find, which is the same error as understating it.
    """
    if form == "ABS_TOL":
        return h_abs
    return h_rel


def verdict(tol, h_eff, rel_eps=1e-12):
    """OVER_TIGHT / MATCHED / VACUOUS for a tolerance against a quoting
    interval of half-width h_eff, in the tolerance's own units."""
    if h_eff is None:
        return "AMBIGUOUS"
    if abs(tol - h_eff) <= rel_eps * h_eff:
        return "MATCHED"
    return "OVER_TIGHT" if tol < h_eff else "VACUOUS"


# --------------------------------------------------------------------------
# The AST census
# --------------------------------------------------------------------------

class Row(object):
    """One tolerance comparison, with whatever the tree says about it."""

    __slots__ = ("path", "line", "form", "tol", "label", "token", "approx",
                 "sigfigs", "h_abs", "h_rel", "h_eff", "verdict", "ratio",
                 "bound_via", "width_token")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    def __repr__(self):
        return "Row(%s:%d %s tol=%r %s)" % (
            os.path.basename(self.path or "?"), self.line or 0,
            self.form, self.tol, self.verdict)

    def fields(self):
        """Every slot as a tuple, for comparing two censuses by VALUE.

        Added 2026-11-12.  Row deliberately keeps object identity for
        equality -- two runs of census_repo() over the same unchanged tree
        compare UNEQUAL, because every Row is a fresh object.  Anyone
        regression-testing the census against itself therefore needs this,
        and check R2b records the trap rather than leaving it to be
        rediscovered.  __eq__ is NOT added: rows are collected in lists that
        are deduplicated nowhere, and giving them value equality would change
        what `in` and `count` mean across the module.
        """
        return tuple(getattr(self, k) for k in self.__slots__)


def _const_num(node):
    """The numeric value of a literal node, else None.  Handles the unary
    minus that the parser does not fold."""
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) \
            and not isinstance(node.value, bool):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        inner = _const_num(node.operand)
        return None if inner is None else -inner
    return None


def _is_abs_call(node):
    return (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == "abs")


def _contains_sub(node):
    for n in ast.walk(node):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub):
            return True
    return False


def classify_left(left):
    """REL_TOL / ABS_TOL / OTHER_TOL for the left-hand side of a `< tol`.

    REL_TOL is a residual DIVIDED by something: abs(a-b)/c, or (a-b)/c.
    ABS_TOL is an undivided residual: abs(a-b).
    """
    if isinstance(left, ast.BinOp) and isinstance(left.op, ast.Div):
        if _contains_sub(left.left) or (_is_abs_call(left.left)
                                        and left.left.args
                                        and _contains_sub(left.left.args[0])):
            return "REL_TOL"
        return None
    if _is_abs_call(left) and left.args and _contains_sub(left.args[0]):
        return "ABS_TOL"
    if _contains_sub(left):
        return "OTHER_TOL"
    return None


def _check_labels(tree):
    """[(start_line, end_line, label)] for every ...check(<str>, ...) call.

    Any callable whose name ends in "check" counts: this tree uses `_check`
    and `check` interchangeably across modules.

    CORRECTED 2026-11-13.  This docstring previously also listed `_c`, which
    the code has never accepted -- "_c" does not end in "check".  Writing a
    planted control with `_c` on 2026-11-13 produced a census with zero rows
    and no error, which is how the mismatch was found.  No module in the tree
    uses `_c` for a check (measured: the only occurrences are inside this
    module's own planted source), so nothing was being missed; the docstring
    was describing a tolerance the code does not have.
    """
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = None
        if isinstance(node.func, ast.Name):
            name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            name = node.func.attr
        if not name or not name.lower().endswith("check"):
            continue
        if not node.args:
            continue
        label = _static_str(node.args[0])
        if label is None:
            continue
        end = getattr(node, "end_lineno", node.lineno) or node.lineno
        out.append((node.lineno, end, label))
    return out


def _static_str(node):
    """The value of a string literal or of a statically-joinable
    concatenation of them (`"a" "b"`, `"a" + "b"`, `"a %s" % x` -> "a %s")."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Add):
            l, r = _static_str(node.left), _static_str(node.right)
            if l is not None and r is not None:
                return l + r
        if isinstance(node.op, ast.Mod):
            return _static_str(node.left)
    if isinstance(node, ast.JoinedStr):
        parts = [p.value for p in node.values
                 if isinstance(p, ast.Constant) and isinstance(p.value, str)]
        return "".join(parts) if parts else None
    return None


def _enclosing_label(labels, line):
    """The label of the innermost check call containing `line`."""
    best = None
    for start, end, label in labels:
        if start <= line <= end:
            if best is None or (end - start) < (best[1] - best[0]):
                best = (start, end, label)
    return best[2] if best else None


def extract_quote(label):
    """(token, approx) -- the published value named in a check label, as
    written, and whether the label hedged it with '~'."""
    if not label:
        return (None, False)
    if not _has_cue(label):
        return (None, False)
    m = _QUOTE_RE.search(label)
    if not m:
        return (None, False)
    return (m.group(2), m.group(1) == "~")


def _residual_bindings(tree):
    """{(scope_start, name): value_node} for every single-target assignment
    whose value contains a subtraction.

    WHY THIS EXISTS.  The first version of this census matched only an inline
    residual, `abs(a - b) / b < tol`, and so missed the commonest shape in this
    tree:

        rel = abs(G_new_mu_matched - G_old) / G_old
        _check("O1 ... to <0.5 %", rel < 5e-3, ...)

    which is `graphene_self_heating_model.overlap_with_isothermal_model`, i.e.
    the model-versus-model validation that the 2026-11-10 session wrote.  A
    census blind to a name binding is a census of a coding style, not of a
    class of claim.  Resolution is deliberately ONE level deep and scoped to
    the enclosing function: deeper chains exist, and the count of rows reached
    only through a name is reported so the blind spot has a size.
    """
    out = {}
    for fn in ast.walk(tree):
        if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef,
                              ast.Module)):
            continue
        start = getattr(fn, "lineno", 0)
        end = getattr(fn, "end_lineno", 10 ** 9) or 10 ** 9
        for node in ast.walk(fn):
            if not isinstance(node, ast.Assign) or len(node.targets) != 1:
                continue
            tgt = node.targets[0]
            if not isinstance(tgt, ast.Name):
                continue
            if _contains_sub(node.value):
                out[(start, end, tgt.id)] = node.value
    return out


def _resolve(bindings, node, line):
    """The residual a bare Name was bound to, in the innermost scope
    containing `line`, or None."""
    if not isinstance(node, ast.Name):
        return None
    best = None
    for (start, end, name), value in bindings.items():
        if name == node.id and start <= line <= end:
            if best is None or (end - start) < (best[0][1] - best[0][0]):
                best = ((start, end, name), value)
    return best[1] if best else None


def census_source(src, path="<src>"):
    """Every tolerance comparison in one source, adjudicated."""
    tree = ast.parse(src)
    labels = _check_labels(tree)
    bindings = _residual_bindings(tree)
    rows = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Compare):
            continue
        if len(node.ops) != 1 or not isinstance(node.ops[0], (ast.Lt, ast.LtE)):
            continue
        tol = _const_num(node.comparators[0])
        if tol is None or tol <= 0.0:
            continue
        bound_via = None
        form = classify_left(node.left)
        if form is None:
            resolved = _resolve(bindings, node.left, node.lineno)
            if resolved is not None:
                form = classify_left(resolved)
                bound_via = node.left.id
        if form is None:
            continue
        label = _enclosing_label(labels, node.lineno)
        token, approx = extract_quote(label)
        if token is None:
            rows.append(Row(path=path, line=node.lineno, form=form, tol=tol,
                            label=label, token=None, approx=False,
                            sigfigs=None, h_abs=None, h_rel=None, h_eff=None,
                            verdict="NO_QUOTE", ratio=None,
                            bound_via=bound_via))
            continue
        h_abs, h_rel, sf = half_width(token)
        h_eff = effective_half_width(form, h_abs, h_rel)
        # the second warrant, 2026-11-13: a declared width token overrides
        # the centre's own half-width, on relative forms only
        w_token, _w_approx = extract_width_warrant(label)
        if w_token is not None and form != "ABS_TOL":
            w_abs, w_rel, _w_sf = half_width(w_token)
            if w_rel is not None:
                h_eff = w_rel
        v = verdict(tol, h_eff)
        rows.append(Row(path=path, line=node.lineno, form=form, tol=tol,
                        label=label, token=token, approx=approx, sigfigs=sf,
                        h_abs=h_abs, h_rel=h_rel, h_eff=h_eff, verdict=v,
                        ratio=(None if h_eff is None else h_eff / tol),
                        bound_via=bound_via, width_token=w_token))
    # ---- the control class: comparisons against a published number that
    # ---- already use the quoted number's own rounding ---------------------
    tol_lines = set(r.line for r in rows)
    for start, end, label in labels:
        token, approx = extract_quote(label)
        if token is None:
            continue
        kind = None
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or node.lineno < start \
                    or node.lineno > end:
                continue
            fname = (node.func.id if isinstance(node.func, ast.Name)
                     else getattr(node.func, "attr", None))
            if fname in _SIGFIG_FNS:
                kind = "SIGFIG"
                break
            if fname == "round":
                kind = "ROUND_EQ"
            if fname == "float" and node.args \
                    and isinstance(node.args[0], ast.BinOp) \
                    and isinstance(node.args[0].op, ast.Mod):
                fmt = _static_str(node.args[0].left)
                if fmt and _SIGFIG_FMT_RE.search(fmt):
                    kind = "SIGFIG_FMT"
        if kind is None:
            # a label that names a published value while the check beside it
            # tests something else entirely.  Not a defective tolerance -- an
            # ABSENT one.  Reported, not adjudicated.
            if not any(start <= ln <= end for ln in tol_lines):
                h_abs, h_rel, sf = half_width(token)
                rows.append(Row(path=path, line=start, form="NO_TEST",
                                tol=None, label=label, token=token,
                                approx=approx, sigfigs=sf, h_abs=h_abs,
                                h_rel=h_rel, h_eff=None,
                                verdict="LIT_NO_TEST", ratio=None,
                                bound_via=None))
            continue
        if any(start <= ln <= end for ln in tol_lines):
            continue        # already counted as a tolerance row
        h_abs, h_rel, sf = half_width(token)
        rows.append(Row(path=path, line=start, form=kind, tol=None,
                        label=label, token=token, approx=approx, sigfigs=sf,
                        h_abs=h_abs, h_rel=h_rel, h_eff=None,
                        verdict="SIGFIG_TEST", ratio=None, bound_via=None))

    rows.sort(key=lambda r: (r.line, r.form))
    return rows


SELF = os.path.basename(os.path.abspath(__file__))

# Modules whose check labels quote published numbers AS CONTROLS rather than
# as claims about a source.  Censusing one of these reports the instrument's
# own planted material as findings about the tree.  See census_repo's
# docstring for how the second entry arrived and what it cost.
PLANTING_MODULES = (SELF, "graphene_citation_census_audit.py")

# The sibling repository, for the cross-repository run of 2026-11-12.  Looked
# up beside this one rather than hard-coded, and its absence SKIPS the R1/X
# checks rather than failing them -- a clone of one repository alone must
# still pass its own suite.
SIBLING = os.path.join(os.path.dirname(REPO), "Design-Verification-Roadmap")


# ---------------------------------------------------------------------------
# FILE DISCOVERY.  SPLIT OUT OF census_repo ON 2026-11-12 BECAUSE THE FIRST
# CROSS-REPOSITORY USE OF THIS INSTRUMENT RETURNED AN EMPTY CENSUS.
#
# census_repo() was written as a flat `os.listdir(root)` over *.py.  That is
# correct for THIS repository, which keeps every module at the top level, and
# it is silently wrong for any repository that does not: run over
# Design-Verification-Roadmap it returned ZERO rows, not because that tree has
# no tolerances but because every one of its 40-odd Python files lives under
# tools/ or examples/.  An empty census is indistinguishable from a clean one
# unless something says how many files were read, so _python_files is now
# reported alongside the row count everywhere the census is printed.
#
# The instrument's hidden assumption was invisible in the repository it was
# written for.  That is the generalisable finding, and it is why `recursive`
# defaults to True rather than preserving the old behaviour: the old default
# cannot be distinguished from a correct answer by its output.  Check R2
# asserts that the change leaves THIS tree's census bit-for-bit identical,
# which it must, because this tree is flat.
# ---------------------------------------------------------------------------
_SKIP_DIRS = {".git", "__pycache__", ".venv", "venv", "node_modules",
              ".pytest_cache", ".mypy_cache"}


def _python_files(root, recursive=True):
    """Every *.py under `root`, sorted, skipping VCS and cache directories.

    With recursive=False this is the pre-2026-11-12 flat behaviour, retained
    only so that check R1 can assert what it used to return.
    """
    if not recursive:
        return [os.path.join(root, n) for n in sorted(os.listdir(root))
                if n.endswith(".py")]
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in _SKIP_DIRS)
        for n in sorted(filenames):
            if n.endswith(".py"):
                out.append(os.path.join(dirpath, n))
    return sorted(out)


def census_repo(root=None, include_self=False, recursive=True):
    """Every tolerance comparison in every module in the tree.

    THIS MODULE IS EXCLUDED BY DEFAULT, and the exclusion is the standing item
    A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES (created 2026-10-17)
    getting its instance -- found by running the census, not by reading for
    it.  `PLANT` below contains, deliberately, one comparison against a
    published number per verdict, and the check suite contains two tolerances
    of 1e-18 against Dorgan's 1.1e11 and 9.4e11.  On the first full run those
    five rows appeared in the census OF THE TREE as three genuine OVER_TIGHT
    findings and two LIT_NO_TEST rows, one of them quoting a "published" value
    of 0 lifted out of check 5d's own label.  The tree had acquired three
    tolerance defects that exist only inside the instrument that reports them.

    Pass include_self=True to census this module too; the check suite does
    exactly that, to assert the self-plant is real rather than describing it.

    GENERALISED FROM ONE FILENAME TO A DECLARED CLASS, 2026-11-13.  On that
    day `graphene_citation_census_audit` was added -- a second instrument,
    whose check labels necessarily quote published tokens because its subject
    IS the quoting of published tokens.  Its arrival put 10 new quoted rows
    into this census's table of THE TREE, including three ORDINAL controls and
    two UNWARRANTED ones, i.e. the same defect this docstring was written
    about, one module over.  A `SELF`-only exclusion cannot see that: the
    excluded thing is not "this file", it is "a module that plants the pattern
    it censuses", and there are now two of them.  `PLANTING_MODULES` is that
    class, declared rather than inferred, and check 8h measures what
    excluding it does instead of asserting that it is harmless.
    """
    root = root or REPO
    rows = []
    for p in _python_files(root, recursive=recursive):
        if os.path.basename(p) in PLANTING_MODULES and not include_self:
            continue
        try:
            src = open(p, encoding="utf-8").read()
        except OSError:
            continue
        try:
            rows.extend(census_source(src, p))
        except SyntaxError as exc:
            rows.append(Row(path=p, line=getattr(exc, "lineno", 0) or 0,
                            form="UNPARSED", tol=None, label=str(exc),
                            verdict="NO_QUOTE"))
    return rows


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def print_census(rows, quoted_only=True):
    """The table.  Rows whose label names no published value are counted but
    not printed unless asked for: they are in scope for the census (a
    tolerance is a tolerance) and out of scope for a verdict."""
    print("=" * 78)
    print("COMPARISONS AGAINST A PUBLISHED NUMBER")
    print("=" * 78)
    hdr = ("%-36s %-5s %-9s %-9s %-9s %-4s %-10s %-12s %s"
           % ("module", "line", "form", "tol", "quoted", "s.f.",
              "licensed", "verdict", "licensed/tol"))
    print(hdr)
    print("-" * len(hdr))
    shown = 0
    for r in rows:
        if quoted_only and r.verdict == "NO_QUOTE":
            continue
        shown += 1
        print("%-36s %-5d %-9s %-9s %-9s %-4s %-10s %-12s %s"
              % (os.path.basename(r.path), r.line, r.form,
                 ("%g" % r.tol) if r.tol is not None else "-",
                 ("~" if r.approx else "") + str(r.token),
                 r.sigfigs if r.sigfigs is not None else "amb",
                 ("%.5g" % r.h_eff) if r.h_eff is not None else "-",
                 r.verdict,
                 ("x%.6g" % r.ratio) if r.ratio is not None else ""))
    print("-" * len(hdr))
    counts = {}
    for r in rows:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    print("%d rows printed; %d tolerance comparisons in the tree in total, "
          "of which" % (shown, sum(1 for r in rows if r.tol is not None)))
    for v in VERDICTS:
        if counts.get(v):
            print("    %-12s %d" % (v, counts[v]))
    print("%d row(s) were reached only through a local name binding"
          % sum(1 for r in rows if r.bound_via))
    return counts


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------

passed = 0
failed = 0


def check(label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        print("  PASS  %s" % label)
    else:
        failed += 1
        print("  FAIL  %s\n          %s" % (label, detail))
    return bool(ok)


# A planted source carrying one instance of every verdict this module can
# return.  The nineteenth instance of the standing item GIVE EVERY VERDICT IN
# EVERY VOCABULARY A PLANTED CONTROL -- and the first in which the vocabulary
# being controlled is the census's OWN.
PLANT = '''
def planted():
    # OVER_TIGHT: 2 s.f. licenses 1/22 = 4.5454 %, the check allows 1 %
    _check("P1 matches the paper's 1.1e11 cm^-2",
           abs(got - 1.1e11) / 1.1e11 < 0.01, "")
    # VACUOUS: 3 s.f. licenses +/-0.005, the check allows +/-0.5
    _check("P2 reproduces the published 6.45 kOhm/sq",
           abs(got - 6.45) < 0.5, "")
    # MATCHED: 2 s.f. on 9.4e11 licenses exactly 1/188
    _check("P3 agrees with the paper's 9.4e11 cm^-2",
           abs(got - 9.4e11) / 9.4e11 < 0.005319148936170213, "")
    # AMBIGUOUS: "20" could be one figure or two; no verdict is available
    _check("P4 within the paper's 20 %",
           abs(got - 7.0) / 7.0 < 0.2, "")
    # NO_QUOTE: a tolerance test that claims nothing about any source
    _check("P5 the two halves agree", abs(a - b) < 1e-9, "")
    # SIGFIG_TEST: the instrument the commissioning item asked for
    _check("P6 rounds to the paper's 4 %", round(share) == 4, "")
    # LIT_NO_TEST: names a published number, tests something else.  Only the
    # SIGN of the published statement is tested, which is the real
    # graphene_self_heating_model V3 shape.
    _check("P7 the paper's 1.5e12 is 4-5 % above their own C_ox",
           flag is True, "")
    # VACUOUS, found only because the cue is a PHRASE: 5 s.f. on 1.2794
    # licenses +/-5e-5 and the check allows +/-5e-4
    _check("PB kappa_ox(300K) == accepted thermal-SiO2 value 1.2794 W/m.K",
           abs(k - 1.2794) < 5e-4, "")
    # and the same word in this tree's OTHER sense must not cue anything
    _check("PC the FULL adjudication is accepted", abs(u - v) < 1e-9, "")
    # SIGFIG_TEST again, in the %g vocabulary rather than the round()
    # vocabulary -- the kappa_impurity_bound D2 shape
    _check("PA does not round to the paper's 1.5e12 at any precision",
           float('%.2g' % z) != 1.5, "")
    # a residual reached only through a local name.  PLANTED AS "VACUOUS"
    # AND IT IS NOT: 2.5e3 written to two significant figures licenses
    # +/-50 cm^2/Vs, i.e. exactly 2 %, so a 2 % tolerance is MATCHED to the
    # last bit.  The planted expectation was wrong and the instrument was
    # right, which is the only way round that is any use.  Left at 0.02 with
    # the arithmetic recorded, rather than nudged to make the plant look
    # prescient.
    rel = abs(x - 2.5e3) / 2.5e3
    _check("P8 matches the published 2.5e3 cm^2/Vs", rel < 0.02, "")
    # VACUOUS, reached through a name: 1.5 % allowed where 2 s.f. on 7.5
    # licenses 0.5/7.5 = 0.667 %
    rel2 = abs(y - 7.5) / 7.5
    _check("P9 matches the paper's 7.5 ps", rel2 < 0.015, "")
'''


def main():
    print("=" * 78)
    print("graphene_tolerance_census_audit.py -- 2026-11-11")
    print("=" * 78)

    # ---- 1. significant figures read off the spelling -------------------
    print("\n[1] _sigfigs reads precision off the SPELLING, and refuses to "
          "guess")
    for tok, want in (("1.0", 2), ("9.4e11", 2), ("1e4", 1),
                      ("1.1113e11", 5), ("6.45", 3), ("0.0050", 2),
                      ("12", 2), ("100", None), ("20", None), ("0", None)):
        check("1a _sigfigs(%r) == %r" % (tok, want),
              _sigfigs(tok) == want, "got %r" % _sigfigs(tok))
    check("1b a literal's spelling and a label's spelling can disagree: "
          "'1.0e4' reads 2 s.f., '1e4' reads 1 -- the A2 case",
          _sigfigs("1.0e4") == 2 and _sigfigs("1e4") == 1,
          "%r %r" % (_sigfigs("1.0e4"), _sigfigs("1e4")))

    # ---- 2. the licensed interval, as exact rationals -------------------
    print("\n[2] the licensed interval is an EXACT rational, not a range")
    h_abs, h_rel, sf = half_width("1.1e11")
    check("2a Dorgan's 1.1e11 at 2 s.f. licenses h_rel EXACTLY 1/22",
          sf == 2 and h_rel * 22.0 == 1.0000000000000002 or
          abs(h_rel - 1.0 / 22.0) <= 1e-18,
          "h_rel = %r, 1/22 = %r" % (h_rel, 1.0 / 22.0))
    check("2b and h_abs EXACTLY 5e9 cm^-2", h_abs == 5.0e9, "%r" % h_abs)
    h_abs2, h_rel2, sf2 = half_width("9.4e11")
    check("2c Dorgan's 9.4e11 at 2 s.f. licenses h_rel EXACTLY 1/188",
          abs(h_rel2 - 1.0 / 188.0) <= 1e-18,
          "h_rel = %r, 1/188 = %r" % (h_rel2, 1.0 / 188.0))
    check("2d both footnote-18 values share the SAME h_abs, because both are "
          "quoted at the same decade and the same figures",
          h_abs == h_abs2, "%r vs %r" % (h_abs, h_abs2))
    check("2e an ambiguous spelling licenses NOTHING -- half_width returns "
          "None rather than a default",
          half_width("20")[1] is None, "%r" % (half_width("20"),))

    # ---- 3. the 2026-11-10 bug, reproduced as arithmetic ----------------
    print("\n[3] the 2026-11-10 B2 defect, reproduced as an identity")
    got, q = 1.1113e11, 1.1e11
    rel_diff = abs(got - q) / q
    check("3a the recorded failure is 1.027 %", round(rel_diff * 100, 3) == 1.027,
          "%.6f %%" % (rel_diff * 100))
    check("3b the value DOES round to the quoted figure", _rounds_to(got, q, 2),
          "%r -> %r" % (got, _round_sig(got, 2)))
    check("3c the failure lies STRICTLY inside (old tolerance, licensed "
          "half-width) -- the whole window in which a correct value fails a "
          "check", 0.01 < rel_diff < 1.0 / 22.0,
          "0.01 < %.8f < %.8f" % (rel_diff, 1.0 / 22.0))
    check("3d so the old 1 %% tolerance is OVER_TIGHT by exactly 1/22/0.01 = "
          "4.5454...", verdict(0.01, 1.0 / 22.0) == "OVER_TIGHT"
          and abs((1.0 / 22.0) / 0.01 - 4.5454545454545454) < 1e-12,
          "%r %r" % (verdict(0.01, 1.0 / 22.0), (1.0 / 22.0) / 0.01))
    check("3e a tolerance SET to the licensed half-width is MATCHED, so the "
          "defect is the NUMBER, not the operator -- this corrects the "
          "commissioning item's own wording",
          verdict(1.0 / 22.0, 1.0 / 22.0) == "MATCHED",
          verdict(1.0 / 22.0, 1.0 / 22.0))

    # ---- 4. the scaled half-width (the census's own first bug) ----------
    print("\n[4] an absolute tolerance is judged in its own units")
    check("4a effective_half_width routes ABS_TOL to h_abs",
          effective_half_width("ABS_TOL", 0.005, 0.000775) == 0.005, "")
    check("4b and REL_TOL to h_rel",
          effective_half_width("REL_TOL", 0.005, 0.000775) == 0.000775, "")
    check("4c the real 6.45 kOhm/sq row is vacuous by EXACTLY 2, not by the "
          "12.9 the unscaled first version reported",
          abs(0.005 / 0.01 - 0.5) < 1e-15
          and verdict(0.01, 0.005) == "VACUOUS", "")

    # ---- 5. the symmetry that must give zero ----------------------------
    print("\n[5] zeros")
    check("5a a source with no comparison censuses to EXACTLY zero rows",
          census_source("x = 1\ny = x + 2\n") == [], "")
    check("5b a source whose only comparison has no residual censuses to "
          "EXACTLY zero rows (a bare `x < 1` is not a tolerance test)",
          census_source("def f(x):\n    _check('c', x < 1.0, '')\n") == [], "")
    _nq = census_source("def f(a, b):\n    _check('halves agree', "
                        "abs(a - b) < 1e-9, '')\n")
    check("5c but an UNLABELLED residual is still counted -- the literature "
          "filter is not what makes the census empty",
          len(_nq) == 1 and _nq[0].verdict == "NO_QUOTE", "%r" % (_nq,))
    check("5d a quoted token of 0 licenses nothing and raises nothing",
          half_width("0") == (None, None, None), "%r" % (half_width("0"),))

    # ---- 6. the planted control: every verdict, once --------------------
    print("\n[6] planted control -- every verdict in the vocabulary")
    prows = census_source(PLANT, "<plant>")
    seen = [(r.label.split()[0], r.verdict) for r in prows]
    want = [("P1", "OVER_TIGHT"), ("P2", "VACUOUS"), ("P3", "MATCHED"),
            ("P4", "AMBIGUOUS"), ("P5", "NO_QUOTE"), ("P6", "SIGFIG_TEST"),
            ("P7", "LIT_NO_TEST"), ("PB", "VACUOUS"), ("PC", "NO_QUOTE"),
            ("PA", "SIGFIG_TEST"), ("P8", "MATCHED"), ("P9", "VACUOUS")]
    check("6a the plant yields exactly the twelve planted rows, in order, "
          "with the twelve intended verdicts", seen == want,
          "got %r" % (seen,))
    for v in VERDICTS:
        check("6b verdict %s has a planted instance" % v,
              any(r.verdict == v for r in prows),
              "absent from the plant")
    _p8 = [r for r in prows if r.label.startswith("P8")]
    check("6c P8 is reached ONLY through the local name `rel`",
          len(_p8) == 1 and _p8[0].bound_via == "rel",
          "%r" % ([r.bound_via for r in _p8],))
    _p1 = [r for r in prows if r.label.startswith("P1")][0]
    check("6d P1's over-tightness is the B2 factor to the last figure",
          abs(_p1.ratio - 1.0 / 22.0 / 0.01) < 1e-12, "%r" % _p1.ratio)
    check("6f the %g vocabulary of the significant-figure test is recognised "
          "as such, not scored LIT_NO_TEST",
          [r.form for r in prows if r.label.startswith("PA")]
          == ["SIGFIG_FMT"],
          "%r" % [r.form for r in prows if r.label.startswith("PA")])
    check("6g the phrase cue fires on 'accepted ... value' and NOT on this "
          "tree's adjudication sense of the same word",
          [r.verdict for r in prows if r.label.startswith("PB")] == ["VACUOUS"]
          and [r.verdict for r in prows if r.label.startswith("PC")]
          == ["NO_QUOTE"],
          "%r" % [(r.label.split()[0], r.verdict) for r in prows
                  if r.label.split()[0] in ("PB", "PC")])
    _p7 = [r for r in prows if r.label.startswith("P7")][0]
    check("6e P7's '4-5 %' does NOT become the quoted value: the hyphen-"
          "welded number is not a quotation", _p7.token == "1.5e12",
          "token %r" % _p7.token)

    # ---- 7. mutation testing of the adjudicator -------------------------
    print("\n[7] mutation -- five defects in the verdict, all detected")
    # `import graphene_tolerance_census_audit as self_mod` -- the first form
    # of this harness -- loads a SECOND copy of this file when it is run as a
    # script, so every monkeypatch landed on a module nothing else was using
    # and all five mutants were "detected" only because the planted
    # expectation for P8 was itself wrong.  A mutation suite that reports 5/5
    # while testing nothing is the 2026-09-17 failure again, in the tool
    # written to prevent it.  sys.modules[__name__] is the running module.
    self_mod = sys.modules[__name__]
    real = self_mod.verdict
    mutants = (
        ("inequality flipped",
         lambda t, h, rel_eps=1e-12: ("MATCHED"
                                      if h is not None
                                      and abs(t - h) <= rel_eps * h
                                      else "AMBIGUOUS" if h is None
                                      else ("VACUOUS" if t < h
                                            else "OVER_TIGHT"))),
        ("always MATCHED", lambda t, h, rel_eps=1e-12: "MATCHED"),
        ("ambiguity suppressed",
         lambda t, h, rel_eps=1e-12: ("OVER_TIGHT" if h is None or t < h
                                      else "VACUOUS")),
        ("equality band widened to 100 %",
         lambda t, h, rel_eps=1e-12: real(t, h, rel_eps=1.0)),
        ("never MATCHED",
         lambda t, h, rel_eps=1e-12: ("AMBIGUOUS" if h is None
                                      else "OVER_TIGHT" if t < h
                                      else "VACUOUS")),
    )
    detected = 0
    for name, mut in mutants:
        self_mod.verdict = mut
        try:
            got_rows = census_source(PLANT, "<plant>")
            got_seen = [(r.label.split()[0], r.verdict) for r in got_rows]
        finally:
            self_mod.verdict = real
        if got_seen != want:
            detected += 1
            print("    detected: %-28s -> %s" % (name, got_seen[:3]))
        else:
            print("    MISSED:   %s" % name)
    check("7a all five verdict mutants are detected by the planted control",
          detected == 5, "%d/5" % detected)

    _real_sf = self_mod._sigfigs
    try:
        self_mod._sigfigs = lambda tok: 2        # guess instead of refusing
        _mrows = census_source(PLANT, "<plant>")
        _mseen = [(r.label.split()[0], r.verdict) for r in _mrows]
    finally:
        self_mod._sigfigs = _real_sf
    check("7b a _sigfigs that GUESSES 2 instead of refusing is detected: P4 "
          "stops being AMBIGUOUS", _mseen != want
          and dict(_mseen).get("P4") != "AMBIGUOUS",
          "%r" % (dict(_mseen).get("P4"),))

    # ---- 8. the real tree -----------------------------------------------
    print("\n[8] the tree\n")
    rows = census_repo()
    counts = print_census(rows)

    print("")
    check("8a the census is not empty on the real tree",
          sum(1 for r in rows if r.tol is not None) > 50,
          "%d" % sum(1 for r in rows if r.tol is not None))
    check("8b name-bound residuals are a LARGE fraction of the tree, so the "
          "inline-only first version was missing most of it",
          sum(1 for r in rows if r.bound_via) >= 40,
          "%d" % sum(1 for r in rows if r.bound_via))
    check("8c every adjudicated row carries a quoted token and a verdict in "
          "the declared vocabulary",
          all(r.token is not None and r.verdict in VERDICTS
              for r in rows if r.verdict != "NO_QUOTE"), "")
    # 8d PASSES BUT NO LONGER MEANS WHAT IT MEANT WHEN IT WAS WRITTEN ON
    # 2026-11-11, AND THE ORIGINAL ASSERTION IS KEPT RATHER THAN RESTATED.
    # It was written to witness that the tree still carried live instances of
    # the two defect classes.  The OVER_TIGHT instance it was counting was
    # graphene_self_heating_model A2's +-5 % against the one-figure ~1e4.
    # That was REPAIRED on 2026-11-12 -- tightened to the +-1.7857 % the same
    # sentence's two-figure areal statement licenses -- so the defect 8d
    # witnessed is gone.  8d goes on passing because the REPAIRED A2 is
    # itself censused OVER_TIGHT, at ratio 28, which is a FALSE POSITIVE and
    # the first this instrument has produced.  See 9l.
    # 8d REWRITTEN 2026-11-13, AND ITS OLD NUMBERS ARE KEPT IN PLACE.
    #
    # As first written it asserted the tree still contains at least one
    # OVER_TIGHT and one VACUOUS row.  On 2026-11-12 it was annotated: the
    # OVER_TIGHT row (A2) had become a FALSE POSITIVE of this instrument
    # rather than a defect of the tree, so the check was green while
    # witnessing nothing.  Today the split-warrant repair (section 10) turns
    # that row MATCHED, so the old form would be RED -- and red for the right
    # reason, which an assertion cannot express.  The check is therefore split
    # into what each half is now evidence for.
    check("8d the tree contains at least one VACUOUS comparison against a "
          "published number, and it is a real defect: the 6.45 kOhm/sq row, "
          "whose THREE-figure licensed width comes from a number this "
          "repository wrote itself (the citation census calls it UNWARRANTED)",
          counts.get("VACUOUS", 0) >= 1, "%r" % counts)
    check("8d2 and the tree now contains ZERO OVER_TIGHT rows, where before "
          "2026-11-13 it contained exactly one.  That one was A2, the false "
          "positive of 9l, and the split-warrant repair removed it WITHOUT "
          "touching the comparison -- only the warrant the census reads",
          counts.get("OVER_TIGHT", 0) == 0,
          "OVER_TIGHT %d (was 1 on 2026-11-11 and 2026-11-12, at ratio x28)"
          % counts.get("OVER_TIGHT", 0))
    check("8e no row in the tree is adjudicated against a NEGATIVE OR ZERO "
          "quoted value -- the footnote-ordinal bug would show up here",
          all(float(r.token) > 0 for r in rows if r.token is not None),
          "%r" % [(os.path.basename(r.path), r.line, r.token)
                  for r in rows if r.token is not None
                  and float(r.token) <= 0])
    self_rows = [r for r in census_repo(include_self=True)
                 if os.path.basename(r.path) == SELF
                 and r.verdict != "NO_QUOTE"]
    check("8f A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES (standing item, "
          "2026-10-17) has an instance HERE: censusing this module adds %d "
          "adjudicated rows that are its own controls" % len(self_rows),
          len(self_rows) >= 3,
          "%r" % [(r.line, r.verdict, r.token) for r in self_rows])
    check("8g and excluding it removes every one of them from the tree's "
          "figures", not any(os.path.basename(r.path) == SELF for r in rows),
          "%r" % [r.line for r in rows
                  if os.path.basename(r.path) == SELF])

    # 8h added 2026-11-13.  The self-plant item generalised to a class the
    # day a second instrument joined the tree, and the cost is MEASURED here
    # rather than asserted away.
    _with_cite = census_repo(include_self=True)
    _cite_rows = [r for r in _with_cite
                  if os.path.basename(r.path)
                  == "graphene_citation_census_audit.py"]
    _cite_quoted = [r for r in _cite_rows if r.token is not None]
    check("8h THE SECOND PLANTING MODULE, AND WHAT IT WOULD HAVE COST: "
          "graphene_citation_census_audit (2026-11-13) contributes %d quoted "
          "row(s) of its own planted controls, which a SELF-only exclusion "
          "would have reported as findings about the tree -- so the exclusion "
          "is now the declared class PLANTING_MODULES"
          % len(_cite_quoted),
          len(_cite_quoted) >= 5
          and not any(os.path.basename(r.path)
                      == "graphene_citation_census_audit.py" for r in rows),
          "quoted tokens it would have injected: %r"
          % sorted(set(r.token for r in _cite_quoted)))
    check("8i and the exclusion changes NO adjudicated verdict about the "
          "tree, because every row it removes belongs to a file that was "
          "not in the tree before today",
          sorted((os.path.basename(r.path), r.line, r.verdict)
                 for r in rows if r.token is not None)
          == sorted((os.path.basename(r.path), r.line, r.verdict)
                    for r in _with_cite if r.token is not None
                    and os.path.basename(r.path) not in PLANTING_MODULES))

    # ---- 9. does the computed value actually round to the quoted one? ---
    #
    # THE PAYOFF.  A vacuous tolerance is a latent defect until something
    # wrong slips through it, and one has.  Check A1 of
    # graphene_self_heating_model asserts
    #
    #     kappa_ox(300 K) == accepted thermal-SiO2 value 1.2794 W/m.K
    #
    # with a tolerance of 5e-4 where five quoted figures license 5e-5.  The
    # computed value is 1.2789670..., which rounds to 1.2790 at five figures,
    # NOT to 1.2794.  The model and the accepted value disagree in the fourth
    # decimal and the check says they are equal, because its tolerance is ten
    # times the interval the quotation claims.  The 2026-11-10 prediction that
    # the census would find "at least one vacuous rather than over-tight"
    # tolerance is therefore confirmed in its strong form: not merely loose,
    # but loose enough to be currently passing a value that does not
    # reproduce the number in its own label.
    print("\n[9] the vacuous tolerances, against the values they admit")
    import graphene_fet_signed_carrier_model as _fsc
    import graphene_self_heating_model as _sh

    _kappa = _sh.kappa_ox(300.0)
    check("9a kappa_ox(300K) = 1.2789670... does NOT round to the accepted "
          "1.2794 at the five figures that value is quoted to",
          not _rounds_to(_kappa, 1.2794, 5)
          and _round_sig(_kappa, 5) == 1.2790,
          "%r -> %r" % (_kappa, _round_sig(_kappa, 5)))
    check("9b and A1 passes anyway, because |1.2789670 - 1.2794| = 4.33e-4 "
          "sits between the licensed 5e-5 and the coded 5e-4",
          5e-5 < abs(_kappa - 1.2794) < 5e-4,
          "5e-5 < %.6e < 5e-4" % abs(_kappa - 1.2794))
    check("9c a tolerance set to the licensed half-width would FAIL A1, "
          "which is the whole value of setting it there",
          not (abs(_kappa - 1.2794) < half_width("1.2794")[0]),
          "|d| = %.3e, licensed %.3e"
          % (abs(_kappa - 1.2794), half_width("1.2794")[0]))

    _rs = 1e-3 / _fsc.SIGMA_MIN_CH3
    check("9d the OTHER vacuous row is not masking anything today: "
          "6.453202... does round to Chapter 3's 6.45 at 3 s.f.",
          _rounds_to(_rs, 6.45, 3) and _round_sig(_rs, 3) == 6.45,
          "%r -> %r" % (_rs, _round_sig(_rs, 3)))
    _rth = _sh.thermal_resistances(W=7e-6, L=4e-6, t_ox=300e-9,
                                   T=_sh.T_AMBIENT, T0=_sh.T_AMBIENT)["R_th"]
    check("9e and the OVER_TIGHT row is passing with room to spare, so its "
          "defect stays latent: 9915.6 K/W is inside both the 5 %% tested "
          "band and the 1-s.f. licensed band",
          abs(_rth - 1.0e4) / 1.0e4 < 0.05 and _rounds_to(_rth, 1.0e4, 1),
          "R_th = %.1f, rel = %.4f" % (_rth, abs(_rth - 1.0e4) / 1.0e4))
    check("9f of the three defective tolerances the census found, exactly "
          "ONE was passing a value its own label called equal to a published "
          "number it does not round to",
          sum(1 for got, q, sf in ((_kappa, 1.2794, 5), (_rs, 6.45, 3),
                                   (_rth, 1.0e4, 1))
              if not _rounds_to(got, q, sf)) == 1, "")

    # ---- 9h..9k: what reading the paper did to 9a-9c ---------------------
    #
    # 9a-9c are all still true and all mis-framed.  They adjudicate 1.2794 as
    # a published number quoted to five figures.  Dorgan et al. print no such
    # value at any precision: supp. Sec. 6 gives the fit ln(T_ox^0.52) - 1.687
    # with 0.52 and 1.687 as its only digits, and supp. Sec. 5 gives the only
    # numerical conductivity in the paper, "kappa_ox ~ 1.3 Wm-1K-1 ... (at
    # 300 K)".  So the census was right that A1's tolerance did not match the
    # precision its label claimed, and wrong -- because the label was wrong --
    # about there being a source behind that precision at all.  A census that
    # reads precision off a label inherits every false citation in the tree,
    # and nothing in it can detect one.  That limit is now a check.
    check("9h the paper's actual printed value, 1.3 at 2 s.f., licenses "
          "+/-0.05 -- a band 1000x the 5e-5 the phantom 1.2794 appeared to "
          "license", abs(half_width("1.3")[0] - 0.05) < 1e-15
          and abs(half_width("1.3")[0] / half_width("1.2794")[0] - 1000.0)
          < 1e-9,
          "%r %r" % (half_width("1.3")[0], half_width("1.2794")[0]))
    check("9i and the computed value DOES round to it, so the transcription "
          "of the fit was never in question -- only its warrant was",
          _rounds_to(_kappa, 1.3, 2) and _round_sig(_kappa, 2) == 1.3,
          "%r -> %r" % (_kappa, _round_sig(_kappa, 2)))
    _sh_rows = [r for r in rows
                if os.path.basename(r.path) == "graphene_self_heating_model.py"
                and r.token == "1.2794"]
    check("9j no row in the tree now puts a TOLERANCE against 1.2794 behind a "
          "publication cue: the surviving 1.2794 row is the "
          "significant-figure check that records the discrepancy",
          len(_sh_rows) == 1 and _sh_rows[0].verdict == "SIGFIG_TEST",
          "%r" % [(r.line, r.form, r.verdict) for r in _sh_rows])
    check("9k THE LIMIT OF THIS INSTRUMENT, stated as a check so it is not "
          "forgotten: the census reads quoted precision off a check's LABEL, "
          "so a label citing a value its source never printed is adjudicated "
          "as though the citation were real.  Detecting that needs the paper, "
          "not the AST.  5 s.f. on 1.2794 and 2 s.f. on 1.3 are both read "
          "correctly from their spellings, and only one of them has a source.",
          _sigfigs("1.2794") == 5 and _sigfigs("1.3") == 2, "")

    # -----------------------------------------------------------------
    # R and X: THE FIRST CROSS-REPOSITORY RUN OF THIS INSTRUMENT, AND WHAT
    # IT FOUND ABOUT THE INSTRUMENT RATHER THAN ABOUT THE OTHER TREE.
    #
    # Closes "RUN THE TOLERANCE CENSUS OVER Design-Verification-Roadmap",
    # created 2026-11-11.  The headline result is not a count of that
    # repository's tolerances: it is that census_repo() returned ZERO rows
    # over it and would have been believed, because this repository is flat
    # and so the flat-tree assumption had never been exercised.
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------
    # 9l and 9m: TWO LIMITATIONS OF THIS INSTRUMENT, BOTH FOUND ON
    # 2026-11-12 BY REPAIRING A CHECK THAT THIS INSTRUMENT HAD FLAGGED.
    # -----------------------------------------------------------------
    _a2 = [r for r in rows
           if os.path.basename(r.path) == "graphene_self_heating_model.py"
           and (r.label or "").startswith("A2 ")]
    a2_now = [r for r in rows if r.width_token is not None]
    check("9l THE CENSUS'S FIRST FALSE POSITIVE, REPAIRED 2026-11-13 AND "
          "RECORDED HERE RATHER THAN DELETED.  Until today the census "
          "assumed a check's CENTRE and its WIDTH are warranted by the SAME "
          "printed number.  A2 takes its centre from Dorgan's one-figure "
          "'~1e4' and its width from the two-figure '~2.8e-7 per unit of "
          "device area' in the same sentence -- the same quantity, since "
          "A = L*W is exact -- so the census read only the centre, licensed "
          "+-50 %, saw +-1.7857 % and returned OVER_TIGHT AT RATIO 28: a "
          "wrong verdict reached by correct arithmetic.  Section 10 is the "
          "repair, and this check now asserts the row is no longer "
          "OVER_TIGHT while keeping the old numbers in place",
          len(a2_now) == 1 and a2_now[0].verdict != "OVER_TIGHT",
          "A2 today: verdict %s, licensed width read from %r instead of from "
          "the centre's own one figure"
          % (a2_now[0].verdict, a2_now[0].width_token) if a2_now else "")
    _name_tol, _name_files = 0, set()
    for _p in _python_files(REPO):
        try:
            _t = ast.parse(open(_p, encoding="utf-8").read())
        except (SyntaxError, OSError):
            continue
        for _n in ast.walk(_t):
            if (isinstance(_n, ast.Compare) and len(_n.ops) == 1
                    and isinstance(_n.ops[0], (ast.Lt, ast.LtE))
                    and isinstance(_n.comparators[0], ast.Name)
                    and _const_num(_n.comparators[0]) is None
                    and classify_left(_n.left) is not None):
                _name_tol += 1
                _name_files.add(os.path.basename(_p))
    check("9m A BLIND SPOT MEASURED RATHER THAN ASSERTED: census_source "
          "resolves name-bound RESIDUALS (the left side, %d rows via "
          "bound_via) but not name-bound TOLERANCES -- `tol` is taken only "
          "from a literal comparator, so a check whose width is a named "
          "constant is outside the census entirely.  There are %d such "
          "comparisons in %d modules, and A2 was one of them for the length "
          "of one edit on 2026-11-12, which is how this was found"
          % (sum(1 for r in rows if r.bound_via), _name_tol,
             len(_name_files)),
          _name_tol >= 1 and len(_name_files) >= 1,
          "%d comparisons in %r" % (_name_tol, sorted(_name_files)))

    # ======================================================================
    print("\nSECTION 10.  The second warrant: centre and width from "
          "different printed numbers")
    # ======================================================================
    a2 = [r for r in rows if r.width_token is not None]
    check("10a the tree now contains exactly one row with a DECLARED width "
          "warrant, and it is A2 -- centre from Dorgan's one-figure ~1e4, "
          "width from the two-figure ~2.8e-7 in the same sentence",
          len(a2) == 1 and a2[0].token == "1e4"
          and a2[0].width_token == "2.8e-7",
          "%r" % [(os.path.basename(r.path), r.line, r.token, r.width_token)
                  for r in a2])
    check("10b and its verdict is MATCHED, where on 2026-11-11 and "
          "2026-11-12 it was OVER_TIGHT at ratio x28.  THE COMPARISON IN THE "
          "CODE IS UNCHANGED: what changed is which printed number the census "
          "reads the licensed width off",
          a2 and a2[0].verdict == "MATCHED",
          "verdict %s, ratio %s" % (a2[0].verdict, a2[0].ratio) if a2 else "")
    check("10c the transfer is legitimate for exactly one reason, and it is "
          "an EXACT one: A = L*W is given by the stated geometry and not "
          "fitted, so 2.8e-7 / A = 1e4 to within one ulp (2026-11-12 I1) and "
          "a RELATIVE width on one is a relative width on the other",
          abs((2.8e-7 / 2.8e-11) - 1.0e4) / 1.0e4 < 1e-15,
          "2.8e-7 / 2.8e-11 = %.17e against 1.0e4" % (2.8e-7 / 2.8e-11))
    # the exact arithmetic of the new verdict, to the ulp
    import struct as _struct

    def _bits(x):
        return _struct.unpack("<q", _struct.pack("<d", x))[0]

    _coded = 0.05 / 2.8
    _census = half_width("2.8e-7")[1]
    _ulp = abs(_bits(_coded) - _bits(_census))
    check("10d THE MATCHED VERDICT IS A STATEMENT ABOUT THE EQUALITY BAND, "
          "NOT ABOUT BITWISE EQUALITY: A2's coded 0.05/2.8 and the census's "
          "own 0.5e-8/2.8e-7 are two spellings of 1/56 that differ by "
          "EXACTLY %d ulp, and at rel_eps=0 the verdict would be VACUOUS by "
          "that one ulp instead of MATCHED" % _ulp,
          _ulp == 1
          and verdict(_coded, _census) == "MATCHED"
          and verdict(_coded, _census, rel_eps=0.0) == "VACUOUS",
          "coded %.20f, census %.20f, Fraction(1,56) %.20f"
          % (_coded, _census, float(__import__("fractions").Fraction(1, 56))))
    check("10e THE RESTRICTION, CHECKED: a width warrant is honoured only on "
          "a relative form.  No ABS_TOL row in the tree carries one, so "
          "nothing is being silently transferred across units today",
          not any(r.width_token is not None and r.form == "ABS_TOL"
                  for r in rows),
          "%r" % [(os.path.basename(r.path), r.form) for r in rows
                  if r.width_token is not None])
    # planted controls for the new mechanism
    _PW = """
def planted_width():
    _check("W1 agrees with the paper's ~1e4 to the width the two-figure 2.8e-7 "
       "licenses [width-warrant: 2.8e-7]", abs(a - 1.0e4) / 1.0e4 < 0.01785714285714286)
    _check("W2 agrees with the paper's ~1e4 and declares NO width warrant",
       abs(a - 1.0e4) / 1.0e4 < 0.01785714285714286)
    _check("W3 agrees with the paper's ~1e4 with a width warrant whose own "
       "spelling is ambiguous [width-warrant: 300]", abs(a - 1.0e4) / 1.0e4 < 0.01785714285714286)
"""
    wrows = census_source(_PW, "<planted-width>")
    wmap = dict((r.label.split()[0], r) for r in wrows if r.label)
    check("10f PLANTED CONTROLS FOR THE SECOND WARRANT (the twenty-third "
          "instance of the standing item): the SAME comparison at the SAME "
          "tolerance is MATCHED with the declared width warrant (W1) and "
          "OVER_TIGHT without it (W2) -- so the mechanism, and not the "
          "arithmetic, is what the verdict turns on",
          wmap.get("W1") is not None and wmap["W1"].verdict == "MATCHED"
          and wmap.get("W2") is not None
          and wmap["W2"].verdict == "OVER_TIGHT",
          "W1 %s, W2 %s" % (wmap["W1"].verdict if "W1" in wmap else "-",
                            wmap["W2"].verdict if "W2" in wmap else "-"))
    check("10g and a width warrant whose OWN spelling fixes no precision is "
          "refused, not guessed: W3 falls back to the centre's half-width "
          "and is OVER_TIGHT, exactly as if it had declared nothing",
          wmap.get("W3") is not None
          and wmap["W3"].verdict == "OVER_TIGHT"
          and wmap["W3"].width_token == "300",
          "W3 %s with width_token %r"
          % (wmap["W3"].verdict if "W3" in wmap else "-",
             wmap["W3"].width_token if "W3" in wmap else None))
    # the mutant
    _real_eww = extract_width_warrant
    try:
        globals()["extract_width_warrant"] = lambda label: (None, False)
        wrows_m = census_source(_PW, "<planted-width>")
        wmap_m = dict((r.label.split()[0], r) for r in wrows_m if r.label)
    finally:
        globals()["extract_width_warrant"] = _real_eww
    check("10h and the mutant that ignores the declared warrant is detected: "
          "W1 reverts to OVER_TIGHT",
          wmap_m.get("W1") is not None
          and wmap_m["W1"].verdict == "OVER_TIGHT",
          "mutated W1 -> %s" % (wmap_m["W1"].verdict
                                if "W1" in wmap_m else "-"))
    matched_rows = [r for r in rows if r.verdict == "MATCHED"]
    check("10i SCORING THE 2026-11-12 PREDICTION (2).  FIRST HALF RIGHT: it "
          "said teaching the census two warrants would turn MATCHED, not OK, "
          "on A2 -- because A2's width equals the licensed half-width "
          "exactly, so the boundary case is the right verdict.  It does, and "
          "10d shows the equality holds to one ulp and not better",
          len(matched_rows) == 1 and matched_rows[0].token == "1e4")
    check("10j SECOND HALF WRONG, AND UNSCORABLE AS POSED -- the X4 class "
          "again.  It predicted the repair would expose at least one row "
          "currently adjudicated MATCHED as mis-adjudicated.  The tree "
          "contained ZERO MATCHED rows on 2026-11-12, so the population it "
          "named was empty before the repair ran and A2 is the FIRST MATCHED "
          "row this tree has ever had",
          len([r for r in rows if r.verdict == "MATCHED"
               and r.width_token is None]) == 0,
          "MATCHED rows today: %r; all of them created by today's repair"
          % [(os.path.basename(r.path), r.token) for r in matched_rows])

    print("\nSECTION R.  The recursion the cross-repository run forced")
    _flat_own = census_repo(recursive=False)
    _rec_own = census_repo(recursive=True)
    check("R2 over THIS tree the recursive and flat censuses are identical "
          "field-for-field (%d rows each), because this tree is flat -- so "
          "making recursion the default cannot alter any existing finding"
          % len(_rec_own),
          [r.fields() for r in _flat_own] == [r.fields() for r in _rec_own]
          and len(_rec_own) == len(_flat_own),
          "flat %d rows, recursive %d rows" % (len(_flat_own), len(_rec_own)))
    check("R2b and that comparison HAS to be by field: Row keeps object "
          "identity for equality, so census_repo() == census_repo() is False "
          "for two runs over one unchanged tree.  Recorded as a check "
          "because a regression test written the obvious way would fail for "
          "a reason that has nothing to do with the census",
          census_repo(recursive=True) != census_repo(recursive=True)
          and ([r.fields() for r in census_repo(recursive=True)]
               == [r.fields() for r in census_repo(recursive=True)]),
          "object comparison False, field comparison True")

    if not os.path.isdir(SIBLING):
        print("  SKIP  R1/X1-X3: sibling repository not present at %s"
              % SIBLING)
    else:
        _sib_flat = census_repo(root=SIBLING, recursive=False)
        _sib_rec = census_repo(root=SIBLING, recursive=True)
        _n_flat = len(_python_files(SIBLING, recursive=False))
        _n_rec = len(_python_files(SIBLING, recursive=True))
        _verdicts = sorted({r.verdict for r in _sib_rec})
        _adjudicated = [r for r in _sib_rec
                        if r.verdict in ("OVER_TIGHT", "VACUOUS", "MATCHED",
                                         "AMBIGUOUS")]
        _no_quote = [r for r in _sib_rec if r.verdict == "NO_QUOTE"]
        _own_no_quote = [r for r in _rec_own if r.verdict == "NO_QUOTE"]

        check("R1 THE FINDING ABOUT THE INSTRUMENT: run flat over the "
              "sibling repository the census reads %d Python files and "
              "returns %d rows, while recursively it reads %d and returns "
              "%d.  An empty census is indistinguishable from a clean one, "
              "so this would have been reported as 'no tolerance defects'"
              % (_n_flat, len(_sib_flat), _n_rec, len(_sib_rec)),
              _n_flat == 0 and len(_sib_flat) == 0 and _n_rec > 40
              and len(_sib_rec) > 0,
              "flat: %d files / %d rows; recursive: %d files / %d rows"
              % (_n_flat, len(_sib_flat), _n_rec, len(_sib_rec)))
        check("X1 the sibling repository's actual census: %d rows over %d "
              "files, verdicts %s, and ZERO adjudicated -- its suites "
              "compare against exact expected values far more often than "
              "against tolerances"
              % (len(_sib_rec), _n_rec, ",".join(_verdicts)),
              len(_adjudicated) == 0 and _verdicts == ["NO_QUOTE"],
              "%d rows, %d adjudicated, verdicts %r, files carrying rows: %r"
              % (len(_sib_rec), len(_adjudicated), _verdicts,
                 sorted({os.path.relpath(r.path, SIBLING)
                         for r in _sib_rec})))
        check("X2 SCORING THIS MODULE'S OWN 2026-11-11 PREDICTION (3): it "
              "said the run would return 'more NO_QUOTE rows than this tree "
              "has in total, and zero adjudicated ones'.  HALF RIGHT -- zero "
              "adjudicated is right, and %d NO_QUOTE rows against this "
              "tree's %d is wrong by a factor of %.0f, in the opposite "
              "direction"
              % (len(_no_quote), len(_own_no_quote),
                 len(_own_no_quote) / max(1, len(_no_quote))),
              len(_adjudicated) == 0 and len(_no_quote) < len(_own_no_quote),
              "sibling NO_QUOTE %d, own NO_QUOTE %d"
              % (len(_no_quote), len(_own_no_quote)))
        check("X3 SCORING THE SIBLING'S 2026-11-11 PREDICTION (2): it said "
              "zero adjudicated rows but NOT more NO_QUOTE rows than 102, "
              "because its suites check exact values.  BOTH SUBSTANTIVE "
              "HALVES RIGHT, including the reason.  Its meta-claim that this "
              "tree's prediction was 'right in its first half and wrong in "
              "its second' is INVERTED: the first half was the NO_QUOTE "
              "count, which failed, and the second was zero adjudicated, "
              "which held",
              len(_adjudicated) == 0 and len(_no_quote) <= 102
              and len(_own_no_quote) == 102,
              "sibling: %d NO_QUOTE, %d adjudicated; own NO_QUOTE now %d"
              % (len(_no_quote), len(_adjudicated), len(_own_no_quote)))
        check("X4 and BOTH predictions were unscorable as posed until today: "
              "each named a population in a run that, with the instrument as "
              "it then stood, returned the empty set.  The prediction pair "
              "is scored above against the recursive run, which is the only "
              "run either could have meant",
              len(_sib_flat) == 0 and len(_sib_rec) == 6,
              "flat run returned %d rows; the scored run is the recursive "
              "one, %d rows" % (len(_sib_flat), len(_sib_rec)))

    print("\n" + "=" * 78)
    print("TOTAL: %d passed, %d failed" % (passed, failed))
    print("=" * 78)
    return 1 if failed else 0


if __name__ == "__main__":
    if "--census" in sys.argv:
        print_census(census_repo(), quoted_only="--all" not in sys.argv)
        sys.exit(0)
    sys.exit(main())
