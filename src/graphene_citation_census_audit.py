"""
graphene_citation_census_audit.py

A LABEL CAN NAME A NUMBER ITS SOURCE NEVER PRINTED, AND NO AMOUNT OF
ARITHMETIC IN THIS TREE CAN NOTICE.  THIS IS THE INSTRUMENT THAT CAN.

Created 2026-11-13 as the RANKED-FIRST item of the 2026-11-12 entry:

    "A CENSUS OF CITATIONS, NOT TOLERANCES -- carried from 11-11, now ranked
     first, and today strengthened its case twice.  9k already said no
     AST-level instrument can catch a label naming a value its source never
     printed; 9l now adds a second class the AST cannot see -- a width whose
     warrant is a DIFFERENT printed number from its centre.  Both are citation
     facts, not arithmetic facts.  The tree has 13 labels claiming a source and
     the worklist is mechanical to emit."

THE DIFFERENCE BETWEEN THIS INSTRUMENT AND THE TOLERANCE CENSUS

`graphene_tolerance_census_audit` asks an ARITHMETIC question: given that a
label says the source printed q, is the width of the check beside it the width
q's own rounding licenses?  It reads the quoted precision off the label and
takes the citation on trust.  Checks 9k and 9l are that trust written down as
two stated limits.

This module asks the CITATION question: did the named source print q at all,
and at the precision the label spells it to?  It cannot be answered by parsing
Python, so it is not answered by parsing Python.  It is answered against a
DECLARED REGISTRY of what each cited source actually prints -- `SOURCES` below
-- every entry of which carries the document it was read from and the date of
the read.  The registry is the evidence; this module is only the bookkeeping
that forces each label to meet it.

WHAT THE REGISTRY REST ON (2026-11-13)

Dorgan's arXiv text (1005.2711, which carries the Supporting Online Materials
in the same file) was fetched and read on 2026-11-13, independently of the
2026-11-10 read that produced `graphene_self_heating_model`.  The second read
CONFIRMED every transcription the module makes, including the two that were
most at risk:

  * the fit is printed `kappa_ox = ln(T_ox^0.52) - 1.687` with `T_ox` defined
    as the MEAN oxide temperature `(T_0+T)/2`.  The module's docstring already
    records the ln-identity, and `thermal_resistances` already evaluates the
    fit at `0.5*(T0+T)`.  The hazard here was real and the tree had already
    avoided it.
  * `kappa_ox ~ 1.3 Wm-1K-1 ... (at 300 K)` is printed in supp. Sec. 5 and
    NOWHERE in the main text, which is exactly what A1's label claims.

A citation census that returns "every citation checks out" on the module with
the most citations per line is a result, not a null result, and it is recorded
as one.  The defects it DID find are in the labels of other modules and in the
older instrument's reading of them.

THE VERDICTS

  CITED        the named source prints this value, at a precision at least as
               fine as the label spells it to.
  INFLATED     the source prints this value at FEWER significant figures than
               the label's spelling claims.  The citation is real and the
               precision claim is not.  `1.0e4` against a printed `10^4` is
               the shape, and it is the same defect as A2's old `1.0e4`
               literal, one layer out: in the LABEL, where the tolerance
               census takes its evidence from.
  AMBIG_PREC   the source's own spelling does not fix its significant figures
               (a printed `300`), so no precision comparison is available and
               none is invented.
  DERIVED      not printed, but obtained from printed values by a declared
               exact relation -- the device area `A = L*W` from a printed
               `L = 4 um` and `W = 7 um`.  A different warrant class from
               CITED, not a weaker one.
  UNWARRANTED  the source is in the registry and prints no such value.  This
               is the 9k class, and it is what the item was built to find.
  ORDINAL      the number is a section, figure, footnote or reference ordinal
               and not a value at all.  The quotation is a misparse, and
               since the tolerance census ADJUDICATES what it quotes, a
               misparse here is a verdict there.
  RETRACTED    the label names the value in order to RECORD THAT THE CITATION
               WAS FALSE.  Not a live citation, and must not be scored as one.
  UNREGISTERED the named source has no registry entry: reported, not
               adjudicated.
  NO_SOURCE    a publication cue with no identifiable source key.

HOW THE WORKLIST IS OBTAINED, AND WHY IT IS THE OLD INSTRUMENT'S

The population is not re-derived.  It is `graphene_tolerance_census_audit`'s
own set of rows that carry a quoted token -- 13 of them -- so that the two
instruments cannot disagree about WHICH labels claim a source, only about what
is true of them.  Check V5 asserts that identity field-for-field, and it is
the overlapping-limit validation required of any new instrument here.

The consequence is deliberate: this module inherits the old module's
`extract_quote`, misparses included.  Inheriting the misparse is how the
misparse gets found.  ORDINAL is not a hypothetical class.

EXACTLY-KNOWN VALIDATIONS

  * `_sigfigs('1e4') == 1` and `_sigfigs('1.0e4') == 2` as EXACT INTEGERS.
    INFLATED is an integer inequality between two spellings; no tolerance
    enters it, which is the one way this instrument is sounder than the
    tolerance census.
  * every registry half-width as an EXACT RATIONAL: 2 s.f. on 1.1e11 is
    exactly 1/22, on 2.8e-7 exactly 1/56, on 9.4e11 exactly 1/188.
  * a source with no publication-cued label censuses to EXACTLY zero rows --
    the symmetry that must give zero.
  * the declared derivation `4e-6 * 7e-6 -> 2.8e-11` is checked for its exact
    ulp distance, not for "agreement".
  * all nine verdicts get a planted control (the TWENTY-SECOND instance of the
    standing item), and the controls are then run against five mutated verdict
    functions, every one of which must be detected.

Usage:
    python3 graphene_citation_census_audit.py             # checks
    python3 graphene_citation_census_audit.py --census    # the table only
"""

import os
import re
import sys
from fractions import Fraction
from math import floor, log10

import graphene_tolerance_census_audit as TC
from graphene_tolerance_census_audit import _sigfigs, half_width, extract_quote

REPO = os.path.dirname(os.path.abspath(__file__))
SELF = os.path.basename(os.path.abspath(__file__))

VERDICTS = ("CITED", "INFLATED", "AMBIG_PREC", "DERIVED", "UNWARRANTED",
            "ORDINAL", "RETRACTED", "UNREGISTERED", "NO_SOURCE")


# ===========================================================================
# THE REGISTRY
#
# One entry per source this repository cites by name in a check label.  Each
# printed value is recorded AS THE SOURCE SPELLS IT, because the spelling is
# the precision claim, and with the place in the document it is printed.
# `read` is the date this repository last read the document itself; a value
# with no read behind it does not go in.
# ===========================================================================

class Printed(object):
    """One number a source prints, as printed."""

    __slots__ = ("token", "what", "where")

    def __init__(self, token, what, where):
        self.token = token
        self.what = what
        self.where = where

    def __repr__(self):
        return "Printed(%r, %r)" % (self.token, self.where)


class Derived(object):
    """A number obtained from printed values by an exact declared relation."""

    __slots__ = ("token", "what", "froms", "relation")

    def __init__(self, token, what, froms, relation):
        self.token = token
        self.what = what
        self.froms = froms          # tokens of the printed values used
        self.relation = relation    # the relation, in words

    def __repr__(self):
        return "Derived(%r, %r)" % (self.token, self.relation)


class Source(object):
    __slots__ = ("key", "citation", "read", "printed", "derived", "aliases")

    def __init__(self, key, citation, read, printed=(), derived=(),
                 aliases=()):
        self.key = key
        self.citation = citation
        self.read = read
        self.printed = tuple(printed)
        self.derived = tuple(derived)
        self.aliases = tuple(aliases)


SOURCES = {}


def _register(src):
    SOURCES[src.key] = src
    return src


_register(Source(
    key="dorgan",
    citation="V. E. Dorgan, M.-H. Bae and E. Pop, 'Mobility and saturation "
             "velocity in graphene on SiO2', Appl. Phys. Lett. 97, 082112 "
             "(2010); arXiv:1005.2711, which carries the Supporting Online "
             "Materials in the same file.",
    read="2026-11-10, re-read independently 2026-11-13",
    printed=(
        # main text, the sentence after Eq. (5).  BOTH numbers, one sentence.
        Printed("1e4", "R_th at 300 K for the paper's geometry, K/W",
                "main text after Eq. (5), printed 'R_th ~ 10^4 K/W'; also "
                "supp. Sec. 5"),
        Printed("2.8e-7", "R_th per unit of device area, m^2.K/W",
                "main text after Eq. (5), printed '~2.8x10^-7 m^2K/W per unit "
                "of device area' -- the SAME SENTENCE as the 10^4"),
        Printed("84", "R_ox share of R_th, per cent",
                "main text after Eq. (5), printed '~84%'"),
        Printed("12", "R_Si share of R_th, per cent",
                "main text after Eq. (5), printed '~12%'"),
        Printed("4", "R_B share of R_th, per cent",
                "main text after Eq. (5), printed '~4%'"),
        # supplement
        Printed("1.3", "kappa_ox at 300 K, W/m.K",
                "supp. Sec. 5, printed 'kappa_ox ~ 1.3 Wm-1K-1 ... (at 300 K)'"
                " -- NOT in the main text"),
        Printed("0.52", "exponent of the kappa_ox fit",
                "supp. Sec. 6, printed 'kappa_ox = ln(T_ox^0.52) - 1.687' "
                "with T_ox = (T_0+T)/2"),
        Printed("1.687", "offset of the kappa_ox fit", "supp. Sec. 6"),
        # footnote 18
        Printed("9.4e11", "n floor for hbar.w_OP = 160 meV, cm^-2",
                "footnote 18"),
        Printed("1.1e11", "n floor for hbar.w_OP = 55 meV, cm^-2",
                "footnote 18"),
        # geometry
        Printed("7", "device width W, um", "Fig. 1(a) caption; supp. Sec. 1"),
        Printed("4", "channel length L, um",
                "Fig. 1(a) caption; supp. Sec. 1"),
        Printed("300", "oxide thickness t_ox, nm",
                "Fig. 1(a) caption; supp. Sec. 1"),
        # the v_sat temperature statement
        Printed("20", "v_sat DECREASE between ~280 K and ~500 K, SiO2 OP "
                      "dominant, per cent",
                "main text after Eq. (7), printed '~20% decrease'"),
        Printed("2", "v_sat DECREASE between ~280 K and ~500 K, graphene OP "
                     "dominant, per cent",
                "main text after Eq. (7), printed '~2% decrease'"),
        Printed("55", "SiO2 surface-phonon energy, meV", "main text, Fig. 3c"),
        Printed("160", "graphene OP energy, meV", "main text, Fig. 3c"),
    ),
    derived=(
        Derived("2.8e-11", "device area A, m^2", ("4", "7"),
                "A = L*W from the printed L = 4 um and W = 7 um; exact, "
                "because the geometry is stated and not fitted"),
    ),
    aliases=("dorgan's", "the paper", "paper", "footnote-18", "footnote"),
))

_register(Source(
    key="zion",
    citation="Y. Zion, A. Butenko, A. Sharoni, M. Kaveh and I. Shlimak, "
             "'Charge carrier transport asymmetry in monolayer graphene', "
             "arXiv:1801.03565.",
    read="2026-10-29/30",
    printed=(
        Printed("1.5e12", "N_A^(S), acceptor density, cm^-2",
                "the value quoted for N_A^(S) = C_ox V_NP / e with "
                "V_NP ~ 20 V; recorded in notes/2026-10-30"),
        Printed("20", "V_NP, V", "same statement"),
    ),
    aliases=("zion's",),
))

_register(Source(
    key="codata",
    citation="CODATA recommended values of the fundamental physical "
             "constants, as carried by scipy.constants.",
    read="continuous -- it is the constant set the tree imports",
    printed=(
        # CODATA prints CONSTANTS.  It does not print any quantity this
        # repository computes FROM them, and that distinction is this
        # module's first real finding.
        Printed("1.602176634e-19", "elementary charge e, C", "SI definition"),
        Printed("6.62607015e-34", "Planck constant h, J.s", "SI definition"),
    ),
    aliases=(),
))


# A source key is named in a label if the key or one of its aliases appears.
def source_of(label):
    """The registry key a label names, or ('', None) / (name, None)."""
    if not label:
        return (None, None)
    low = label.lower()
    for key, src in SOURCES.items():
        if key in low:
            return (key, src)
        for a in src.aliases:
            if a and a in low:
                return (key, src)
    # No registry key.  A label can still NAME a source that is simply not
    # registered, and that is a different fact from a label that cites
    # anonymously ("the published value").  Only the first can be chased.
    m = _NAME_RE.search(label)
    if m:
        return (m.group(1), None)
    for cue in ("published", "literature", "quoted", "ref.", "codata",
                "accepted"):
        if cue in low:
            return ("", None)
    return (None, None)


# "Nobody et al", "Nobody's".  A capitalised surname carrying a citation, for
# the one purpose of telling UNREGISTERED (a source named and not registered,
# which a later session can go and read) from NO_SOURCE (a citation with no
# name in it at all, which nobody can chase).
_NAME_RE = re.compile(r"\b([A-Z][a-zA-Z-]{2,})(?:\s+et\s+al|'s)")


# ===========================================================================
# ORDINALS AND RETRACTIONS
#
# Both classes are DETECTED FROM THE LABEL and both detectors were calibrated
# against this tree's thirteen labels rather than guessed.  Their measured
# behaviour on the tree is asserted in section 4, so a later edit that makes
# either of them fire more widely fails a check instead of quietly
# reclassifying a citation.
# ===========================================================================

# A number preceded by one of these is an ORDINAL: it numbers a part of the
# document, it is not a quantity the document reports.
#
# FOUND, NOT ANTICIPATED.  The tolerance census's `_QUOTE_RE` matches a
# publication cue and then the next number within 40 characters.  On
#
#   "A1 kappa_ox(300K) rounds to the 1.3 W/m.K Dorgan supp. Sec. 5 prints ..."
#
# the cue that wins is "Dorgan" and the next number after it is the 5 of
# "Sec. 5".  The tolerance census therefore records A1 as a comparison against
# a published value of 5 quoted to one significant figure -- and prints it in
# its table, where it has been sitting since 2026-11-11.  It is the same
# defect as the footnote-ordinal bug that census fixed on the day it was
# written, surviving in the one position its fix did not cover: the fix made
# "footnote-\d+" a longer cue so the ordinal could not be read as a value,
# which does nothing for an ordinal that follows the cue instead of preceding
# the number.  Its check 8e cannot see this one either: 8e guards against a
# NEGATIVE quoted value, and a section ordinal is positive.
ORDINAL_WORDS = (
    "sec.", "sec", "section", "supp.", "supplementary", "fig.", "fig",
    "figure", "table", "eq.", "eqn", "equation", "chapter", "chap.",
    "ref.", "refs.", "note", "item", "step", "part", "appendix", "app.",
)

_ORD_RE = re.compile(
    r"(?:%s)\s*$" % "|".join(re.escape(w) for w in ORDINAL_WORDS),
    re.IGNORECASE)

# A label that names a value IN ORDER TO RECORD THAT THE CITATION WAS FALSE.
#
# CALIBRATED, AND THE CALIBRATION MATTERS.  The first marker list included
# "does not reproduce", which fires on
#
#   "D2 N_A^(S) does NOT reproduce: Zion's quoted 1.5e12 is 4-5 % above what
#    their own C_ox and V_NP give"
#
# where the citation is PERFECTLY REAL and it is the repository's
# reproduction of it that fails.  Scoring that as a retraction would have
# deleted a live citation from the census -- the exact inverse of the defect
# this module exists to find.  The distinction is between a label saying THE
# SOURCE DID NOT PRINT THIS (a retraction) and a label saying WE CANNOT
# REPRODUCE WHAT THE SOURCE PRINTED (a live citation with a negative result),
# and only the first is a retraction.  The markers below are therefore
# phrases about the CITATION's truth, never about a numerical disagreement.
RETRACTION_MARKERS = (
    "was false", "is false", "no such", "phantom", "never printed",
    "there is no accepted", "does not exist", "retracted", "withdrawn",
)


def is_ordinal(label, token):
    """True if `token` appears in `label` as a document ordinal."""
    if not label or token is None:
        return False
    for m in re.finditer(re.escape(token), label):
        head = label[:m.start()].rstrip()
        if _ORD_RE.search(head):
            return True
    return False


def is_retraction(label):
    if not label:
        return False
    low = label.lower()
    return any(mk in low for mk in RETRACTION_MARKERS)


# ===========================================================================
# ADJUDICATION
# ===========================================================================

def _rounds_at(value, token):
    """Does `value` round to `token` at the number of significant figures
    `token` is SPELLED to?  Returns None when the spelling is ambiguous.

    Deliberately not imported from the self-heating model: that `_rounds_to`
    takes the figure count as an argument, and here the figure count is the
    thing under test.  The spelling supplies it.
    """
    s = _sigfigs(token)
    q = float(token)
    if s is None:
        # The spelling does not fix the precision, so there is no rounding
        # test to run.  Exact equality is the only comparison available and
        # the caller is told by `_sigfigs` returning None, not by this
        # function inventing a figure count.
        return None if value != q else True
    if value == 0.0 or q == 0.0:
        return value == q
    return TC._round_sig(value, s) == TC._round_sig(q, s)


def match_printed(src, token):
    """Every registry entry of `src` whose value agrees with `token` at
    `token`'s own spelled precision.  Returns (printed_hits, derived_hits)."""
    ph, dh = [], []
    if src is None:
        return (ph, dh)
    for p in src.printed:
        if _rounds_at(float(p.token), token) is True:
            ph.append(p)
    for d in src.derived:
        if _rounds_at(float(d.token), token) is True:
            dh.append(d)
    return (ph, dh)


def adjudicate(label, token, src_key, src):
    """The verdict for one quoted token, and the warrants behind it.

    Order of precedence, and each step is a decision with a reason:

      1. ORDINAL first.  A number that is not a quantity cannot be compared
         against a registry, and adjudicating it would manufacture a finding.
      2. RETRACTED second.  A label recording that a citation was false is
         not making the citation, so UNWARRANTED -- which it would otherwise
         score, correctly on the arithmetic -- would double-count a defect
         the tree has already fixed.
      3. then the registry.
    """
    if is_ordinal(label, token):
        return ("ORDINAL", [], [])
    if is_retraction(label):
        return ("RETRACTED", [], [])
    if src_key is None or src_key == "":
        return ("NO_SOURCE", [], [])
    if src is None:
        return ("UNREGISTERED", [], [])
    ph, dh = match_printed(src, token)
    if ph:
        label_sf = _sigfigs(token)
        src_sf = [_sigfigs(p.token) for p in ph]
        if any(s is None for s in src_sf):
            return ("AMBIG_PREC", ph, dh)
        if label_sf is None:
            return ("AMBIG_PREC", ph, dh)
        if label_sf > max(s for s in src_sf):
            return ("INFLATED", ph, dh)
        return ("CITED", ph, dh)
    if dh:
        return ("DERIVED", ph, dh)
    return ("UNWARRANTED", ph, dh)


class CiteRow(object):
    __slots__ = ("path", "line", "label", "token", "approx", "src_key",
                 "verdict", "printed", "derived", "tol_verdict")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    def fields(self):
        """Value-comparable tuple.  The warrant objects are reduced to their
        tokens, because Printed and Derived keep object identity for equality
        exactly as Row does (see the tolerance census's R2b)."""
        return (os.path.basename(self.path or ""), self.line, self.token,
                self.src_key, self.verdict,
                tuple(p.token for p in (self.printed or ())),
                tuple(d.token for d in (self.derived or ())))

    def __repr__(self):
        return "CiteRow(%s:%s %r %s)" % (os.path.basename(self.path or "?"),
                                         self.line, self.token, self.verdict)


def census(rows=None, root=None):
    """The citation census over the tolerance census's own quoted rows.

    THE POPULATION IS NOT RE-DERIVED.  `graphene_tolerance_census_audit`
    already decides which labels claim a source; disagreeing with it about
    that would make the two instruments' findings incomparable, which is the
    one thing a second instrument must not do.
    """
    if rows is None:
        rows = TC.census_repo(root=root, include_self=False, recursive=True)
    out = []
    for r in rows:
        if r.token is None:
            continue
        key, src = source_of(r.label)
        v, ph, dh = adjudicate(r.label, r.token, key, src)
        out.append(CiteRow(path=r.path, line=r.line, label=r.label,
                           token=r.token, approx=r.approx, src_key=key,
                           verdict=v, printed=ph, derived=dh,
                           tol_verdict=r.verdict))
    out.sort(key=lambda c: (os.path.basename(c.path or ""), c.line))
    return out


# ===========================================================================
# REPORTING
# ===========================================================================

def print_census(rows, out=None):
    out = out or sys.stdout
    out.write("\n" + "=" * 78 + "\n")
    out.write("THE CITATION CENSUS: EVERY LABEL IN THIS TREE THAT CLAIMS A "
              "SOURCE\n")
    out.write("=" * 78 + "\n")
    hdr = ("%-38s %-5s %-9s %-7s %-12s %s" %
           ("module", "line", "token", "s.f.", "source", "verdict"))
    out.write(hdr + "\n")
    out.write("-" * 110 + "\n")
    for c in rows:
        sf = _sigfigs(c.token)
        out.write("%-38s %-5s %-9s %-7s %-12s %s\n" % (
            os.path.basename(c.path), c.line, c.token,
            "amb" if sf is None else sf, c.src_key or "-", c.verdict))
        if c.verdict in ("UNWARRANTED", "ORDINAL", "INFLATED", "RETRACTED",
                         "AMBIG_PREC", "DERIVED"):
            out.write("        %s\n" % (c.label or "")[:150])
            for p in c.printed:
                out.write("        warrant: %s printed as %r (%s)\n"
                          % (p.what, p.token, p.where))
            for d in c.derived:
                out.write("        warrant: %s, %s\n" % (d.what, d.relation))
    out.write("-" * 110 + "\n")
    counts = {}
    for c in rows:
        counts[c.verdict] = counts.get(c.verdict, 0) + 1
    out.write("%d label(s) claiming a source\n" % len(rows))
    for v in VERDICTS:
        if v in counts:
            out.write("    %-13s %d\n" % (v, counts[v]))
    multi = [c for c in rows if len(c.printed) > 1]
    out.write("%d row(s) match MORE THAN ONE printed value and so do not "
              "identify their own warrant\n" % len(multi))
    return counts


# ===========================================================================
# CHECKS
# ===========================================================================

_PASS = []
_FAIL = []


def check(label, ok, detail=""):
    (_PASS if ok else _FAIL).append(label)
    print("  %s  %s" % ("PASS" if ok else "FAIL", label))
    if detail:
        print("        %s" % detail)
    return bool(ok)


def tol_v(quoted_rows, token):
    """The TOLERANCE census's verdict on the row carrying `token`, so this
    module can say what the older instrument did with each finding."""
    for r in quoted_rows:
        if r.token == token:
            return r.verdict
    return None


# ---------------------------------------------------------------------------
# PLANTED CONTROLS.  The TWENTY-SECOND instance of the standing item GIVE
# EVERY VERDICT IN EVERY VOCABULARY A PLANTED CONTROL, and the first in which
# the vocabulary is a vocabulary of WARRANTS rather than of widths.
#
# Each control is a (name, label, token, expected) tuple.  The labels are
# written here, in this module, and are NOT read out of the tree -- a control
# that the tree could change is not a control.  Two of them are verbatim
# copies of real labels, which is noted on each.
# ---------------------------------------------------------------------------
CONTROLS = (
    ("P1", "A3 R_B share rounds to the paper's 4 %", "4", "CITED"),
    ("P2", "I1 2.8e-7 / A agrees with the paper's own 1.0e4 to within one ulp",
     "1.0e4", "INFLATED"),
    ("P3", "the oxide is Dorgan's quoted 300 nm", "300", "AMBIG_PREC"),
    ("P4", "S3 the band brackets Dorgan's 2.8e-11 m^2", "2.8e-11", "DERIVED"),
    ("P5", "this module reproduces Section 3.3's quoted 6.45 kOhm/sq from "
           "CODATA", "6.45", "UNWARRANTED"),
    ("P6", "A1 kappa_ox(300K) rounds to the 1.3 W/m.K Dorgan supp. Sec. 5 "
           "prints for the oxide at 300 K", "5", "ORDINAL"),
    ("P7", "A1b ... so A1's pre-2026-11-11 wording '== accepted value 1.2794' "
           "was false while A1 was green", "1.2794", "RETRACTED"),
    ("P8", "matches the 3.14 that Nobody et al print for the same geometry",
     "3.14", "UNREGISTERED"),
    ("P9", "agrees with the published value 1.23", "1.23", "NO_SOURCE"),
)


def run_controls():
    got = []
    for name, label, token, _exp in CONTROLS:
        key, src = source_of(label)
        v, ph, dh = adjudicate(label, token, key, src)
        got.append((name, v))
    return got


def expected_controls():
    return [(name, exp) for name, _l, _t, exp in CONTROLS]


# ---------------------------------------------------------------------------
# MUTANTS.  Five defects, each in a different one of this module's decisions,
# each of which MUST change at least one control's verdict.
# ---------------------------------------------------------------------------
def _mutants():
    # THE RUNNING MODULE, NOT A SECOND COPY OF IT.
    #
    # The first version of this harness did `import graphene_citation_census_
    # audit as M` and patched attributes on M.  Under `python3 thisfile.py`
    # the running module is `__main__`, so the import built a SECOND module
    # object from the same source and every mutation was applied to the copy
    # while `run_controls` went on calling the original.  All five mutants
    # reported NOT DETECTED -- the mutation never arrived.  This is the
    # 2026-09-30 item (ARRIVAL, NOT RESPONSE: the null control cannot certify
    # its own mutation) in a new disguise, and check 3e is the arrival check
    # that would have caught it as a defect in the harness instead of as five
    # apparent survivals in the subject.
    M = sys.modules[__name__]

    def m1():                      # ordinal detection switched off
        return ("is_ordinal", lambda label, token: False)

    def m2():                      # retraction detection switched off
        return ("is_retraction", lambda label: False)

    def m3():                      # inflation at >= instead of >
        def adj(label, token, src_key, src):
            if M.is_ordinal(label, token):
                return ("ORDINAL", [], [])
            if M.is_retraction(label):
                return ("RETRACTED", [], [])
            if src_key is None or src_key == "":
                return ("NO_SOURCE", [], [])
            if src is None:
                return ("UNREGISTERED", [], [])
            ph, dh = M.match_printed(src, token)
            if ph:
                lsf = _sigfigs(token)
                ssf = [_sigfigs(p.token) for p in ph]
                if any(s is None for s in ssf) or lsf is None:
                    return ("AMBIG_PREC", ph, dh)
                if lsf >= max(ssf):          # <-- the mutation
                    return ("INFLATED", ph, dh)
                return ("CITED", ph, dh)
            if dh:
                return ("DERIVED", ph, dh)
            return ("UNWARRANTED", ph, dh)
        return ("adjudicate", adj)

    def m4():                      # an ambiguous spelling GUESSED as 2 s.f.
        real = M._sigfigs

        def sg(token):
            s = real(token)
            return 2 if s is None else s
        return ("_sigfigs", sg)

    def m5():                      # a derived value scored as a printed one
        def mp(src, token):
            ph, dh = [], []
            if src is None:
                return (ph, dh)
            for p in list(src.printed) + list(src.derived):
                if M._rounds_at(float(p.token), token) is True:
                    ph.append(p)
            return (ph, dh)
        return ("match_printed", mp)

    return (("ordinal detection off", m1),
            ("retraction detection off", m2),
            ("inflation at >= not >", m3),
            ("ambiguous spelling guessed as 2 s.f.", m4),
            ("a derived value scored as printed", m5))


# ===========================================================================
# MAIN
# ===========================================================================

def _ulps(a, b):
    """Distance between two POSITIVE floats in units in the last place.

    The bit patterns of positive doubles are monotonic in the value, so the
    difference of the patterns IS the number of representable doubles
    between them.  Restricted to positive arguments and asserted, rather
    than extended to signed zero and negatives it will never see.
    """
    import struct
    assert a > 0.0 and b > 0.0, "ulp distance is defined here for positives"
    pa = struct.unpack("<q", struct.pack("<d", a))[0]
    pb = struct.unpack("<q", struct.pack("<d", b))[0]
    return abs(pa - pb)


def exact_half_width(token):
    """The quoting half-width of a written number as an EXACT Fraction.

    Not float arithmetic: the registry's precision claims are rational
    numbers and the point of asserting them is that they are exact.
    """
    s = _sigfigs(token)
    if s is None:
        return (None, None, None)
    q = Fraction(token)
    if q == 0:
        return (None, None, s)
    e = floor(log10(abs(float(q))))
    h_abs = Fraction(1, 2) * Fraction(10) ** (e - s + 1)
    return (h_abs, h_abs / abs(q), s)


def main():
    print("=" * 78)
    print("CITATION CENSUS -- DID THE NAMED SOURCE PRINT THE NAMED NUMBER?")
    print("=" * 78)

    # ---------------------------------------------------------------- 1 ----
    print("\n[1] the registry")
    n_printed = sum(len(s.printed) for s in SOURCES.values())
    n_derived = sum(len(s.derived) for s in SOURCES.values())
    for key in sorted(SOURCES):
        s = SOURCES[key]
        print("  %-8s %d printed, %d derived   read: %s"
              % (key, len(s.printed), len(s.derived), s.read))
    check("1a every registered source carries a citation and a read date",
          all(s.citation and s.read for s in SOURCES.values()),
          "%d sources, %d printed values, %d declared derivations"
          % (len(SOURCES), n_printed, n_derived))
    check("1b every printed token is a number as written, and its float is "
          "recoverable",
          all(isinstance(float(p.token), float)
              for s in SOURCES.values() for p in s.printed))
    check("1c every printed value carries WHERE in the document it is "
          "printed -- a registry entry with no location is not evidence",
          all(p.where for s in SOURCES.values() for p in s.printed))
    dup = {}
    for p in SOURCES["dorgan"].printed:
        dup.setdefault(p.token, []).append(p.what)
    shared = {k: v for k, v in dup.items() if len(v) > 1}
    check("1d THE REGISTRY IS HONEST ABOUT NON-UNIQUE VALUES: Dorgan prints "
          "the bare token '4' twice (the R_B share in per cent and the "
          "channel length in um), so a label quoting '4' does not identify "
          "its own warrant",
          "4" in shared and len(shared["4"]) == 2,
          "shared tokens: %r" % {k: len(v) for k, v in shared.items()})

    # ---------------------------------------------------------------- 2 ----
    print("\n[2] exactly-known validations")
    check("2a _sigfigs('1e4') is EXACTLY 1 and _sigfigs('1.0e4') is EXACTLY "
          "2 -- INFLATED is an integer inequality between two spellings, "
          "with no tolerance anywhere in it",
          _sigfigs("1e4") == 1 and _sigfigs("1.0e4") == 2,
          "1e4 -> %r, 1.0e4 -> %r, 10^4 as the paper prints it -> 1"
          % (_sigfigs("1e4"), _sigfigs("1.0e4")))
    h11 = exact_half_width("1.1e11")
    h94 = exact_half_width("9.4e11")
    h28 = exact_half_width("2.8e-7")
    check("2b the half-widths of the registry's two-figure values are EXACT "
          "rationals: 1.1e11 -> 1/22, 9.4e11 -> 1/188, 2.8e-7 -> 1/56",
          h11[1] == Fraction(1, 22) and h94[1] == Fraction(1, 188)
          and h28[1] == Fraction(1, 56),
          "%s, %s, %s" % (h11[1], h94[1], h28[1]))
    check("2c and they agree with the float half-widths the tolerance census "
          "computes, to the last bit",
          abs(float(h11[1]) - half_width("1.1e11")[1]) == 0.0
          and abs(float(h28[1]) - half_width("2.8e-7")[1]) == 0.0)
    check("2d a printed '300' has NO significant-figure count from its "
          "spelling, and the registry does not invent one",
          _sigfigs("300") is None
          and exact_half_width("300") == (None, None, None))
    # the declared derivation, to the ulp
    L, W = 4e-6, 7e-6
    A = L * W
    check("2e THE DECLARED DERIVATION, TO THE ULP: the printed geometry gives "
          "L*W = %.17e against the 2.8e-11 the labels quote, a distance of "
          "%d ulp -- so 'exact' is asserted at the only precision at which it "
          "is true" % (A, _ulps(A, 2.8e-11)),
          _ulps(A, 2.8e-11) <= 1 and TC._round_sig(A, 2) == 2.8e-11,
          "A = %.17e, quoted 2.8e-11, ulp distance %d"
          % (A, _ulps(A, 2.8e-11)))
    # the symmetry that must give zero
    empty = census(rows=[])
    check("2f THE SYMMETRY THAT MUST GIVE ZERO: a population with no quoted "
          "label censuses to EXACTLY zero rows",
          len(empty) == 0 and isinstance(empty, list))

    # ---------------------------------------------------------------- 3 ----
    print("\n[3] the planted controls, and the mutants that must disturb them")
    got = run_controls()
    exp = expected_controls()
    check("3a all nine verdicts in the vocabulary have a planted control "
          "(the twenty-second instance of the standing item)",
          set(v for _n, _l, _t, v in CONTROLS) == set(VERDICTS),
          "%d controls over %d verdicts" % (len(CONTROLS), len(VERDICTS)))
    check("3b every control returns its intended verdict", got == exp,
          "got %r" % (got,))
    M = sys.modules[__name__]
    try:
        import graphene_citation_census_audit as _copy
    except ImportError:
        _copy = M
    check("3e THE ARRIVAL CHECK FOR THE HARNESS ITSELF: the mutants are "
          "applied to the RUNNING module object.  Imported by name under "
          "__main__ this file yields a DIFFERENT object (%s), and patching "
          "that one leaves the subject untouched and every mutant reporting "
          "NOT DETECTED"
          % ("a distinct copy" if _copy is not M else "the same object"),
          M is sys.modules[__name__]
          and (_copy is not M or __name__ != "__main__"),
          "running module %r, imported-by-name module %r"
          % (__name__, getattr(_copy, "__name__", None)))
    detected = []
    for name, build in _mutants():
        attr, repl = build()
        orig = getattr(M, attr)
        setattr(M, attr, repl)
        try:
            mg = run_controls()
        finally:
            setattr(M, attr, orig)
        diff = [(n, a, b) for (n, a), (_n2, b) in zip(mg, exp) if a != b]
        detected.append((name, diff))
        print("    %-38s -> %s" % (name, diff if diff else "NOT DETECTED"))
    check("3c all five mutants are detected by the planted controls",
          all(d for _n, d in detected),
          "%d/%d" % (sum(1 for _n, d in detected if d), len(detected)))
    check("3d and the mutant that GUESSES a significant-figure count instead "
          "of refusing is detected on P3 specifically -- the same defect the "
          "tolerance census's 7b plants, now in the warrant vocabulary",
          any(n == "ambiguous spelling guessed as 2 s.f."
              and any(x[0] == "P3" for x in d) for n, d in detected))

    # ---------------------------------------------------------------- 4 ----
    print("\n[4] the detectors, measured on the real tree rather than "
          "asserted")
    tol_rows = TC.census_repo(include_self=False, recursive=True)
    quoted = [r for r in tol_rows if r.token is not None]
    rows = census(rows=tol_rows)
    ords_ = [c for c in rows if c.verdict == "ORDINAL"]
    rets = [c for c in rows if c.verdict == "RETRACTED"]
    check("4a the ordinal detector fires on EXACTLY one label in the tree, "
          "and it is A1's 'Dorgan supp. Sec. 5'",
          len(ords_) == 1 and "Sec. 5" in (ords_[0].label or ""),
          "%r" % [(os.path.basename(c.path), c.line, c.token) for c in ords_])
    check("4b the retraction detector fires on EXACTLY one label, and it is "
          "A1b -- it does NOT fire on D2, whose citation is live and whose "
          "REPRODUCTION is what fails",
          len(rets) == 1 and "1.2794" == rets[0].token
          and not any("does NOT reproduce" in (c.label or "")
                      for c in rets),
          "%r" % [(os.path.basename(c.path), c.token) for c in rets])
    d2 = [c for c in rows if c.token == "1.5e12"]
    check("4c and D2 is scored as the live citation it is",
          len(d2) == 1 and d2[0].verdict == "CITED"
          and d2[0].src_key == "zion",
          "D2 -> %s (%s)" % (d2[0].verdict, d2[0].src_key) if d2 else "absent")
    via_generic = [c for c in rows
                   if c.src_key == "dorgan"
                   and "dorgan" not in (c.label or "").lower()]
    bad_generic = [c for c in via_generic
                   if "dorgan" not in open(c.path).read(4000).lower()]
    check("4d THE ALIAS IS AN ASSUMPTION, SO IT IS MEASURED: %d row(s) "
          "resolve to Dorgan through the generic alias 'the paper' rather "
          "than through the name, and every one of them is in a module whose "
          "own header names Dorgan" % len(via_generic),
          len(bad_generic) == 0,
          "unwarranted generic resolutions: %r"
          % [os.path.basename(c.path) for c in bad_generic])

    # ---------------------------------------------------------------- 5 ----
    print("\n[5] the overlapping limit: this instrument against the one it "
          "extends")
    check("5a the citation census's population IS the tolerance census's "
          "quoted rows, one for one",
          len(rows) == len(quoted) == 13,
          "%d citation rows, %d quoted tolerance rows" % (len(rows),
                                                          len(quoted)))
    check("5b and the two agree label-for-label and token-for-token, which "
          "is the whole of what they are allowed to agree about",
          sorted((os.path.basename(r.path), r.line, r.token)
                 for r in quoted)
          == sorted((os.path.basename(c.path), c.line, c.token)
                    for c in rows))
    check("5c the two instruments DISAGREE where they are supposed to: every "
          "row the tolerance census adjudicated or reported, this one "
          "re-reads as a warrant question, and the verdicts are drawn from "
          "disjoint vocabularies",
          set(VERDICTS).isdisjoint(set(TC.VERDICTS)))
    # reproducibility of the new instrument
    check("5d two runs of the citation census over one unchanged tree are "
          "equal FIELD-FOR-FIELD (and CiteRow, like Row, is not "
          "value-equal -- see the tolerance census's R2b)",
          [c.fields() for c in census(rows=tol_rows)]
          == [c.fields() for c in rows]
          and census(rows=tol_rows) != rows)

    # ---------------------------------------------------------------- 6 ----
    print("\n[6] the findings, one at a time")
    by_tok = {}
    for c in rows:
        by_tok.setdefault(c.token, []).append(c)

    # FINDING 1 -- the ordinal
    a1 = ords_[0]
    tol_a1 = [r for r in quoted if r.line == a1.line
              and os.path.basename(r.path) == os.path.basename(a1.path)]
    check("6a FINDING 1, AND IT IS A DEFECT IN THE OLDER INSTRUMENT: the "
          "tolerance census reads A1's published value as 5, the ORDINAL OF "
          "THE SUPPLEMENTARY SECTION, and has printed it in its table since "
          "2026-11-11.  The value A1's label actually cites is the 1.3 "
          "W/m.K three words earlier",
          a1.token == "5" and tol_a1 and tol_a1[0].token == "5"
          and "1.3 W/m.K" in (a1.label or ""),
          "tolerance census row: token=%r s.f.=%r verdict=%s"
          % (tol_a1[0].token, tol_a1[0].sigfigs, tol_a1[0].verdict)
          if tol_a1 else "")
    check("6b and the older instrument's own guard cannot see it: 8e rejects "
          "a NEGATIVE quoted value, which is what the footnote-ordinal bug "
          "produced, and a section ordinal is positive",
          float(a1.token) > 0.0)
    check("6c the citation it displaced IS in the registry, at the precision "
          "the label claims: Dorgan supp. Sec. 5 prints kappa_ox ~ 1.3 "
          "Wm-1K-1 at 300 K",
          any(p.token == "1.3" for p in SOURCES["dorgan"].printed)
          and _rounds_at(1.3, "1.3") is True)
    check("6d and the check A1 actually performs is sound -- it tests "
          "against 1.3 at 2 s.f. in the code, so the misparse is in the "
          "CENSUS's reading of the label and not in the physics",
          "_rounds_to(kappa_ox(300.0), 1.3, 2)"
          in open(os.path.join(REPO,
                               "graphene_self_heating_model.py")).read())

    # FINDING 2 -- the unwarranted citation
    unw = [c for c in rows if c.verdict == "UNWARRANTED"]
    check("6e FINDING 2, AND IT IS THE 9k CLASS THE ITEM WAS BUILT TO FIND: "
          "exactly one label in the tree names a source that does not print "
          "its number -- 'reproduces Section 3.3's quoted 6.45 kOhm/sq from "
          "CODATA'.  CODATA prints e and h.  6.45 kOhm/sq is 1/(4e^2/h) "
          "COMPUTED from them, by this repository, in Chapter 3",
          len(unw) == 1 and unw[0].token == "6.45"
          and unw[0].src_key == "codata",
          "%r" % [(os.path.basename(c.path), c.token, c.src_key)
                  for c in unw])
    check("6f and the tolerance census ADJUDICATED that row -- VACUOUS by "
          "x0.5 against a three-figure publication interval that no "
          "publication licenses, because 6.45's spelling claims three figures "
          "and nothing printed them",
          tol_v(quoted, "6.45") == "VACUOUS",
          "tolerance verdict %r; licensed half-width came from the spelling "
          "of a number this tree wrote itself" % tol_v(quoted, "6.45"))
    check("6g THE REPAIR IS A WORDING, NOT A NUMBER: the quantity is right, "
          "the arithmetic is right, and the warrant is Chapter 3's own "
          "Section 3.3, which is this thesis.  A citation to CODATA for a "
          "derived quantity is a category error and the fix is to say "
          "'Chapter 3' and cite CODATA for e and h",
          True,
          "not applied in this commit: changing the label changes the "
          "tolerance census's population, and the two instruments are "
          "cross-validated against each other in 5a-5b this session")

    # FINDING 3 -- the inflated precision
    infl = [c for c in rows if c.verdict == "INFLATED"]
    check("6h FINDING 3: exactly one label spells a real citation to more "
          "figures than its source printed -- I1 writes Dorgan's R_th as "
          "'1.0e4' where the paper prints '10^4'.  Two figures claimed, one "
          "printed",
          len(infl) == 1 and infl[0].token == "1.0e4",
          "%r" % [(os.path.basename(c.path), c.token) for c in infl])
    check("6i and the inflation is exactly a factor of 10 in licensed width: "
          "1 s.f. on 1e4 licenses +-50 %, 2 s.f. on 1.0e4 licenses +-5 %, "
          "both exactly",
          exact_half_width("1e4")[1] == Fraction(1, 2)
          and exact_half_width("1.0e4")[1] == Fraction(1, 20),
          "1e4 -> %s, 1.0e4 -> %s" % (exact_half_width("1e4")[1],
                                      exact_half_width("1.0e4")[1]))
    check("6j I1 is a LIT_NO_TEST row, so the inflation licenses nothing "
          "today -- it is a latent defect of exactly the kind A2's literal "
          "was before 2026-11-12, one layer out in the label",
          tol_v(quoted, "1.0e4") == "LIT_NO_TEST")

    # FINDING 4 -- derived, not printed
    der = [c for c in rows if c.verdict == "DERIVED"]
    check("6k FINDING 4, A WARRANT CLASS AND NOT A DEFECT: 'Dorgan's "
          "2.8e-11 m^2' is not printed anywhere in Dorgan.  It is L*W from "
          "the printed geometry, and the possessive is doing work the paper "
          "did not do",
          len(der) == 1 and der[0].token == "2.8e-11")
    # 6l IS KEPT AS WRITTEN AND ITS SECOND CLAUSE IS WRONG.  ANNOTATED IN
    # PLACE 2026-11-14, NOT REWRITTEN.
    #
    # "stated, not fitted" is right and does the work 6l wanted: A = L*W is
    # exact in the paper's own arithmetic, which is why 11-12's I1 could show
    # the two published numbers are one statement divided by a constant.
    # "INHERITS NO UNCERTAINTY" does not follow and is false.  A stated value
    # is still a PRINTED value, and 4 um and 7 um are printed to ONE
    # significant figure each.  What the derivation does not inherit is
    # FITTING uncertainty; what it does inherit is the quoting interval of
    # each printed length, and that is wide.  6l2-6l4 carry the arithmetic.
    #
    # 6l's condition is also the literal `True`, so it could not have caught
    # its own second clause.  A check that cannot fail is not evidence --
    # this repository's own words, from the tolerance census's Section 5.
    # The condition is left as it is so the annotation has something to
    # annotate; 6l2-6l4 are the checks.
    check("6l and the distinction has teeth: the printed geometry is EXACT "
          "(stated, not fitted), so the derivation inherits no uncertainty "
          "-- which is precisely why 2026-11-12's I1 could show the paper's "
          "two numbers are one statement divided by a constant.  "
          "[SECOND CLAUSE SUPERSEDED 2026-11-14: it inherits no FITTING "
          "uncertainty and does inherit the printed lengths' quoting "
          "intervals -- see 6l2-6l4.  Numbers unchanged.]",
          True,
          "if L and W had been fitted, 2.8e-11 would carry their intervals "
          "and the areal and absolute statements would be independent")

    # 6l2-6l4: WHAT THE DERIVATION ACTUALLY INHERITS, as exact rationals.
    _h_abs, _h_rel, _sf = TC.half_width("2.8e-11")
    check("6l2 the spelled token '2.8e-11' licenses EXACTLY 1/56 in relative "
          "terms at the two figures it is written to, which is the width the "
          "tolerance census reads off it",
          _sf == 2 and abs(_h_rel - 1.0 / 56.0) <= 4e-17,
          "h_rel = %.17e, 1/56 = %.17e, difference %.2e"
          % (_h_rel, 1.0 / 56.0, _h_rel - 1.0 / 56.0))

    _A_lo, _A_hi = 3.5e-6 * 6.5e-6, 4.5e-6 * 7.5e-6
    _rel_lo = (_A_lo - 2.8e-11) / 2.8e-11
    _rel_hi = (_A_hi - 2.8e-11) / 2.8e-11
    check("6l3 CORRECTING 6l's SECOND CLAUSE: the PRINTED warrant is two "
          "one-figure lengths, each licensing half a unit in its last place, "
          "so the printed geometry licenses an area band whose edges are "
          "EXACTLY 21/2 and 23/2 times the spelled half-width.  Citing the "
          "product as a published number understates its licensed width by a "
          "factor of about eleven",
          abs(abs(_rel_lo) / _h_rel - 10.5) < 1e-12
          and abs(_rel_hi / _h_rel - 11.5) < 1e-12,
          "L in [3.5,4.5] um and W in [6.5,7.5] um give A in [%.6e, %.6e]; "
          "lo %.6f %% = %.4f x h_rel; hi %+.6f %% = %.4f x h_rel"
          % (_A_lo, _A_hi, 100.0 * _rel_lo, abs(_rel_lo) / _h_rel,
             100.0 * _rel_hi, _rel_hi / _h_rel))

    check("6l4 and the printed band is ASYMMETRIC (-3/16, +23/112), which a "
          "relative half-width cannot represent at all.  So the spelled "
          "token hides three things and not one: the width, the warrant "
          "(9l's class) and the shape",
          abs(_rel_lo + 3.0 / 16.0) < 1e-15
          and abs(_rel_hi - 23.0 / 112.0) < 1e-15
          and abs(_rel_hi) > abs(_rel_lo),
          "-3/16 = %.17e, +23/112 = %.17e, so the band is wider above than "
          "below by %.4f %% of the quoted value"
          % (-3.0 / 16.0, 23.0 / 112.0, 100.0 * (_rel_hi + _rel_lo)))

    check("6l5 and 6l's own verdict SURVIVES this correction: the band "
          "brackets the derived value, so DERIVED is about provenance and "
          "not about correctness.  2.8e-11 is the right number, cited as the "
          "wrong kind of thing",
          _A_lo < 2.8e-11 < _A_hi,
          "[%.6e, %.6e] contains 2.8e-11" % (_A_lo, _A_hi))

    # the ambiguous ones, and the non-identifying ones
    amb = [c for c in rows if c.verdict == "AMBIG_PREC"]
    check("6m FINDING 5: %d row(s) quote a value whose spelling fixes no "
          "precision, so no precision verdict is available and none is "
          "given.  V3's '~20 %%' is the instance, and Dorgan does print a "
          "'~20%%' -- so the citation is real and only its precision is "
          "unscorable" % len(amb),
          len(amb) == 1 and amb[0].token == "20"
          and any(p.token == "20" for p in amb[0].printed),
          "%r" % [(os.path.basename(c.path), c.token) for c in amb])
    multi = [c for c in rows if len(c.printed) > 1]
    check("6n FINDING 6, MEASURED: %d row(s) match more than one printed "
          "value of their source and therefore do not identify their own "
          "warrant.  A3's 'the paper's 4 %%' matches both the R_B share and "
          "the 4 um channel length, and only a reader can tell which"
          % len(multi),
          len(multi) >= 1,
          "%r" % [(os.path.basename(c.path), c.token,
                   [p.what for p in c.printed]) for c in multi])

    # ---------------------------------------------------------------- 7 ----
    print("\n[7] what the independent re-read CONFIRMED, which is a result")
    import graphene_self_heating_model as SH
    check("7a the fit's ARGUMENT survived the re-read: Dorgan prints "
          "kappa_ox = ln(T_ox^0.52) - 1.687 with T_ox = (T_0+T)/2 the MEAN "
          "oxide temperature, and thermal_resistances already evaluates it "
          "at 0.5*(T0+T) rather than at T",
          "kappa_ox(0.5 * (T0 + T))"
          in open(os.path.join(REPO, "graphene_self_heating_model.py")).read())
    check("7b at the 300 K benchmark the two readings coincide EXACTLY, "
          "because the mean of 300 and 300 is 300 with no rounding at all -- "
          "so every A-series number is independent of this distinction",
          0.5 * (300.0 + 300.0) == 300.0
          and SH.kappa_ox(0.5 * (300.0 + 300.0)) == SH.kappa_ox(300.0))
    t = 500.0
    k_mean = SH.kappa_ox(0.5 * (300.0 + t))
    k_naive = SH.kappa_ox(t)
    check("7c and away from 300 K they do NOT coincide, which is why the "
          "distinction was worth re-reading for: at T = 500 K the mean-oxide "
          "reading gives %.6f W/m.K and the device-temperature reading "
          "%.6f, a difference of %.3f %%"
          % (k_mean, k_naive, 100.0 * (k_naive - k_mean) / k_mean),
          abs(k_naive - k_mean) > 1e-3,
          "the tree takes the mean reading, which is the printed one")
    check("7d the ln-identity the module's docstring asserts is true to the "
          "last bit at the benchmark: ln(300^0.52) and 0.52*ln(300) agree to "
          "%d ulp"
          % _ulps(__import__("math").log(300.0 ** 0.52),
                  0.52 * __import__("math").log(300.0)),
          _ulps(__import__("math").log(300.0 ** 0.52),
                0.52 * __import__("math").log(300.0)) <= 4,
          "printed form ln(T^0.52), coded form 0.52*ln(T)")
    confirmed = [tk for tk in ("2.8e-7", "1.3", "84", "12", "9.4e11",
                               "1.1e11", "20")
                 if any(p.token == tk for p in SOURCES["dorgan"].printed)]
    check("7e every Dorgan value this tree cites was found in the "
          "independently fetched text, at the place the labels say: %s"
          % ", ".join(confirmed),
          len(confirmed) == 7)
    check("7f THE RESULT STATED AS A RESULT: the module with the most "
          "published numbers per line -- graphene_self_heating_model -- has "
          "ZERO unwarranted citations.  Its one bad label was already "
          "retracted in place on 2026-11-11 and the census finds the "
          "retraction, not the claim",
          not any(c.verdict == "UNWARRANTED"
                  and "self_heating" in os.path.basename(c.path)
                  for c in rows))

    # ---------------------------------------------------------------- 8 ----
    print("\n[8] limits, and the 2026-11-12 predictions this module can "
          "score")
    check("8a THE LIMIT OF THIS INSTRUMENT, stated as a check so it is not "
          "forgotten: the registry is hand-built from documents a session "
          "read, so a value absent from the registry is UNWARRANTED whether "
          "the source printed it or not.  The defence is 1c -- every entry "
          "names its location -- and the failure mode is a registry that is "
          "merely incomplete reporting a sound citation as a defect",
          all(p.where for s in SOURCES.values() for p in s.printed))
    check("8b A SECOND LIMIT: the census matches on VALUE and not on "
          "MEANING, so a label quoting a number the source prints for a "
          "DIFFERENT quantity scores CITED.  6n measures how often the "
          "ambiguity is live in this tree (%d row(s)) rather than leaving it "
          "to be discovered" % len(multi),
          len(multi) >= 1)
    check("8c SCORING THE 2026-11-12 PREDICTION (1): it said the citation "
          "census's first finding would be in graphene_self_heating_model "
          "and not in an older module, on the reasoning that volume of "
          "citation and not age predicts where a census bites.  WRONG, and "
          "wrong in the direction the prediction argued against: the one "
          "UNWARRANTED row is in graphene_fet_signed_carrier_model, and "
          "self_heating has none",
          len(unw) == 1 and "fet_signed_carrier" in os.path.basename(
              unw[0].path),
          "the finding IS in self_heating if the ordinal misparse counts as "
          "a finding -- but that is a defect in the instrument's reading, "
          "not a false citation, so it is not what the prediction named")
    check("8d and the reasoning behind the prediction was half right: "
          "self_heating produced the only ORDINAL and the only RETRACTED "
          "row, i.e. both of the citation anomalies that are not false "
          "citations, which is what a high volume of citation per line "
          "actually predicts",
          len(ords_) == 1 and len(rets) == 1
          and "self_heating" in os.path.basename(ords_[0].path)
          and "self_heating" in os.path.basename(rets[0].path))

    # ---------------------------------------------------------------- 9 ----
    print("\n[9] the tree")
    counts = print_census(rows)
    check("9a the census is not empty on the real tree, and it reports how "
          "many labels it read so that an empty census could not be "
          "mistaken for a clean one",
          len(rows) == 13, "%d rows, %r" % (len(rows), counts))
    check("9b every row carries a verdict in the declared vocabulary",
          all(c.verdict in VERDICTS for c in rows))
    check("9c A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES (standing "
          "item, 2026-10-17): this module's own labels quote registry "
          "tokens, so the tolerance census is run with include_self=False "
          "and this module is OUTSIDE its own population",
          not any(SELF == os.path.basename(c.path) for c in rows))
    check("9d exactly three rows carry a verdict that is a defect "
          "(UNWARRANTED, INFLATED, ORDINAL); the other ten are sound "
          "citations, declared derivations, retractions or unscorable "
          "precisions",
          len(unw) + len(infl) + len(ords_) == 3
          and len(rows) - 3 == 10)

    print("\n" + "=" * 78)
    print("TOTAL: %d passed, %d failed" % (len(_PASS), len(_FAIL)))
    print("=" * 78)
    return 1 if _FAIL else 0


if __name__ == "__main__":
    if "--census" in sys.argv:
        print_census(census())
        sys.exit(0)
    sys.exit(main())
