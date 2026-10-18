"""
graphene_standing_list_audit.py

A STANDING ITEM THAT ASSERTS A STATE OF THE TREE IS A CLAIM, AND NOTHING HERE
HAS EVER CHECKED ONE.

Created 2026-10-18, after two of this repository's own standing items were
found, in ten minutes of reading, to assert the opposite of what the tree
contains.  Both were in the "Not yet covered (candidates for future runs)"
list of the 2026-10-17 entry, i.e. in the one artefact every session reads to
decide what to do:

  (1) "CHAPTERS 2 AND 3 ARE STILL UNDRAFTED -- their computational results
      have been complete for weeks.  This is the largest remaining block of
      pure writing and has now gone three sessions without being touched."

      thesis_draft/02-electronic-structure.md was drafted 2026-09-28 (commit
      5bbd5f3, eleven sections) and thesis_draft/03-transport-and-optical-
      properties.md on 2026-09-29 (c1d5abc, ten sections).  Both have been in
      the tree for nineteen days.  AND THE LIST KNEW: from 2026-10-01 to
      2026-10-05 this item was carried, correctly, in the amended form
      "Chapters 2 and 3 ARE DRAFTED but Section 4.8.2's form of error is not
      audited for elsewhere".  It then vanished, and returned on 2026-10-17
      in its PRE-09-28 wording.  The gap is measured in Section 3 rather than
      stated here, so that this docstring cannot go stale the way its subject
      did.

      So this is not a stale item carried forward past its repair.  It is a
      RETIRED WORDING REVIVED -- and revived into first-rank contention: the
      2026-10-17 entry ends it "Worth ranking against item 16 deliberately
      rather than by default", so the one item offered as a competitor to the
      top measured item was describing writing committed three weeks earlier.

  (2) Item 16's precondition (ii), recorded 2026-10-17 as "Section 4.7's prose
      and Section 4.10's annotation -- is still untouched and is still
      writing, not plumbing".

      Commit 7604c05, 2026-10-13, is titled "thesis: Section 4.6.12, and the
      two sentences in Sections 4.7 and 4.10 that the measurement makes false,
      annotated in place".  Both annotations are in the file.  Section 4.7's
      names the per-metal numbers (-1.2583 % Cu, +20.8352 % Ni, +39.2281 % Au,
      +47.6925 % Pd), states that "Only Pd gives a small, plausible positive
      residual" is "true of the committed model and false of the corrected
      one", and says it is "left standing with this annotation rather than
      rewritten, PENDING THE ITEM 16 FLIP".  Section 4.10's amendment says the
      9.6 % margin "stands as computed" and is a margin against input error
      that a known modelling error already exceeds fivefold.

      Precondition (ii) is therefore not outstanding writing.  What is
      outstanding is the writing that can only be done WITH the flip -- when
      "pending the item 16 flip" itself goes stale and the committed table
      becomes the superseded column.  The 2026-10-17 instruction, "precondition
      (ii) and the flip in the same session", is right; its reason is wrong,
      and the difference matters, because a session that believes (ii) is
      missing can discharge it by writing a sentence about a state of the code
      that does not exist yet.

WHY AN INSTRUMENT AND NOT TWO CORRECTIONS

Design-Verification-Roadmap built tools/open_item_staleness_audit.py on
2026-10-13 for exactly this failure, and on 2026-10-17 that audit retired an
item of its own for the identical reason: "it names work finished before the
item started being carried".  That repository caught its instance on the same
day this one revived its instance.  The difference between the two days is one
module, and this is it.

AND THE LIST IS NOT A CLOSED SYSTEM, WHICH IS THE PART AN AUDIT CANNOT SEE

Both of this session's findings share a cause that no check in this file can
reach: the standing list has an EXTERNAL INPUT.  The scheduled task that runs
these sessions carries a frozen description of repository state dated
2026-09-17, including the sentence "Chapters 2-3 (electronic structure,
transport/optical) are still undrafted despite their computational results
being complete -- the largest remaining block of pure writing".  That is the
2026-10-17 bullet almost verbatim, eleven entries after the list had retired
it.  A session that reads its own brief as current state will revive retired
items, and will do so with the authority of an instruction rather than the
tentativeness of a memory.

That cannot be checked from inside the repository -- the brief is not in the
tree -- so it is recorded here as a NOT_MECHANISABLE row with its reason,
rather than silently omitted.  The checkable half is the half this file does:
ANY standing item that asserts a state of the tree is re-verified against the
tree before it is allowed to rank anything.

CLASSES
  CONFIRMED         the item's assertion about the tree is true
  CONTRADICTED      the tree says otherwise -- the item must not rank
  REVIVED           an earlier entry carried this item in a form the tree
                    confirms, and a later entry returned to a form it denies
  NOT_MECHANISABLE  the item makes no checkable assertion about the tree, or
                    its subject is outside the tree; reason REQUIRED
"""

import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, 'AUTOMATION_LOG.md')


# ---------------------------------------------------------------------------
# Reading the standing list out of the log
# ---------------------------------------------------------------------------

ENTRY_RE = re.compile(r'^## (20\d\d-\d\d-\d\d)', re.M)
NYC_RE = re.compile(r'\*\*Not yet covered[^\n]*\n(.*?)(?=\n\*\*Commits this run|\Z)',
                    re.S)


def entries(text):
    """[(date, body)] for every dated entry, oldest first."""
    marks = [(m.group(1), m.start()) for m in ENTRY_RE.finditer(text)]
    out = []
    for i, (d, s) in enumerate(marks):
        e = marks[i + 1][1] if i + 1 < len(marks) else len(text)
        out.append((d, text[s:e]))
    return out


def bullets(body):
    """Top-level '- ' bullets of an entry's Not-yet-covered list, each
    flattened to one line.  Continuation lines are indented, so a top-level
    bullet is '- ' at column 0."""
    m = NYC_RE.search(body)
    if not m:
        return []
    out, cur = [], None
    for line in m.group(1).split('\n'):
        if line.startswith('- '):
            if cur is not None:
                out.append(' '.join(cur.split()))
            cur = line[2:]
        elif cur is not None:
            if line.strip() == '':
                continue
            if line.startswith(('  ', '\t')):
                cur += ' ' + line.strip()
            else:
                out.append(' '.join(cur.split()))
                cur = None
    if cur is not None:
        out.append(' '.join(cur.split()))
    return out


# ---------------------------------------------------------------------------
# The baseline: assertions about the tree, each with a predicate
# ---------------------------------------------------------------------------
# Recorded by (id, trigger, expectation, predicate) with NO line numbers and
# no entry dates in the trigger -- a baseline that goes stale on an unrelated
# edit teaches the operator to ignore it (the 2026-10-17 rule).
#
# `expect_true` is what the ITEM asserts.  The predicate returns what the TREE
# says.  CONTRADICTED is expect_true != tree.

def _nonblank(path):
    with open(path, encoding='utf-8') as fh:
        return [ln for ln in fh.read().split('\n') if ln.strip()]


def _sections(path):
    with open(path, encoding='utf-8') as fh:
        return len(re.findall(r'^## ', fh.read(), re.M))


def p_chapters_23_undrafted(root):
    """Tree-side truth of 'Chapters 2 and 3 are undrafted'."""
    detail = []
    drafted = True
    for name in ('02-electronic-structure.md',
                 '03-transport-and-optical-properties.md'):
        path = os.path.join(root, 'thesis_draft', name)
        if not os.path.exists(path):
            detail.append('%s MISSING' % name)
            drafted = False
            continue
        n, s = len(_nonblank(path)), _sections(path)
        detail.append('%s %d non-blank lines, %d sections' % (name, n, s))
        if n < 200 or s < 5:
            drafted = False
    return (not drafted), '; '.join(detail)


# THIS FUNCTION'S FIRST DRAFT FAILED ITS OWN PLANT (P5), AND THE REASON IS THE
# SUBJECT OF THIS MODULE.  It tested `phrase in open(path).read()` for a
# sentence that the markdown file WRAPS -- the Section 4.10 amendment reads
#     > ... the fragility was realised by
#     > the model and not by its inputs.*
# so the byte string is absent while the claim is present -- and for numbers
# the file writes with a Unicode minus (U+2212) where the baseline had an
# ASCII hyphen.  The checker therefore reported "Section 4.10 is untouched"
# about a file that carries the amendment, i.e. it CONFIRMED the very item it
# was written to contradict, and it did so by answering a weaker question than
# the one asked: "are these exact bytes here", not "is this claim here".
#
# That is the same defect as 2026-10-17's textual grep for the rebinding route,
# and as the two docstrings that passed an audit by spelling while being the
# thing it audited for.  A verdict must be gated on the mechanism it reports
# (Design-Verification-Roadmap, 2026-10-16).  The fix is to normalise before
# matching; the FAIL is kept on the record in the transcript committed with
# this module.

_MINUSES = {'\u2212': '-', '\u2013': '-', '\u2014': '-', '\u00a0': ' '}


def _normalise(text):
    """Markdown-insensitive form: drop blockquote markers and list bullets,
    fold Unicode dashes and non-breaking spaces to ASCII, collapse all
    whitespace.  A claim that survives re-wrapping must survive this."""
    for k, v in _MINUSES.items():
        text = text.replace(k, v)
    text = re.sub(r'(?m)^[ \t]*(?:>[ \t]*)+', ' ', text)
    text = re.sub(r'(?m)^[ \t]*[-*][ \t]+', ' ', text)
    return re.sub(r'\s+', ' ', text)


def p_item16_ii_untouched(root):
    """Tree-side truth of 'Section 4.7's prose and Section 4.10's annotation
    are still untouched'.  Untouched means NEITHER annotation is present."""
    path = os.path.join(root, 'thesis_draft', '04-graphene-fet-device-physics.md')
    if not os.path.exists(path):
        return True, '04-graphene-fet-device-physics.md MISSING'
    txt = _normalise(open(path, encoding='utf-8').read())
    a47 = 'true of the committed model and false of the corrected one' in txt
    a410 = 'the fragility was realised by the model and not by its inputs' in txt
    nums = all(_normalise(s) in txt for s in ('+47.6925 %', '-203.276'))
    detail = ('4.7 annotation present=%s; 4.10 amendment present=%s; '
              'per-metal numbers present=%s' % (a47, a410, nums))
    return (not (a47 and a410)), detail


# RETIRED rows are the standing half, and they are why this suite does not go
# vacuous the moment its finding is acted on.  Once an item is retired from the
# list, "absent" is no longer a clean result to be read as a pass by default:
# the row asserts BOTH that the wording stays gone AND that the tree-side fact
# which justified retiring it still holds.  So the suite fails if the wording
# comes back (the failure mode this module exists for -- it has already
# happened once) and equally if the retirement turns out to have been
# premature, e.g. if the chapters were later truncated below the threshold.
#
# This is the R4/R5 shape of 2026-10-17: a count can be satisfied by the wrong
# rows, and a check that only counts what is present reads a cleared list and a
# broken parser identically.  Without it, this file would be clean on the day
# the finding landed and silent forever after -- which is precisely how the
# Name census came to be read as a clean bill of health for five sessions.

BASELINE = [
    # id, trigger, what the item asserts, tree predicate, gloss, retired_on
    ('A1', re.compile(r'chapters? 2 and 3 are still undrafted', re.I),
     True, p_chapters_23_undrafted,
     'the item asserts both chapters are undrafted', ''),
    ('A2', re.compile(r'\(ii\)[^|]*?is still untouched', re.I),
     True, p_item16_ii_untouched,
     "the item asserts Section 4.7's prose and Section 4.10's annotation are untouched",
     ''),
]

# Set to the empty string to re-arm a row as live (i.e. to put the item back on
# the list deliberately).  RETIRED_ON carries a DATE rather than a bool so the
# transcript says when the judgement was made and by which entry.


# Items whose subject is outside the tree.  REASON REQUIRED -- an exemption
# with no reason is a fault in its own right (the 2026-10-01 rule).
NOT_MECHANISABLE = [
    ('N1', 'the scheduled task brief is not in the tree',
     'The brief that drives these sessions carries a frozen 2026-09-17 '
     'description of repository state and is the probable source of the '
     'revived A1 wording.  Nothing in this repository can read it, so the '
     'check that would catch the cause cannot live here.  The compensating '
     'rule is procedural and is recorded in Section 7.9: a claim about tree '
     'state that arrives from outside the tree is re-verified before it ranks.'),
]


def audit(root, logtext):
    """-> (rows, bullets_seen).  rows: (id, cls, item_assert, tree_says, detail)"""
    ents = entries(logtext)
    latest = ents[-1] if ents else ('<none>', '')
    seen = bullets(latest[1])
    rows = []
    for bid, trig, expect_true, pred, what, retired_on in BASELINE:
        hit = [b for b in seen if trig.search(b)]
        tree, detail = pred(root)
        if retired_on:
            if hit:
                # The wording is back on the list after being retired.  That
                # is the whole subject of this module, so it is a FAILURE and
                # not a report.
                rows.append((bid, 'REVIVED', expect_true, tree,
                             'retired %s and present again in the %s list; %s'
                             % (retired_on, latest[0], detail)))
            elif expect_true == tree:
                # Absent, but the tree now agrees with the retired assertion:
                # the retirement has become wrong.
                rows.append((bid, 'RETIREMENT-STALE', expect_true, tree,
                             'retired %s on the grounds that the tree said '
                             'otherwise, and the tree no longer does; %s'
                             % (retired_on, detail)))
            else:
                rows.append((bid, 'RETIRED', expect_true, tree,
                             'retired %s; the tree still contradicts the '
                             'retired wording, so the retirement stands; %s'
                             % (retired_on, detail)))
            continue
        if not hit:
            rows.append((bid, 'ABSENT', expect_true, tree,
                         'no bullet in the %s list matches; %s' % (latest[0], detail)))
            continue
        cls = 'CONFIRMED' if (expect_true == tree) else 'CONTRADICTED'
        rows.append((bid, cls, expect_true, tree, '%s; %s' % (what, detail)))
    return rows, seen


# A row is a fault iff it is in one of these classes.  Named once, so the
# transcript and the exit code cannot disagree about what counts as a failure.
FAULT_CLASSES = ('CONTRADICTED', 'REVIVED', 'RETIREMENT-STALE', 'ABSENT')


def revivals(logtext, trig, antitrig):
    """Entries whose bullet matches `trig`, and entries whose bullet matches
    `antitrig` (the retired/corrected wording).  A trig entry strictly LATER
    than the last antitrig entry is a REVIVAL."""
    t, a = [], []
    for d, body in entries(logtext):
        bs = bullets(body)
        if any(trig.search(b) for b in bs):
            t.append(d)
        if any(antitrig.search(b) for b in bs):
            a.append(d)
    revived = [d for d in t if a and d > a[-1]]
    return t, a, revived


# ---------------------------------------------------------------------------
# Plants: the checker measured against exactly known answers
# ---------------------------------------------------------------------------

PLANT_LOG = """\
## 2026-01-01 - a plant

**Not yet covered (candidates for future runs):**

- **CHAPTERS 2 AND 3 ARE STILL UNDRAFTED** - planted; the tree says otherwise.
- **FLIP THE DEFAULT** - (ii) - Section 4.7's prose and Section 4.10's
  annotation - is still untouched and is still writing, not plumbing.
- **SOMETHING WITH NO TREE ASSERTION AT ALL** - carried.

**Commits this run:** 0
"""

# A list carrying A1's wording and NOT A2's.  Used by P3b to check that the
# rows are INDEPENDENT -- that one row going ABSENT does not change another
# row's verdict, and that a partial list is not silently treated as a whole
# one.  Written as a constant with a reader rather than left as scaffolding:
# an unread module-level constant is a provenance claim with no code behind it
# (graphene_dead_name_sweep, 2026-10-01), and this module would have shipped
# one.
PLANT_LOG_PARTIAL = """\
## 2026-01-01 - a plant carrying A1's wording only

**Not yet covered (candidates for future runs):**

- **CHAPTERS 2 AND 3 ARE STILL UNDRAFTED** - planted; A2's wording is absent.

**Commits this run:** 0
"""

PLANT_EMPTY = """\
## 2026-01-01 - a plant with no list at all

Nothing here.
"""


def main():
    global BASELINE
    passed = failed = 0

    def check(label, ok, detail=''):
        nonlocal passed, failed
        if ok:
            passed += 1
            print('  [PASS] %s' % label)
        else:
            failed += 1
            print('  [FAIL] %s' % label)
        if detail:
            for ln in str(detail).split('\n'):
                print('         %s' % ln)

    print('=' * 78)
    print('STANDING-LIST AUDIT -- do the items that rank the work')
    print('tell the truth about the tree?')
    print('2026-10-18.  Two items in the 2026-10-17 list assert the opposite')
    print('of what the tree contains, and one of them was revived from a')
    print('wording this list had already retired.')
    print('=' * 78)

    logtext = open(LOG, encoding='utf-8').read()

    # -- Section 1: plants, against exactly known answers -------------------
    print('\n' + '-' * 78)
    print('1.  POSITIVE CONTROL -- the checker on plants whose answers are known')
    print('-' * 78)

    # A plant tree in which chapters 2 and 3 really ARE undrafted, and in
    # which 4.7/4.10 really ARE untouched: both items must come back
    # CONFIRMED.  Same checker, same baseline, opposite verdicts -- which is
    # what distinguishes a checker from a constant.
    with tempfile.TemporaryDirectory() as td:
        os.makedirs(os.path.join(td, 'thesis_draft'))
        for nm in ('02-electronic-structure.md',
                   '03-transport-and-optical-properties.md'):
            with open(os.path.join(td, 'thesis_draft', nm), 'w') as fh:
                fh.write('# stub\n')
        with open(os.path.join(td, 'thesis_draft',
                               '04-graphene-fet-device-physics.md'), 'w') as fh:
            fh.write('## 4.7 nothing annotated here\n')
        _live = [(b[0], b[1], b[2], b[3], b[4], '') for b in BASELINE]
        _saved, BASELINE = BASELINE, _live
        try:
            rows, seen = audit(td, PLANT_LOG)
        finally:
            BASELINE = _saved
        got = {r[0]: r[1] for r in rows}
        check('P1 empty-tree plant, rows LIVE: A1 CONFIRMED (chapters really absent)',
              got.get('A1') == 'CONFIRMED', 'got %s' % got.get('A1'))
        check('P2 empty-tree plant, rows LIVE: A2 CONFIRMED (annotations really absent)',
              got.get('A2') == 'CONFIRMED', 'got %s' % got.get('A2'))
        check('P3 plant list parsed to exactly 3 bullets',
              len(seen) == 3, 'got %d: %r' % (len(seen), [b[:40] for b in seen]))
        
        # P3b ROW INDEPENDENCE, on a list carrying A1's wording and not A2's.
        _sv, BASELINE = BASELINE, _live
        try:
            _rows, _seen = audit(HERE, PLANT_LOG_PARTIAL)
        finally:
            BASELINE = _sv
        _g = {r[0]: r[1] for r in _rows}
        check('P3b rows are independent: a list with only A1 gives '
              'A1 CONTRADICTED and A2 ABSENT',
              _g.get('A1') == 'CONTRADICTED' and _g.get('A2') == 'ABSENT',
              'A1=%s A2=%s bullets=%d' % (_g.get('A1'), _g.get('A2'), len(_seen)))

    # The real tree with the plant list.  Run once with the retirements
    # CLEARED -- this is the verdict the rows carried on 2026-10-18, before
    # they were retired, and the transcript committed with this module in its
    # red state is that run.  Keeping it as a plant means the finding cannot
    # evaporate when the list is corrected (the 2026-10-17 vacuity rule).
    live = [(b[0], b[1], b[2], b[3], b[4], '') for b in BASELINE]
    saved, BASELINE = BASELINE, live
    try:
        rows, _ = audit(HERE, PLANT_LOG)
    finally:
        BASELINE = saved
    got = {r[0]: r[1] for r in rows}
    check('P4 real tree, plant list, rows LIVE: A1 CONTRADICTED',
          got.get('A1') == 'CONTRADICTED', 'got %s' % got.get('A1'))
    check('P5 real tree, plant list, rows LIVE: A2 CONTRADICTED',
          got.get('A2') == 'CONTRADICTED', 'got %s' % got.get('A2'))

    # And once AS RETIRED against the same plant list: the retired wording is
    # present, so both must come back REVIVED.  This is the check that would
    # have caught 2026-10-17 on the day it happened.
    #
    # Both plants force the retirement state explicitly rather than reading
    # whatever the baseline currently carries, so this whole section is
    # invariant under retiring or re-arming a row.  A plant whose expected
    # answer moves when the configuration moves is not a control.
    retd = [(b[0], b[1], b[2], b[3], b[4], b[5] or '2026-10-18')
            for b in saved]
    BASELINE = retd
    try:
        rows, _ = audit(HERE, PLANT_LOG)
    finally:
        BASELINE = saved
    got = {r[0]: r[1] for r in rows}
    check('P5b retired rows + plant list carrying the wording -> REVIVED',
          got.get('A1') == 'REVIVED' and got.get('A2') == 'REVIVED',
          'A1=%s A2=%s' % (got.get('A1'), got.get('A2')))
    check('P5c REVIVED is a fault class, so a revival cannot pass quietly',
          'REVIVED' in FAULT_CLASSES)

    # VACUITY.  A list the parser cannot find must not read as a clean bill of
    # health -- that is exactly how the inert-knob census came to be read as
    # one for five sessions (2026-10-17).  ABSENT is not CONFIRMED.
    saved, BASELINE = BASELINE, live
    try:
        rows, seen = audit(HERE, PLANT_EMPTY)
    finally:
        BASELINE = saved
    got = {r[0]: r[1] for r in rows}
    check('P6 VACUITY, rows LIVE: no list found -> every row ABSENT, '
          'none CONFIRMED',
          seen == [] and set(got.values()) == {'ABSENT'},
          'bullets=%d classes=%r' % (len(seen), sorted(set(got.values()))))
    check('P6b ABSENT is a fault class, so an unparseable list cannot read '
          'as a clean run', 'ABSENT' in FAULT_CLASSES)

    # The retired case on the same empty list is RETIRED and that is correct --
    # but only because the tree predicate still runs.  Assert that it did: a
    # RETIRED row that skipped its predicate would be a bare "the wording is
    # gone", which is the report this module refuses to accept as a check.
    BASELINE = retd
    try:
        rows, _ = audit(HERE, PLANT_EMPTY)
    finally:
        BASELINE = saved
    treevals = {r[0]: r[3] for r in rows}
    check('P6c retired rows still evaluate the tree predicate on an empty list',
          treevals.get('A1') is False and treevals.get('A2') is False,
          'A1 tree=%r A2 tree=%r' % (treevals.get('A1'), treevals.get('A2')))

    # DISJOINTNESS (the Design-Verification-Roadmap 2026-10-17 rule): the
    # mechanisable and not-mechanisable sets must not overlap, so that
    # loosening either cannot move a contradiction into the benign column.
    mech = {b[0] for b in BASELINE}
    nmech = {n[0] for n in NOT_MECHANISABLE}
    check('P7 baseline ids and NOT_MECHANISABLE ids are disjoint',
          not (mech & nmech), 'overlap %r' % sorted(mech & nmech))
    check('P8 every NOT_MECHANISABLE row carries a reason of >= 40 chars',
          all(len(n[2]) >= 40 for n in NOT_MECHANISABLE),
          '%d rows' % len(NOT_MECHANISABLE))

    # -- Section 2: the real list -------------------------------------------
    print('\n' + '-' * 78)
    print('2.  THE REAL STANDING LIST, AS OF THE LATEST ENTRY')
    print('-' * 78)

    ents = entries(logtext)
    rows, seen = audit(HERE, logtext)
    print('  latest entry: %s' % ents[-1][0])
    print('  top-level bullets in its Not-yet-covered list: %d' % len(seen))
    print('  baseline assertions checked: %d' % len(rows))
    print()
    print('  %-4s %-14s %-9s %-9s' % ('id', 'class', 'item says', 'tree says'))
    for bid, cls, ia, ts, detail in rows:
        print('  %-4s %-14s %-9s %-9s' % (bid, cls, ia, ts))
        for ln in re.findall(r'.{1,66}(?:\s|$)', detail):
            if ln.strip():
                print('       %s' % ln.strip())
    print()
    for nid, subj, reason in NOT_MECHANISABLE:
        print('  %-4s NOT_MECHANISABLE  %s' % (nid, subj))
        for ln in re.findall(r'.{1,66}(?:\s|$)', reason):
            if ln.strip():
                print('       %s' % ln.strip())

    # NON-VACUITY on the real list: if the parser silently found nothing, the
    # rows above are an artefact and the two FAILs below would read as PASSes.
    print()
    check('R1 the latest list parsed to at least 10 top-level bullets',
          len(seen) >= 10, 'got %d' % len(seen))
    check('R2 no baseline row is in a fault class %r' % (FAULT_CLASSES,),
          all(r[1] not in FAULT_CLASSES for r in rows),
          'faults: %r' % [(r[0], r[1]) for r in rows if r[1] in FAULT_CLASSES])
    check('R3 every retired row is still justified by the tree',
          all(r[1] != 'RETIREMENT-STALE' for r in rows),
          'stale retirements: %r'
          % [r[0] for r in rows if r[1] == 'RETIREMENT-STALE'])
    check('R4 no retired wording has returned to the list',
          all(r[1] != 'REVIVED' for r in rows),
          'revived: %r' % [r[0] for r in rows if r[1] == 'REVIVED'])

    # -- Section 3: the revival ---------------------------------------------
    print('\n' + '-' * 78)
    print('3.  A1 IS NOT A STALE ITEM.  IT IS A RETIRED WORDING REVIVED')
    print('-' * 78)

    trig = re.compile(r'chapters? 2 ?(?:and|-|–) ?3 (?:are still|remain) undrafted',
                      re.I)
    anti = re.compile(r'chapters? 2 and 3 are drafted', re.I)
    t, a, revived = revivals(logtext, trig, anti)
    print('  entries asserting UNDRAFTED : %s' % ', '.join(t))
    print('  entries asserting DRAFTED   : %s' % ', '.join(a))
    print('  revived (undrafted AFTER the last drafted) : %s'
          % (', '.join(revived) or 'none'))
    print()
    alldates = [d for d, _ in entries(logtext)]
    gap = [d for d in alldates if a and d > a[-1] and d < (revived[0] if revived else 'ZZ')]
    print('  entries between the last corrected one and the revival : %d%s'
          % (len(gap), (' (%s..%s)' % (gap[0], gap[-1])) if gap else ''))
    print()
    print('  Chapter 2 was drafted 2026-09-28 and Chapter 3 2026-09-29, so')
    print('  every UNDRAFTED entry up to 09-27 was true when written.  The')
    print('  list then tracked the drafting correctly for %d entries' % len(a))
    print('  (%s..%s), carried nothing for %d, and the entry after those is'
          % (a[0], a[-1], len(gap)))
    print('  the revival.  An item absent long enough stops being remembered')
    print('  as retired and starts being available to be re-derived.')

    check('  V1 the list did at one point carry the corrected wording',
          len(a) >= 1, 'entries: %r' % a)
    check('  V2 at least one entry revived the retired wording',
          len(revived) >= 1, 'revived: %r' % revived)
    check('  V3 the revival is strictly later than the last corrected entry',
          bool(revived) and bool(a) and min(revived) > max(a),
          'first revival %s vs last corrected %s'
          % (revived[0] if revived else '-', a[-1] if a else '-'))

    # -- Section 4: git corroboration for A2 --------------------------------
    print('\n' + '-' * 78)
    print('4.  A2 CORROBORATED FROM THE HISTORY, NOT ONLY FROM THE FILE')
    print('-' * 78)
    try:
        # `git -S` matches raw bytes, so the phrase chosen here is one the
        # file does NOT wrap.  A phrase that wraps would come back empty and
        # read as "never committed" -- the P5 defect, one layer out.
        out = subprocess.run(
            ['git', '-C', HERE, 'log', '--format=%h %ad %s', '--date=short',
             '-S', 'true of the committed model and false of the corrected one',
             '--', 'thesis_draft/04-graphene-fet-device-physics.md'],
            capture_output=True, text=True, timeout=60)
        hist = out.stdout.strip()
    except Exception as exc:                                  # pragma: no cover
        hist = ''
        print('  git unavailable: %r' % (exc,))
    if hist:
        for ln in hist.split('\n'):
            print('  %s' % ln)
        check('  G1 the Section 4.7 annotation was committed on 2026-10-13',
              '2026-10-13' in hist, hist.split('\n')[0])
        check('  G2 and its commit message names BOTH sections',
              'Sections 4.7 and 4.10' in hist, hist.split('\n')[0])
    else:
        print('  (no history available; the in-file evidence in Section 2')
        print('   stands on its own and is the stronger of the two)')

    print('\n' + '=' * 78)
    print('TOTAL: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
