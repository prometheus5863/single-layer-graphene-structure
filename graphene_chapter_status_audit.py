"""
graphene_chapter_status_audit.py

A GUARD AGAINST A CORRECTION THAT DID NOT HOLD.

On 2026-10-18 this repository audited its own standing list and found, as item
A1, that the list's claim

    "CHAPTERS 2 AND 3 ARE STILL UNDRAFTED -- their computational results have
     been complete for weeks.  This is the largest remaining block of pure
     writing"

was false: Chapter 2 was drafted 2026-09-28 (`5bbd5f3`, eleven sections) and
Chapter 3 on 2026-09-29 (`c1d5abc`, ten sections).  The 10-18 entry's own title
records the mechanism -- "one of them came from outside the repository".

THE CLAIM THEN CAME BACK.  It is absent from 10-18 through 10-25 and present
again in the "Not yet covered" list of 2026-10-26, 2026-10-27 and 2026-10-28:
three consecutive sessions re-asserting a claim this repository had already
refuted in writing, against files in its own tree.

That is a different kind of defect from the ones this repository usually finds.
The others are claims nothing checked.  This one WAS checked, and the check did
not survive, because the correction lived in a log entry and log entries are
not read by anything.  The authority for the claim is outside the tree, so
every session that reconstructs the list from outside re-imports it, and no
amount of correcting it in prose prevents the next re-import.

A prose correction cannot fix a one-way channel.  A committed check can, which
is what this module is: it measures each chapter's drafted state FROM THE FILE,
and refuses a "Not yet covered" list that describes a measured-drafted chapter
as undrafted.  Same shape as the repair pattern named on 2026-10-28 -- a
registry entry carries its reason, and a check confirms the reason is still
true of the source -- applied for the first time to a claim whose source is a
markdown file rather than a Python name.

Two planted controls, per the standing requirement that every verdict have
one: Section 3 plants a fabricated "Chapter 5 ... still undrafted" sentence and
requires the detector to FLAG it, and Section 4 writes a genuine stub chapter
under TMPDIR and requires the classifier to call it UNDRAFTED.  Without both,
a detector that never fires and a classifier that calls everything drafted
would each read as green.
"""

import os
import re
import sys
import tempfile

LOG = 'AUTOMATION_LOG.md'
CHAPTER_DIR = 'thesis_draft'

# A chapter counts as DRAFTED when it has real sections and real prose.  The
# thresholds are set well below every drafted chapter in the tree and well
# above a stub; Section 1 prints the margin on both sides so that a chapter
# drifting toward either threshold is visible rather than silently reclassified.
#
# THE FIRST VERSION OF THIS METRIC WAS WRONG AND IS RECORDED RATHER THAN
# QUIETLY REPLACED.  It counted LINES of at least 40 characters and required 5
# sections, and on its first run it classified `01-introduction.md` UNDRAFTED:
# 3 sections, 97 such lines.  Chapter 1 is 5340 words -- MORE than Chapters 2,
# 3 and 5 -- written in 120 unwrapped lines averaging ~314 characters, and it
# has three `##` sections because an introduction legitimately has three.  So
# the metric was measuring the author's line-wrapping and the section count was
# calibrated on the long chapters.  Words outside code fences are
# wrapping-independent, and Section 5 checks that property by re-wrapping a
# real chapter and requiring the measurement not to move.
MIN_SECTIONS = 3
MIN_PROSE_WORDS = 1500

# The registry.  Each entry carries the claim and the reason the claim is made,
# and Section 2 confirms the reason is still true of the file.
CLAIMED = {
    '01-introduction.md': ('DRAFTED', 'drafted before the log begins'),
    '02-electronic-structure.md': ('DRAFTED', 'drafted 2026-09-28, 5bbd5f3, eleven sections'),
    '03-transport-and-optical-properties.md': ('DRAFTED', 'drafted 2026-09-29, c1d5abc, ten sections'),
    '04-graphene-fet-device-physics.md': ('DRAFTED', 'the chapter the device work lands in'),
    '05-graphene-interconnects.md': ('DRAFTED', 'drafted with Chapter 4'),
    '06-graphene-photodetectors.md': ('DRAFTED', 'drafted through Section 6.8'),
    '07-discussion-and-outlook.md': ('DRAFTED', 'drafted; carries the standing item list'),
}

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
# Section 1: measure a chapter
# ---------------------------------------------------------------------------

def measure_text(text):
    """(n_sections, n_prose_words) for chapter text, code fences excluded."""
    n_sections = 0
    n_words = 0
    in_fence = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith('```'):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if stripped.startswith('## '):
            n_sections += 1
        if not stripped.startswith('#'):
            n_words += len(stripped.split())
    return n_sections, n_words


def measure(path):
    with open(path, encoding='utf-8') as fh:
        return measure_text(fh.read())


def classify(path):
    n_sections, n_words = measure(path)
    state = ('DRAFTED' if (n_sections >= MIN_SECTIONS and n_words >= MIN_PROSE_WORDS)
             else 'UNDRAFTED')
    return state, n_sections, n_words


# ---------------------------------------------------------------------------
# Section 2: the "Not yet covered" detector
# ---------------------------------------------------------------------------

CH_WORD = r'(?:chapters?|chs?\.?)'
UNDRAFTED_WORD = r'(?:undrafted|not (?:yet )?drafted|still to be (?:drafted|written)|unwritten)'

def latest_entry(log_text):
    """The text of the most recent '## 2026-..-..' entry in the log."""
    starts = [m.start() for m in re.finditer(r'(?m)^## 20\d\d-\d\d-\d\d', log_text)]
    if not starts:
        return ''
    return log_text[starts[-1]:]


def chapter_numbers_in(fragment):
    """Chapter numbers named in a fragment of text, e.g. 'Chapters 2-3' -> {2,3}."""
    nums = set()
    for m in re.finditer(r'(\d+)\s*(?:[-‒–—]|and|,|&)\s*(\d+)', fragment):
        a, b = int(m.group(1)), int(m.group(2))
        if 1 <= a <= 12 and 1 <= b <= 12 and a <= b:
            nums.update(range(a, b + 1))
    for m in re.finditer(r'\b(\d{1,2})\b', fragment):
        n = int(m.group(1))
        if 1 <= n <= 12:
            nums.add(n)
    return nums


def undrafted_claims(entry_text):
    """
    Chapter numbers an entry claims are undrafted.

    A sentence qualifies when a chapter word and an undrafted word appear within
    the same sentence, in either order -- the 10-26..10-28 wording puts the
    chapter first ("Chapters 2-3 are still undrafted") and the 10-18 quotation
    of it does too, but a future rewording must not slip past.
    """
    found = {}
    for sentence in re.split(r'(?<=[.;\n])\s+', entry_text):
        low = sentence.lower()
        if re.search(CH_WORD, low) and re.search(UNDRAFTED_WORD, low):
            # Skip a sentence that is explicitly REPORTING the false claim.
            if re.search(r'(?:was false|refuted|had been drafted|are drafted|'
                         r'is false|no longer|stale|regress)', low):
                continue
            m = re.search(CH_WORD + r'[^.;\n]*', low)
            for n in chapter_numbers_in(m.group(0) if m else low):
                found.setdefault(n, sentence.strip()[:140])
    return found


def chapter_file_for(number):
    prefix = '%02d-' % number
    for name in CLAIMED:
        if name.startswith(prefix):
            return name
    return None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('CHAPTER STATUS AUDIT -- does the standing list agree with the tree?')
    print('=' * 78)

    print('\nSECTION 1.  Measured state of every chapter, against its registry claim')
    print('  %-44s %-10s %-9s %-7s %-8s'
          % ('file', 'claimed', 'measured', '## secs', 'words'))
    measured = {}
    for name, (claim, reason) in sorted(CLAIMED.items()):
        path = os.path.join(CHAPTER_DIR, name)
        if not os.path.exists(path):
            check('registry names a file that exists: %s' % name, False,
                  'missing from %s/' % CHAPTER_DIR)
            continue
        state, n_sec, n_prose = classify(path)
        measured[name] = state
        print('  %-44s %-10s %-9s %-7d %-8d'
              % (name, claim, state, n_sec, n_prose))
        check('registry claim matches the file: %s' % name, state == claim,
              'claimed %s because "%s"; measured %s (%d sections, %d prose words; '
              'thresholds %d / %d)'
              % (claim, reason, state, n_sec, n_prose, MIN_SECTIONS, MIN_PROSE_WORDS))

    print('\n  Margin to the DRAFTED thresholds, so a drift is visible:')
    for name in sorted(measured):
        n_sec, n_prose = measure(os.path.join(CHAPTER_DIR, name))
        print('    %-44s sections +%-4d words +%-6d'
              % (name, n_sec - MIN_SECTIONS, n_prose - MIN_PROSE_WORDS))

    print('\nSECTION 2.  The live check: does the latest log entry contradict the tree?')
    log_text = open(LOG, encoding='utf-8').read()
    entry = latest_entry(log_text)
    header = entry.split('\n', 1)[0] if entry else '(none)'
    print('  latest entry: %s' % header)
    claims = undrafted_claims(entry)
    if not claims:
        print('  no chapter is described as undrafted in that entry.')
    contradictions = []
    for number, sentence in sorted(claims.items()):
        name = chapter_file_for(number)
        state = measured.get(name)
        verdict = 'CONTRADICTS TREE' if state == 'DRAFTED' else 'consistent'
        print('  Chapter %-2d -> %-44s %-9s %s'
              % (number, name or '(no such chapter)', state or 'n/a', verdict))
        print('      "%s"' % sentence)
        if state == 'DRAFTED':
            contradictions.append(number)
    check('the latest log entry describes no measured-DRAFTED chapter as undrafted',
          not contradictions,
          'contradicted chapters: %s.  A prose correction was already made on '
          '2026-10-18 (item A1) and did not hold; this check is the committed '
          'form of it.' % (sorted(contradictions) if contradictions else 'none'))

    print('\nSECTION 3.  Planted control: the detector must FIRE on a fabricated claim')
    planted = (entry if entry else '## 2026-01-01 stub\n') + \
        '\n- **Chapter 5 is still undrafted** and is the largest remaining block.\n'
    planted_claims = undrafted_claims(planted)
    check('detector flags a planted "Chapter 5 ... still undrafted"',
          5 in planted_claims,
          'planted-claim chapters detected: %s' % sorted(planted_claims))
    check('planted control does not perturb the real verdict',
          set(planted_claims) - {5} == set(claims),
          'real %s, planted-minus-5 %s'
          % (sorted(claims), sorted(set(planted_claims) - {5})))

    print('\nSECTION 4.  Planted control: the classifier must call a stub UNDRAFTED')
    tmpdir = tempfile.mkdtemp(prefix='chapstatus_', dir=os.environ.get('TMPDIR', '/tmp'))
    stub = os.path.join(tmpdir, '08-stub.md')
    with open(stub, 'w', encoding='utf-8') as fh:
        fh.write('# Chapter 8: Stub\n\n## 8.1 To be written\n\nTBD.\n')
    stub_state, s_sec, s_prose = classify(stub)
    check('a stub chapter classifies as UNDRAFTED', stub_state == 'UNDRAFTED',
          '%d sections, %d prose words -> %s' % (s_sec, s_prose, stub_state))
    print('         (written under TMPDIR, not in the tree, per the standing rule')
    print('          that a module must not plant its pattern where other audits read)')

    print('\nSECTION 5.  The metric must not measure the line-wrapping')
    print('  The first version of this metric counted long LINES and called')
    print('  Chapter 1 undrafted because it is written unwrapped.  Re-wrap a real')
    print('  chapter two ways and require the measurement not to move.')
    victim = os.path.join(CHAPTER_DIR, '01-introduction.md')
    original = open(victim, encoding='utf-8').read()
    base = measure_text(original)
    # (a) collapse every paragraph onto one line; (b) hard-wrap at 60 columns.
    collapsed = re.sub(r'(?<!\n)\n(?![\n#`\-|*0-9])', ' ', original)
    import textwrap
    rewrapped_lines = []
    for line in original.splitlines():
        if line.strip().startswith(('#', '`', '-', '|', '*')) or not line.strip():
            rewrapped_lines.append(line)
        else:
            # break_on_hyphens and break_long_words must BOTH be off: with the
            # defaults, textwrap splits "device-physics" across lines and the
            # variant then has one word more than the original.  The first run
            # of this check failed by exactly +1 word for that reason -- a
            # defect in the re-wrapping harness, not in the metric, and fixed
            # here rather than absorbed into a tolerance.
            rewrapped_lines.extend(
                textwrap.wrap(line, 60, break_on_hyphens=False,
                              break_long_words=False) or [''])
    rewrapped = '\n'.join(rewrapped_lines)
    for label, variant in (('paragraphs collapsed to one line', collapsed),
                           ('hard-wrapped at 60 columns', rewrapped)):
        got = measure_text(variant)
        check('word count is wrapping-independent (%s)' % label,
              got[1] == base[1],
              'original %d sections / %d words, variant %d / %d'
              % (base[0], base[1], got[0], got[1]))

    print('\n' + '=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
