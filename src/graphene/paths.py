"""Where a script's figures and transcripts go.

One helper, ``results_path(chapter, name)``, replaces the bare relative
filenames every module used to pass to ``savefig``/``open``.  It does not
change the working directory.

Routing rule
------------
* A module that lives in this repository's ``src/`` directory writes to
  ``results/figures/<chapter>/`` (images) or ``results/transcripts/<chapter>/``
  (everything else), creating the directory if needed.
* A module that lives anywhere else -- in practice a COPY of a module that a
  mutation harness or the pristine-transcript audit has placed in a temporary
  directory -- gets the bare ``name`` back, i.e. it writes into its own
  working directory exactly as it did before the reorganisation.  This keeps
  mutated or historical runs from overwriting the committed results.
* ``GRAPHENE_RESULTS_ROOT``, if set, replaces ``<repo>/results`` for canonical
  modules (used by CI to regenerate into a scratch tree before comparing).

The caller is identified from the calling frame's ``__file__``; pass
``caller=__file__`` explicitly to override.
"""
import os
import re
import sys

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.dirname(_PKG_DIR)
REPO_ROOT = os.path.dirname(SRC_DIR)

_IMAGE_EXT = {".png", ".pdf", ".svg", ".jpg", ".jpeg", ".gif"}


def results_root():
    return os.environ.get("GRAPHENE_RESULTS_ROOT") or os.path.join(REPO_ROOT, "results")


def _is_canonical(caller):
    if not caller:
        return False
    d = os.path.dirname(os.path.realpath(caller))
    return d == os.path.realpath(SRC_DIR)


def results_path(chapter, name, caller=None):
    """Return the path a module should write output ``name`` to.

    ``chapter`` is a short label such as ``'ch4'`` or ``'legacy'``.
    """
    if caller is None:
        caller = sys._getframe(1).f_globals.get("__file__")
    if not _is_canonical(caller):
        return name
    kind = "figures" if os.path.splitext(name)[1].lower() in _IMAGE_EXT else "transcripts"
    d = os.path.join(results_root(), kind, chapter)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def committed_path(chapter, name):
    """Location of a committed output, for audits that compare against it."""
    kind = "figures" if os.path.splitext(name)[1].lower() in _IMAGE_EXT else "transcripts"
    return os.path.join(REPO_ROOT, "results", kind, chapter, name)


# --------------------------------------------------------------------------
# Locating files that lived at the repository root before the reorganisation
# (added 2026-11-16, merge of main @ e4fbc26).  Modules written on main name
# their inputs by root-relative bare names -- 'AUTOMATION_LOG.md',
# 'thesis_draft/04-...md', 'foo_output.txt', 'graphene_x.py'.  These helpers
# map such a name to where the file now lives *inside the tree rooted at
# ``root``*, so the same call works on this checkout, on a pristine clone of
# it, and on a flat pre-reorg clone (where the bare name is found first).
# --------------------------------------------------------------------------
_RENAMED_DIRS = {"thesis_draft": "thesis"}


def locate(root, name):
    """Absolute path of ``name`` (root-relative, flat-layout spelling) in ``root``."""
    flat = os.path.join(root, name)
    if os.path.exists(flat):
        return flat
    head, _, rest = name.replace("\\", "/").partition("/")
    if head in _RENAMED_DIRS:
        return os.path.join(root, _RENAMED_DIRS[head], rest)
    if name.endswith(".py") and os.path.exists(os.path.join(root, "src", name)):
        return os.path.join(root, "src", name)
    base = os.path.basename(name)
    for kind in ("transcripts", "figures"):
        d = os.path.join(root, "results", kind)
        if os.path.isdir(d):
            for ch in sorted(os.listdir(d)):
                p = os.path.join(d, ch, base)
                if os.path.exists(p):
                    return p
    return flat


def rel_locate(root, name):
    """``locate`` as a path relative to ``root`` with '/' separators (for ``git show``)."""
    return os.path.relpath(locate(root, name), root).replace(os.sep, "/")


def repo_root_of(src_dir):
    """The repository root for a module directory: its parent when the
    directory is a ``src/`` with ``graphene/`` beside the modules, else itself."""
    if os.path.basename(os.path.normpath(src_dir)) == "src" and os.path.isdir(os.path.join(src_dir, "graphene")):
        return os.path.dirname(os.path.normpath(src_dir))
    return src_dir


_RECORD_A = "appendix-A-verification-record.md"
_RECORD_B = "appendix-B-failure-taxonomy.md"


def chapter_text(path):
    """A thesis chapter as the audits that read its prose must see it: the body
    plus that chapter's part of the verification record.

    The 2026-11 restructure moved each chapter's audit annotations out of the
    body into Appendix A (section A.<n+1> for chapter n) and Chapter 7's
    Sections 7.5-7.7 into Appendix B. Nothing was deleted, so a check of "is
    this claim in the chapter" must read both. A file with no NN- prefix, or a
    directory without the appendices (a planted stub, a pre-reorg clone), reads
    as the file alone."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    d, name = os.path.split(path)
    m = re.match(r"(\d\d)-", name)
    if not m:
        return text
    n = int(m.group(1))
    rec = os.path.join(d, _RECORD_A)
    if os.path.exists(rec):
        with open(rec, encoding="utf-8") as fh:
            sec = re.search(r"^## A\.%d .*?(?=^## A\.|\Z)" % (n + 1), fh.read(), re.M | re.S)
        if sec:
            text += "\n" + sec.group(0).split("\n", 1)[-1]
    rec = os.path.join(d, _RECORD_B)
    if n == 7 and os.path.exists(rec):
        with open(rec, encoding="utf-8") as fh:
            text += "\n" + fh.read().split("\n", 1)[-1]
    return text
