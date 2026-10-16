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
