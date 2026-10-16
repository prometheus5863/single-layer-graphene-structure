"""Support package for the single-layer graphene device models.

The physical models themselves are still flat modules in ``src/`` (see
docs/INVENTORY.md, "Layout decision"); this package currently holds only
infrastructure that must not be scanned by the repository's own audits.
"""
from .paths import results_path, committed_path  # noqa: F401
