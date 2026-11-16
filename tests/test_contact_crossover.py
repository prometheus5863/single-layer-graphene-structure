"""Exact-value tests for the contact-doping crossover (correction of record, 2026-11-16).

The contact doping is measured from the n/p crossover W0 = 5.4 eV, not from
graphene's own 4.5 eV (thesis §4.7, Appendix A corrections table). These are
tests against values the model must reproduce EXACTLY, not plausible ranges.

    python tests/test_contact_crossover.py        # or: python -m pytest tests/
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
os.environ.setdefault("MPLBACKEND", "Agg")

import graphene_contact_doping_model as cd  # noqa: E402

N_BULK = 2.0e16  # m^-2, the on-state channel density §4.7 uses


def test_crossover_is_the_literature_value():
    assert cd.W_CROSSOVER == 5.4
    assert cd.W_GRAPHENE == 4.5  # graphene's own work function is unchanged (Chapter 6 uses it)


def test_metal_at_the_crossover_leaves_graphene_undoped():
    # W = W0 -> dWF = 0 -> dV = 0 -> n_contact = 0 exactly (not merely small).
    assert cd.contact_edge_carrier_density(cd.W_CROSSOVER) == 0.0


def test_sign_follows_the_crossover():
    for d in (0.05, 0.3, 0.75):
        assert cd.contact_edge_carrier_density(cd.W_CROSSOVER - d) > 0  # below W0: n-type
        assert cd.contact_edge_carrier_density(cd.W_CROSSOVER + d) < 0  # above W0: p-type


def test_physisorbed_metals_dope_with_their_literature_sign():
    # Cu (4.65 eV) n-dopes and Pt (5.65 eV) p-dopes graphene, as Giovannetti et al. find.
    assert cd.contact_edge_carrier_density(cd.METAL_WORK_FUNCTIONS["Cu"]) > 0
    assert cd.contact_edge_carrier_density(cd.METAL_WORK_FUNCTIONS["Pt"]) < 0


def test_no_junction_when_contact_matches_the_channel():
    # The model's exact zero of R_extra: n_contact == n_bulk gives a flat profile.
    R, _, n_x = cd.junction_extra_resistance_from_ncontact(N_BULK, N_BULK)
    assert R == 0.0
    assert (n_x == N_BULK).all()


def test_neutral_contact_resistance_is_bounded_by_the_puddle_floor():
    # At W = W0 the contact pins the edge density at zero, so the profile rises
    # from 0 to n_bulk and R_extra is positive -- but finite, because the
    # puddle-regularised conductivity never falls below sigma_min. It cannot
    # exceed the resistance of the whole junction length at sigma_min.
    R, _, _, n_c = cd.junction_extra_resistance(cd.W_CROSSOVER, N_BULK)
    assert n_c == 0.0
    sigma_min = float(cd._sheet_conductivity_puddle_regularized(0.0))
    bound_ohm_um = 1e-6 / sigma_min * 1e6  # L_junction = 1 um, per um of width
    assert 0.0 < R < bound_ohm_um


def test_model_scope_is_declared():
    assert set(cd.PHYSISORBED_METALS) == {"Cu", "Au", "Pt"}
    for m in ("Ni", "Pd", "Ti"):
        assert not cd.in_model(m)
        assert cd.OUTSIDE_MODEL_METALS[m] == "chemisorbed"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print("PASS", t.__name__)
        except AssertionError as exc:
            failed += 1
            print("FAIL", t.__name__, exc)
    print("%d/%d passed" % (len(tests) - failed, len(tests)))
    sys.exit(1 if failed else 0)
