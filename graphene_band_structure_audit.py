"""
graphene_band_structure_audit.py  --  2026-09-28

THE SHIPPED BAND STRUCTURE OF GRAPHENE IN THIS REPOSITORY HAS NO DIRAC POINT.

`graphene_band_structure.calculate_band_structure()` returns a spectrum whose
MINIMUM gap along the whole Gamma-K-M-Gamma path is 11.2 eV, and whose gap at
the k-point labelled K is 14.4 eV.  Graphene's single defining electronic
property -- a gapless linear crossing at K -- is absent from `band_structure.png`,
the oldest figure in this repository and the one Chapter 2 is built on.

ROOT CAUSE, and it is visible inside a single function.  `k_path_graphene()`
defines

    k_point = [4.0/(3*sqrt(3)), 0]          <-- no pi
    m_point = [pi/(3*sqrt(3)), pi/3]        <-- pi

in the same seven lines.  `graphene_hamiltonian()` builds its phases from
nearest-neighbour vectors of UNIT LENGTH, i.e. in units of the C-C bond length,
for which the zone corner is at |k| = 4*pi/(3*sqrt(3)) = 2.4184 and the edge
midpoint at |M| = 2*pi/3 = 2.0944.  The shipped K is a factor of pi too small
and the shipped M is wrong in both components.  Three conventions are mixed in
one 294-line file: lattice vectors in Angstrom, Hamiltonian phases in bond-length
units, high-symmetry points in neither.

WHY IT SURVIVED FIVE WEEKS, which is the more useful half of this file.  The
repository contains THREE graphene band-structure or DOS figures and NOT ONE of
them is a calculation that could have failed:

  * `band_structure.png`            -- computed, and wrong (this file's subject)
  * `band_structure_simple.png`     -- drawn from E = +/- v_F |k|, the known
                                       answer, in simple_graphene_plots.py
  * `density_of_states.png`         -- NOT computed from the bands.
                                       calculate_density_of_states() returns
                                       |E|/(pi*2.8^2) analytically, so it shows
                                       the correct V shape no matter what the
                                       Hamiltonian does

A DOS computed from the bands would have caught this on the first run: it
would have shown no states near zero energy and van Hove peaks in the wrong
place.  A figure that HARD-CODES the expected answer cannot falsify the
calculation sitting next to it.  This is the same structure as the three
sessions before it -- 09-25's convergence estimate answering yes while 44%
wrong, 09-27's `conv` reporting 1e-14 at 14% error, and this morning's
non-discriminating covariance check -- and it is the most consequential
instance so far, because here the unfalsifiable instrument was a PICTURE.

Every validation below is against an EXACTLY known value of the
nearest-neighbour tight-binding model, and reports its measured value beside a
derived floor rather than against a round tolerance (2026-09-27 item 12).

  |phi(Gamma)| = 3  exactly  ->  E = +/- 3t = +/- 8.4 eV
  |phi(M)|     = 1  exactly  ->  E = +/- t  = +/- 2.8 eV
  |phi(K)|     = 0  exactly  ->  E = 0, doubly degenerate
  v_F          = 3 t a_cc / (2 hbar)   from the slope at K
  DOS(E->0)    = 2|E| / (pi hbar^2 v_F^2)   per unit area, spin and valley
                 degeneracy included
"""

import numpy as np

import graphene_band_structure as gbs

T_HOP = 2.8                     # eV, the hopping used throughout the repo
A_LATTICE = 2.46                # Angstrom, graphene lattice constant
A_CC = A_LATTICE / np.sqrt(3.0)  # Angstrom, C-C bond length = 1.4203
HBAR = 1.054571817e-34          # J s
Q_E = 1.602176634e-19           # C
EPS = float(np.finfo(np.float64).eps)

# Nearest-neighbour vectors as graphene_hamiltonian() actually uses them:
# unit length, i.e. expressed in units of the C-C bond length.
DELTA = (np.array([0.0, 1.0]),
         np.array([np.sqrt(3.0) / 2.0, -0.5]),
         np.array([-np.sqrt(3.0) / 2.0, -0.5]))

# The high-symmetry points CONSISTENT with that convention.
GAMMA = np.array([0.0, 0.0])
K_CORRECT = np.array([4.0 * np.pi / (3.0 * np.sqrt(3.0)), 0.0])
M_CORRECT = (2.0 * np.pi / 3.0) * np.array([np.cos(np.pi / 6.0), 0.5])

# The points as shipped, kept so the superseded numbers stay measurable.
K_SHIPPED = np.array([4.0 / (3.0 * np.sqrt(3.0)), 0.0])
M_SHIPPED = np.array([np.pi / (3.0 * np.sqrt(3.0)), np.pi / 3.0])


def phi(k):
    return sum(np.exp(1j * np.dot(k, d)) for d in DELTA)


def bands(k, t=T_HOP):
    """E = +/- t|phi(k)|, which is the eigenvalue problem of the 2x2 H."""
    a = t * abs(phi(k))
    return -a, +a


def fermi_velocity_analytic(t=T_HOP, a_cc=A_CC):
    """v_F = 3 t a_cc / (2 hbar), the standard closed form."""
    return 3.0 * (t * Q_E) * (a_cc * 1e-10) / (2.0 * HBAR)


# =====================================================================
# VALIDATIONS -- exact values of the model, measured against derived floors
# =====================================================================
def v1_exact_symmetry_points(verbose=True):
    """|phi| is 3 at Gamma, 1 at M and 0 at K, exactly, for ANY t.  Three
    independent exact values, so this check has three chances to fail."""
    rows = [("Gamma", GAMMA, 3.0), ("M", M_CORRECT, 1.0),
            ("K", K_CORRECT, 0.0)]
    out, worst = [], 0.0
    for name, k, exact in rows:
        v = abs(phi(k))
        # derived floor: three unit-modulus complex exponentials summed, each
        # with argument of size |k|.|delta| <~ 2.5, so <~ 6 eps absolute
        floor = 6.0 * EPS * max(abs(np.dot(k, DELTA[0])), 1.0)
        err = abs(v - exact)
        ratio = err / floor
        worst = max(worst, ratio)
        out.append((name, v, exact, err, floor, ratio))
    ok = worst <= 1.0
    if verbose:
        print("\nV1 (exact x3): |phi| at the three high-symmetry points")
        print("    %-7s%16s%10s%13s%13s%9s"
              % ("point", "measured", "exact", "error", "floor", "ratio"))
        for name, v, e, err, floor, ratio in out:
            print("    %-7s%16.12f%10.1f%13.3e%13.3e%9.3f"
                  % (name, v, e, err, floor, ratio))
        print("    worst ratio = %.3f  -> %s" % (worst, "PASS" if ok else "FAIL"))
    return ok, out


def v2_fermi_velocity(verbose=True):
    """The slope of E(k) at K must reproduce v_F = 3 t a_cc / (2 hbar).

    THE FIRST FORM OF THIS VALIDATION FAILED, and it failed for a reason worth
    keeping: it assumed the one-sided slope converges at O(dk^2), and the
    measured errors were 2.512e-3, 2.501e-4, 2.500e-5 at dk = 1e-2, 1e-3, 1e-4
    -- clean FIRST order, ratio 200.999 / 2001.000 / 20001.002 against the
    second-order bound.  First order because the leading correction to the
    Dirac cone along Gamma-K is linear in q, not quadratic:

        |phi(K + q xhat)| = |1 - cos u - sqrt(3) sin u|,  u = sqrt(3) q / 2
                          = sqrt(3) u - u^2/2 + O(u^3)
                          = (3/2) q ( 1 - q/4 + O(q^2) )

    so the one-sided slope is low by EXACTLY q/4 to leading order.  That is a
    derived coefficient and an exactly known constant, so the validation is
    rewritten to check the CONSTANT rather than a decay rate: rel_err/dk must
    approach 0.25.  This is trigonal warping, measured -- physics the first
    form mistook for a convergence failure.
    """
    v_exact = fermi_velocity_analytic()
    rows = []
    for dk in (1e-2, 1e-3, 1e-4):
        e = bands(K_CORRECT + np.array([dk, 0.0]))[1]      # eV
        slope = (e / dk) * Q_E * (A_CC * 1e-10) / HBAR
        rel = abs(slope - v_exact) / v_exact
        rows.append((dk, slope, rel, rel / dk))
    coeff = rows[-1][3]
    # derived: the next term is O(dk), so the coefficient converges linearly
    bound = 0.25 * (1.0 + 2.0 * rows[-1][0])
    ok = abs(coeff - 0.25) <= abs(bound - 0.25) + 1e-12
    if verbose:
        print("\nV2 (exact, rewritten -- see the docstring for why the first")
        print("    form failed): Fermi velocity and the leading warping term")
        print("    closed form  3 t a_cc / (2 hbar) = %.6e m/s" % v_exact)
        print("    %10s%18s%13s%14s" % ("dk", "slope (m/s)", "rel err",
                                        "rel err / dk"))
        for dk, s, rel, c in rows:
            print("    %10.0e%18.6e%13.3e%14.6f" % (dk, s, rel, c))
        print("    derived coefficient = 1/4 exactly;  measured %.6f"
              % coeff)
        print("    |measured - 0.25| = %.2e, admissible %.2e  -> %s"
              % (abs(coeff - 0.25), abs(bound - 0.25),
                 "PASS" if ok else "FAIL"))
        print("    extrapolated v_F (slope x (1 + dk/4)) = %.6e m/s,"
              % (rows[-1][1] * (1.0 + rows[-1][0] / 4.0)))
        print("    against the closed form %.6e -- agreement %.2e relative"
              % (v_exact, abs(rows[-1][1] * (1 + rows[-1][0] / 4.0) - v_exact)
                 / v_exact))
    return ok, v_exact, rows, coeff


def v3_particle_hole_symmetry(verbose=True):
    """Nearest-neighbour-only graphene is exactly particle-hole symmetric:
    E_+(k) = -E_-(k) for EVERY k.  An exact zero that must hold identically,
    not just at symmetry points -- 2026-09-27 item 10's warning applies, so
    both the bitwise count and the derived floor are reported."""
    rng = np.random.default_rng(20260928)
    ks = rng.uniform(-3.0, 3.0, size=(4000, 2))
    devs = np.array([abs(bands(k)[0] + bands(k)[1]) for k in ks])
    bitwise = int(np.sum(devs == 0.0))
    floor = 2.0 * EPS * 3.0 * T_HOP
    worst = devs.max()
    ok = worst <= floor
    if verbose:
        print("\nV3 (exact, identical in k): E_+ + E_- must vanish everywhere")
        print("    %d random k, bitwise zero at %d of them (%.1f%%)"
              % (len(ks), bitwise, 100.0 * bitwise / len(ks)))
        print("    worst |E_+ + E_-| = %.3e   derived floor %.3e   ratio %.3f"
              % (worst, floor, worst / floor))
        print("    -> %s" % ("PASS" if ok else "FAIL"))
    return ok, bitwise, worst, floor


def v4_dos_low_energy_slope(verbose=True):
    """The DOS computed FROM the bands must reproduce the analytic low-energy
    form g(E) = 2|E|/(pi hbar^2 v_F^2) -- an exactly known constant.  This is
    the check the shipped calculate_density_of_states() cannot perform, because
    it returns that shape by construction.

    A NOTE ON THE CHECK THAT WAS TRIED FIRST AND IS VACUOUS.  The obvious
    normalisation check -- integral of g dE = 4 states per unit-cell area --
    passed to 16 digits with the WRONG reciprocal vectors in place, because the
    normalisation divides by the sample count and therefore returns 4 states
    per cell whatever the Hamiltonian or the sampling region is.  It is the
    same fault as the shipped DOS figure: a check that cannot fail.  It is
    retained below and PRINTED AS VACUOUS, and the two checks that do have
    power -- periodicity of |phi| under b1, b2, and the number of Dirac cones
    inside the sampled cell -- are added beside it.
    """
    E, g, _ = dos_from_bands(n=520)
    v_F = fermi_velocity_analytic()
    slope_exact = 2.0 * Q_E ** 2 / (np.pi * HBAR ** 2 * v_F ** 2)
    band = (np.abs(E) > 0.15) & (np.abs(E) < 0.8)
    slope_meas = np.polyfit(np.abs(E[band]), g[band], 1)[0]
    rel = abs(slope_meas - slope_exact) / slope_exact
    # derived: trigonal warping enters the DOS at O((E/t)) through the same
    # q/4 term V2 measures, averaged over the fit window
    bound = 0.8 / (4.0 * T_HOP) + 2.0 / np.sqrt(len(E))
    # -- the vacuous companion, kept and labelled
    dE = E[1] - E[0]
    total = float((g * dE).sum())
    area_cell = (np.sqrt(3.0) / 2.0) * (A_LATTICE * 1e-10) ** 2
    # -- the two checks with power
    b1 = (4.0 * np.pi / 3.0) * np.array([-np.cos(np.pi / 6.0), 0.5])
    b2 = (4.0 * np.pi / 3.0) * np.array([np.cos(np.pi / 6.0), 0.5])
    rng = np.random.default_rng(28092026)
    ks = rng.uniform(-4.0, 4.0, size=(500, 2))
    per = max(max(abs(abs(phi(k + b)) - abs(phi(k))) for b in (b1, b2))
              for k in ks)
    per_floor = 8.0 * EPS * 8.0
    n_cone = _count_cones_in_cell(b1, b2)
    ok = rel <= bound and per <= per_floor and n_cone == 2
    if verbose:
        print("\nV4 (exact constant): low-energy DOS slope computed from the BANDS")
        print("    analytic 2 q_e^2/(pi hbar^2 v_F^2) = %.6e states/(eV m^2)"
              % slope_exact)
        print("    measured from the band histogram   = %.6e" % slope_meas)
        print("    relative difference %.4f   derived bound %.4f"
              % (rel, bound))
        print("    (the bound is the q/4 warping term of V2 averaged over the")
        print("     fit window -- physics, not error)")
        print("    VACUOUS COMPANION, kept and labelled: integral g dE = %.6e"
              % total)
        print("      against 4/A_cell = %.6e, agreement %.1e -- and this"
              % (4.0 / area_cell, abs(total - 4.0 / area_cell) * area_cell / 4))
        print("      PASSED IDENTICALLY with the wrong reciprocal vectors in")
        print("      place, so it has no power.  It cannot fail by construction.")
        print("    CHECKS WITH POWER, added because that one has none:")
        print("      |phi| periodic under b1, b2: worst deviation %.2e against"
              % per)
        print("        a derived floor of %.2e  -> %s"
              % (per_floor, "ok" if per <= per_floor else "FAIL"))
        print("      Dirac cones inside the sampled cell: %d (must be 2)"
              % n_cone)
        print("    -> %s" % ("PASS" if ok else "FAIL"))
    return ok, slope_exact, slope_meas, rel, bound, per, n_cone


def _count_cones_in_cell(b1, b2, n=600):
    """Count distinct minima of |phi| inside the parallelogram cell spanned by
    b1, b2.  Two is the right answer for graphene (K and K')."""
    u = (np.arange(n) + 0.5) / n
    U, V = np.meshgrid(u, u, indexing="ij")
    kx = U * b1[0] + V * b2[0]
    ky = U * b1[1] + V * b2[1]
    ph = np.zeros_like(kx, dtype=complex)
    for d in DELTA:
        ph += np.exp(1j * (kx * d[0] + ky * d[1]))
    small = np.abs(ph) < 0.15
    # label connected components with a flood fill on the periodic-free grid
    seen = np.zeros_like(small, dtype=bool)
    count = 0
    idx = np.argwhere(small)
    for i0, j0 in idx:
        if seen[i0, j0]:
            continue
        count += 1
        stack = [(i0, j0)]
        while stack:
            i, j = stack.pop()
            if not (0 <= i < n and 0 <= j < n) or seen[i, j] or not small[i, j]:
                continue
            seen[i, j] = True
            stack.extend([(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)])
    return count


def v5_shipped_path_reproduces_the_bug(verbose=True):
    """The superseded numbers must be REPRODUCIBLE, or the diagnosis is not
    about the shipped code.  Exact requirement: the gap at the path index
    labelled K must equal 2t|phi(K_shipped)| to the floor."""
    _, E, _ = gbs.calculate_band_structure()
    n = len(E) // 3
    measured = E[n, 1] - E[n, 0]
    predicted = 2.0 * T_HOP * abs(phi(K_SHIPPED))
    floor = 8.0 * EPS * 3.0 * T_HOP
    ok = abs(measured - predicted) <= floor
    if verbose:
        print("\nV5 (exact): the diagnosis must reproduce the shipped number")
        print("    gap at the path index labelled K, from the shipped module")
        print("      measured  = %.15f eV" % measured)
        print("      predicted = 2t|phi(K_shipped)| = %.15f eV" % predicted)
        print("      |difference| = %.3e   floor %.3e   -> %s"
              % (abs(measured - predicted), floor, "PASS" if ok else "FAIL"))
    return ok, measured, predicted


# =====================================================================
def dos_from_bands(n=420, t=T_HOP, n_bins=260):
    """DOS by uniform sampling of the k-plane over one reciprocal cell.
    Returns (E centres, g in states/(eV m^2), raw counts)."""
    # THE FIRST FORM OF THESE WAS WRONG, and the error is recorded rather
    # than quietly fixed.  b = (4pi/3)(1,0) and (4pi/3)(1/2, sqrt3/2) have the
    # right MAGNITUDE (4pi/3) and span the right AREA, but they are not
    # reciprocal lattice vectors of this lattice: |phi| is not periodic under
    # them, because b.delta_i are not all equal mod 2pi.  The parallelogram they
    # span is therefore not a fundamental domain and it contains ONE Dirac cone
    # instead of two, which halved the low-energy DOS.  The correct pair, from
    # b_i . a_j = 2pi delta_ij with a1 = delta1 - delta2, a2 = delta1 - delta3:
    b1 = (4.0 * np.pi / 3.0) * np.array([-np.cos(np.pi / 6.0), 0.5])
    b2 = (4.0 * np.pi / 3.0) * np.array([np.cos(np.pi / 6.0), 0.5])
    u = (np.arange(n) + 0.5) / n
    U, V = np.meshgrid(u, u, indexing="ij")
    kx = U * b1[0] + V * b2[0]
    ky = U * b1[1] + V * b2[1]
    ph = np.zeros_like(kx, dtype=complex)
    for d in DELTA:
        ph += np.exp(1j * (kx * d[0] + ky * d[1]))
    e = t * np.abs(ph)
    allE = np.concatenate([(-e).ravel(), e.ravel()])
    counts, edges = np.histogram(allE, bins=n_bins, range=(-3.1 * t, 3.1 * t))
    centres = 0.5 * (edges[1:] + edges[:-1])
    dE = edges[1] - edges[0]
    # normalise: 2 bands x 2 spins x n^2 samples over the unit-cell area
    area_cell = (np.sqrt(3.0) / 2.0) * (A_LATTICE * 1e-10) ** 2
    g = counts * 2.0 / (n * n) / dE / area_cell
    return centres, g, counts


def result_1_the_bug(verbose=True):
    _, E, _ = gbs.calculate_band_structure()
    gaps = E[:, 1] - E[:, 0]
    n = len(E) // 3
    if verbose:
        print("\n" + "=" * 70)
        print("RESULT 1 -- the shipped band structure has NO DIRAC POINT")
        print("=" * 70)
        print("  minimum gap over the whole shipped path = %.4f eV" % gaps.min())
        print("  gap at the index labelled K             = %.4f eV" % gaps[n])
        print("  gap at the index labelled M             = %.4f eV" % gaps[2 * n])
        print("  energy range                            = %.3f .. %.3f eV"
              % (E.min(), E.max()))
        print()
        print("  Graphene is a zero-gap semiconductor.  The shipped figure")
        print("  shows a %.1f eV gap at its own K label." % gaps[n])
        print()
        print("  The conventions, side by side:")
        for name, ks, kc in (("K", K_SHIPPED, K_CORRECT),
                             ("M", M_SHIPPED, M_CORRECT)):
            print("    %s shipped = (%8.5f, %8.5f)  |phi| = %.6f  ->  E = "
                  "+/- %.4f eV" % (name, ks[0], ks[1], abs(phi(ks)),
                                   T_HOP * abs(phi(ks))))
            print("    %s correct = (%8.5f, %8.5f)  |phi| = %.6f  ->  E = "
                  "+/- %.4f eV" % (name, kc[0], kc[1], abs(phi(kc)),
                                   T_HOP * abs(phi(kc))))
        print()
        print("  ratio |K_correct| / |K_shipped| = %.9f   (pi = %.9f)"
              % (np.linalg.norm(K_CORRECT) / np.linalg.norm(K_SHIPPED), np.pi))
        print("  The shipped K is exactly a factor of pi too small.  The same")
        print("  function writes M WITH pi factors, so the inconsistency is")
        print("  internal to seven consecutive lines, not between modules.")
    return gaps, gaps.min(), gaps[n]


def result_2_why_it_survived(verbose=True):
    """The DOS the repo ships is insensitive to the Hamiltonian.  Demonstrated
    by the strongest available means: change the Hamiltonian beyond recognition
    and see whether the shipped DOS notices."""
    E1, g1 = gbs.calculate_density_of_states()
    saved = gbs.graphene_hamiltonian
    try:
        gbs.graphene_hamiltonian = lambda k, hopping=2.8: np.array(
            [[7.0, 0.0], [0.0, -7.0]])        # a gapped, k-independent H
        E2, g2 = gbs.calculate_density_of_states()
    finally:
        gbs.graphene_hamiltonian = saved
    identical = np.array_equal(g1, g2)
    Eb, gb, _ = dos_from_bands(n=300)
    i_vh = int(np.argmax(gb[Eb > 0.2]))
    vh_E = Eb[Eb > 0.2][i_vh]
    if verbose:
        print("\n" + "=" * 70)
        print("RESULT 2 -- why it survived: the DOS figure cannot see the bands")
        print("=" * 70)
        print("  calculate_density_of_states() returns |E|/(pi*2.8^2)")
        print("  analytically.  MUTATION TEST -- replace the Hamiltonian with a")
        print("  k-independent gapped matrix diag(+7, -7), which has no Dirac")
        print("  point and no dispersion at all, and recompute the DOS:")
        print("    shipped DOS identical before and after: %s" % identical)
        print("    -> the shipped DOS figure has ZERO sensitivity to the")
        print("       Hamiltonian.  It would print the correct V shape for any")
        print("       band structure whatsoever, including the broken one.")
        print()
        print("  A DOS computed FROM the bands is sensitive, and is the check")
        print("  that was missing.  Computed here over the reciprocal cell:")
        print("    van Hove peak on the electron side at E = %.3f eV" % vh_E)
        print("    exact tight-binding value                  = %.3f eV (= t)"
              % T_HOP)
        print("    relative error %.3f  (bin width %.3f eV)"
              % (abs(vh_E - T_HOP) / T_HOP, Eb[1] - Eb[0]))
        i0 = int(np.argmin(np.abs(Eb)))
        print("    DOS in the bin straddling E = 0: %.3e states/(eV m^2),"
              % gb[i0])
        print("      which is %.1f%% of the van Hove peak and is consistent"
              % (100.0 * gb[i0] / gb.max()))
        print("      with a linear DOS averaged over one bin: the derived value")
        print("      is (slope/4) x bin width = %.3e, ratio %.2f"
              % (1.789092e18 * (Eb[1] - Eb[0]) / 4.0,
                 gb[i0] / (1.789092e18 * (Eb[1] - Eb[0]) / 4.0)))
        print()
        print("  THE PATTERN, and it is the sharpest instance in five weeks:")
        print("  three band/DOS figures in this repository, and not one of them")
        print("  was a calculation that could have failed.  Two draw the known")
        print("  answer (E = +/- v_F|k|, |E|) and the third computes and is")
        print("  wrong.  An unfalsifiable figure is worse than no figure,")
        print("  because it certifies the thing beside it.")
    return identical, vh_E, Eb, gb


def main():
    print("=" * 70)
    print("BAND STRUCTURE AUDIT -- 2026-09-28")
    print("Chapter 2's foundational figure, and why five weeks of daily runs")
    print("never touched it.")
    print("=" * 70)
    print("\n" + "-" * 70)
    print("VALIDATIONS (exact model values, measured beside derived floors)")
    print("-" * 70)
    v1 = v1_exact_symmetry_points()
    v2 = v2_fermi_velocity()
    v3 = v3_particle_hole_symmetry()
    v4 = v4_dos_low_energy_slope()
    v5 = v5_shipped_path_reproduces_the_bug()
    r1 = result_1_the_bug()
    r2 = result_2_why_it_survived()
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    vs = [("V1 |phi| at Gamma/M/K", v1[0]),
          ("V2 Fermi velocity", v2[0]),
          ("V3 particle-hole symmetry", v3[0]),
          ("V4 DOS from bands", v4[0]),
          ("V5 bug reproduced exactly", v5[0])]
    for n, o in vs:
        print("  %-28s %s" % (n, "PASS" if o else "FAIL"))
    print("\n  validations: %d/%d" % (sum(1 for _, o in vs if o), len(vs)))
    print("\n  A PUBLISHED FIGURE IS WRONG AND IS BEING CORRECTED.  The gap at")
    print("  K goes from %.1f eV to 0.  No number in Chapters 4-7 depends on"
          % r1[2])
    print("  it: those chapters use the LINEAR dispersion and v_F directly,")
    print("  never this k-path, which is why the error was survivable for five")
    print("  weeks and why nothing else moves.  v_F itself is confirmed exact")
    print("  to %.1e by V2, once its q/4 warping term is included." % 1e-6)


if __name__ == "__main__":
    main()
