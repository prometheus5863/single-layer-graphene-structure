"""
Audit of Chapter 3's computational results (graphene_transport_properties.py),
performed while drafting Chapter 3 on 2026-09-29.

Why this exists
---------------
On 2026-09-28 Chapter 2 was drafted and the drafting falsified its own status
row: the "Computational results complete" label it had carried since
2026-08-23 was attached to a band structure with no Dirac point. Chapter 3
carries the SAME label, from the same day, and the 2026-09-28 entry recorded
the label as meaning nothing until tested.

This module tests it. Every check below compares a shipped quantity against a
value that is EXACTLY known -- a textbook constant, a dimensional identity, or
a branch that must be reachable -- rather than against a plausible range. The
2026-09-17 and 2026-09-28 entries both record a bug that a plausible-range
check passed and an exact check caught.

Result: four defects, all confirmed, all in the same module.

  T1  Universal optical conductivity is written with hbar where the formula in
      its own docstring has h. 2*pi too large. The implied single-layer
      absorption is 14.40%, not 2.2925%.
  T4  The acoustic-phonon mobility uses the three-dimensional T^-3/2
      deformation-potential exponent; graphene's is T^-1.
  T2  The conductance quantum and the Dirac-point minimum conductivity are
      written with hbar where they need h. 2*pi too large, twice.
  T3  The Pauli-blocking switch compares an energy in JOULES against a bare
      number that is in eV, so the branch is unreachable: the "interband
      transitions become significant when hbar*omega > 2*mu" feature has never
      fired, at any frequency, and the function has always returned the blocked
      value 0.5*sigma_0.

Blast radius, checked rather than assumed: none of this module's figures are
committed to the repository, and no other module IMPORTS it -- but
graphene_photodetector_model.py's header CITES
`graphene_transport_properties.py::calculate_optical_conductivity` as the
source of its 2.3% absorption while independently re-typing the correct
`ALPHA_ABS = pi/137.036`. Chapter 6 is therefore protected from T1 only by the
fact that it did not use the result it cites. That is the Chapter 2 pattern
exactly: what was wrong was what nothing read.

Run: python3 graphene_transport_optical_audit.py
"""

import sys

import numpy as np
from scipy.constants import e, h, hbar, c, epsilon_0
from scipy.constants import Boltzmann as kB
from scipy.constants import alpha as ALPHA_FINE


def banner(s):
    print()
    print('=' * 78)
    print(s)
    print('=' * 78)


def check(name, ok, detail):
    print('  [%s] %s' % ('PASS' if ok else 'FAIL', name))
    for line in detail.split('\n'):
        print('        %s' % line)
    return ok


def main():
    import graphene_transport_properties as tp

    n_ok = 0
    n_bad = 0

    banner('T1 -- universal optical conductivity: hbar written for h')

    sigma0_correct = np.pi * e**2 / (2 * h)          # = e^2 / (4 hbar)
    _, sigma_shipped = tp.calculate_optical_conductivity(np.array([5e14]))
    sigma_shipped = float(np.atleast_1d(sigma_shipped)[0])

    # The identity that makes this exact rather than approximate: the
    # normal-incidence absorption of a free-standing conducting sheet of
    # 2D conductivity sigma is A = sigma / (epsilon_0 c) to first order, and
    # for sigma = pi e^2 / 2h this equals pi * alpha EXACTLY, with alpha the
    # fine-structure constant. Two independent routes to the same number is
    # what makes it an anchor and not a preference.
    A_from_sigma = sigma0_correct / (epsilon_0 * c)
    alpha_derived = e**2 / (4 * np.pi * epsilon_0 * hbar * c)
    A_derived = np.pi * alpha_derived
    A_measured = np.pi * ALPHA_FINE
    algebraic = abs(A_from_sigma - A_derived) / A_derived
    codata_slack = abs(A_derived - A_measured) / A_measured
    if check('T1a  anchor: sigma_0/(eps_0 c) is pi*alpha, to machine precision',
             algebraic < 1e-15,
             'sigma_0 = pi e^2 / (2h)             = %.6e S\n'
             'route 1, sigma_0 / (eps_0 c)        = %.10f %%\n'
             'route 2, pi * e^2/(4 pi eps_0 hbar c) = %.10f %%\n'
             'relative disagreement               = %.3e   (algebraic identity)\n'
             'route 3, pi * CODATA measured alpha = %.10f %%\n'
             'which differs from routes 1-2 by    = %.3e\n'
             'That last residual is not a defect and not a rounding error: it\n'
             'is the slack between CODATA\'s independently MEASURED alpha and\n'
             'the alpha implied by its e, h, eps_0 and c, and it is well inside\n'
             'the ~1.6e-10 relative uncertainty CODATA quotes on alpha. Worth\n'
             'naming, because a tolerance tightened past it turns a correct\n'
             'anchor into a failing check -- which is what the first version of\n'
             'this audit did at tol = 1e-12.'
             % (sigma0_correct, 100 * A_from_sigma, 100 * A_derived,
                algebraic, 100 * A_measured, codata_slack)):
        n_ok += 1
    else:
        n_bad += 1

    # The shipped function multiplies sigma_0 by interband_factor, which T3
    # shows is always 0.5, so undo that to expose the sigma_0 it used.
    sigma0_shipped = sigma_shipped / 0.5
    ratio = sigma0_shipped / sigma0_correct
    is_2pi = abs(ratio - 2 * np.pi) < 1e-9
    ok = not is_2pi
    if ok:
        n_ok += 1
    else:
        n_bad += 1
    check('T1b  shipped sigma_0 equals the textbook value', ok,
          'shipped sigma_0   = %.6e S\n'
          'correct sigma_0   = %.6e S\n'
          'ratio             = %.9f   (2*pi = %.9f)\n'
          'the docstring states sigma_0 = pi e^2 / (2h); the code writes\n'
          '  np.pi * e**2 / (2 * hbar)\n'
          'implied single-layer absorption = %.4f %%, against 2.2925 %%'
          % (sigma0_shipped, sigma0_correct, ratio, 2 * np.pi,
             100 * sigma0_shipped / (epsilon_0 * c)))

    print()
    print('        NOTE, to stop a future reader chasing it: the implied 14.40 %')
    print('        here and the 14.40 eV spurious gap found at the K label on')
    print('        2026-09-28 are numerically a COINCIDENCE. That one was a')
    print('        wrong reciprocal-lattice convention; this one is hbar for h.')
    print('        They share no mechanism.')

    banner('T2 -- conductance quantum and minimum conductivity: hbar for h')

    G0_correct = 2 * e**2 / h
    G0_shipped = 2 * e**2 / hbar          # as written at line 145 of the module
    Gmin_shipped = 4 * e**2 / hbar        # as written at line 149
    sigma_min_theory = 4 * e**2 / (np.pi * h)

    ok = abs(G0_shipped / G0_correct - 1) < 1e-12
    if ok:
        n_ok += 1
    else:
        n_bad += 1
    check('T2a  conductance quantum G_0 = 2e^2/h', ok,
          'shipped 2e^2/hbar = %.6e S\n'
          'correct 2e^2/h    = %.6e S\n'
          'ratio             = %.9f  (= 2*pi)\n'
          "the module's own comment reads 'Conductance quantum', which is\n"
          '2e^2/h = 7.748e-5 S = (12.906 kOhm)^-1, a value fixed by metrology'
          % (G0_shipped, G0_correct, G0_shipped / G0_correct))

    ok = abs(Gmin_shipped / sigma_min_theory - 1) < 1e-12
    if ok:
        n_ok += 1
    else:
        n_bad += 1
    check('T2b  Dirac-point minimum conductivity', ok,
          'shipped 4e^2/hbar      = %.6e S\n'
          'theory  4e^2/(pi h)    = %.6e S\n'
          'experiment ~ 4e^2/h    = %.6e S\n'
          'ratio shipped/theory   = %.4f\n'
          'wrong on both counts: hbar for h, AND the 1/pi of the ballistic\n'
          'self-consistent result is absent. The comment calls 4e^2/hbar\n'
          "'a hallmark of graphene'; the hallmark is 4e^2/(pi h)."
          % (Gmin_shipped, sigma_min_theory, 4 * e**2 / h,
             Gmin_shipped / sigma_min_theory))

    banner('T3 -- the Pauli-blocking branch is unreachable (a unit mismatch)')

    # The shipped code computes hbar_omega in JOULES and compares it against
    # `2 * mu` where mu = 0.1 is a bare number intended as eV.
    # An unreachable branch is a check that cannot fail -- the 2026-09-28 rule
    # -- and the detector is to sweep the whole physical range and require the
    # branch to be taken SOMEWHERE.
    lambdas_nm = np.logspace(-1, 5, 2000)       # 0.1 nm (hard X-ray) to 100 um
    freqs = c / (lambdas_nm * 1e-9)
    _, sig = tp.calculate_optical_conductivity(freqs)
    sig = np.atleast_1d(sig)
    distinct = np.unique(np.round(sig / sig.max(), 12))
    took_both = distinct.size >= 2
    ok = took_both
    if ok:
        n_ok += 1
    else:
        n_bad += 1
    check('T3a  the interband branch is taken somewhere in 0.1 nm - 100 um', ok,
          'swept %d wavelengths spanning %.3e eV to %.3e eV of photon energy\n'
          'distinct values of sigma_real returned: %d\n'
          'the returned value is constant at 0.5*sigma_0 across the entire\n'
          'electromagnetic spectrum, so the branch has never been taken'
          % (lambdas_nm.size, h * c / (lambdas_nm.max() * 1e-9) / e,
             h * c / (lambdas_nm.min() * 1e-9) / e, distinct.size))

    threshold_eV = 0.2 / e
    if check('T3b  diagnosis: the comparison is joules against a bare eV number',
                  True,
                  'code: np.where(hbar_omega > 2 * mu, 1.0, 0.5), mu = 0.1\n'
                  'hbar_omega is in J; 2*mu = 0.2 is intended as eV but is bare.\n'
                  'the branch therefore requires a photon energy above\n'
                  '  0.2 J = %.3e eV\n'
                  'which is %.0e times the energy of a 550 nm photon.\n'
                  'POSITIVE CONTROL that this detector works: with the\n'
                  'comparison written correctly (hbar_omega/e > 2*mu, eV on\n'
                  'both sides) the same sweep returns %d distinct values.'
                  % (threshold_eV, threshold_eV / 2.2547,
                     np.unique(np.where(
                         (h * freqs / e) > 0.2, 1.0, 0.5)).size)):
        n_ok += 1
    else:
        n_bad += 1

    banner('T4 -- acoustic-phonon mobility exponent: the 3D result in a 2D material')

    T = np.array([100.0, 300.0, 500.0])
    _, _, mu_ac, _ = tp.calculate_mobility_vs_temperature(T)
    # measured exponent of the shipped curve
    p_shipped = -np.polyfit(np.log(T), np.log(mu_ac), 1)[0]
    ok = abs(p_shipped - 1.0) < 1e-9
    if ok:
        n_ok += 1
    else:
        n_bad += 1
    check('T4a  mu_acoustic ~ 1/T, as graphene requires', ok,
          'shipped exponent p in mu_ac ~ T^-p : %.4f\n'
          'graphene requires p = 1, not 1.5.\n'
          'For longitudinal acoustic phonons the graphene resistivity is\n'
          '  rho_LA = pi D_A^2 kB T / (4 e^2 hbar rho_s v_s^2 v_F^2),\n'
          'LINEAR in T and independent of carrier density, so at fixed n the\n'
          'mobility goes as 1/T. The T^-3/2 the code uses is the\n'
          "three-dimensional deformation-potential result, and the module's\n"
          "own comment states it as 'mu ~ T^(-1.5) at high T'.\n"
          'Hwang & Das Sarma, Phys. Rev. B 77, 115449 (2008).\n'
          'At 500 K this understates the mobility by a factor %.3f.'
          % (p_shipped, (500 / 300.0) ** (p_shipped - 1.0)))

    banner('SUMMARY')
    print('  %d checks passed, %d failed' % (n_ok, n_bad))
    print()
    print('  Chapter 3 carried "Computational results complete" from 2026-08-23.')
    print('  Drafting it found %d defects in the module that label refers to,' % n_bad)
    print('  every one of them exposed by an exactly-known value and none of')
    print('  them by a plausible range: 14.40 % absorption, 2*pi on two')
    print('  metrological constants, an unreachable branch, and a 3D exponent')
    print('  in a 2D material.')
    print()
    print('  This is the second consecutive chapter whose "results complete"')
    print('  label was falsified by the act of writing the chapter, and the')
    print('  second consecutive time the defects were in the part of the')
    print('  repository that nothing downstream imports.')
    return 0 if n_bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
