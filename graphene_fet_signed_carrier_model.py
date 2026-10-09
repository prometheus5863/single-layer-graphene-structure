"""
graphene_fet_signed_carrier_model.py

A SIGNED-CARRIER FET MODEL.  Chapter 7 Section 7.9 item 9, created 2026-10-07,
promoted 2026-10-08 from sibling to *prerequisite* of item 12 (re-derive
Chapter 4 on the exact charge relation), and untouched until today.

------------------------------------------------------------------------------
WHY THIS BLOCKS ITEM 12
------------------------------------------------------------------------------
`graphene_fet_model.carrier_density()` returns

    n = sqrt(n_eff**2 + n_puddle**2)                                      (D2)

which is a MAGNITUDE.  The sign is computed -- `n_eff = C_ox*dV/e * C_q/(C_q +
C_ox)` is odd in `dV` because `C_q` is even -- and then discarded by the
quadrature floor.  Everything downstream takes `abs(n)` (`sheet_conductivity`),
so for the committed conductivity numbers the discard is harmless.  It is NOT
harmless for item 12.  The exact charge relation declared on 2026-10-08 is

    V_g - V_dirac - V_ch = e*n/C_ox + E_F(n)/e                            (4.28)

and both terms on the right are ODD in n.  Solving it on a magnitude, which is
what `graphene_potential_declaration.n_quasi_fermi()` does (it opens with
`dV = np.abs(dV)`), yields the right |n| but leaves (4.28) with no branch:
there is no object in the code that knows whether the channel at a given point
is n-type or p-type.  Re-deriving Chapter 4 on (4.28) *first* would therefore
rebuild the branch-sign flaw inside a new equation, which is precisely the
reason 10-08 ordered item 9 ahead of item 12.

This module supplies the missing branch, and it does so without changing a
single committed number.  That claim is checked, not asserted: Section 1's
decomposition reproduces (D2) as an identity, and Section 6 re-runs the
shipped transfer characteristic unchanged.

------------------------------------------------------------------------------
THE DECOMPOSITION
------------------------------------------------------------------------------
Near the Dirac point a disordered graphene sheet is not "almost empty"; it is
an electron-hole PUDDLE LANDSCAPE carrying both species at once [Martin et al.,
Nature Physics 4, 144 (2008); Adam et al., PNAS 104, 18392 (2007)].  Write the
net (electrostatic) density as signed and split it into species:

    n_net = n_e - n_h          (set by electrostatics, SIGNED, odd in dV)
    S     = sqrt(n_net**2 + n_puddle**2)
    n_e   = (S + n_net)/2,     n_h = (S - n_net)/2

Three things follow, and all three are exact rather than fitted:

  (i)   n_e + n_h == S == exactly what (D2) returns.  The shipped magnitude is
        revealed to be the TOTAL CONDUCTING DENSITY, not the net density and
        not a per-species density.  Section 7 shows that identification is what
        makes the ten-session-old Section 4.9 floor item decidable at all.
  (ii)  n_e * n_h == n_puddle**2 / 4, a mass-action law with no free parameter.
        The quadrature floor of (D2), chosen in August for its smoothness, is
        algebraically equivalent to a constant electron-hole product -- the
        same closure the puddle literature writes down from physics.
  (iii) n_e(+dV) == n_h(-dV) BITWISE.  Graphene's ambipolar symmetry becomes a
        machine-checkable identity instead of a figure one looks at.

(iii) is this session's exactly-known-value check, in the sense 09-28 through
10-08 established: a symmetry that must return zero, not a plausible range.

------------------------------------------------------------------------------
WHAT THIS MODULE ADDS THAT THE MAGNITUDE MODEL CANNOT EXPRESS
------------------------------------------------------------------------------
Section 4   a BRANCHED root of (4.28), odd by construction, so item 12 has a
            signed charge relation to be re-derived on.
Section 5   the internal charge-neutrality point: at finite V_ds the channel
            splits into a p region and an n region, and the shipped model
            represented that p-n junction as a dip in |n|.
Section 7   the Section 4.9 floor item, open since 2026-09-29 and untouched
            for ELEVEN sessions, answered -- with a factor the signed reading
            is what pins down.

Run:  python graphene_fet_signed_carrier_model.py
Writes: fet_signed_carrier_model.png, fet_signed_carrier_output.txt (via tee)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import graphene_diffusion_current_model as dd

E_CHARGE = 1.602176634e-19
H_PLANCK = 6.62607015e-34
A_FERMI = dd.A_FERMI          # V_F = E_F/e = A_FERMI*sqrt(n), from the repo's v_F

# Chapter 3 Section 3.3's MEASURED Dirac-point minimum conductivity, 4 q_e^2/h.
# Quoted there as 1.5496e-04 S and 6.45 kOhm/sq.  Not a new parameter: it is
# this thesis's own Chapter 3 number, recomputed here from CODATA so the
# comparison in Section 7 cannot drift against a re-typed literal.
SIGMA_MIN_CH3 = 4.0 * E_CHARGE ** 2 / H_PLANCK

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
# Section 1: the signed net density and the species split
# ---------------------------------------------------------------------------

def n_net(V_g, V_ch=0.0):
    """
    SIGNED net sheet density [1/m^2] from the repository's own series-factor
    form.  Byte-for-byte the expression inside carrier_density() with the
    quadrature floor and the magnitude NOT applied.  Positive = n-type.
    """
    dV = np.asarray(V_g, dtype=float) - V_ch - gfet.V_dirac
    n_es = gfet.C_ox * dV / E_CHARGE
    C_q = gfet.quantum_capacitance(dV, T=gfet.T)
    return n_es * C_q / (C_q + gfet.C_ox)


def n_net_of_dV(dV, n_puddle=None):
    """
    The same signed net density as a function of the LOCAL OVERDRIVE dV
    directly, with no round trip through V_dirac.

    This entry point exists because of Section 2's finding: writing the bias as
    V_dirac +/- x and letting n_net() subtract V_dirac back is not an exact
    operation in binary floating point, and it manufactures a violation of
    graphene's ambipolar symmetry that the model itself does not contain.  Every
    sweep in this repository is built as a linspace over V_g, so any BITWISE
    symmetry claim has to say which parameterisation it is bitwise in.
    """
    dV = np.asarray(dV, dtype=float)
    n_es = gfet.C_ox * dV / E_CHARGE
    C_q = gfet.quantum_capacitance(dV, T=gfet.T)
    return n_es * C_q / (C_q + gfet.C_ox)


def split_species(n_net_val, n_puddle=None):
    """
    Electron and hole sheet densities [1/m^2] from the signed net density,
    closed by the puddle landscape.  Both are >= 0 by construction.
    """
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    n_net_val = np.asarray(n_net_val, dtype=float)
    S = np.sqrt(n_net_val ** 2 + n_puddle ** 2)
    return 0.5 * (S + n_net_val), 0.5 * (S - n_net_val)


def sheet_conductivity_two_carrier(n_e, n_h, mu_e=None, mu_h=None):
    """
    Two-carrier Drude conductivity, sigma = e*(mu_e*n_e + mu_h*n_h) [S/sq].

    With mu_e == mu_h == gfet.mu this is e*mu*(n_e + n_h) and therefore
    IDENTICAL to the shipped sheet_conductivity(carrier_density(...)); the
    electron-hole mobility asymmetry that would separate them is a real effect
    in contact-doped devices and is NOT introduced here, because this module's
    contract is to add a branch without moving a number.  Asymmetric mobility
    is logged as a follow-on item instead.
    """
    mu_e = gfet.mu if mu_e is None else mu_e
    mu_h = gfet.mu if mu_h is None else mu_h
    return E_CHARGE * (mu_e * np.asarray(n_e) + mu_h * np.asarray(n_h))


# ---------------------------------------------------------------------------
# Section 4: the branched root of the exact charge relation (4.28)
# ---------------------------------------------------------------------------

def n_net_exact(dV):
    """
    SIGNED root of (4.28): e*n/C_ox + A_FERMI*sign(n)*sqrt(|n|) = dV.

    Both left-hand terms are odd in n, so the relation is odd and its root is
    odd in dV.  With x = sqrt(|n|), a = e/C_ox, b = A_FERMI the magnitude
    solves a*x^2 + b*x - |dV| = 0 on BOTH branches; the branch is then carried
    by sign(dV).  The 2026-10-08 solver stopped at the magnitude; this function
    is the same closed form with the branch restored, which is the whole of
    what item 12 was waiting on.
    """
    dV = np.asarray(dV, dtype=float)
    a = E_CHARGE / gfet.C_ox
    b = A_FERMI
    x = (-b + np.sqrt(b * b + 4.0 * a * np.abs(dV))) / (2.0 * a)
    return np.sign(dV) * x * x


def exact_relation_residual(n_signed):
    """Left-hand side of (4.28) evaluated on a SIGNED n, in volts."""
    n_signed = np.asarray(n_signed, dtype=float)
    return (E_CHARGE * n_signed / gfet.C_ox
            + A_FERMI * np.sign(n_signed) * np.sqrt(np.abs(n_signed)))


# ---------------------------------------------------------------------------
# Section 5: where the channel changes branch
# ---------------------------------------------------------------------------

def channel_branch_census(V_g, Vds, n_segments=50):
    """
    Branch structure of the channel on the SAME integration grid
    transfer_characteristic() uses: V_ch uniform from 0 (source) to Vds (drain).

    Returns (frac_n, frac_p, V_ch_cnp) where V_ch_cnp is the channel potential
    at which the local branch flips, or None if the channel is unipolar.
    """
    V_profile = np.linspace(0.0, Vds, n_segments)
    nn = n_net(V_g, V_profile)
    frac_n = float(np.mean(nn > 0.0))
    frac_p = float(np.mean(nn < 0.0))
    V_cnp = V_g - gfet.V_dirac
    inside = (0.0 < V_cnp) and (V_cnp < Vds)
    return frac_n, frac_p, (V_cnp if inside else None)


def ambipolar_window(Vds):
    """
    The V_g interval over which the charge-neutrality point sits INSIDE the
    channel, i.e. the channel is simultaneously p-type and n-type.
    dV = V_g - V_ch - V_dirac vanishes for some V_ch in (0, Vds) exactly when
    V_dirac < V_g < V_dirac + Vds.  Exact, not swept.
    """
    return gfet.V_dirac, gfet.V_dirac + Vds


# ---------------------------------------------------------------------------
# Section 7: the Chapter 4 / Chapter 3 floor adjudication
# ---------------------------------------------------------------------------

def sigma_min_ch4(n_puddle=None):
    """
    Chapter 4's Dirac-point sheet conductivity [S/sq].  At dV = 0 the split
    gives n_e = n_h = n_puddle/2, so sigma_min = e*mu*n_puddle EXACTLY -- and
    it is the signed decomposition that makes clear n_puddle is the TOTAL
    (n_e + n_h), which is the step that turns the comparison below from a
    factor-of-2 convention argument into a number.
    """
    if n_puddle is None:
        n_puddle = gfet.n_puddle
    n_e, n_h = split_species(0.0, n_puddle=n_puddle)
    return float(sheet_conductivity_two_carrier(n_e, n_h))


def n_puddle_matching_ch3():
    """The total puddle density that reproduces Chapter 3's 4 q_e^2/h floor."""
    return SIGMA_MIN_CH3 / (E_CHARGE * gfet.mu)


def onoff_at(Vds, Vg_on, n_puddle):
    """
    (I_on, I_off, ratio) from the SHIPPED transfer characteristic with
    gfet.n_puddle temporarily set to n_puddle.  I_off is taken at V_g =
    V_dirac, the current minimum of the ambipolar curve.
    """
    saved = gfet.n_puddle
    try:
        gfet.n_puddle = n_puddle
        I = gfet.transfer_characteristic(np.array([Vg_on, gfet.V_dirac]), Vds=Vds)
    finally:
        gfet.n_puddle = saved
    I_on, I_off = float(I[0]), float(I[1])
    return I_on, I_off, I_on / I_off


def intrinsic_onoff_at(Vds, Vg_on, n_puddle):
    """
    The same ratio with the contact resistance removed, so the attenuation
    between the channel and the terminals can be quoted rather than inferred.
    """
    saved_np = gfet.n_puddle
    saved_rc = gfet.Rc_total
    try:
        gfet.n_puddle = n_puddle
        gfet.Rc_total = 0.0
        I = gfet.transfer_characteristic(np.array([Vg_on, gfet.V_dirac]), Vds=Vds)
    finally:
        gfet.n_puddle = saved_np
        gfet.Rc_total = saved_rc
    I_on, I_off = float(I[0]), float(I[1])
    return I_on, I_off, I_on / I_off


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(Vds_demo=0.5):
    dV = np.linspace(-2.0, 2.0, 1201)
    nn = n_net(gfet.V_dirac + dV, 0.0)
    n_e, n_h = split_species(nn)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))

    ax = axes[0]
    ax.semilogy(dV, n_e * 1e-4, linewidth=2, label=r'$n_e$ (electrons)')
    ax.semilogy(dV, n_h * 1e-4, linewidth=2, label=r'$n_h$ (holes)')
    ax.semilogy(dV, np.abs(nn) * 1e-4, '--', color='k', linewidth=1.4,
                label=r'$|n_{net}| = |n_e - n_h|$')
    ax.semilogy(dV, (n_e + n_h) * 1e-4, ':', color='crimson', linewidth=2.2,
                label=r'$n_e + n_h$ = shipped $n(V_g)$')
    ax.axhline(gfet.n_puddle * 0.5 * 1e-4, color='gray', linestyle='-.',
               linewidth=1, label=r'$n_{puddle}/2$')
    ax.set_xlabel(r'$V_g - V_{Dirac}$ (V)')
    ax.set_ylabel(r'sheet density (cm$^{-2}$)')
    ax.set_title('1. The species the magnitude model hid\n'
                 r'$n_e n_h = n_{puddle}^2/4$ exactly; $n_e(+\Delta V)=n_h(-\Delta V)$ bitwise')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, alpha=0.3, which='both')

    ax = axes[1]
    Vg = np.linspace(gfet.V_dirac - 0.6, gfet.V_dirac + Vds_demo + 0.6, 400)
    fn = np.array([channel_branch_census(v, Vds_demo)[0] for v in Vg])
    fp = np.array([channel_branch_census(v, Vds_demo)[1] for v in Vg])
    ax.plot(Vg, fn * 100, linewidth=2, label='channel fraction n-type')
    ax.plot(Vg, fp * 100, linewidth=2, label='channel fraction p-type')
    lo, hi = ambipolar_window(Vds_demo)
    ax.axvspan(lo, hi, color='gold', alpha=0.25,
               label='ambipolar: p-n junction INSIDE channel')
    ax.set_xlabel(r'$V_g$ (V)')
    ax.set_ylabel('fraction of channel (%)')
    ax.set_title('2. The internal p-n junction\n'
                 r'$V_{ds}$ = %.2f V; window $V_{Dirac} < V_g < V_{Dirac}+V_{ds}$'
                 % Vds_demo)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    np_scan = np.logspace(np.log10(5e14), np.log10(5e16), 400)
    R_sheet = 1.0 / np.array([sigma_min_ch4(p) for p in np_scan])
    ax.loglog(np_scan * 1e-4, R_sheet * 1e-3, linewidth=2,
              label=r'Ch. 4 floor $1/(e\mu n_{puddle})$')
    ax.axhline(1.0 / SIGMA_MIN_CH3 * 1e-3, color='crimson', linestyle='--',
               linewidth=2, label=r'Ch. 3 $\S$3.3 measured $4q_e^2/h$ = 6.45 k$\Omega$/sq')
    ax.axvline(gfet.n_puddle * 1e-4, color='gray', linestyle=':', linewidth=2,
               label=r'shipped $n_{puddle}$ = 5$\times$10$^{11}$ cm$^{-2}$')
    ax.axvline(n_puddle_matching_ch3() * 1e-4, color='seagreen', linestyle='-.',
               linewidth=2, label=r'$n_{puddle}$ matching Ch. 3')
    ax.set_xlabel(r'$n_{puddle} = n_e + n_h$ at the Dirac point (cm$^{-2}$)')
    ax.set_ylabel(r'Dirac-point sheet resistance (k$\Omega$/sq)')
    ax.set_title('3. Two chapters, two floors\n'
                 'the eleven-session-old $\\S$4.9 item, decidable once $n_{puddle}$ is named')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('fet_signed_carrier_model.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('A SIGNED-CARRIER FET MODEL -- Chapter 7 Section 7.9 item 9')
    print('the prerequisite item 12 has been waiting on since 2026-10-08')
    print('=' * 78)
    print()
    print('Device constants read FROM graphene_fet_model (not re-typed here):')
    print('  C_ox      = %.6e F/m^2' % gfet.C_ox)
    print('  n_puddle  = %.6e 1/m^2   (%.3e 1/cm^2)'
          % (gfet.n_puddle, gfet.n_puddle * 1e-4))
    print('  mu        = %.4f m^2/(V.s)' % gfet.mu)
    print('  V_dirac   = %.4f V' % gfet.V_dirac)
    print('  A_FERMI   = %.6e V.m          (from graphene_diffusion_current_model)' % A_FERMI)
    print()

    # -- Section 1 ---------------------------------------------------------
    print('-' * 78)
    print('SECTION 1.  The decomposition reproduces (D2) as an identity')
    print('-' * 78)
    dV_grid = np.linspace(-2.0, 2.0, 801)
    Vg_grid = gfet.V_dirac + dV_grid
    nn = n_net(Vg_grid, 0.0)
    n_e, n_h = split_species(nn)
    shipped = gfet.carrier_density(Vg_grid, 0.0)

    tot = n_e + n_h
    rel_tot = np.max(np.abs(tot - shipped) / shipped)
    bitwise_tot = bool(np.array_equal(tot, shipped))
    check('n_e + n_h reproduces carrier_density() over 801 points',
          rel_tot < 1e-15,
          'max relative difference %.3e ; bitwise identical: %s'
          % (rel_tot, bitwise_tot))

    net_back = n_e - n_h
    S = np.sqrt(nn ** 2 + gfet.n_puddle ** 2)
    abs_net = np.abs(net_back - nn)
    ulp_net = float(np.max(abs_net / np.spacing(S)))
    rel_vs_S = float(np.max(abs_net / S))
    rel_vs_n = float(np.max(abs_net / np.maximum(np.abs(nn), 1.0)))
    check('n_e - n_h recovers the net density to <= 1 ulp OF THE TOTAL',
          ulp_net <= 1.0,
          'max %.3e 1/m^2 = %.2f ulp of S ; %.3e relative to S'
          % (float(np.max(abs_net)), ulp_net, rel_vs_S))
    print('         AND A FINDING, not a rounding note: measured against the')
    print('         NET density instead of the total the same residual is')
    print('         %.3e -- a factor %.0f larger -- because n_net -> 0 at' % (rel_vs_n, rel_vs_n / rel_vs_S))
    print('         neutrality while S stays at the puddle floor.  Writing')
    print('         this check as a ratio and choosing its denominator is')
    print('         what decides PASS from FAIL, which is Section 7.9\'s')
    print('         10-06 ratio-criterion item in a fresh instance.')
    print('         DESIGN RULE: never reconstruct the net density from the')
    print('         species pair; call n_net() / n_net_of_dV() directly.')

    prod = n_e * n_h
    target = (gfet.n_puddle ** 2) / 4.0
    rel_prod = np.max(np.abs(prod - target) / target)
    check('MASS ACTION: n_e*n_h == n_puddle^2/4 with no free parameter',
          rel_prod < 1e-13,
          'target %.6e ; max relative residual %.3e' % (target, rel_prod))

    check('both species non-negative everywhere',
          bool(np.all(n_e >= 0.0) and np.all(n_h >= 0.0)),
          'min n_e %.6e , min n_h %.6e' % (n_e.min(), n_h.min()))

    # -- Section 2: the exactly-known-value checks -------------------------
    print()
    print('-' * 78)
    print('SECTION 2.  Validation against exactly known values')
    print('            (symmetries that must return ZERO, not plausible ranges)')
    print('-' * 78)

    # (a) the symmetry stated in the variable the model is a function OF
    a_e, a_h = split_species(n_net_of_dV(dV_grid))
    b_e, b_h = split_species(n_net_of_dV(-dV_grid))
    check('AMBIPOLAR SYMMETRY n_e(+dV) == n_h(-dV), BITWISE in dV, 801 points',
          bool(np.array_equal(a_e, b_h) and np.array_equal(a_h, b_e)),
          'max absolute difference %.3e 1/m^2 -- exact zero required and met'
          % float(max(np.max(np.abs(a_e - b_h)), np.max(np.abs(a_h - b_e)))))

    # (b) the SAME symmetry, routed through V_dirac the way every sweep in
    #     this repository routes it.  This is the session's falsified check.
    n_e_p, n_h_p = split_species(n_net(gfet.V_dirac + dV_grid, 0.0))
    n_e_m, n_h_m = split_species(n_net(gfet.V_dirac - dV_grid, 0.0))
    mirror_resid = float(np.max(np.abs(n_e_p - n_h_m)))
    bitwise_vg = bool(np.array_equal(n_e_p, n_h_m))
    dv_p = (gfet.V_dirac + dV_grid) - gfet.V_dirac
    dv_m = (gfet.V_dirac - dV_grid) - gfet.V_dirac
    dv_odd = bool(np.array_equal(dv_m, -dv_p))
    # The V_g route is ACCOUNTED FOR rather than thresholded: feed the
    # reconstructed arguments into the dV-level function and require bitwise
    # agreement.  If that holds, the whole residual above is the argument's.
    pred_e, pred_h = split_species(n_net_of_dV(dv_p))
    pred_e_m, pred_h_m = split_species(n_net_of_dV(dv_m))
    accounted = bool(np.array_equal(pred_e, n_e_p)
                     and np.array_equal(pred_h_m, n_h_m))
    check('the V_g-route residual is FULLY ACCOUNTED FOR by its argument',
          accounted,
          'feeding the reconstructed dV into the dV-level model reproduces '
          'the V_g-route species BITWISE, so the %.1f /m^2 below is the '
          'argument\'s error and none of it is the model\'s' % mirror_resid)
    print('         THIS CHECK FIRST RAN AS "BITWISE" AND FAILED, AT %.1f /m^2.'
          % mirror_resid)
    print('         The symmetry is NOT broken -- (a) above is bitwise.  What')
    print('         is not exact is the ARGUMENT: (V_dirac - x) - V_dirac is')
    print('         not -((V_dirac + x) - V_dirac) in binary floating point.')
    print('         Measured directly: that reconstruction is odd bitwise? %s,'
          % dv_odd)
    print('         worst |dV_- + dV_+| = %.3e V, at x = %.1f V, which the'
          % (float(np.max(np.abs(dv_m + dv_p))),
             float(dV_grid[int(np.argmax(np.abs(dv_m + dv_p)))])))
    print('         C_ox/e lever of %.3e /(V.m^2) turns into the %.1f /m^2'
          % (gfet.C_ox / E_CHARGE, mirror_resid))
    print('         above.  A symmetry check can fail because of how its')
    print('         argument was BUILT rather than because the symmetry is')
    print('         absent -- the 09-30 arrival lesson, in the input instead')
    print('         of in the mutation.')

    n0 = n_net(gfet.V_dirac, 0.0)
    e0, h0 = split_species(n0)
    check('dV = 0 gives n_net EXACTLY zero',
          float(n0) == 0.0, 'n_net(0) = %r' % float(n0))
    check('dV = 0 gives n_e == n_h BITWISE (the Dirac point is neutral)',
          float(e0) == float(h0),
          'n_e = n_h = %.6e 1/m^2 = n_puddle/2 = %.6e'
          % (float(e0), gfet.n_puddle / 2.0))

    nq = n_net_exact(dV_grid)
    nq_m = n_net_exact(-dV_grid)
    check('(4.28) signed root is ODD, BITWISE: n(-dV) == -n(+dV)',
          bool(np.array_equal(nq_m, -nq)),
          'max absolute difference %.3e 1/m^2' % float(np.max(np.abs(nq_m + nq))))

    resid = exact_relation_residual(nq) - dV_grid
    scale = np.maximum(np.abs(dV_grid), 1e-3)
    rel_resid = float(np.max(np.abs(resid) / scale))
    check('(4.28) is SATISFIED by its signed root (algebraic identity)',
          rel_resid < 1e-14,
          'max relative residual %.3e over dV in [-2, 2] V' % rel_resid)

    sigma_two = sheet_conductivity_two_carrier(n_e, n_h)
    sigma_shipped = gfet.sheet_conductivity(shipped)
    ulps_sig = float(np.max(np.abs(sigma_two - sigma_shipped)
                            / np.spacing(sigma_shipped)))
    rel_sig = float(np.max(np.abs(sigma_two - sigma_shipped) / sigma_shipped))
    # Criterion: the precision Chapter 4 actually commits.  Its transcripts
    # quote four significant figures; 1e-12 is six orders BELOW that and still
    # eleven orders above float noise, so it is a statement about the chapter
    # rather than a number tuned until the test went green.
    COMMITTED_PRECISION = 1e-12
    check('two-carrier sigma == shipped sigma at mu_e = mu_h, to the precision '
          'Chapter 4 commits (1e-12 rel)',
          rel_sig < COMMITTED_PRECISION,
          'max %.3e S/sq = %.3e relative = %.2f ulp ; bitwise: %s'
          % (float(np.max(np.abs(sigma_two - sigma_shipped))), rel_sig,
             ulps_sig, bool(np.array_equal(sigma_two, sigma_shipped))))
    print('         THIS CHECK WAS FIRST WRITTEN AS "BITWISE" AND FAILED at')
    print('         %.2f ulp, and it was NOT fixed by widening a ulp budget' % ulps_sig)
    print('         until it passed -- 09-29\'s alpha anchor established that')
    print('         the repair for a failing exactness check is to state the')
    print('         identity being checked and name the slack.  The identity:')
    print('         0.5*(S+n) + 0.5*(S-n) == S algebraically but reassociates')
    print('         the sum; e*mu* then adds two roundings.  The slack: %.3e'
          % rel_sig)
    print('         relative, against four significant figures in the chapter.')
    print('         10-02 said byte identity is not correctness.  Today adds')
    print('         the converse: ABSENCE of byte identity is not a defect,')
    print('         and three ulp of summation order is not a physics claim.')

    # -- Section 3: the branch census along the channel --------------------
    print()
    print('-' * 78)
    print('SECTION 3.  The internal p-n junction the magnitude model could not')
    print('            represent')
    print('-' * 78)
    print('  The shipped transfer_characteristic() averages the LOCAL channel')
    print('  resistance over V_ch from 0 to V_ds.  Where dV = V_g - V_ch -')
    print('  V_dirac changes sign the channel changes TYPE, and |n| merely')
    print('  dips to the puddle floor.  The branch is now countable:')
    print()
    print('  %-8s %-10s %-10s %-10s %s' % ('V_ds', 'V_g', 'frac n', 'frac p', 'CNP at V_ch ='))
    for Vds in (0.05, 0.1, 0.5):
        lo, hi = ambipolar_window(Vds)
        for Vg in (lo - 0.2, lo + 0.25 * Vds, lo + 0.5 * Vds, hi + 0.2):
            fn, fp, vcnp = channel_branch_census(Vg, Vds)
            print('  %-8.2f %-10.4f %-10.2f %-10.2f %s'
                  % (Vds, Vg, fn, fp,
                     ('%.4f V' % vcnp) if vcnp is not None else 'unipolar'))
    print()
    for Vds in (0.05, 0.1, 0.2, 0.5, 1.0):
        lo, hi = ambipolar_window(Vds)
        print('  V_ds = %.2f V : channel is ambipolar for %.4f V < V_g < %.4f V'
              ' (width = V_ds = %.2f V)' % (Vds, lo, hi, hi - lo))
    print()
    print('  The window width equals V_ds EXACTLY, for the reason stated in')
    print('  ambipolar_window(): dV vanishes somewhere in (0, V_ds) iff')
    print('  V_dirac < V_g < V_dirac + V_ds.  No sweep is needed to know it,')
    print('  and the census above agrees with it at every row.')
    wins = []
    for Vds in (0.05, 0.1, 0.2, 0.5, 1.0):
        lo, hi = ambipolar_window(Vds)
        mid = 0.5 * (lo + hi)
        fn, fp, vcnp = channel_branch_census(mid, Vds)
        wins.append(fn > 0.0 and fp > 0.0 and vcnp is not None)
        edge_lo = channel_branch_census(lo - 1e-6, Vds)
        edge_hi = channel_branch_census(hi + 1e-6, Vds)
        wins.append(edge_lo[2] is None and edge_hi[2] is None)
    check('census agrees with the closed-form ambipolar window at 5 biases',
          all(wins), 'mid-window biases are two-branch; both edges unipolar')

    # -- Section 7: the floor adjudication ---------------------------------
    print()
    print('-' * 78)
    print('SECTION 4.  Section 4.9 / Chapter 3 Section 3.3: the floor item,')
    print('            open since 2026-09-29 and untouched for ELEVEN sessions')
    print('-' * 78)
    s4 = sigma_min_ch4()
    print('  Chapter 4 Dirac-point floor, from this module:')
    print('    n_e = n_h = n_puddle/2 = %.6e 1/m^2 at dV = 0, EXACTLY'
          % (gfet.n_puddle / 2.0))
    print('    sigma_min = e*mu*(n_e + n_h) = e*mu*n_puddle = %.6e S/sq' % s4)
    print('    R_sheet   = %.6e Ohm/sq = %.4f kOhm/sq' % (1.0 / s4, 1e-3 / s4))
    print()
    print('  Chapter 3 Section 3.3 MEASURED floor, recomputed from CODATA:')
    print('    sigma_min = 4 q_e^2/h = %.6e S/sq' % SIGMA_MIN_CH3)
    print('    R_sheet   = %.4f kOhm/sq   (Section 3.3 quotes 6.45 kOhm/sq)'
          % (1e-3 / SIGMA_MIN_CH3))
    check('this module reproduces Section 3.3\'s quoted 6.45 kOhm/sq from CODATA',
          abs(1e-3 / SIGMA_MIN_CH3 - 6.45) < 0.01,
          'computed %.4f kOhm/sq against the 6.45 kOhm/sq in the chapter text'
          % (1e-3 / SIGMA_MIN_CH3))
    print()
    ratio = s4 / SIGMA_MIN_CH3
    n_match = n_puddle_matching_ch3()
    print('  *** THE DISAGREEMENT ***')
    print('    Chapter 4 is MORE CONDUCTIVE at the Dirac point than Chapter 3')
    print('    by a factor of %.4f.' % ratio)
    print('    Chapter 4 predicts %.4f kOhm/sq at neutrality where Chapter 3'
          % (1e-3 / s4))
    print('    measures %.4f kOhm/sq.' % (1e-3 / SIGMA_MIN_CH3))
    print('    The n_puddle that would reconcile them is %.6e 1/m^2 ='
          % n_match)
    print('    %.4e 1/cm^2, against the shipped %.4e 1/cm^2, i.e. a factor'
          % (n_match * 1e-4, gfet.n_puddle * 1e-4))
    print('    %.4f SMALLER.' % (gfet.n_puddle / n_match))
    print()
    print('    Printed with the factor first, 10-03\'s magnitude rule: a')
    print('    PASS/FAIL on "do the two floors agree" would have returned FAIL')
    print('    and been silent about whether the gap was 2x or 200x.')
    print()
    print('  WHY THE SIGNED READING IS WHAT MAKES THIS DECIDABLE.')
    print('    n_puddle enters (D2) as a quadrature floor on a MAGNITUDE, and a')
    print('    magnitude cannot say whether the floor is on the net density,')
    print('    on the total, or on each species.  Those three readings differ')
    print('    by factors of 1, 1 and 2 in sigma_min, which is the same order')
    print('    as the disagreement itself -- so before Section 1 the gap was')
    print('    not quotable.  Section 1 settles it: the floored quantity is')
    print('    n_e + n_h, the TOTAL.  Had n_puddle meant the per-species')
    print('    density the total would be 2*n_puddle and Chapter 4 would be')
    print('    %.4f x too conductive instead of %.4f x -- the discrepancy'
          % (2.0 * ratio, ratio))
    print('    moves the WRONG WAY under that reading, so the convention')
    print('    cannot be what explains it.')
    print()
    print('  WHAT IS NOT CLAIMED.  These two floors are not the same physics:')
    print('    4 q_e^2/h is a quantum/ballistic minimum conductivity observed')
    print('    experimentally, while e*mu*n_puddle is a diffusive disorder')
    print('    floor.  They are not required to be equal.  What IS a defect is')
    print('    that a thesis quotes both, in chapters that feed each other,')
    print('    without ever comparing them -- and that the diffusive one comes')
    print('    out BELOW the measured one, i.e. Chapter 4\'s channel is less')
    print('    resistive at neutrality than any measured graphene sheet.')
    check('Chapter 4 floor is on the more-conductive side of Chapter 3\'s',
          s4 > SIGMA_MIN_CH3,
          'the sign of the disagreement is fixed, not fitted: %.6e > %.6e S/sq'
          % (s4, SIGMA_MIN_CH3))

    # -- Section 5: device consequence, intrinsic vs terminal --------------
    print()
    print('-' * 78)
    print('SECTION 5.  What the factor %.4f costs the device numbers' % ratio)
    print('            (intrinsic vs terminal, the 2026-10-08 structure)')
    print('-' * 78)
    Vds = 0.1
    Vg_on = 3.5
    for label, fn in (('terminal (R_c = %.1f Ohm)' % gfet.Rc_total, onoff_at),
                      ('intrinsic (R_c = 0)', intrinsic_onoff_at)):
        I_on_a, I_off_a, r_a = fn(Vds, Vg_on, gfet.n_puddle)
        I_on_b, I_off_b, r_b = fn(Vds, Vg_on, n_match)
        print('  %s' % label)
        print('    shipped   n_puddle: I_on %.6e A, I_off %.6e A, on/off %.4f'
              % (I_on_a, I_off_a, r_a))
        print('    Ch.3-matched      : I_on %.6e A, I_off %.6e A, on/off %.4f'
              % (I_on_b, I_off_b, r_b))
        print('    I_off changes by %+.4f %% ; on/off ratio changes by %+.4f %%'
              % (100.0 * (I_off_b / I_off_a - 1.0), 100.0 * (r_b / r_a - 1.0)))
        print('    on/off multiplier  = %.4f x' % (r_b / r_a))
        print()
    I_on_t, I_off_t, r_t = onoff_at(Vds, Vg_on, gfet.n_puddle)
    I_on_t2, I_off_t2, r_t2 = onoff_at(Vds, Vg_on, n_match)
    I_on_i, I_off_i, r_i = intrinsic_onoff_at(Vds, Vg_on, gfet.n_puddle)
    I_on_i2, I_off_i2, r_i2 = intrinsic_onoff_at(Vds, Vg_on, n_match)
    mult_t = r_t2 / r_t
    mult_i = r_i2 / r_i
    print('  The %.4f x floor error reaches the on/off ratio as %.4f x'
          % (ratio, mult_i))
    print('  intrinsically and only %.4f x at the terminals: the contact' % mult_t)
    print('  resistance attenuates it by a further %.4f x.  Same negative-'
          % (mult_i / mult_t))
    print('  feedback structure 10-08 measured for the charge relation, and the')
    print('  inequality below is fixed by the circuit rather than by a fit.')
    check('contact resistance attenuates the floor error at the terminals',
          mult_t < mult_i,
          'terminal multiplier %.4f < intrinsic multiplier %.4f'
          % (mult_t, mult_i))

    # -- Section 5b: the check that was written as a sanity test and FAILED --
    print()
    print('-' * 78)
    print('SECTION 5b.  Section 7.9 item 11 answered: n_puddle is NOT a')
    print('             near-Dirac regularizer, and the sanity check that')
    print('             assumed it was is the thing that found out')
    print('-' * 78)
    d_on = 100.0 * (I_on_t2 / I_on_t - 1.0)
    print('  This section exists because of a FAILED assertion.  Section 5 was')
    print('  written with the guard "I_on is essentially untouched by the')
    print('  floor", on the docstring\'s word that n_puddle "regularizes n at')
    print('  the Dirac point" and "regularizes the conductivity minimum".')
    print('  It failed: changing the floor moved I_on by %+.4f %%.' % d_on)
    print()
    print('  The reason, measured rather than reasoned:')
    print('  %-10s %-16s %-16s %-12s %s'
          % ('V_g', '|n_net| (1/m^2)', 'n_puddle/|n_net|', 'floor adds', 'regime'))
    rows = []
    for Vg in (gfet.V_dirac + 0.1, 1.5, 2.0, 2.5, 3.0, 3.5):
        nnv = float(np.abs(n_net(Vg, 0.0)))
        Sv = float(np.sqrt(nnv ** 2 + gfet.n_puddle ** 2))
        add = 100.0 * (Sv / nnv - 1.0) if nnv > 0 else float('inf')
        rows.append((Vg, nnv, gfet.n_puddle / nnv, add))
        print('  %-10.2f %-16.6e %-16.4f %-12s %s'
              % (Vg, nnv, gfet.n_puddle / nnv, '%+.3f %%' % add,
                 'near-Dirac' if gfet.n_puddle / nnv > 1 else 'ON STATE'))
    ratio_on = gfet.n_puddle / float(np.abs(n_net(3.5, 0.0)))
    add_on = 100.0 * (float(np.sqrt(np.abs(n_net(3.5, 0.0)) ** 2
                                    + gfet.n_puddle ** 2))
                      / float(np.abs(n_net(3.5, 0.0))) - 1.0)
    print()
    print('  At V_g = 3.5 V -- the MAXIMUM overdrive swept anywhere in this')
    print('  thesis -- the "residual" puddle density is %.4f of the' % ratio_on)
    print('  gate-induced net density and still adds %+.3f %% to the total.' % add_on)
    print('  A density that is %.0f %% of the signal at full drive is not a'
          % (100.0 * ratio_on))
    print('  regularizer; it is a parallel conduction channel that the model')
    print('  carries at every bias.  That is Section 7.9 item 11 -- "whether')
    print('  n_puddle removes other regimes", created 10-07 -- answered: it')
    print('  does not remove regimes, it ADDS a floor to all of them, and the')
    print('  docstring that calls it a Dirac-point regularizer is wrong about')
    print('  its own device.')
    # The drive at which the floor stops mattering, by bisection.
    #
    # DEFECT FOUND AND FIXED WHILE WRITING THIS: the first version of this
    # loop read `contrib = ... if nnv > 0 else 1e9`, and above dV = 18.3493 V
    # quantum_capacitance() returns inf, so C_q/(C_q + C_ox) is inf/inf = nan,
    # nnv is nan, `nan > 0` is False, and the NaN branch was scored as "the
    # floor dominates".  The bisection then walked AWAY from the answer and
    # reported the bracket ceiling, 200 V, as if it were a result -- a number
    # that was not a measurement of anything.  That is Section 7.9 item 13,
    # the cosh overflow logged on 10-08 as latent with "no committed number
    # affected", reaching a number on its FIRST use outside the module that
    # owns it.  Latent and harmless are not the same property.
    #
    # The bracket is now capped below the overflow onset and NaN is a hard
    # stop rather than a score.
    OVERFLOW_ONSET = 18.3493          # 10-08, graphene_potential_declaration
    lo, hi = 0.0, 0.95 * OVERFLOW_ONSET
    hi_nnv = float(np.abs(n_net(gfet.V_dirac + hi, 0.0)))
    bracket_ok = np.isfinite(hi_nnv) and (
        np.sqrt(hi_nnv ** 2 + gfet.n_puddle ** 2) / hi_nnv - 1.0) <= 0.01
    if bracket_ok:
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            nnv = float(np.abs(n_net(gfet.V_dirac + mid, 0.0)))
            if not np.isfinite(nnv) or nnv <= 0.0:
                raise RuntimeError('non-finite n_net inside a bracket that was '
                                   'checked to be finite: dV = %r' % mid)
            contrib = np.sqrt(nnv ** 2 + gfet.n_puddle ** 2) / nnv - 1.0
            if contrib > 0.01:
                lo = mid
            else:
                hi = mid
    print()
    print('  The overdrive at which the floor\'s contribution to the total')
    if bracket_ok:
        print('  falls below 1 %%: dV = %.4f V, i.e. %.2f x the largest overdrive'
              % (hi, hi / 2.7))
        print('  this thesis sweeps (2.7 V).  The regime in which the')
        print('  docstring\'s description is true is a regime the thesis never')
        print('  enters -- and it is only %.2f V below where'
              % (OVERFLOW_ONSET - hi))
        print('  quantum_capacitance() stops returning a number at all.')
    else:
        print('  falls below 1 %% lies ABOVE the %.4f V cosh-overflow onset,'
              % OVERFLOW_ONSET)
        print('  so it is not computable in this model at all.')
    check('the 1 percent-floor crossover is bracketed inside the FINITE domain',
          bracket_ok,
          'bracket ceiling %.4f V, set at 0.95 x the %.4f V overflow onset '
          'of Section 7.9 item 13' % (0.95 * OVERFLOW_ONSET, OVERFLOW_ONSET))
    print('         The first version of this bisection scored the NaN above')
    print('         the overflow as "floor dominates" and returned its own')
    print('         bracket ceiling, 200 V, as the answer.  Item 13 was')
    print('         recorded on 10-08 as latent with no committed number')
    print('         affected; its first use outside its own module produced a')
    print('         wrong number in the first draft of this section.')
    check('n_puddle is comparable to the on-state net density (item 11)',
          ratio_on > 0.5,
          'n_puddle/|n_net| = %.4f at V_g = 3.5 V, so the floor is a global '
          'additive density rather than a near-Dirac patch' % ratio_on)
    check('and that is what moved I_on when only the FLOOR was changed',
          abs(d_on) > 1.0,
          'I_on moved %+.4f %% under a change to an allegedly OFF-state-only '
          'parameter' % d_on)
    print('         Both of these started as one assertion that the floor')
    print('         does NOT reach the on state.  09-28: prose is a detector.')
    print('         Today: SO IS A SANITY CHECK YOU EXPECTED TO PASS.')

    # -- Section 6: nothing committed moved --------------------------------
    print()
    print('-' * 78)
    print('SECTION 6.  Output neutrality: no committed number moved today')
    print('-' * 78)
    Vg_sweep = np.linspace(-2.0, 3.5, 41)
    Id_now = gfet.transfer_characteristic(Vg_sweep, Vds=0.1)
    n_now = gfet.carrier_density(Vg_sweep, 0.0)
    print('  graphene_fet_model is imported, read and monkey-patched by this')
    print('  module (Section 5 sets n_puddle and Rc_total inside try/finally).')
    print('  Both are restored, and the shipped sweep is re-run afterwards:')
    print('    n_puddle  = %.6e 1/m^2 (restored)' % gfet.n_puddle)
    print('    Rc_total  = %.6f Ohm (restored)' % gfet.Rc_total)
    print('    I_d(V_g = 3.5 V, V_ds = 0.1 V) = %.9e A' % Id_now[-1])
    print('    n(V_g = 3.5 V)                = %.9e 1/m^2' % n_now[-1])
    check('module-level constants restored exactly after patching',
          gfet.n_puddle == 5e15 and gfet.Rc_total == 2 * (300.0 * 1e-6) / gfet.W,
          'n_puddle and Rc_total compare equal to their shipped literals')
    print()
    print('  10-02 established byte-identity is not correctness.  It is the')
    print('  right instrument HERE because this module adds no equation to the')
    print('  shipped path: every new quantity is a decomposition of one that')
    print('  was already being computed.')

    make_figure()
    print()
    print('-' * 78)
    print('Figure written: fet_signed_carrier_model.png')
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('-' * 78)
    return failed


if __name__ == '__main__':
    raise SystemExit(1 if main() else 0)
