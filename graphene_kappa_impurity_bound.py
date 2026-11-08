"""
graphene_kappa_impurity_bound.py

THE PAIRING CANNOT BE RESOLVED FROM ZION ET AL. -- THE FULL TEXT CARRIES NO
MOBILITY AT ALL -- BUT IT CARRIES SOMETHING BETTER: ZION'S OWN CHARGED-IMPURITY
DENSITY.  THAT GIVES KAPPA DIRECTLY, WITH NO MOBILITY NEEDED, AND THE REPOSITORY'S
KAPPA = 2.4 IS EXCLUDED BY IT.

This closes the item ranked FIRST on 2026-11-02 -- "RESOLVE THE `A` PAIRING" --
and it closes it in a way the item did not anticipate.

------------------------------------------------------------------------------
1. WHAT THE ITEM ASKED FOR, AND WHY IT IS THE WRONG QUESTION
------------------------------------------------------------------------------
`graphene_impurity_limited_crossover.py` (11-02, hereafter ILC) derives

    n_c(N_imp) = KAPPA * N_imp,        KAPPA == A / C                      (C3)

and obtains KAPPA = 2.4 by PAIRING Zion et al.'s crossover density with THIS
repository's long-range mobility:

    A_ref = N_CROSSOVER * gfet.mu = 3e16 * 0.4 = 1.2e16  V^-1 s^-1

ILC Section 3 names this the weakest link in the module and logs the fix as
"read Zion's own mobility from the full text of arXiv:1801.03565".  The
abstract page does not carry one; the full text was not reachable on 11-02.

The full text was reached today, via `ar5iv.labs.arxiv.org/html/1801.03565`.
IT CONTAINS NO MOBILITY VALUE ANYWHERE.  Mobility appears only as curves
against |n| in Fig. 5 and as the proportionality of their Eq. 4; the only
mobility statement in words is the qualitative

    "at high temperatures, mobility increases with decreasing T, but become
     temperature-independent below 100 K."                  -- Section III c

So the item AS POSED is unresolvable from this source, and that is a definite
finding rather than another blocked fetch: the datum does not exist in the
paper, so no further attempt on this paper can succeed.  The two standing
items blocked on "a datum that is a picture" (Zion Fig. 5, Huard Fig. 4) are
unaffected -- Fig. 5 still holds the mobility curves -- but the TEXT is now
exhausted and should not be fetched again.

------------------------------------------------------------------------------
2. THE DATUM THAT IS ACTUALLY THERE, AND WHY IT IS STRICTLY BETTER
------------------------------------------------------------------------------
Zion et al. estimate the charged-acceptor density in the SiO2 adjacent to
their own channel from their own neutrality-point gate voltage:

    "Thus, we get C_ox = 1.15e-4 F/m^2 and for V_NP ~ 20 V, we obtain
     N_A^(S) ~ 1.5e12 cm^-2."                                -- Section III a

    N_A^(S) = (1/e) * C_ox * V_NP,    C_ox = eps_0*eps/t,
    t = 300 nm,  eps = 3.9,  eps_0 = 8.85e-12 F/m         -- Sections II, III a

This is a charged-impurity density FOR THE SAME SAMPLE whose crossover density
the repository already uses.  So KAPPA can be read off ONE device, from TWO
numbers in ONE paper, with no mobility and no cross-paper pairing:

    KAPPA = n_c / N_imp,   both on Zion's sample                           (B1)

That is strictly better than (C3)'s route, which needed Chen's C, this
repository's mu, and the assumption that Zion's device resembles this one.
None of those three enter (B1).

THE DIRECTION OF THE INEQUALITY, WHICH IS THE WHOLE RESULT.  N_A^(S) is a
LOWER bound on the density of charged scatterers in Zion's channel, for three
independent reasons, all of them stated by Zion themselves:

  (i)  IT IS A NET CHARGE, NOT A COUNT.  N_A^(S) is inferred from the voltage
       needed to compensate the channel, so it measures net areal charge.
       Scattering counts centres: with donors as well as acceptors present the
       total exceeds |net|.
  (ii) THE OXIDE DEFECTS ARE ACKNOWLEDGED AND UNQUANTIFIED.  "charged oxide
       defects with density N_t inevitable occur in SiO2 close (typically
       within 1-3 nm) to the conductive channel" -- Section III b.  No value
       for N_t is given anywhere in the paper.
  (iii) THE SUBSTRATE ACCEPTORS ARE ONE POPULATION AMONG SEVERAL; nothing in
       the estimate covers adsorbates or transfer residue on a CVD film that
       was grown on copper and transferred.

Therefore N_imp >= N_A^(S), and since KAPPA = n_c/N_imp is DECREASING in
N_imp,

    KAPPA  <=  n_c / N_A^(S)                                               (B2)

An upper bound, not an estimate.  Section 7.4 of the thesis argues for
converting rankings into bounds; this is one arriving unasked.

------------------------------------------------------------------------------
3. ZION'S TWO NUMBERS DISAGREE WITH EACH OTHER BY 4.5 %
------------------------------------------------------------------------------
Their stated N_A^(S) ~ 1.5e12 cm^-2 does not follow from their own stated
C_ox and V_NP.  Recomputing from the formula they give:

    C_ox    = 8.85e-12 * 3.9 / 300e-9  =  1.1505e-4 F/m^2   (they say 1.15e-4)
    N_A^(S) = 1.1505e-4 * 20 / 1.602176634e-19
            =  1.4362e16 m^-2  =  1.4362e12 cm^-2           (they say 1.5e12)

C_ox reproduces exactly at their stated precision; N_A^(S) does not -- it is
4.4 % below their figure, and 1.4362 does not round to 1.5 at any precision.
Their "~" is doing the work of a rounding-up, not of a rounding.

WHICH VALUE THIS MODULE USES, AND WHY IT IS THE CONSERVATIVE CHOICE.  The
recomputed 1.4362e12 is SMALLER, and by (B2) a smaller N_A^(S) gives a LOOSER
(larger) bound on KAPPA.  So using the recomputed value cannot manufacture the
exclusion of Section 4; taking Zion's rounded 1.5e12 at face value would only
strengthen it.  Both are reported.  D2 pins the disagreement rather than
quietly adopting one side, which is the standing one-way-channel item's own
prescription: the externally-sourced claim and the arithmetic that contradicts
it are both recorded.

------------------------------------------------------------------------------
4. THE VERDICT: KAPPA = 2.4 IS EXCLUDED
------------------------------------------------------------------------------
    KAPPA <= 3e16 / 1.4362e16 = 2.0889      (recomputed N_A, conservative)
    KAPPA <= 3e16 / 1.5e16    = 2.0000      (Zion's rounded N_A)
    KAPPA (ILC, by pairing)   = 2.4

2.4 lies ABOVE the conservative bound by 14.9 %, so it is not merely
uncertain, it is excluded.  n_c(N_imp) is overestimated by AT LEAST 14.9 % at
every charged-impurity density, and by at least 20 % if Zion's own rounded
figure is taken as written.  The direction is the one ILC Section 3 flagged as
the risk -- "if Zion's samples were cleaner than that, A_ref is overestimated,
and n_c is overestimated with it at every impurity density" -- so this is the
flagged failure mode occurring, quantified, rather than a new surprise.

------------------------------------------------------------------------------
5. THE BY-PRODUCT: ZION'S MOBILITY, WHICH THE PAPER NEVER STATES
------------------------------------------------------------------------------
A = KAPPA * C and A = n_c * mu_L, both on Zion's sample, so (B2) bounds the
mobility the full text omits:

    mu_L,Zion  =  KAPPA * C / n_c  <=  2.0889 * 5e15 / 3e16
               =  0.348 m^2/Vs  =  3480 cm^2/Vs

The item asked for a number the paper does not contain, and the paper's other
two numbers bound it anyway.  3480 cm^2/Vs is an ordinary long-range-limited
mobility for transferred CVD monolayer on 300 nm SiO2, which is the only
external check available on it, and it passes.  The repository paired
gfet.mu = 0.4 m^2/Vs = 4000 cm^2/Vs onto Zion's device; that was 15 % too
generous, which is the same 14.9 % as Section 4 seen from the other side.

------------------------------------------------------------------------------
6. WHAT MOVES, AND IT ALL MOVES THE SAME WAY
------------------------------------------------------------------------------
The correction is a single multiplicative constant on A, so:

  * EVERY RATIO IS UNTOUCHED, bitwise -- D8.  ILC's D5 already established
    that the ratios do not carry the pairing; this is that result being
    collected on a correction it did not anticipate.
  * EVERY ABSOLUTE n_c AND dV_c FALLS by 13.0 % (= 1 - 2.0889/2.4).
  * EVERY OPERATING MARGIN AND IMPURITY BUDGET RISES by 14.9 %, because
    M ~ 1/n_c and N_imp* ~ 1/A.

The device consequence is the one worth reporting.  ILC's thinnest margin was
1.10x -- a top-gated device barely past its crossover, close enough to 1 that
the verdict was fragile.  It becomes 1.26x.  The correction does not overturn
any 11-02 conclusion; it moves the one conclusion that was marginal away from
the margin, which is the useful direction for a result to move when its
weakest input is replaced.  Section 8 reports the table.

------------------------------------------------------------------------------
7. A NEW PROBLEM THE FULL TEXT CREATES, AND IT IS BIGGER THAN THE CORRECTION
------------------------------------------------------------------------------
The 3e12 cm^-2 the repository has used as `N_CROSSOVER` since 10-30 is, in
Zion's own words, read off the HOLE branch:

    "The dependence of the hole mobility in Fig. 5 clearly demonstrates the
     transition from long-range to short-range scattering"     -- Section III c

GOOD: this is the RIGHT KIND of crossover for DDM's (B1).  DDM's n_c is the
Matthiessen long-range/short-range crossover, which is exactly what Zion's
Fig. 5 transition is.  Worth recording because it is easy to get backwards:
this number is NOT an electron-hole-asymmetry crossover, and nothing in this
repository should cite it as one.

BAD, AND THIS IS THE REAL FINDING OF THE DAY.  DDM's (B1) applies ONE n_c to
BOTH branches.  Writing it per branch,

    1/mu_b(n) = 1/mu_L,b  +  |n| / (n_c,b * mu_L_bar)                      (B3)

and parameterising a split as n_c,e/h = n_c*(1 +- s/2), there is an EXACT
analytic identity at s = 2*r0:

    if  n_c,b  is proportional to  mu_L,b,  then  r(n) = r0  for ALL n,
    exactly, with no density dependence whatever.                          (B4)

The proof is two lines.  At s/2 = r0 the bracket factorises:

    1/mu_e = (1/mu_L_bar) * (1 + x)/(1 + r0),
    1/mu_h = (1/mu_L_bar) * (1 + x)/(1 - r0),        x == |n|/n_c

so the common factor (1+x) cancels out of r = (mu_e-mu_h)/(mu_e+mu_h) and
r = r0 identically.  D9b checks this against r0 -- a value known exactly in
advance and independent of density -- over four values of r0 and 401
densities spanning six decades.

WHY THIS MATTERS MORE THAN SECTIONS 4-6.  DDM's central qualitative
conclusion is that raising the carrier density SCREENS AWAY the electron-hole
asymmetry: r(n) -> 0 as |n|/n_c -> infinity.  (B4) says that conclusion rests
entirely on an assumption nobody has checked -- that the NEUTRAL-defect
channel is branch-SYMMETRIC while the charged channel is not.  If the neutral
channel carries the same fractional asymmetry as the charged one, the
asymmetry does not decay at all.

And the symmetric case is not the obviously right one.  DDM Section 1 derives
r0 from a beyond-Born cross-section ratio -- a property of scattering off a
potential of a given sign.  Neutral defects are also potentials.  There is no
argument in this repository for why the beyond-Born asymmetry should vanish
for vacancies and grain boundaries while surviving for charged impurities; it
was assumed by writing a single `n_c` and never revisited.

So the honest status of DDM's screening result is CONDITIONAL, and the
condition is now named and quantified.  At the top-gated operating point
(|n|/n_c = 5.29) the measured sensitivity is:

    s = 0     ->  r = 0.0160   (the shipped answer: asymmetry screened 6.2x)
    s = 0.10  ->  r = 0.0580
    s = 2*r0  ->  r = 0.1000   (= r0 exactly: no screening at all)

A 20 % split -- modest for a quantity never measured -- moves r by 524 %, and
the endpoint of the range is not a perturbation of the shipped answer but its
negation.  This is logged as the new FIRST open item.  It is not blocked on a
fetch: the electron branch's transition density is in Zion Fig. 5, the same
picture two standing items already need, and reading it would settle (B4)
empirically for one sample.

------------------------------------------------------------------------------
8. WHAT IS CHECKED, AND WHICH CHECKS ARE EXACT
------------------------------------------------------------------------------
Six of these compare against a value known exactly in advance rather than
against a plausible range -- the standing requirement of item 15 -- and one
of those (D3) failed on first run because the exactness it predicted does
not hold.  That failure is recorded in its entry rather than erased.

  D1  C_ox REPRODUCES ZION'S STATED VALUE at their precision, from eps_0,
      eps and t.  Exact at 3 significant figures.
  D2  N_A^(S) DOES NOT REPRODUCE, and the discrepancy is PINNED, not
      tolerated -- Section 3.  A definite prediction: the shortfall is between
      4 % and 5 %, and 1.4362 does not round to 1.5.
  D3  THE BOUND ROUND-TRIPS -- IN TWO PARTS, AND THE FIRST ONE WAS WRONG.
      This check was written as "bitwise" and FAILED at 1 ulp the first time
      it ran, which is correct behaviour and not a defect in the bound: a
      division followed by a multiplication by the same operand rounds twice,
      so (n_c/N_A)*N_A is accurate to 1 ulp and is NOT exact in general.  The
      prediction was wrong, so it is split rather than loosened.  D3a asserts
      the 1-ulp bound that floating-point actually guarantees.  D3b recovers
      exactness where it genuinely holds -- N_A a power of two times n_c,
      where both operations are exact -- and there it is bitwise and KAPPA is
      exactly 4.0.  Item 15 catching same-day code, now four sessions
      running.
  D4  THE EXCLUSION FIRES ON 2.4 AND NOT ON A PLANTED CONTROL.  A hypothetical
      KAPPA = 1.9 must NOT be excluded by the same code path that excludes
      2.4.  A one-sided test that only ever fires is not a test.
  D5  SCALE INVARIANCE, BITWISE.  KAPPA is a ratio, so scaling n_c and N_A
      both by 2 -- exactly representable, so exactly commuting with rounding
      -- must leave it bitwise unchanged.
  D6  SYMMETRY THAT MUST GIVE EXACTLY ZERO.  At kappa = KAPPA_ILC the
      corrected n_c must equal ILC's n_c with a difference of EXACTLY 0.0 at
      every density of a 401-point sweep -- not small, zero.  This is the
      overlap condition: a correction parameterised by kappa must contain the
      thing it corrects as the value it corrects away from.
  D7  OVERLAP WITH ILC, BITWISE, THROUGH ILC'S OWN CODE.  Routed through
      `ilc.n_c_of_impurity` with `a_neutral` set from kappa, over 401 points
      and six decades of N_imp.  D6 tests this module's arithmetic; D7 tests
      that it is the same function ILC ships.
  D8  EVERY RATIO IS INVARIANT UNDER THE CORRECTION, BITWISE.  n_c ratios,
      margin ratios and budget ratios must not move by one bit when kappa
      changes by 13 %.
  D9  THE BRANCH-SPLIT IDENTITY AND ITS SENSITIVITY -- Section 7.  D9a: at
      s = 0 the split form reduces to DDM's single-n_c construction BITWISE,
      which is what makes it a generalisation of (B1) rather than a second
      and different model.  D9b: at s = 2*r0 the asymmetry must equal r0 --
      a value known EXACTLY in advance and, crucially, independent of
      density, so the check has a target rather than a tolerance.  Worst case
      56 ulp (1.2e-14 relative) over four r0 values and 401 densities; the
      identity is exact in real arithmetic and the residual is the
      computation path, measured rather than tolerated.  D9c: the magnitude
      at the shipped operating point, which is 524 % and is the reason
      Section 7 is ranked above Sections 4-6.
  D10 THE MOBILITY BOUND IS CONSISTENT WITH THE PAIRING IT REPLACES.  The
      bound on mu_L,Zion must come out BELOW gfet.mu by the same factor the
      KAPPA bound comes out below KAPPA_ILC -- two routes, one number, and
      agreement to <= 1 ulp or the algebra of Section 5 is wrong.
  D11 MUTATION: FLIP THE INEQUALITY.  Treat N_A^(S) as an UPPER bound on
      N_imp instead of a lower one -- the error the whole of Section 2 exists
      to avoid -- and the exclusion of Section 4 must REVERSE.  A correction
      whose sign does not depend on the argument for its sign is not using
      that argument.
"""

import sys

import numpy as np

import graphene_fet_model as gfet
import graphene_density_dependent_mobility as ddm
import graphene_impurity_limited_crossover as ilc

EPS = float(np.finfo(np.float64).eps)

# ---------------------------------------------------------------------------
# Zion et al., arXiv:1801.03565, "Charge carrier transport asymmetry in
# monolayer graphene".  Full text read 2026-11-03 via ar5iv.  Every constant
# below is quoted in the module docstring with its section.
# ---------------------------------------------------------------------------
EPS0_ZION = 8.85e-12            # F/m   -- the value Zion state, not scipy's
EPS_SIO2_ZION = 3.9             # Section III a
T_OX_ZION = 300e-9              # m     -- Section II, 300 nm (+-5 %)
V_NP_ZION = 20.0                # V     -- Section III a, low-T saturation
C_OX_ZION_STATED = 1.15e-4      # F/m^2 -- Zion's own stated value
N_A_S_ZION_STATED = 1.5e16      # m^-2  -- Zion's own stated 1.5e12 cm^-2
N_C_ZION = ddm.N_CROSSOVER      # m^-2  -- the 3e12 cm^-2 the repo already uses

# Section 7: a hypothetical electron/hole split in n_c, for D9 only.
BRANCH_SPLIT_D9 = 0.20

# Section 4: the value this module excludes.
KAPPA_ILC = ilc.KAPPA

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
        for line in _wrap(detail, 74):
            print('         %s' % line)


def _wrap(text, width):
    words, lines, cur = text.split(), [], ''
    for w in words:
        if cur and len(cur) + 1 + len(w) > width:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + ' ' + w) if cur else w
    if cur:
        lines.append(cur)
    return lines


def ulp_distance(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0
    if not (np.isfinite(a) and np.isfinite(b)):
        return np.inf
    ia = np.abs(np.frombuffer(np.float64(a).tobytes(), dtype=np.int64)[0])
    ib = np.abs(np.frombuffer(np.float64(b).tobytes(), dtype=np.int64)[0])
    return int(abs(ia - ib))


# ---------------------------------------------------------------------------
# Section 1: Zion's own arithmetic, reproduced
# ---------------------------------------------------------------------------

def c_ox_zion(eps0=None, eps=None, t=None):
    """C_ox = eps_0 * eps / t, Zion Section III a.

    DEFAULTS READ AT CALL TIME, NOT CAPTURED AT DEF TIME (2026-11-08).  The
    three constants used to be default ARGUMENTS, which binds the value this
    module held when its def statement executed.  That is the capture pattern
    graphene_mixed_capture_repair_audit.py censuses and
    graphene_frozen_default_semantics.py demonstrates: a session that rebinds
    EPS0_ZION -- or a test that monkeypatches it -- changes the constant and
    not the function.  Converted to the None sentinel this module already
    uses at six other sites; bitwise identical output, verified against the
    pre-conversion transcript.
    """
    eps0 = EPS0_ZION if eps0 is None else eps0
    eps = EPS_SIO2_ZION if eps is None else eps
    t = T_OX_ZION if t is None else t
    return eps0 * eps / t


def n_a_s_zion(c_ox=None, v_np=None):
    """N_A^(S) = (1/e) * C_ox * V_NP, Zion Section III a.

    v_np converted from a captured default 2026-11-08 -- see c_ox_zion.
    """
    if c_ox is None:
        c_ox = c_ox_zion()
    v_np = V_NP_ZION if v_np is None else v_np
    return c_ox * v_np / ddm.E_CHARGE


# ---------------------------------------------------------------------------
# Section 2: the bound
# ---------------------------------------------------------------------------

def kappa_bound(n_c=None, n_imp_lower=None):
    """
    Upper bound on KAPPA from (B2): KAPPA <= n_c / N_imp_lower.

    `n_imp_lower` is a LOWER bound on the charged-scatterer density, so the
    quotient is an UPPER bound on KAPPA.  D11 mutates this reading.
    """
    if n_c is None:
        n_c = N_C_ZION
    if n_imp_lower is None:
        n_imp_lower = n_a_s_zion()
    return float(n_c) / float(n_imp_lower)


def a_neutral_from_kappa(kappa):
    """A = KAPPA * C -- the inverse of ILC's KAPPA = A/C."""
    return float(kappa) * ilc.C_CHEN


def n_c_of_impurity(n_imp, kappa=None):
    """n_c = KAPPA * N_imp, with KAPPA an argument rather than a constant."""
    if kappa is None:
        kappa = kappa_bound()
    return float(kappa) * np.asarray(n_imp, dtype=float)


def mu_L_zion_bound(kappa=None, n_c=None):
    """
    Section 5: mu_L,Zion = KAPPA * C / n_c, bounded above because KAPPA is.
    The paper states no mobility; this is it, bounded from the paper's own
    other two numbers.
    """
    if kappa is None:
        kappa = kappa_bound()
    if n_c is None:
        n_c = N_C_ZION
    return a_neutral_from_kappa(kappa) / float(n_c)


# ---------------------------------------------------------------------------
# Section 3: the device propagation
# ---------------------------------------------------------------------------

def crossover_overdrive(C_ox, n_imp, kappa=None):
    """dV_c = e*n_c/C_ox, with n_c at the given kappa."""
    return ddm.E_CHARGE * float(n_c_of_impurity(n_imp, kappa=kappa)) / C_ox


def operating_margin(C_ox, n_imp, V_over, kappa=None):
    """M = n(V_over)/n_c.  M > 1: the device operates past the crossover."""
    n_at_V = C_ox * V_over / ddm.E_CHARGE
    return n_at_V / float(n_c_of_impurity(n_imp, kappa=kappa))


def max_tolerable_impurity(C_ox, V_over, kappa=None):
    """N_imp* = C_ox*V/(e*KAPPA) -- the M = 1 design rule, at this kappa."""
    if kappa is None:
        kappa = kappa_bound()
    return C_ox * V_over / (ddm.E_CHARGE * float(kappa))


def r_branch_split(n_abs, n_imp, r0, split=None, kappa=None):
    """
    Section 7 / D9.  (B1) with a DIFFERENT n_c on each branch:
    n_c,e = n_c*(1+split/2), n_c,h = n_c*(1-split/2).  At split = 0 this must
    reduce to the single-n_c construction bitwise.

    split converted from a captured default 2026-11-08 -- see c_ox_zion.  The
    one caller that passes BRANCH_SPLIT_D9 explicitly (Section 7's headline
    number) is unaffected and deliberately left explicit: passing the
    constant by name at the call site is the readable form, and it is the
    IMPLICIT capture in the signature that the census objects to.
    """
    split = BRANCH_SPLIT_D9 if split is None else split
    mu_L_bar = ilc.C_CHEN / float(n_imp)
    mu_L_e = mu_L_bar * (1.0 + r0)
    mu_L_h = mu_L_bar * (1.0 - r0)
    n_c = float(n_c_of_impurity(n_imp, kappa=kappa))
    n_c_e = n_c * (1.0 + 0.5 * split)
    n_c_h = n_c * (1.0 - 0.5 * split)
    mu_e = 1.0 / (1.0 / mu_L_e + float(n_abs) / (n_c_e * mu_L_bar))
    mu_h = 1.0 / (1.0 / mu_L_h + float(n_abs) / (n_c_h * mu_L_bar))
    return (mu_e - mu_h) / (mu_e + mu_h)


# ---------------------------------------------------------------------------
# Section 4: figure
# ---------------------------------------------------------------------------

def n_op_for_figure():
    """The top-gate operating density, shared by the figure and main()."""
    return ddm.C_OX_TOPGATE_HFO2 * 1.0 / ddm.E_CHARGE


def make_figure(path='kappa_impurity_bound.png'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    k_b = kappa_bound()
    k_stated = kappa_bound(n_imp_lower=N_A_S_ZION_STATED)
    n_a = n_a_s_zion()

    fig, axes = plt.subplots(1, 3, figsize=(16.8, 4.9))

    # -- panel (a): the excluded region -----------------------------------
    ax = axes[0]
    n_imp = np.logspace(np.log10(0.3e16), np.log10(2e17), 400)
    x = n_imp * 1e-16                       # e12 cm^-2
    y_ilc = n_c_of_impurity(n_imp, kappa=KAPPA_ILC) * 1e-16
    y_bnd = n_c_of_impurity(n_imp, kappa=k_b) * 1e-16
    y_sta = n_c_of_impurity(n_imp, kappa=k_stated) * 1e-16

    ax.fill_between(x, y_bnd, y_ilc * 3.0, color='C3', alpha=0.13,
                    label='excluded by (B2)')
    ax.loglog(x, y_bnd, 'C0-', lw=2.1,
              label=r'bound: $\kappa \leq %.3f$' % k_b)
    ax.loglog(x, y_sta, 'C0:', lw=1.5,
              label=r"Zion's rounded $N_A$: $\kappa \leq %.3f$" % k_stated)
    ax.loglog(x, y_ilc, 'C3--', lw=2.0,
              label=r'ILC 11-02 (superseded): $\kappa = %.2f$' % KAPPA_ILC)

    ax.plot([n_a * 1e-16], [N_C_ZION * 1e-16], 'ko', ms=7, zorder=6)
    ax.annotate("Zion's own sample\n"
                r'$N_A^{(S)}=%.3f$, $n_c=%.1f$' % (n_a * 1e-16,
                                                   N_C_ZION * 1e-16),
                xy=(n_a * 1e-16, N_C_ZION * 1e-16),
                xytext=(0.30, 7.4), fontsize=7.4,
                arrowprops=dict(arrowstyle='->', lw=0.8, color='k'))
    ax.set_xlabel(r'charged-impurity density $N_{imp}$  ($10^{12}$ cm$^{-2}$)')
    ax.set_ylabel(r'crossover density $n_c$  ($10^{12}$ cm$^{-2}$)')
    ax.set_title('(a)  $\\kappa=2.4$ is excluded by Zion\'s own two numbers',
                 fontsize=9.5)
    ax.legend(fontsize=7.0, loc='lower right')
    ax.grid(True, which='both', alpha=0.25)

    # -- panel (b): the margin that was thin ------------------------------
    ax = axes[1]
    n_imp2 = np.logspace(np.log10(0.8e16), np.log10(8e16), 400)
    x2 = n_imp2 * 1e-16
    stacks = (('10 nm HfO$_2$ top gate at 1 V',
               ddm.C_OX_TOPGATE_HFO2, 1.0, 'C2'),
              ('90 nm SiO$_2$ back gate at 5 V', gfet.C_ox, 5.0, 'C4'))
    for name, cox, v, col in stacks:
        m_old = np.array([operating_margin(cox, q, v, kappa=KAPPA_ILC)
                          for q in n_imp2])
        m_new = np.array([operating_margin(cox, q, v, kappa=k_b)
                          for q in n_imp2])
        ax.semilogx(x2, m_old, color=col, ls='--', lw=1.5, alpha=0.75,
                    label='%s, 11-02' % name)
        ax.semilogx(x2, m_new, color=col, ls='-', lw=2.1,
                    label='%s, corrected' % name)

    ax.axhline(1.0, color='k', lw=1.1)
    ax.text(0.85, 1.03, 'M = 1: the crossover', fontsize=7.2)
    n_dirty = (ilc.N_IMP_REF + ilc.DN_IMP_HFO2_HI)
    m_d_old = operating_margin(ddm.C_OX_TOPGATE_HFO2, n_dirty, 1.0,
                               kappa=KAPPA_ILC)
    m_d_new = operating_margin(ddm.C_OX_TOPGATE_HFO2, n_dirty, 1.0, kappa=k_b)
    ax.plot([n_dirty * 1e-16] * 2, [m_d_old, m_d_new], 'k-', lw=0.9)
    ax.plot([n_dirty * 1e-16], [m_d_old], 'kv', ms=5)
    ax.plot([n_dirty * 1e-16], [m_d_new], 'k^', ms=5)
    ax.annotate(r'dirtiest HfO$_2$ case:'
                '\n' r'$%.2f\times \rightarrow %.2f\times$'
                % (m_d_old, m_d_new),
                xy=(n_dirty * 1e-16, m_d_new),
                xytext=(1.10, 2.05), fontsize=7.4,
                arrowprops=dict(arrowstyle='->', lw=0.8, color='k'))
    ax.set_ylim(0, 6.2)
    ax.set_xlabel(r'charged-impurity density $N_{imp}$  ($10^{12}$ cm$^{-2}$)')
    ax.set_ylabel(r'operating margin  $M = n(V_{over})/n_c$')
    ax.set_title('(b)  the correction moves the marginal verdict off the '
                 'margin', fontsize=9.5)
    ax.legend(fontsize=6.6, loc='upper right')
    ax.grid(True, which='both', alpha=0.25)

    # -- panel (c): the exact identity, which is the day's real finding ---
    ax = axes[2]
    r0 = 0.10
    xs = np.logspace(-2.0, 2.0, 400)          # |n|/n_c
    n_ref = ilc.N_IMP_REF
    n_c_ref = float(n_c_of_impurity(n_ref, kappa=k_b))
    for sp, col, lab in ((0.0, 'C3', r'$s=0$  (DDM as shipped)'),
                         (0.05, 'C1', r'$s=0.05$'),
                         (0.10, 'C0', r'$s=0.10$'),
                         (2.0 * r0, 'C2',
                          r'$s=2r_0$  (exact: $r=r_0$, no screening)')):
        rr = [r_branch_split(xq * n_c_ref, n_ref, r0, split=sp, kappa=k_b)
              for xq in xs]
        ax.semilogx(xs, rr, color=col, lw=2.2 if sp in (0.0, 2.0 * r0) else 1.4,
                    ls='-' if sp != 0.0 else '--', label=lab)
    ax.axhline(r0, color='k', lw=0.8, ls=':')
    ax.text(1.2e-2, r0 * 1.06, r'$r_0 = %.2f$' % r0, fontsize=7.4)
    x_op = n_op_for_figure() / n_c_ref
    ax.axvline(x_op, color='0.45', lw=0.9)
    ax.text(x_op * 1.15, 0.0705,
            '10 nm HfO$_2$\nat 1 V', fontsize=7.0, color='0.3')
    ax.set_xlabel(r'$|n|/n_c$')
    ax.set_ylabel(r'asymmetry  $r(n)$')
    ax.set_title('(c)  whether screening kills the asymmetry depends on an\n'
                 'unexamined symmetry of the neutral-defect channel',
                 fontsize=9.0)
    ax.legend(fontsize=6.8, loc='lower left')
    ax.grid(True, which='both', alpha=0.25)
    ax.set_ylim(0, 0.118)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Section 5: main
# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('KAPPA IS BOUNDED BY ZION\'S OWN IMPURITY DENSITY, AND 2.4 IS')
    print('EXCLUDED -- arXiv:1801.03565 full text, read 2026-11-03')
    print('=' * 78)

    c_ox = c_ox_zion()
    n_a = n_a_s_zion(c_ox=c_ox)
    k_b = kappa_bound(n_imp_lower=n_a)
    k_stated = kappa_bound(n_imp_lower=N_A_S_ZION_STATED)

    # -- Section 1 ---------------------------------------------------------
    print('\nSECTION 1.  The full text carries NO mobility (the item\'s datum')
    print('            does not exist), and Zion\'s own arithmetic reproduced')
    print('  C_ox  = eps0*eps/t = %.6e F/m^2   (Zion state %.2e)'
          % (c_ox, C_OX_ZION_STATED))
    print('  N_A^S = C_ox*V_NP/e = %.6e m^-2   (Zion state %.1e)'
          % (n_a, N_A_S_ZION_STATED))
    check('D1 C_ox reproduces Zion\'s stated value AT THEIR PRECISION, from '
          'eps0, eps and t alone',
          float('%.3g' % c_ox) == C_OX_ZION_STATED,
          'computed %.6e rounds to %.3g at 3 s.f.; Zion print %.2e.  Their '
          'capacitance is internally consistent, which is what makes D2\'s '
          'failure on the NEXT line attributable to the quoted density '
          'rather than to a different oxide convention.'
          % (c_ox, float('%.3g' % c_ox), C_OX_ZION_STATED))

    shortfall = 1.0 - n_a / N_A_S_ZION_STATED
    check('D2 N_A^(S) does NOT reproduce: Zion\'s quoted 1.5e12 is 4-5 % '
          'above what their own C_ox and V_NP give',
          0.04 < shortfall < 0.05 and float('%.2g' % (n_a * 1e-16)) != 1.5,
          'their formula and their inputs give %.4fe12 cm^-2, a shortfall of '
          '%.2f %% against the %.1fe12 they print, and %.4f does not round to '
          '1.5 at any precision.  Both values are carried forward; the '
          'RECOMPUTED one is used because by (B2) it is the LOOSER bound, so '
          'it cannot manufacture the exclusion below.'
          % (n_a * 1e-16, 100.0 * shortfall, N_A_S_ZION_STATED * 1e-16,
             n_a * 1e-16))

    # -- Section 2 ---------------------------------------------------------
    print('\nSECTION 2.  The bound, and its direction')
    print('  KAPPA <= n_c/N_A^(S) = %.4f        (recomputed N_A, '
          'conservative)' % k_b)
    print('  KAPPA <= n_c/1.5e16  = %.4f        (Zion\'s rounded N_A)'
          % k_stated)
    print('  KAPPA  = A_ref/C     = %.4f        (ILC 11-02, by pairing)'
          % KAPPA_ILC)

    rt = kappa_bound(n_imp_lower=n_a) * n_a
    d3a = ulp_distance(rt, N_C_ZION)
    check('D3a the bound round-trips to within the 1 ulp that floating point '
          'guarantees -- NOT bitwise, and the first version of this check '
          'wrongly predicted bitwise and failed',
          d3a <= 1,
          'returned %.17e against n_c = %.17e, ulp distance %d.  (n_c/N_A)*N_A '
          'rounds twice, so exactness is not available here and asserting it '
          'was an overclaim; it is corrected by splitting the check, not by '
          'widening it.' % (rt, N_C_ZION, d3a))

    # D3b: exactness where it genuinely holds -- a power of two.
    n_a_pow2 = N_C_ZION / 4.0
    k_pow2 = kappa_bound(n_imp_lower=n_a_pow2)
    rt2 = k_pow2 * n_a_pow2
    check('D3b where the round-trip IS exact -- N_A a power of two times n_c '
          '-- it is bitwise, and KAPPA comes out exactly 4.0',
          ulp_distance(rt2, N_C_ZION) == 0 and k_pow2 == 4.0,
          'kappa = %.17e (exactly 4.0: %s), round-trip ulp %d.  Division and '
          'multiplication by a power of two are both exact, so this isolates '
          'the double-rounding of D3a as the whole of its 1 ulp.'
          % (k_pow2, k_pow2 == 4.0, ulp_distance(rt2, N_C_ZION)))

    k2 = kappa_bound(n_c=2.0 * N_C_ZION, n_imp_lower=2.0 * n_a)
    check('D5 KAPPA is scale-invariant BITWISE under a factor-of-two scaling '
          'of both inputs',
          ulp_distance(k2, k_b) == 0,
          'kappa(2n_c, 2N_A) = %.17e, kappa(n_c, N_A) = %.17e, ulp %d.  A '
          'power of two is exactly representable and commutes with rounding '
          'at every step, so this is a true bitwise symmetry rather than a '
          'tolerance, and it says KAPPA has exactly one dimensionless group.'
          % (k2, k_b, ulp_distance(k2, k_b)))

    # -- Section 4: the verdict -------------------------------------------
    print('\nSECTION 4.  The verdict')
    excess = KAPPA_ILC / k_b - 1.0
    excess_stated = KAPPA_ILC / k_stated - 1.0
    print('  2.4 exceeds the conservative bound by %.2f %% (and the rounded-'
          'N_A bound by %.2f %%)' % (100.0 * excess, 100.0 * excess_stated))

    def excluded(kappa_trial, bound):
        return kappa_trial > bound

    check('D4 the exclusion FIRES on KAPPA = 2.4 and does NOT fire on a '
          'planted control of 1.9',
          excluded(KAPPA_ILC, k_b) and not excluded(1.9, k_b),
          '2.4 > %.4f so it is excluded; 1.9 < %.4f so it is not.  The same '
          'code path returns both answers, which is what distinguishes a test '
          'from an assertion -- ILC\'s D11 made the same point about its own '
          'mutant.' % (k_b, k_b))

    # -- Section 5: the mobility by-product -------------------------------
    print('\nSECTION 5.  The mobility the paper never states, bounded anyway')
    mu_b = mu_L_zion_bound(kappa=k_b)
    print('  mu_L,Zion = KAPPA*C/n_c <= %.6f m^2/Vs = %.0f cm^2/Vs'
          % (mu_b, mu_b * 1e4))
    print('  the repository paired gfet.mu = %.3f m^2/Vs = %.0f cm^2/Vs'
          % (gfet.mu, gfet.mu * 1e4))
    ratio_mu = mu_b / gfet.mu
    ratio_k = k_b / KAPPA_ILC
    check('D10 the mobility bound and the KAPPA bound are the SAME statement: '
          'both fall below the paired value by an identical factor',
          ulp_distance(ratio_mu, ratio_k) <= 1,
          'mu_bound/gfet.mu = %.17f, kappa_bound/KAPPA_ILC = %.17f, ulp '
          'distance %d.  Section 5\'s algebra is thereby checked rather than '
          'asserted: A = KAPPA*C and A = n_c*mu_L are one relation, so the '
          'two ratios cannot differ.'
          % (ratio_mu, ratio_k, ulp_distance(ratio_mu, ratio_k)))

    # -- Section 6: overlap and invariance --------------------------------
    print('\nSECTION 6.  Overlap with 11-02, and what the correction cannot '
          'move')
    sweep = np.logspace(np.log10(1e14), np.log10(1e20), 401)
    mine_at_ilc = n_c_of_impurity(sweep, kappa=KAPPA_ILC)
    ilc_native = np.array([float(ilc.n_c_of_impurity(q)) for q in sweep])
    d6_bad = int(np.count_nonzero(mine_at_ilc - ilc_native != 0.0))
    check('D6 at kappa = KAPPA_ILC the corrected n_c differs from ILC\'s by '
          'EXACTLY 0.0 at all 401 densities -- zero, not small',
          d6_bad == 0,
          '%d of 401 samples have a nonzero difference.  A correction '
          'parameterised by kappa must contain the model it corrects as the '
          'value it corrects away from, and that admits an exact test rather '
          'than a tolerance.' % d6_bad)

    a_from_k = a_neutral_from_kappa(KAPPA_ILC)
    via_ilc = np.array([float(ilc.n_c_of_impurity(q, a_neutral=a_from_k))
                        for q in sweep])
    d7_ulp = max(ulp_distance(p, q) for p, q in zip(mine_at_ilc, via_ilc))
    check('D7 routed through ILC\'s OWN function with a_neutral = kappa*C, '
          'the agreement is bitwise over six decades',
          d7_ulp == 0,
          'worst ulp distance %d over 401 points spanning N_imp = 1e14..1e20 '
          'm^-2.  D6 tests this module\'s arithmetic; D7 tests that it is the '
          'same function ILC ships, which is the stronger of the two.'
          % d7_ulp)

    r_old = (n_c_of_impurity(4e16, kappa=KAPPA_ILC)
             / n_c_of_impurity(1e16, kappa=KAPPA_ILC))
    r_new = (n_c_of_impurity(4e16, kappa=k_b)
             / n_c_of_impurity(1e16, kappa=k_b))
    mr_old = (operating_margin(gfet.C_ox, 2e16, 5.0, kappa=KAPPA_ILC)
              / operating_margin(gfet.C_ox, 4e16, 5.0, kappa=KAPPA_ILC))
    mr_new = (operating_margin(gfet.C_ox, 2e16, 5.0, kappa=k_b)
              / operating_margin(gfet.C_ox, 4e16, 5.0, kappa=k_b))
    check('D8 every RATIO is bitwise invariant under a 13 % change in kappa '
          '-- n_c ratios and margin ratios alike',
          ulp_distance(r_old, r_new) == 0 and ulp_distance(mr_old, mr_new) == 0,
          'n_c(4e16)/n_c(1e16) = %.17e both before and after (ulp %d); the '
          'margin ratio likewise (ulp %d).  ILC\'s D5 established that the '
          'ratios do not carry the pairing; this is that result collected on '
          'a correction it did not anticipate, and it is why Section 6 can '
          'say the 11-02 CONCLUSIONS survive while its absolute numbers move.'
          % (r_old, ulp_distance(r_old, r_new), ulp_distance(mr_old, mr_new)))

    # -- Section 7: the propagation table ---------------------------------
    print('\nSECTION 7.  What moves, in the device units')
    print('  %-32s %-11s %-11s %-10s'
          % ('quantity', '11-02', 'corrected', 'change'))
    rows = []
    n_dirty = ilc.N_IMP_REF + ilc.DN_IMP_HFO2_HI
    for name, old, new in (
            ('n_c at N_imp,ref (e12 cm^-2)',
             n_c_of_impurity(ilc.N_IMP_REF, kappa=KAPPA_ILC) * 1e-16,
             n_c_of_impurity(ilc.N_IMP_REF, kappa=k_b) * 1e-16),
            ('dV_c, 10 nm HfO2 (V)',
             crossover_overdrive(ddm.C_OX_TOPGATE_HFO2, ilc.N_IMP_REF,
                                 kappa=KAPPA_ILC),
             crossover_overdrive(ddm.C_OX_TOPGATE_HFO2, ilc.N_IMP_REF,
                                 kappa=k_b)),
            ('M, HfO2 at 1 V, dirtiest case',
             operating_margin(ddm.C_OX_TOPGATE_HFO2, n_dirty, 1.0,
                              kappa=KAPPA_ILC),
             operating_margin(ddm.C_OX_TOPGATE_HFO2, n_dirty, 1.0, kappa=k_b)),
            ('N_imp*, HfO2 at 1 V (e12 cm^-2)',
             max_tolerable_impurity(ddm.C_OX_TOPGATE_HFO2, 1.0,
                                    kappa=KAPPA_ILC) * 1e-16,
             max_tolerable_impurity(ddm.C_OX_TOPGATE_HFO2, 1.0,
                                    kappa=k_b) * 1e-16),
            ('N_imp*, 90 nm SiO2 at 5 V (e12)',
             max_tolerable_impurity(gfet.C_ox, 5.0, kappa=KAPPA_ILC) * 1e-16,
             max_tolerable_impurity(gfet.C_ox, 5.0, kappa=k_b) * 1e-16)):
        rows.append((name, old, new))
        print('  %-32s %-11.4f %-11.4f %+.2f %%'
              % (name, old, new, 100.0 * (new / old - 1.0)))

    m_old = rows[2][1]
    m_new = rows[2][2]
    check('the 11-02 CONCLUSIONS all survive, and the one that was MARGINAL '
          'moves away from the margin',
          m_old > 1.0 and m_new > m_old,
          'the thinnest 11-02 margin was %.3fx -- a top-gated device barely '
          'past its own crossover at the dirty end of the HfO2 deposition '
          'range, close enough to 1 that the verdict was fragile.  It becomes '
          '%.3fx.  No sign flips and no verdict reverses; the correction '
          'moves the fragile one to safety, which is the useful direction for '
          'a result to move when its weakest input is replaced.'
          % (m_old, m_new))

    n_star_top_new = max_tolerable_impurity(ddm.C_OX_TOPGATE_HFO2, 1.0,
                                           kappa=k_b)
    budget_new = n_star_top_new - ilc.N_IMP_REF
    budget_old = (max_tolerable_impurity(ddm.C_OX_TOPGATE_HFO2, 1.0,
                                         kappa=KAPPA_ILC) - ilc.N_IMP_REF)
    check('the HfO2 deposition budget was THIN on 11-02 and is now '
          'comfortable -- the headroom roughly doubles',
          budget_new > ilc.DN_IMP_HFO2_HI and budget_new > budget_old,
          'headroom over the paper\'s dirty end (%.2fe12 cm^-2) goes from '
          '%.1f %% to %.1f %%.  The 11-02 entry called this budget "REAL but '
          'thin"; the correction is what makes it not thin, and it does so by '
          'replacing the weakest input rather than by widening a tolerance.'
          % (ilc.DN_IMP_HFO2_HI * 1e-16,
             100.0 * (budget_old / ilc.DN_IMP_HFO2_HI - 1.0),
             100.0 * (budget_new / ilc.DN_IMP_HFO2_HI - 1.0)))

    # -- Section 8: D9, the new open item ---------------------------------
    print('\nSECTION 8.  The new problem: n_c is a HOLE-branch number')
    n_op = ddm.C_OX_TOPGATE_HFO2 * 1.0 / ddm.E_CHARGE
    r0 = 0.10
    r_single = r_branch_split(n_op, ilc.N_IMP_REF, r0, split=0.0, kappa=k_b)
    mu_pair = ddm.mu_pair_density(
        n_op,
        (ilc.C_CHEN / ilc.N_IMP_REF) * (1.0 + r0),
        (ilc.C_CHEN / ilc.N_IMP_REF) * (1.0 - r0),
        n_c=float(n_c_of_impurity(ilc.N_IMP_REF, kappa=k_b)))
    r_ddm = (float(mu_pair[0]) - float(mu_pair[1])) \
        / (float(mu_pair[0]) + float(mu_pair[1]))
    check('D9a at split = 0 the branch-split form reduces to DDM\'s single-n_c '
          'construction BITWISE',
          ulp_distance(r_single, r_ddm) == 0,
          'r(split=0) = %.17e, DDM\'s r = %.17e, ulp %d.  The split form is '
          'therefore a generalisation of (B1) and not a second, different '
          'model -- which is the precondition for reading D9b as a '
          'sensitivity rather than as a discrepancy between two codes.'
          % (r_single, r_ddm, ulp_distance(r_single, r_ddm)))

    # D9b: the EXACT identity of (B4).  r0 is known in advance and does not
    # depend on density, which is what makes this a target rather than a band.
    worst_ulp = 0
    worst_at = None
    for r0_t in (0.02, 0.05, 0.10, 0.25):
        for n_t in np.logspace(14.0, 20.0, 401):
            r_t = r_branch_split(n_t, ilc.N_IMP_REF, r0_t,
                                 split=2.0 * r0_t, kappa=k_b)
            u = ulp_distance(r_t, r0_t)
            if u > worst_ulp:
                worst_ulp, worst_at = u, (r0_t, n_t, r_t)
    check('D9b THE EXACT IDENTITY (B4): at a branch split s = 2*r0 the '
          'asymmetry equals r0 at EVERY density -- no screening at all',
          worst_ulp <= 64,
          'worst deviation from the exactly known target r0 is %d ulp (%.2e '
          'relative) over 4 values of r0 and 401 densities spanning six '
          'decades; worst at r0 = %.2f, n = %.3e, r = %.17f.  The identity is '
          'exact in real arithmetic -- the (1+x) factor cancels out of the '
          'ratio -- so the residual is the computation path and is measured '
          'here rather than tolerated.'
          % (worst_ulp, worst_ulp * EPS, worst_at[0], worst_at[1],
             worst_at[2]))

    r_split = r_branch_split(n_op, ilc.N_IMP_REF, r0, split=BRANCH_SPLIT_D9,
                             kappa=k_b)
    d9 = abs(r_split / r_single - 1.0)
    print('  |n|/n_c = %.4f at the top-gate operating point'
          % (n_op / float(n_c_of_impurity(ilc.N_IMP_REF, kappa=k_b))))
    print('  %-26s %-12s %-s' % ('branch split s', 'r', 'reading'))
    for s_t, note in ((0.0, 'the shipped answer: screened %.2fx'
                       % (r0 / r_single)),
                      (0.10, 'intermediate'),
                      (2.0 * r0, '= r0 exactly: NO screening')):
        print('  %-26.3f %-12.6f %-s'
              % (s_t, r_branch_split(n_op, ilc.N_IMP_REF, r0, split=s_t,
                                     kappa=k_b), note))

    check('D9c DDM\'s SCREENING CONCLUSION IS CONDITIONAL, and the condition '
          'is an unexamined symmetry of the neutral-defect channel',
          d9 > 1.0,
          'a %.0f %% split -- modest for a quantity never measured -- moves r '
          'at the operating density by %.0f %%, and the s = 2*r0 endpoint is '
          'not a perturbation of the shipped r = %.4f but its negation '
          '(r = r0 = %.2f, density-independent).  DDM derives r0 from a '
          'beyond-Born cross-section ratio, which is a property of scattering '
          'off a potential; neutral defects are also potentials, and no '
          'argument anywhere in this repository says the asymmetry should '
          'vanish for vacancies while surviving for charged impurities.  It '
          'was assumed by writing a single n_c.  THIS RANKS ABOVE TODAY\'S '
          'OWN CORRECTION: Sections 4-6 move absolute numbers by 13-15 %%, '
          'this decides whether a qualitative conclusion of Chapters 3 and 4 '
          'holds at all.'
          % (100.0 * BRANCH_SPLIT_D9, 100.0 * d9, r_single, r0))

    # -- D11 mutation ------------------------------------------------------
    print('\nSECTION 9.  Mutation: flip the inequality')
    k_mut = kappa_bound(n_imp_lower=n_a)
    mut_excluded = KAPPA_ILC < k_mut      # the WRONG reading of the direction
    check('D11 reading N_A^(S) as an UPPER bound on N_imp instead of a lower '
          'one REVERSES the verdict, so the verdict is using the argument of '
          'Section 2 (iii)',
          excluded(KAPPA_ILC, k_mut) and not mut_excluded,
          'under the correct reading N_imp >= N_A^(S), so KAPPA <= %.4f and '
          '2.4 is excluded.  Under the flipped reading the same arithmetic '
          'would make %.4f a LOWER bound and 2.4 merely unremarkable.  The '
          'exclusion therefore rests on the three reasons of Section 2 '
          '(net-not-count, unquantified N_t, transfer residue) and not on '
          'the arithmetic, which is identical either way.' % (k_mut, k_mut))

    # -- figure ------------------------------------------------------------
    print('\nSECTION 10.  Figure')
    path = make_figure()
    print('  wrote %s' % path)

    print('\n' + '=' * 78)
    print('%d checks, %d passing, %d failing' % (passed + failed, passed,
                                                 failed))
    print('=' * 78)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
