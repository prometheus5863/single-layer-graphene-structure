"""
graphene_areal_thermal_resistance.py

THE RANKED-FIRST ITEM CLOSES AND ITS STATED JUSTIFICATION IS FALSE.  At
Dorgan's fixed geometry the per-unit-area check is NOT an independent
constraint on R_th and A: it is check A2 with the SAME CENTRE and a 2.8x
TIGHTER tolerance, and the two residuals are bitwise identical.  The
independence the item claimed exists, but it lives in the AREA SCALING, not
in the single-geometry comparison -- and measured there it says the paper's
areal figure of merit is transferable only within a factor 1.8 in area.

This closes "CHECK R_th PER UNIT AREA AGAINST THE PAPER'S ~2.8e-7 m^2.K/W",
created and ranked FIRST on 2026-11-11, and in passing it DECIDES the item
ranked third the same day ("REPAIR A2's OVER-TIGHT TOLERANCE, OR DECIDE NOT
TO") -- see Section 7.  It also corrects a significant-figure statement made
twice yesterday, in the note and in the log, whose conclusion was right and
whose stated reason was wrong (Section 6).

------------------------------------------------------------------------------
1. THE TWO NUMBERS, AND THAT THEY ARE ONE STATEMENT
------------------------------------------------------------------------------
Dorgan, Bae & Pop, APL 97, 082112 (2010), main text, ONE sentence:

    "At 300 K for our geometry R_th ~ 10^4 K/W, or ~2.8x10^-7 m^2 K/W per
     unit of device area"

with W = 7 um, L = 4 um, t_ox = 300 nm.  The paper therefore prints the SAME
physical fact twice, at two precisions:

    1e4     K/W        -- one significant figure   -> licenses +-50 %
    2.8e-7  m^2.K/W    -- two significant figures  -> licenses +-1.79 %

A = L*W = 2.8e-11 m^2 is EXACT, given by the stated geometry and not fitted.
So the second statement is the first divided by an exact constant, and

    2.8e-7 / 2.8e-11 == 1.0e4     exactly, in binary floating point

which Section 3 asserts rather than describes.  This is why the per-area
comparison cannot be independent of A2 at this geometry: the two tests have
the same centre, and only their licensed widths differ.

------------------------------------------------------------------------------
2. THE TRAP, WHICH IS THE WHOLE REASON THE ITEM WAS WORTH PREDICTING
------------------------------------------------------------------------------
"Per unit area" reads as a division.  It is a MULTIPLICATION:

    R_th scales as 1/A  (R_B = 1/(hA), R_ox = t_ox/(k_ox A))

so the area-normalised, geometry-transferable figure of merit is R_th * A,
with units K/W * m^2 = m^2.K/W -- which is the unit the paper prints.  The
quotient R_th / A has units K/(W.m^2), which is NOT what the paper prints,
and numerically

    R_th / A = 3.54e14          -- twenty-one decades from 2.8e-7
    R_th * A = 2.776e-7         -- rounds to 2.8e-7 at two figures

The 2026-11-11 log predicted a session would compute the quotient first and
have a number "seven decades off".  It is TWENTY-ONE decades off; the
prediction was right about the trap and wrong about its size, because it
compared the quotient to R_th's own 1e4 rather than to the areal 2.8e-7.
Checks P4 and P5 plant the quotient form as a control: a check suite that
could not tell the two apart would pass both.

------------------------------------------------------------------------------
3. WHY THE INDEPENDENCE CLAIM IS FALSE (THE EXACT ZERO)
------------------------------------------------------------------------------
The 11-11 entry's reason for ranking this first was:

    "It is an INDEPENDENT constraint: A2 tests R_th alone and therefore tests
     the area and the conductivity only through their product, while this
     tests the ratio.  Of the two it is the one that can fail for a reason A2
     cannot see."

Every clause of that is wrong, and the module asserts so:

  (i)  A2 does NOT test R_th "alone" in any sense that excludes A.  R_th is
       COMPUTED from A by Eq. (5); there is no path by which A2 sees R_th
       without seeing A.
  (ii) The areal quantity is not a "ratio" of the two, it is their product.
  (iii) Because 2.8e-7/2.8e-11 is 1.0e4, the relative residual of R_th*A
       against 2.8e-7 IS the relative residual of R_th against 1.0e4.
       Nothing can fail one and pass the other, so neither "can fail for a
       reason the other cannot see".

ONE CORRECTION MADE WITHIN THIS SESSION, RECORDED RATHER THAN TIDIED AWAY.
I1 and I2 were first written as EXACT equalities -- `2.8e-7/A == 1.0e4` and
`residual difference == 0.0` -- on the reasoning above, and BOTH FAILED on
the first run.  4e-6*7e-6 is 2.7999999999999997e-11, one ulp below the
literal 2.8e-11, so the exact statement is about real numbers and the machine
is working in binary.  The surviving claim is sharper than a tolerance and
weaker than an identity: the two residuals differ by 1.9e-16, one ulp of 1.0,
while two genuinely independent constraints would differ at the 1e-2 level
the residuals themselves occupy -- fourteen decades larger.  P1 now carries
the representability fact that makes this so.

The standing rule since 2026-09-17, when a plausible-range check passed a
real bug, is that every new instrument carries at least one comparison
against an exactly known value.  I3 is that comparison here and it IS an
exact zero; I2 is not, and saying so is the point of this paragraph.

------------------------------------------------------------------------------
4. THE ALGEBRAIC IDENTITY (THE SECOND EXACT ZERO)
------------------------------------------------------------------------------
Multiplying Eq. (5) through by A:

    R_th * A  ==  1/h  +  t_ox/k_ox  +  sqrt(A)/(2 k_Si)
                  \_______________/     \_____________/
                   A-INDEPENDENT          grows as sqrt(A)

This is the content of the whole exercise.  Two of the three terms are
area-independent -- so for them the paper's areal number IS transferable --
and the silicon spreading term is not, because R_Si goes as 1/sqrt(A) rather
than 1/A.  Check I3 requires the identity to hold to machine precision.

------------------------------------------------------------------------------
5. WHERE THE INDEPENDENCE ACTUALLY LIVES, AND WHAT IT MEASURES
------------------------------------------------------------------------------
Section 3 kills the claim at ONE geometry.  Across geometries the areal
quantity carries real information that A2 cannot: A2 at a different area
tests a different number, while R_th*A should be nearly invariant.  Nearly.
areal_decomposition() and areal_band() measure how nearly, and the answer is
a bound on how far the paper's number may be carried:

  - the A-independent floor is 1/h + t_ox/k_ox, reached as A -> 0 (check S1,
    an exact limit);
  - the spreading term contributes ~11.9 % of Dorgan's own 2.776e-7;
  - and the areal figure of merit leaves the paper's own two-figure band
    after only a factor ~1.8 in area.

That last number is the module's new physical result.  "Per unit of device
area" invites reuse of 2.8e-7 as a scale-free areal figure of merit; it is
scale-free only within a factor of about 1.8 in A, and this repository's own
FET (W = 1 um, L = 200 nm, A = 2e-13 m^2) is 140x smaller than Dorgan's --
far outside that band.  Section 5's check S4 quantifies the error that
reusing 2.8e-7 there would make, and it is NOT small.

------------------------------------------------------------------------------
6. A SIGNIFICANT-FIGURE STATEMENT CORRECTED, NOT DELETED
------------------------------------------------------------------------------
Both yesterday's note and yesterday's log say of 2.8e-7:

    "One significant figure licenses +-1.8 %."

2.8e-7 has TWO significant figures.  Two figures license +-1.79 %, which is
the quoted tolerance; ONE figure would license +-16.7 %.  So the number
carried forward was correct and the reason given for it was not -- the
tolerance is right by two-figure arithmetic while being attributed to
one-figure arithmetic that would give a width 9.3x wider.  Checks F1-F3 pin
all three numbers so this cannot drift back.

This matters beyond bookkeeping: the same sentence is what makes A2's +-5 %
look over-tight (Section 7), and under the stated one-figure reading it would
be the per-area check that was over-tight by 9x.

------------------------------------------------------------------------------
7. THE RANKED-THIRD ITEM IS DECIDED: A2 IS NOT OVER-TIGHT; ITS WARRANT WAS
------------------------------------------------------------------------------
The 11-11 list ranked third "REPAIR A2's OVER-TIGHT TOLERANCE, OR DECIDE NOT
TO", calling +-5 % against a one-figure ~1e4 "a tighter claim than the source
supports", with widening to +-50 % the alternative and "nearly vacuous" the
objection.  Sections 1 and 3 dissolve the dilemma:

    the paper states the SAME fact to TWO figures elsewhere in the SAME
    sentence, and that statement licenses +-1.79 % on the identical
    constraint.

So +-5 % is not a tighter claim than the source supports.  It is LOOSER than
the source supports, by 2.8x.  A2's defect was never its width but its
CITATION: it quoted the one-figure half of a sentence whose other half is
two-figure, and then used a width that neither half licenses.  The repair is
therefore to re-warrant A2 and TIGHTEN it to the licensed +-1.79 %, which it
passes with 0.84 % of residual -- and NOT to widen it to +-50 %.  That is
done in graphene_self_heating_model.py in the same commit range as this
module, with A2's original wording and number preserved in place.

This is the third time a census finding in this tree has turned out to be
about a check's warrant rather than its arithmetic (A1 on 11-11 was the
second).  The pattern is now explicit on the standing list.

------------------------------------------------------------------------------
WHAT THIS MODULE IS NOT
------------------------------------------------------------------------------
  - It does not re-derive Eq. (5) or revisit the L/t_ox = 2.22 validity
    violation; R_ox is still over-estimated and every R_th here is still the
    upper bound graphene_self_heating_model.py documents it to be.  The areal
    identity and the band are properties OF that bound.
  - The factor-1.8 band is computed at T = T0 = 300 K with t_ox held at
    Dorgan's 300 nm.  A device of different oxide thickness moves the
    A-independent floor, which moves the band; that sweep is a standing item,
    not a result here.
  - sqrt(A) in R_Si presumes the paper's square-ish spreading geometry.  At
    extreme aspect ratio neither R_Si nor the band is reliable, and nothing
    here tests that.
"""

import sys

import numpy as np

import graphene_fet_model as gfet
import graphene_self_heating_model as sh

# Dorgan's own geometry, from the main-text sentence of Section 1.
W_DORGAN = 7e-6
L_DORGAN = 4e-6
T_OX_DORGAN = 300e-9

# The two published numbers and the figures they are printed to.
R_TH_PAPER = 1.0e4           # K/W,       1 significant figure
R_TH_PAPER_SIGFIGS = 1
AREAL_PAPER = 2.8e-7         # m^2.K/W,   2 significant figures
AREAL_PAPER_SIGFIGS = 2

_CHECKS = []


def check(name, ok, detail):
    _CHECKS.append((name, bool(ok), detail))
    return bool(ok)


def relative_halfwidth(quoted, sigfigs):
    """+-half-width, relative, licensed by `quoted` printed to `sigfigs`."""
    mantissa = abs(quoted) / 10.0 ** np.floor(np.log10(abs(quoted)))
    return 0.5 * 10.0 ** (1 - sigfigs) / mantissa


def areal_thermal_resistance(W=W_DORGAN, L=L_DORGAN, t_ox=T_OX_DORGAN,
                             T=sh.T_AMBIENT, T0=sh.T_AMBIENT,
                             h=sh.H_BOUNDARY):
    """R_th * A in m^2.K/W -- the paper's 'per unit of device area' quantity.

    A PRODUCT, not a quotient.  See Section 2; P4/P5 plant the quotient.
    """
    t = sh.thermal_resistances(W=W, L=L, t_ox=t_ox, T=T, T0=T0, h=h)
    return t["R_th"] * t["A"]


def areal_decomposition(W=W_DORGAN, L=L_DORGAN, t_ox=T_OX_DORGAN,
                        T=sh.T_AMBIENT, T0=sh.T_AMBIENT, h=sh.H_BOUNDARY):
    """The three terms of Section 4's identity, in m^2.K/W, as a dict.

    `floor` is the A-independent sum and `spreading` the sqrt(A) term, so
    floor + spreading == R_th*A identically (check I3).
    """
    A = L * W
    k_ox = sh.kappa_ox(0.5 * (T0 + T))
    boundary = 1.0 / h
    oxide = t_ox / k_ox
    spreading = np.sqrt(A) / (2.0 * sh.kappa_Si(T0))
    return {"boundary": boundary, "oxide": oxide, "spreading": spreading,
            "floor": boundary + oxide, "total": boundary + oxide + spreading,
            "A": A, "k_ox": k_ox}


def areal_band(t_ox=T_OX_DORGAN, T=sh.T_AMBIENT, T0=sh.T_AMBIENT,
               h=sh.H_BOUNDARY, quoted=AREAL_PAPER,
               sigfigs=AREAL_PAPER_SIGFIGS):
    """The range of AREAS over which R_th*A stays inside the paper's own
    two-figure band around `quoted`.  Returns a dict in m^2.

    Inverts Section 4's identity for sqrt(A).  `None` for a bound the floor
    already violates: if the A-independent floor lies above the band's top no
    area satisfies it, and if the floor lies above the band's bottom then
    every area down to zero satisfies the lower edge.
    """
    hw = relative_halfwidth(quoted, sigfigs)
    lo_val, hi_val = quoted * (1.0 - hw), quoted * (1.0 + hw)
    d = areal_decomposition(W=1.0, L=1.0, t_ox=t_ox, T=T, T0=T0, h=h)
    floor = d["floor"]
    two_kSi = 2.0 * sh.kappa_Si(T0)

    def area_for(val):
        s = (val - floor) * two_kSi
        return None if s <= 0.0 else s * s

    return {"halfwidth": hw, "value_lo": lo_val, "value_hi": hi_val,
            "floor": floor, "A_lo": area_for(lo_val), "A_hi": area_for(hi_val)}


# ---------------------------------------------------------------------------
# CHECKS.  Published two-figure comparisons, exact algebraic zeros, exact
# limits, and planted controls.  No plausible ranges anywhere.
# ---------------------------------------------------------------------------
def verify_published_comparison():
    """P: the comparison the item asked for, and the quotient form planted."""
    t = sh.thermal_resistances(W=W_DORGAN, L=L_DORGAN, t_ox=T_OX_DORGAN,
                               T=sh.T_AMBIENT, T0=sh.T_AMBIENT)
    A, R = t["A"], t["R_th"]
    # DEFECT FOUND BY MUTANT M1 BEFORE COMMIT, AND RECORDED RATHER THAN
    # QUIETLY FIXED.  This line first read `product = R * A`, computing the
    # headline quantity inline from sh.thermal_resistances() and never
    # calling this module's own areal_thermal_resistance().  M1 replaces that
    # function with the quotient form and was reported DETECTED-ELSEWHERE:
    # P2 and P3 could not see the injected trap at all, because the function
    # every other module would import was not on the path the headline check
    # tested.  P2 and P3 now go through the public function, which is what
    # makes M1 a test of them.  This is the same defect class as the 10-xx
    # "installed surface" items -- a check that validates an expression
    # rather than the thing that ships.
    product, quotient = areal_thermal_resistance(
        W=W_DORGAN, L=L_DORGAN, t_ox=T_OX_DORGAN), R / A
    hw = relative_halfwidth(AREAL_PAPER, AREAL_PAPER_SIGFIGS)

    # P1 WAS WRITTEN AS `A == 4e-6*7e-6`, WHICH IS A TAUTOLOGY -- A is
    # computed as L*W, so that check could not fail and tested nothing.  It
    # is replaced by the fact it was reaching for, which is NOT a tautology
    # and which drives I1 and I2 below: L*W is not the literal 2.8e-11.
    check("P1 A = L*W is NOT the literal 2.8e-11 the paper's geometry names "
          "-- 4e-6*7e-6 lands one ulp below it -- so every 'exact' identity "
          "between the paper's two numbers holds only to the "
          "representability of 2.8e-11 (this is what I1/I2 then measure)",
          A != 2.8e-11 and abs(A - 2.8e-11) <= np.spacing(2.8e-11),
          "L*W = %.17e, literal 2.8e-11 = %.17e, difference %.3e = %.1f ulp"
          % (A, 2.8e-11, A - 2.8e-11,
             abs(A - 2.8e-11) / np.spacing(2.8e-11)))
    check("P2 R_th*A rounds to the paper's 2.8e-7 m^2.K/W at the 2 s.f. it "
          "is printed to -- the comparison the 11-11 item ranked first",
          sh._rounds_to(product, AREAL_PAPER, AREAL_PAPER_SIGFIGS),
          "areal_thermal_resistance() = %.6e -> %.1e at 2 s.f. "
          "(paper: 2.8e-7)"
          % (product, float(sh._round_sig(product, AREAL_PAPER_SIGFIGS))))
    check("P3 and it also sits inside the +-%.4f %% the two printed figures "
          "license, with room to spare" % (100.0 * hw),
          abs(product - AREAL_PAPER) / AREAL_PAPER < hw,
          "residual %.4f %% of %.1e, licensed +-%.4f %%"
          % (100.0 * (product - AREAL_PAPER) / AREAL_PAPER, AREAL_PAPER,
             100.0 * hw))
    # -- the planted quotient ---------------------------------------------
    check("P4 PLANTED CONTROL: the quotient R_th/A does NOT round to 2.8e-7 "
          "at any precision, so P2 is testing the product and not merely "
          "'some area-combined number'",
          not sh._rounds_to(quotient, AREAL_PAPER, AREAL_PAPER_SIGFIGS),
          "R_th/A = %.6e m^-2.K/W, which is %.1f decades above 2.8e-7"
          % (quotient, np.log10(quotient / AREAL_PAPER)))
    check("P5 and the quotient is 21 decades off, not the 7 the 11-11 "
          "prediction stated -- the prediction compared R_th/A against "
          "R_th's own 1e4 rather than against the areal 2.8e-7",
          round(np.log10(quotient / AREAL_PAPER)) == 21
          and round(np.log10(quotient / R_TH_PAPER)) == 11,
          "log10(quotient/2.8e-7) = %.3f and log10(quotient/1e4) = %.3f; "
          "the predicted '7 decades' matches neither reference, so the "
          "prediction named the right trap and no correct size for it"
          % (np.log10(quotient / AREAL_PAPER),
             np.log10(quotient / R_TH_PAPER)))
    check("P6 the product carries the unit the paper prints and the quotient "
          "does not: scaling W and L by s scales A by s^2, R_th by ~1/s^2, "
          "so R_th*A is invariant to leading order while R_th/A moves by s^4",
          0.11 < (areal_thermal_resistance(W=2 * W_DORGAN, L=2 * L_DORGAN)
                  / product - 1.0) < 0.13,
          "doubling both edges moves R_th*A by %.2f %% (the sqrt(A) term "
          "only) while R_th/A falls by a factor %.1f"
          % (100.0 * (areal_thermal_resistance(W=2 * W_DORGAN,
                                               L=2 * L_DORGAN) / product - 1),
             quotient / (sh.thermal_resistances(W=2 * W_DORGAN,
                                                L=2 * L_DORGAN,
                                                t_ox=T_OX_DORGAN)["R_th"]
                         / (4 * A))))
    return product, quotient, A, R


def verify_independence_is_false(product, A, R):
    """I: the exact zeros that refute the 11-11 independence claim."""
    r_areal = (product - AREAL_PAPER) / AREAL_PAPER
    r_total = (R - R_TH_PAPER) / R_TH_PAPER

    # I1 AND I2 WERE FIRST WRITTEN AS EXACT EQUALITIES ("== 1.0e4", "== 0.0")
    # AND BOTH ARE FALSE, caught by running them before committing.  The
    # algebra is exact; the arithmetic is not, because 2.8e-7, 2.8e-11 and
    # 4e-6*7e-6 are none of them exactly representable (see P1).  The claim
    # that survives is stronger than a tolerance and weaker than an identity:
    # the two agree to the last place of 1.0, while two genuinely independent
    # constraints would differ at the 1e-2 level the residuals occupy.
    check("I1 2.8e-7 / A agrees with the paper's own 1.0e4 to within one ulp "
          "but is NOT exactly equal to it (%.17e), so the two published "
          "numbers are one statement divided by a constant that binary "
          "floating point cannot hold exactly"
          % (AREAL_PAPER / A),
          AREAL_PAPER / A != R_TH_PAPER
          and abs(AREAL_PAPER / A - R_TH_PAPER) <= np.spacing(R_TH_PAPER),
          "2.8e-7/A = %.17e, 1.0e4 = %.17e, difference %.3e = %.1f ulp"
          % (AREAL_PAPER / A, R_TH_PAPER, AREAL_PAPER / A - R_TH_PAPER,
             abs(AREAL_PAPER / A - R_TH_PAPER) / np.spacing(R_TH_PAPER)))
    check("I2 the relative residual of R_th*A against 2.8e-7 and that of "
          "R_th against 1e4 agree to 14 significant figures, their "
          "difference being one ulp of 1.0 and nothing more -- the per-area "
          "check is NOT independent of A2, it IS A2 with a different "
          "licensed width (refutes the 11-11 ranking's stated reason)",
          abs(r_areal - r_total) <= 4.0 * np.spacing(1.0)
          and abs(r_areal - r_total) < 1e-10 * abs(r_total),
          "areal residual %.17e, total residual %.17e, difference %.3e "
          "(= %.1f ulp of 1.0); independent constraints would differ at "
          "~%.1e, which is %.0f decades larger"
          % (r_areal, r_total, r_areal - r_total,
             abs(r_areal - r_total) / np.spacing(1.0), abs(r_total),
             np.log10(abs(r_total) / abs(r_areal - r_total))))
    check("I2b in consequence NOTHING can fail one and pass the other at "
          "equal width, which is the precise denial of 'it can fail for a "
          "reason A2 cannot see'",
          (abs(r_areal) < 0.05) == (abs(r_total) < 0.05),
          "both %.4f %% at the same centre; only the licensed widths differ "
          "(+-1.79 %% vs the +-5 %% A2 coded)" % (100.0 * r_areal))

    d = areal_decomposition()
    check("I3 EXACT ZERO: R_th*A == 1/h + t_ox/k_ox + sqrt(A)/(2 k_Si) to "
          "machine precision -- Section 4's identity, which is what makes "
          "two of the three terms area-independent",
          abs(d["total"] - product) <= 4.0 * np.spacing(product),
          "identity %.17e vs R_th*A %.17e, residual %.1e (<= 4 ulp = %.1e)"
          % (d["total"], product, abs(d["total"] - product),
             4.0 * np.spacing(product)))
    check("I4 the decomposition's three terms sum to the total with no "
          "fourth term hiding, and each is positive",
          abs(d["boundary"] + d["oxide"] + d["spreading"] - d["total"]) == 0.0
          and min(d["boundary"], d["oxide"], d["spreading"]) > 0.0,
          "boundary %.4e + oxide %.4e + spreading %.4e"
          % (d["boundary"], d["oxide"], d["spreading"]))
    return d


def verify_significant_figures():
    """F: the 11-11 'one significant figure licenses +-1.8 %' correction."""
    hw2 = relative_halfwidth(AREAL_PAPER, 2)
    hw1 = relative_halfwidth(AREAL_PAPER, 1)
    hw1_rounded = relative_halfwidth(sh._round_sig(AREAL_PAPER, 1), 1)

    check("F1 2.8e-7 printed to TWO figures licenses +-1.7857 %, which is "
          "the +-1.8 % both the 11-11 note and the 11-11 log quote",
          abs(hw2 * 2.8 - 0.05) < 1e-17 and round(100.0 * hw2, 1) == 1.8,
          "two-figure half-width = %.6f %% -> 1.8 %% at 2 s.f."
          % (100.0 * hw2))
    check("F2 but 2.8e-7 HAS two significant figures, not the one those two "
          "places state -- the quoted tolerance was right and its stated "
          "reason wrong",
          sh._round_sig(AREAL_PAPER, 2) == AREAL_PAPER
          and sh._round_sig(AREAL_PAPER, 1) != AREAL_PAPER,
          "2.8e-7 is unchanged by rounding to 2 s.f. and becomes %r at 1 "
          "s.f., so it carries 2" % sh._round_sig(AREAL_PAPER, 1))
    check("F3 a genuine ONE-figure reading of the same datum would license "
          ">= 9x the width used: %.1f %% on the printed mantissa, %.1f %% "
          "after rounding to 3e-7 -- so the stated reason, had it been "
          "followed, would have made this check nearly vacuous"
          % (100.0 * hw1, 100.0 * hw1_rounded),
          hw1 / hw2 >= 9.0 and hw1_rounded / hw2 >= 9.0,
          "one-figure/two-figure width ratio = %.2f (printed mantissa) and "
          "%.2f (rounded to 3e-7)" % (hw1 / hw2, hw1_rounded / hw2))
    check("F4 and R_th's own ~1e4, at ONE printed figure, licenses +-50 % -- "
          "which is why A2's coded +-5 % was called over-tight on 11-11 and "
          "why Section 7 reverses that: the SAME sentence licenses +-1.79 %",
          relative_halfwidth(R_TH_PAPER, R_TH_PAPER_SIGFIGS) == 0.5,
          "one-figure half-width on 1e4 = %.1f %%, two-figure half-width on "
          "2.8e-7 = %.4f %%, ratio %.1f"
          % (100.0 * relative_halfwidth(R_TH_PAPER, 1), 100.0 * hw2,
             0.5 / hw2))
    return hw2


def verify_area_scaling(d, hw2):
    """S: the exact limit, the measured band, and the repo's own device."""
    check("S1 EXACT LIMIT: as A -> 0 the areal figure of merit tends to the "
          "A-independent floor 1/h + t_ox/k_ox, because only R_Si carries a "
          "sqrt(A)",
          abs(areal_thermal_resistance(W=1e-30, L=1e-30) - d["floor"])
          <= 4.0 * np.spacing(d["floor"]),
          "limit %.17e vs floor %.17e"
          % (areal_thermal_resistance(W=1e-30, L=1e-30), d["floor"]))
    check("S2 the floor is NOT inside the paper's two-figure band, so the "
          "areal number cannot be treated as scale-free: the spreading term "
          "contributes %.2f %% of Dorgan's own value and is the only part "
          "that carries the geometry"
          % (100.0 * d["spreading"] / d["total"]),
          abs(d["floor"] - AREAL_PAPER) / AREAL_PAPER > hw2,
          "floor %.4e is %.2f %% below 2.8e-7, licensed +-%.2f %%"
          % (d["floor"], 100.0 * (AREAL_PAPER - d["floor"]) / AREAL_PAPER,
             100.0 * hw2))

    b = areal_band()
    span = b["A_hi"] / b["A_lo"]
    check("S3 the band in AREA over which R_th*A stays inside the paper's "
          "own two-figure band spans a factor %.3f, and it brackets "
          "Dorgan's 2.8e-11 m^2" % span,
          b["A_lo"] < 2.8e-11 < b["A_hi"] and 1.5 < span < 2.5,
          "A in [%.4e, %.4e] m^2, i.e. [%.3fx, %.3fx] of Dorgan's area; "
          "factor %.4f wide"
          % (b["A_lo"], b["A_hi"], b["A_lo"] / 2.8e-11, b["A_hi"] / 2.8e-11,
             span))

    A_repo = gfet.W * gfet.L
    areal_repo = areal_thermal_resistance(W=gfet.W, L=gfet.L,
                                          t_ox=T_OX_DORGAN)
    err = (AREAL_PAPER - areal_repo) / areal_repo
    check("S4 this repository's own FET (W = 1 um, L = 200 nm, A = %.1e m^2) "
          "is %.0fx smaller in area than Dorgan's and therefore FAR outside "
          "the S3 band: reusing 2.8e-7 as a scale-free areal figure of merit "
          "there would over-state R_th by %.1f %%"
          % (A_repo, 2.8e-11 / A_repo, 100.0 * err),
          A_repo < b["A_lo"] and err > 0.10,
          "areal value at the repo geometry = %.4e m^2.K/W (floor %.4e, "
          "spreading only %.3f %% of it), vs the paper's 2.8e-7 -- an "
          "over-statement of %.2f %%"
          % (areal_repo, d["floor"],
             100.0 * (areal_repo - d["floor"]) / areal_repo, 100.0 * err))
    check("S5 and the error in S4 is bounded by the floor, not unbounded: no "
          "geometry can push the areal figure of merit below 1/h + "
          "t_ox/k_ox, so the over-statement from reusing 2.8e-7 can never "
          "exceed 2.8e-7/floor - 1",
          err < AREAL_PAPER / d["floor"] - 1.0 + 1e-12,
          "worst-case over-statement = %.2f %% (at A -> 0); at the repo "
          "geometry %.2f %% of that worst case is already realised"
          % (100.0 * (AREAL_PAPER / d["floor"] - 1.0),
             100.0 * err / (AREAL_PAPER / d["floor"] - 1.0)))
    return b, areal_repo, err


# ---------------------------------------------------------------------------
# MUTATION.  Inject a defect, re-run the whole suite, require it to fail --
# and record WHICH check caught it against which check was predicted to.
# The three-valued reporting is borrowed from the sibling repository's
# mutation harnesses: DETECTED-AS-PREDICTED, DETECTED-ELSEWHERE and
# UNDETECTED are different results and only the last is a defect in the
# suite, but the middle one means the predicted check is not doing the work
# its name claims.
# ---------------------------------------------------------------------------
def _run_suite():
    """Run every check into a fresh list and return {name: ok}."""
    global _CHECKS
    saved, _CHECKS = _CHECKS, []
    try:
        product, quotient, A, R = verify_published_comparison()
        d = verify_independence_is_false(product, A, R)
        hw2 = verify_significant_figures()
        verify_area_scaling(d, hw2)
        out = {n: ok for n, ok, _ in _CHECKS}
    except Exception as exc:                       # a mutant may crash
        out = {"CRASH": False, "exception": repr(exc)}
    finally:
        _CHECKS = saved
    return out


def _mutants():
    g = globals()
    base_areal = areal_thermal_resistance
    base_decomp = areal_decomposition
    base_hw = relative_halfwidth

    def quotient_not_product(**kw):
        t = sh.thermal_resistances(W=kw.get("W", W_DORGAN),
                                   L=kw.get("L", L_DORGAN),
                                   t_ox=kw.get("t_ox", T_OX_DORGAN))
        return t["R_th"] / t["A"]

    def drop_spreading(**kw):
        d = base_decomp(**kw)
        d["spreading"] = 0.0
        d["total"] = d["floor"]
        return d

    def drop_boundary(**kw):
        d = base_decomp(**kw)
        d["boundary"] = 0.0
        d["floor"] = d["oxide"]
        d["total"] = d["oxide"] + d["spreading"]
        return d

    def linear_not_sqrt(**kw):
        d = base_decomp(**kw)
        d["spreading"] = d["A"] / (2.0 * sh.kappa_Si(sh.T_AMBIENT))
        d["total"] = d["floor"] + d["spreading"]
        return d

    def negative_spreading(**kw):
        d = base_decomp(**kw)
        d["spreading"] = -d["spreading"]
        d["total"] = d["floor"] + d["spreading"]
        return d

    def halfwidth_without_mantissa(quoted, sigfigs):
        return 0.5 * 10.0 ** (1 - sigfigs)         # the decade form

    return [
        ("M1 areal figure of merit computed as R_th/A instead of R_th*A "
         "(Section 2's trap, injected)",
         lambda: g.__setitem__("areal_thermal_resistance",
                               quotient_not_product),
         lambda: g.__setitem__("areal_thermal_resistance", base_areal),
         ("P2", "P3")),
        ("M2 spreading term dropped from the identity",
         lambda: g.__setitem__("areal_decomposition", drop_spreading),
         lambda: g.__setitem__("areal_decomposition", base_decomp),
         ("I3",)),
        ("M3 boundary term 1/h dropped from the identity",
         lambda: g.__setitem__("areal_decomposition", drop_boundary),
         lambda: g.__setitem__("areal_decomposition", base_decomp),
         ("I3",)),
        ("M4 R_Si written with A instead of sqrt(A) -- the error that would "
         "make the areal figure of merit look scale-free",
         lambda: g.__setitem__("areal_decomposition", linear_not_sqrt),
         lambda: g.__setitem__("areal_decomposition", base_decomp),
         ("I3",)),
        ("M5 spreading term given the wrong sign",
         lambda: g.__setitem__("areal_decomposition", negative_spreading),
         lambda: g.__setitem__("areal_decomposition", base_decomp),
         ("I4",)),
        ("M6 relative_halfwidth loses its mantissa division, i.e. the exact "
         "defect Section 6 corrects, injected into the machinery",
         lambda: g.__setitem__("relative_halfwidth",
                               halfwidth_without_mantissa),
         lambda: g.__setitem__("relative_halfwidth", base_hw),
         ("F1",)),
        ("M7 the paper's number read at ONE figure as 3e-7, which is what "
         "the 11-11 sentence literally said",
         lambda: g.__setitem__("AREAL_PAPER", 3e-7),
         lambda: g.__setitem__("AREAL_PAPER", 2.8e-7),
         ("P2",)),
    ]


def run_mutations():
    """Returns (rows, n_detected, n_as_predicted)."""
    base = _run_suite()
    assert all(base.values()), "mutation run started from a failing suite"
    rows = []
    for name, apply_mut, undo, predicted in _mutants():
        apply_mut()
        try:
            got = _run_suite()
        finally:
            undo()
        failing = sorted(n for n, ok in got.items() if ok is False)
        hit = [f for f in failing
               if any(f.startswith(p + " ") or f == p for p in predicted)]
        verdict = ("UNDETECTED" if not failing
                   else "DETECTED-AS-PREDICTED" if hit
                   else "DETECTED-ELSEWHERE")
        rows.append({"mutant": name, "predicted": predicted,
                     "verdict": verdict,
                     "n_failing": len(failing),
                     "first_failing": failing[0] if failing else None})
    n_det = sum(1 for r in rows if r["verdict"] != "UNDETECTED")
    n_pred = sum(1 for r in rows if r["verdict"] == "DETECTED-AS-PREDICTED")
    return rows, n_det, n_pred


# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
def make_figure(fname="areal_thermal_resistance.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = areal_decomposition()
    b = areal_band()
    hw = b["halfwidth"]
    A_grid = np.logspace(-14.0, -9.0, 400)
    areal = np.array([areal_thermal_resistance(W=np.sqrt(a), L=np.sqrt(a),
                                               t_ox=T_OX_DORGAN)
                      for a in A_grid])
    R_grid = np.array([sh.thermal_resistances(W=np.sqrt(a), L=np.sqrt(a),
                                              t_ox=T_OX_DORGAN)["R_th"]
                       for a in A_grid])
    A_dorgan, A_repo = 2.8e-11, gfet.W * gfet.L

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.2))

    ax.semilogx(A_grid, areal * 1e7, lw=2.0, color="#1f4e79",
                label=r"$R_{th}\,A$ (Dorgan Eq. 5 $\times A$)")
    ax.axhline(d["floor"] * 1e7, ls="--", lw=1.4, color="#777777",
               label=r"area-independent floor $1/h + t_{ox}/\kappa_{ox}$")
    ax.axhspan(AREAL_PAPER * (1 - hw) * 1e7, AREAL_PAPER * (1 + hw) * 1e7,
               color="#2e8b57", alpha=0.18,
               label=r"paper's $2.8\times10^{-7}\pm1.79\%$ (2 s.f.)")
    ax.axhline(AREAL_PAPER * 1e7, color="#2e8b57", lw=1.2)
    for a, lbl, c in ((b["A_lo"], None, "#2e8b57"),
                      (b["A_hi"], None, "#2e8b57")):
        ax.axvline(a, ls=":", lw=1.2, color=c)
    ax.plot([A_dorgan], [areal_thermal_resistance() * 1e7], "o", ms=9,
            color="#c00000", zorder=5,
            label=r"Dorgan $A=2.8\times10^{-11}\,$m$^2$: 2.776")
    ax.plot([A_repo], [areal_thermal_resistance(W=gfet.W, L=gfet.L,
                                                t_ox=T_OX_DORGAN) * 1e7],
            "s", ms=9, color="#7030a0", zorder=5,
            label=r"this repo's FET $A=2\times10^{-13}\,$m$^2$: 2.474")
    ax.set_xlabel(r"device area $A = LW$  (m$^2$)")
    ax.set_ylabel(r"$R_{th}A$  ($10^{-7}$ m$^2\!\cdot\!$K/W)")
    ax.set_title("The areal figure of merit is not scale-free:\n"
                 "the paper's own 2 s.f. band spans a factor %.2f in area"
                 % (b["A_hi"] / b["A_lo"]), fontsize=11)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.3)
    ax.annotate("factor %.2f" % (b["A_hi"] / b["A_lo"]),
                xy=(np.sqrt(b["A_lo"] * b["A_hi"]), 2.60), fontsize=8,
                ha="center", color="#2e8b57")

    ax2.loglog(A_grid, R_grid, lw=2.0, color="#1f4e79",
               label=r"$R_{th}\propto 1/A$ (two of three terms)")
    ax2.axhline(R_TH_PAPER, color="#2e8b57", lw=1.2,
                label=r"paper's $\sim10^4$ K/W (1 s.f.)")
    ax2.plot([A_dorgan], [sh.thermal_resistances(W=7e-6, L=4e-6,
                                                 t_ox=T_OX_DORGAN)["R_th"]],
             "o", ms=9, color="#c00000", zorder=5)
    ax2.plot([A_repo], [sh.thermal_resistances(W=gfet.W, L=gfet.L,
                                               t_ox=T_OX_DORGAN)["R_th"]],
             "s", ms=9, color="#7030a0", zorder=5)
    ax2.set_xlabel(r"device area $A = LW$  (m$^2$)")
    ax2.set_ylabel(r"$R_{th}$  (K/W)")
    ax2.set_title("Why normalising is worth doing at all:\n"
                  r"$R_{th}$ moves %.0f$\times$ over the same range that "
                  r"$R_{th}A$ moves %.0f %%"
                  % (R_grid[0] / R_grid[-1],
                     100.0 * (areal[-1] / areal[0] - 1.0)), fontsize=11)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    return fname


# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("AREAL THERMAL RESISTANCE: R_th PER UNIT AREA vs DORGAN'S 2.8e-7")
    print("Closes the item ranked FIRST on 2026-11-11, and refutes the")
    print("reason it was ranked first.")
    print("=" * 78)

    product, quotient, A, R = verify_published_comparison()
    print("\nSECTION 1-2.  The published comparison and the planted quotient")
    print("  A        = %.6e m^2   (L*W, one ulp below the literal 2.8e-11)"
          % A)
    print("  R_th     = %.4f K/W" % R)
    print("  R_th * A = %.6e m^2.K/W   <- the paper's quantity" % product)
    print("  R_th / A = %.6e m^-2.K/W  <- the trap, %.1f decades away"
          % (quotient, np.log10(quotient / AREAL_PAPER)))

    d = verify_independence_is_false(product, A, R)
    print("\nSECTION 3-4.  Independence refuted, and the identity")
    print("  relative residual of R_th*A vs 2.8e-7 : %.6f %%"
          % (100.0 * (product - AREAL_PAPER) / AREAL_PAPER))
    print("  relative residual of R_th  vs 1.0e4   : %.6f %%"
          % (100.0 * (R - R_TH_PAPER) / R_TH_PAPER))
    print("  they differ by %.2e, one ulp of 1.0 -- the same constraint"
          % abs(((product - AREAL_PAPER) / AREAL_PAPER)
                - ((R - R_TH_PAPER) / R_TH_PAPER)))
    print("  identity  R_th*A = 1/h + t_ox/k_ox + sqrt(A)/(2 k_Si):")
    print("    1/h            = %.6e m^2.K/W  (%5.2f %%)  area-independent"
          % (d["boundary"], 100.0 * d["boundary"] / d["total"]))
    print("    t_ox/k_ox      = %.6e m^2.K/W  (%5.2f %%)  area-independent"
          % (d["oxide"], 100.0 * d["oxide"] / d["total"]))
    print("    sqrt(A)/(2kSi) = %.6e m^2.K/W  (%5.2f %%)  GEOMETRY"
          % (d["spreading"], 100.0 * d["spreading"] / d["total"]))
    print("    sum            = %.6e m^2.K/W  (residual %.1e)"
          % (d["total"], abs(d["total"] - product)))

    hw2 = verify_significant_figures()
    print("\nSECTION 6.  Significant figures")
    print("  2.8e-7 at 2 s.f. licenses +-%.4f %%   <- the +-1.8 %% quoted"
          % (100.0 * hw2))
    print("  2.8e-7 read at 1 s.f. would license +-%.2f %% (or +-%.2f %% as "
          "3e-7)" % (100.0 * relative_halfwidth(AREAL_PAPER, 1),
                     100.0 * relative_halfwidth(3e-7, 1)))
    print("  1e4 at 1 s.f. licenses +-%.1f %%"
          % (100.0 * relative_halfwidth(R_TH_PAPER, 1)))
    print("  => the 11-11 note and log said 'one significant figure "
          "licenses")
    print("     +-1.8 %'.  The width is right; the reason is not.")

    b, areal_repo, err = verify_area_scaling(d, hw2)
    print("\nSECTION 5.  Where the independence actually lives")
    print("  A-independent floor      = %.6e m^2.K/W" % d["floor"])
    print("  paper's 2 s.f. band      = [%.4e, %.4e] m^2.K/W"
          % (b["value_lo"], b["value_hi"]))
    print("  => admissible AREA range = [%.4e, %.4e] m^2, a factor %.3f"
          % (b["A_lo"], b["A_hi"], b["A_hi"] / b["A_lo"]))
    print("     (Dorgan's area sits at %.3fx the low edge)"
          % (2.8e-11 / b["A_lo"]))
    print("  this repo's FET: A = %.2e m^2, areal = %.5e m^2.K/W"
          % (gfet.W * gfet.L, areal_repo))
    print("  reusing 2.8e-7 there over-states R_th by %.2f %% "
          "(worst case %.2f %%)"
          % (100.0 * err, 100.0 * (AREAL_PAPER / d["floor"] - 1.0)))

    print("\nSECTION 7.  The ranked-third item, DECIDED")
    print("  A2's coded +-5 % is not over-tight -- it is 2.8x LOOSER than")
    print("  the same sentence licenses.  Tighten to +-%.4f %%, do NOT widen"
          % (100.0 * hw2))
    print("  to the one-figure +-50 %%.  Residual to spare: %.4f %%."
          % (100.0 * hw2 - abs(100.0 * (R - R_TH_PAPER) / R_TH_PAPER)))

    print("\nSECTION 8.  Mutation")
    rows, n_det, n_pred = run_mutations()
    for r in rows:
        print("  %-22s %-2d failing  predicted %-10s  %s"
              % (r["verdict"], r["n_failing"], "/".join(r["predicted"]),
                 r["mutant"][:46]))
    print("  %d/%d mutants detected, %d of them by the predicted check"
          % (n_det, len(rows), n_pred))

    print("\nSECTION 9.  Figure")
    print("  wrote %s" % make_figure())

    print("\n" + "=" * 78)
    passed = sum(1 for _, ok, _ in _CHECKS if ok)
    failed = sum(1 for _, ok, _ in _CHECKS if not ok)
    for name, ok, detail in _CHECKS:
        print("  [%s] %s" % ("PASS" if ok else "FAIL", name))
        print("         %s" % detail)
    print("=" * 78)
    print("%d checks, %d passing, %d failing; %d/%d mutants detected"
          % (passed + failed, passed, failed, n_det, len(rows)))
    print("=" * 78)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
