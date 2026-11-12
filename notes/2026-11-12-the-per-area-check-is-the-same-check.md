# The per-area check is the same check, and the tolerance was too loose

**2026-11-12.** The item ranked first on 2026-11-11 closes. It passes, and the
reason it was ranked first was wrong. Along the way the item ranked third is
decided against its own framing, a significant-figure claim made twice
yesterday is corrected, the tolerance census is made recursive because its
first cross-repository run returned the empty set, and the census produces its
first false positive — on the very check this session repaired.

---

## 1. The item, and the sentence it rests on

The 2026-11-11 standing list ranked first:

> **CHECK `R_th` PER UNIT AREA AGAINST THE PAPER'S `~2.8e-7 m^2.K/W`** —
> created today, **ranked first**, cheap, and the one item today's literature
> read put on the table with a number already in hand. One significant figure
> licenses ±1.8 %. It is an *independent* constraint: A2 tests `R_th` alone and
> therefore tests the area and the conductivity only through their product,
> while this tests the ratio. Of the two it is the one that can fail for a
> reason A2 cannot see.

The datum is one sentence of Dorgan, Bae & Pop:

> "At 300 K for our geometry `R_th ≈ 10^4 K/W`, or `~2.8×10^-7 m^2 K/W` per
> unit of device area"

with `W = 7 μm`, `L = 4 μm`, `t_ox = 300 nm`.

## 2. The result

`graphene_areal_thermal_resistance.py`, 20 checks, 20 passing, 7/7 mutants
detected and all 7 by the predicted check.

    R_th     = 9915.5790 K/W
    A = L*W  = 2.8e-11 m^2
    R_th * A = 2.776362e-07 m^2.K/W   ->  2.8e-7 at 2 s.f.   PASSES
    R_th / A = 3.541278e+14 m^-2.K/W  ->  21.1 decades away

The residual is −0.844 % against the ±1.79 % the two printed figures license.

**The trap was real and its predicted size was wrong.** Yesterday's prediction
(1) said a session would compute `R_th / A` first and "have a number seven
decades off". It is 21.1 decades off against the areal number, and 10.5 against
`R_th`'s own `1e4`. The prediction named the right trap and no correct size for
it; it had compared the quotient to the wrong one of the paper's two numbers.
Checks P4 and P5 keep the quotient as a planted control.

## 3. Why the ranking's reason is false

`A = L*W` is given by the stated geometry, not fitted. So

    2.8e-7 / A  =  1.0e4   to within one ulp,

which makes the paper's two numbers one statement divided by a constant. The
consequence is measured, not argued:

    relative residual of R_th*A against 2.8e-7 : -8.44209690618646577e-03
    relative residual of R_th   against 1.0e4  : -8.44209690618627148e-03
    difference                                 : -1.94e-16  (0.9 ulp of 1.0)

Two genuinely independent constraints would differ at the `8.4e-3` the
residuals themselves occupy — fourteen decades larger. Every clause of the
ranking's reason fails:

- A2 does not test `R_th` "alone" in any sense that excludes `A`; `R_th` is
  *computed* from `A` by Eq. (5).
- The areal quantity is not the "ratio" of `R_th` and `A`, it is their product.
- Nothing can fail one and pass the other at equal width, so neither "can fail
  for a reason the other cannot see".

**The per-area check is A2, with the same centre and a 2.8× tighter licensed
width.** That is still worth having. It is not what the item said it was.

### An error of mine, inside the same section

I1 and I2 were first written as *exact* equalities — `2.8e-7/A == 1.0e4` and
`residual difference == 0.0` — and both **failed on the first run**. `4e-6*7e-6`
is `2.7999999999999997e-11`, one ulp below the literal `2.8e-11`: the identity
is exact in algebra and 1-ulp in binary. The surviving claim is sharper than a
tolerance and weaker than an identity, and the module says so in place rather
than quietly asserting the weaker thing. P1, which had been written as
`A == L*W` — a tautology that could not fail — now carries the representability
fact instead.

## 4. Where the independence actually lives

Multiplying Eq. (5) through by `A`:

    R_th * A  ==  1/h  +  t_ox/k_ox  +  sqrt(A)/(2 k_Si)
                  \_________________/   \_______________/
                    area-independent      grows as sqrt(A)

I3 requires this to machine precision and the residual is exactly `0.0`. Two of
the three terms carry no geometry at all, because `R_B` and `R_ox` go as `1/A`
while `R_Si` goes as `1/sqrt(A)`. Numerically, at Dorgan's geometry:

| term | value (m²·K/W) | share | geometry? |
|---|---|---|---|
| `1/h` | 1.000000e-08 | 3.60 % | no |
| `t_ox/k_ox` | 2.345643e-07 | 84.49 % | no |
| `sqrt(A)/(2 k_Si)` | 3.307189e-08 | 11.91 % | **yes** |
| total | 2.776362e-07 | | |

So the areal figure of merit is scale-free only to the extent the spreading
term is small, and it is 11.91 % of the value.

### The new number: a factor 1.765 in area

Inverting the identity for the edges of the paper's own two-figure band gives
the range of **areas** over which `R_th*A` stays inside it:

    A in [2.3714e-11, 4.1857e-11] m^2  =  [0.847x, 1.495x] of Dorgan's area
    a factor 1.765 wide

"Per unit of device area" invites reuse of `2.8e-7` as a transferable areal
figure of merit. It is transferable within a factor of about 1.8 in area and no
further. **This repository's own FET is 140× smaller** (`W = 1 μm`,
`L = 200 nm`, `A = 2e-13 m²`), far outside the band, with an areal value of
`2.47359e-07` m²·K/W — only 1.13 % of which is spreading. Reusing `2.8e-7`
there would over-state `R_th` by **13.20 %**, against a hard worst case of
14.49 % set by the floor as `A → 0` (S4, S5). The bound is directional and
finite, which is worth more than the point value.

## 5. The ranked-third item, decided against its own framing

> **REPAIR A2's OVER-TIGHT TOLERANCE, OR DECIDE NOT TO** — ranked third and
> genuinely a *decision*: ±5 % against a value the paper prints to one figure
> is a tighter claim than the source supports, but it is also the only thing
> currently pinning `R_th` to its decade. Widening it to the licensed ±50 %
> would make it nearly vacuous.

The dilemma was false. The paper prints the same fact to **two** figures in the
**same sentence**, and that half licenses ±1.7857 % on the identical
constraint. So ±5 % was never a tighter claim than the source supports — it was
**2.8× looser**. A2's defect was its *citation*, not its width: it quoted the
one-figure half of a two-figure sentence and then used a width neither half
licenses.

A2 is therefore **tightened to ±1.7857 %, not widened to ±50 %**, and passes
with 0.94 % of the licensed width to spare. A2a keeps the old ±5 % comparison,
which also passes, so the record shows the repair changed the warrant and not
the verdict. A2b is the per-area comparison itself. A2c asserts A2 and A2b are
*one* constraint rather than leaving a reader to infer independence — which is
exactly the inference the 11-11 list made wrongly.

This is the **third** time a census finding in this tree has turned out to be
about a check's warrant rather than its arithmetic.

## 6. "One significant figure licenses ±1.8 %" — right width, wrong reason

Both yesterday's note (§9) and yesterday's log said this of `2.8e-7`. It has
**two** significant figures. Two figures license ±1.7857 %, which *is* the
quoted number; a genuine one-figure reading would license 16.7–17.9 %, nine to
ten times wider. The width carried forward was correct and the arithmetic
offered for it was not — and had the stated reason been followed the check would
have been nearly vacuous. F1–F4 pin all three widths.

## 7. The census: an empty set, and a false positive

The 11-11 list also created **RUN THE TOLERANCE CENSUS OVER
`Design-Verification-Roadmap`**. Run as it stood, `census_repo(root=...)` read
**0 Python files and returned 0 rows** over that repository — not because it has
no tolerances, but because `census_repo` was a flat `os.listdir` and all 47 of
that tree's Python files live under `tools/` or `examples/`. **An empty census
is indistinguishable from a clean one.** The instrument's hidden assumption was
invisible in the repository it was written for; that is the generalisable
finding, and it is why `recursive` now defaults to `True` rather than preserving
the old behaviour.

Recursively: **6 rows over 47 files, all `NO_QUOTE`, zero adjudicated.**

### Scoring both repositories' predictions

| prediction | claim | verdict |
|---|---|---|
| this tree, 11-11 (3) | more `NO_QUOTE` rows than this tree's 102; zero adjudicated | **half right** — zero adjudicated holds; 6 vs 102 is wrong by 17× *in the opposite direction* |
| sibling, 11-11 (2) | zero adjudicated; *not* more than 102, because its suites check exact values | **both substantive halves right**, reason included |
| sibling, 11-11 (2) meta | "the sibling's prediction is right in its first half and wrong in its second" | **inverted** — the first half was the count, which failed |

And both were **unscorable as posed** until today: each named a population in a
run that, with the instrument as it then stood, returned the empty set.

### Two limitations, found by repairing a check the instrument had flagged

- **9l, the census's first false positive.** It assumes a check's *centre* and
  its *width* are warranted by the same printed number. Repaired A2 takes its
  centre from the one-figure `~1e4` and its width from the two-figure
  `~2.8e-7`, so the census licenses ±50 %, sees ±1.7857 %, and returns
  `OVER_TIGHT` at ratio 28. Wrong verdict, right arithmetic. This is a second
  stated limit alongside 9k (a label citing a value its source never printed).
- **9m, a measured blind spot.** Name-bound *residuals* are resolved through
  `bound_via` (40+ rows); name-bound *tolerances* are not, because `tol` is read
  only from a literal comparator. **7 such comparisons in 7 modules.** A2 was
  the eighth for the length of one edit today — binding its half-width to the
  expression `0.05/2.8` made it `LIT_NO_TEST`, the census's own worst category.
  Repairing A2's warrant had briefly made it unauditable.

Check 8d, written yesterday to witness that the tree still carried live
instances of both defect classes, is annotated rather than restated: it still
passes, but the genuine `OVER_TIGHT` it counted was repaired today and the row
it now counts is 9l's false positive.

## 8. One more guard that earned itself immediately

A2's half-width is now a literal, with `0.05/2.8` in a comment, because the
expression form is what the census cannot see. A2d was added to stop the literal
drifting from its derivation — and **caught a one-ulp error in that literal on
its first run**: `0.017857142857142856` against the true
`0.01785714285714286`. The failure is kept in the comment as the evidence that
A2d does anything.

## 9. Counts

New module 20/20, 7/7 mutants. `self_heating` 46 → 47 checks (A2 re-warranted,
A2a/A2b/A2c/A2d added), 0 failing, mutation suite still 10/10. Census 66 → 69
checks, 0 failing. Three errors of mine caught by running rather than reading
(two exact-equality claims, one tautological check, one unauditable tolerance,
one one-ulp literal — five, counting properly). One ranked-first item closed and
its justification refuted; one ranked-third item decided; one cross-repository
item closed; two instrument limitations recorded.

## Sources

- Dorgan, Bae & Pop, "Mobility and Saturation Velocity in Graphene on SiO2",
  *Appl. Phys. Lett.* **97**, 082112 (2010) —
  [arXiv:1005.2711](https://arxiv.org/pdf/1005.2711) ·
  [publisher PDF incl. SOM](https://poplab.stanford.edu/pdfs/Dorgan-GrapheneVsat-apl10.pdf).
  Main text for the `R_th ≈ 10^4 K/W` / `~2.8×10^-7 m^2 K/W` sentence and the
  geometry; supp. Sec. 6 for `kappa_ox` and `kappa_Si`; supp. Sec. 5 for
  `kappa_ox ≈ 1.3 W/m·K`. **No new source was fetched today** — the whole
  session ran on the datum yesterday's read already had in hand, which is what
  made the item cheap.
- `notes/2026-11-11-a-census-of-comparisons-against-published-numbers.md` §9,
  which created the item and carries the significant-figure error §6 corrects.
- `notes/2026-11-10-self-heating-the-thermal-network-and-two-exact-benchmarks.md`
  for the thermal network and the two exact benchmarks A2 belongs to.
