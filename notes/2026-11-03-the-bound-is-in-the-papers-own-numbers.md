# The bound was in the paper's own numbers, and the mobility it was asked for is not in the paper at all

2026-11-03. Literature record for `graphene_kappa_impurity_bound.py`.

## 1. What was fetched, and what was there

**Source.** Zion, Haran, Sharoni, Shaul, Naveh *et al.*, "Charge carrier
transport asymmetry in monolayer graphene", arXiv:1801.03565. Full text read
today via `https://ar5iv.labs.arxiv.org/html/1801.03565`. The abstract page
was the only part reachable on 2026-10-30 and 2026-11-02; ar5iv served the
whole paper today without incident.

**The item this was meant to discharge** was ranked first on 11-02:

> RESOLVE THE `A` PAIRING ... Blocked on **a datum in the full text of
> arXiv:1801.03565** — the abstract page does not carry a mobility.

**The datum is not in the full text either.** There is no numerical mobility
anywhere in the paper. Mobility appears only as curves against `|n|` in Fig. 5
and as the proportionality of their Eq. 4,

> μ = e·v_F²·τ(E_F)/E_F ∝ (1/√n)·τ(E_F)

with the only mobility statement in words being qualitative:

> "at high temperatures, mobility increases with decreasing T, but become
> temperature-independent below 100 K." — Section III c

So the item **as posed is unresolvable from this source**, and that is a
definite finding rather than another blocked fetch. The text is exhausted.
It should not be fetched again, and the 11-02 note that "a copy of either PDF
in the connected folder unblocks four items at once" is now wrong about this
one: a PDF of Zion would not help, because the number is not in the prose. It
would help with Fig. 5, which is a different blocker (§5 below).

**Fetch log, honestly.** One search, two fetches of the same ar5iv URL (the
second to pull verbatim wording for the quotes below). No failures, no
CAPTCHA. `arxiv.org` and `ar5iv.labs.arxiv.org` both behaved; the standing
note that PubMed/PMC return reCAPTCHA in this environment was not tested and
is neither confirmed nor refuted today.

## 2. What *is* there, and why it is strictly better than what was asked for

Zion estimate the charged-acceptor density in the SiO₂ next to their own
channel, from their own neutrality-point gate voltage:

> "Thus, we get C_ox = 1.15×10⁻⁴ F/m² and for V_NP ≈ 20 V, we obtain
> N_A^(S) ≈ 1.5×10¹² cm⁻²." — Section III a, point (i)

with `N_A^(S) = (1/e)·C_ox·V_NP`, `C_ox = ε₀ε/t`, and from Sections II and
III a: `t = 300 nm` (±5 %), `ε = 3.9`, `ε₀ = 8.85×10⁻¹² F/m`.

This is a **charged-impurity density for the same sample** whose crossover
density this repository has used as `N_CROSSOVER` since 10-30. So

&nbsp;&nbsp;&nbsp;&nbsp;`κ = n_c / N_imp`

is readable off **one** device, from **two** numbers, in **one** paper. It
needs no mobility, no Chen scaling constant, and no assumption that Zion's
device resembles this repository's — which is the whole of what made the
11-02 pairing the weakest link in that module.

## 3. The paper's own two numbers disagree by 4.3 %

Recomputing from the formula and the inputs Zion print:

| quantity | Zion print | recomputed from their own inputs |
|---|---|---|
| `C_ox = ε₀ε/t` | `1.15e-4 F/m²` | `1.1505e-4 F/m²` ✓ |
| `N_A^(S) = C_ox·V_NP/e` | `1.5e12 cm⁻²` | **`1.4362e12 cm⁻²`** ✗ |

The capacitance reproduces exactly at their stated precision, which is what
makes the second line attributable to the quoted density rather than to a
different oxide convention. `1.4362` does not round to `1.5` at any
precision; their "≈" is doing the work of a rounding-**up**.

**Which value the module uses.** The recomputed one — because it is
*smaller*, and a smaller `N_A^(S)` gives a **looser** bound on `κ` (§4). The
conservative choice is therefore also the one that cannot manufacture the
result. Both are reported, and the disagreement is pinned by a check (D2)
rather than silently resolved. This is the standing *one-way channel* item's
own prescription being followed: the externally-sourced claim and the
arithmetic that contradicts it are both recorded, and neither is quietly
dropped.

## 4. The direction of the inequality is the whole result

`N_A^(S)` is a **lower** bound on the density of charged scatterers, for
three independent reasons, all of them Zion's own:

1. **It is a net charge, not a count.** It is inferred from the voltage
   needed to compensate the channel. Scattering counts centres; with donors
   as well as acceptors present, the total exceeds the net.
2. **The oxide defects are acknowledged and unquantified.** "charged oxide
   defects with density N_t inevitable occur in SiO2 close (typically within
   1-3 nm) to the conductive channel" — Section III b. **No value for `N_t`
   is given anywhere in the paper.**
3. **One population among several.** Nothing in the estimate covers
   adsorbates or transfer residue on a CVD film grown on copper and
   transferred.

So `N_imp ≥ N_A^(S)`, and since `κ = n_c/N_imp` decreases in `N_imp`,

&nbsp;&nbsp;&nbsp;&nbsp;**`κ ≤ n_c/N_A^(S) = 2.089`** (recomputed `N_A`), or `≤ 2.000` at Zion's rounded figure.

The 11-02 value `κ = 2.4` lies **above** the conservative bound by 14.9 %. It
is not uncertain, it is excluded. The direction is the one the 11-02 module
had itself flagged as its risk — "if Zion's samples were cleaner than that,
`A_ref` is overestimated" — so this is a flagged failure mode occurring and
being quantified, not a new surprise.

The module's D11 mutates the inequality: read `N_A^(S)` as an *upper* bound
and the verdict reverses. The exclusion therefore rests on the three reasons
above and not on the arithmetic, which is identical either way.

**By-product.** `A = κ·C` and `A = n_c·μ_L` on Zion's sample, so the mobility
the paper never prints is bounded by its other two numbers:

&nbsp;&nbsp;&nbsp;&nbsp;`μ_{L,Zion} ≤ 2.089 × 5e15 / 3e16 = 0.348 m²/Vs = 3481 cm²/Vs`

3481 cm²/Vs is an ordinary long-range-limited mobility for transferred CVD
monolayer on 300 nm SiO₂, which is the only external check available on it,
and it passes. The repository had paired `gfet.mu = 0.4 m²/Vs = 4000 cm²/Vs`
onto Zion's device; that was 15 % too generous — the same 14.9 % from the
other side.

## 5. The finding that outranks the correction

Zion read their `3 × 10¹² cm⁻²` off the **hole** branch:

> "The dependence of the hole mobility in Fig. 5 clearly demonstrates the
> transition from long-range to short-range scattering" — Section III c

Two readings follow, and they point opposite ways.

**It is the right *kind* of crossover.** The repository's `n_c` is the
Matthiessen long-range/short-range crossover, which is exactly what Zion's
Fig. 5 transition is. Worth recording because it is easy to get backwards:
this number is **not** an electron–hole-asymmetry crossover, and nothing here
should cite it as one. (The paper does not state an asymmetry crossover
density at all; its abstract says only that "a significant asymmetry was
observed at low density of charge carriers".)

**But the repository applies it to both branches.** Written per branch,

&nbsp;&nbsp;&nbsp;&nbsp;`1/μ_b(n) = 1/μ_{L,b} + |n|/(n_{c,b}·μ̄_L)`, &nbsp;&nbsp; `n_{c,e/h} = n_c(1 ± s/2)`

there is an exact identity at `s = 2r₀`. The bracket factorises,

&nbsp;&nbsp;&nbsp;&nbsp;`1/μ_e = (1/μ̄_L)(1+x)/(1+r₀)`, &nbsp;&nbsp; `1/μ_h = (1/μ̄_L)(1+x)/(1−r₀)`, &nbsp;&nbsp; `x ≡ |n|/n_c`

the common `(1+x)` cancels out of `r = (μ_e−μ_h)/(μ_e+μ_h)`, and **`r(n) = r₀`
at every density, exactly, with no screening whatever.**

| branch split `s` | `r` at `|n|/n_c = 5.29` | reading |
|---|---|---|
| `0` | `0.0160` | the shipped answer: `r₀` screened 6.2× |
| `0.05` | `0.0370` | |
| `0.10` | `0.0580` | |
| `2r₀ = 0.20` | `0.1000` | `= r₀` exactly: **no screening at all** |

So the central qualitative conclusion of the 10-30 module — that density
screens the electron–hole asymmetry away — rests entirely on the
**neutral-defect channel being branch-symmetric while the charged channel is
not**. That was never argued. It was made by writing a single `n_c`.

And the symmetric case is not the safe default. `r₀` is derived from a
beyond-Born cross-section ratio, which is a property of scattering off a
potential of a given sign; vacancies and grain boundaries are also
potentials. There is no argument anywhere in this repository for the
asymmetry vanishing for neutral defects while surviving for charged ones.

**What would settle it, and it is one figure.** The electron branch's
transition density is in Zion Fig. 5, next to the hole one. Reading both
gives `s` for one sample directly. That is the **same picture** the standing
*lower bound on ℓ* and *`Δ_c` form for the chemisorbed metals* items are
already blocked on, so one figure read discharges three items — which makes
it the cheapest high-value unblocking action now available, displacing the
Zion full-text fetch that this note closes.

## 6. One check was wrong, and it failed the same day it was written

D3 was written as "the bound round-trips **bitwise**" and failed at 1 ulp on
first run. The failure is correct and the check was not: `(n_c/N_A)·N_A`
rounds twice, so exactness is not available in general. It was **split, not
loosened** — D3a asserts the 1 ulp floating point actually guarantees, D3b
recovers bitwise exactness where it genuinely holds (`N_A` a power of two
times `n_c`, where `κ` comes out exactly `4.0` and the round-trip is exact),
which isolates the double-rounding as the whole of D3a's 1 ulp.

This is the standing *audit every tolerance* item catching code written the
same day for the **fourth** consecutive session. The pattern across those
four is now consistent enough to state: the checks that fail on first run are
the ones predicting **exactness**, and they fail by one or two ulp in a
direction that turns out to be a theorem about rounding rather than a defect
in the model. That is an argument for writing them, not against.

## 7. Sources

- Zion *et al.*, "Charge carrier transport asymmetry in monolayer graphene",
  arXiv:1801.03565 — <https://arxiv.org/abs/1801.03565>, full text read via
  <https://ar5iv.labs.arxiv.org/html/1801.03565>. Sections II, III a, III b,
  III c and Fig. 5 as cited above.
- Chen, Jang, Adam, Fuhrer, Williams, Ishigami, "Charged impurity scattering
  in graphene", Nature Physics **4**, 377 (2008), arXiv:0708.2408 — the
  `μ·n_imp = C = 5×10¹⁵ V⁻¹s⁻¹` scaling law, used unchanged from 11-02.
- This repository: `graphene_density_dependent_mobility.py` (10-30),
  `graphene_impurity_limited_crossover.py` (11-02),
  `graphene_kappa_impurity_bound.py` (today),
  `notes/2026-11-02-the-crossover-is-an-impurity-density.md`.
