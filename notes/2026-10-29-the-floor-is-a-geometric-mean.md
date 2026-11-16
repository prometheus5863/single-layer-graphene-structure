# The floor is a geometric mean, and the Dirac point is not a doping measurement

2026-10-29. Session subject: the **ASYMMETRIC ELECTRON/HOLE MOBILITY** item,
created 2026-10-09 and carried untouched for twenty sessions. It is the only
item on the standing list that was deferred *in writing inside the source* —
`graphene_fet_signed_carrier_model.sheet_conductivity_two_carrier` added the
`mu_e`/`mu_h` hooks, declined to use them, and said so — so the item came with
its own interface already designed. Closing it took one new module and two
in-place annotations, and it produced a result that is a correction of record
rather than an addition.

---

## 1. Why the two branches are not equivalent, and why the literature disagrees

Within the Born approximation, Coulomb scattering of a Dirac fermion from a
charged impurity goes as the *square* of the impurity charge, so the impurity's
sign drops out and electrons and holes scatter identically. Any measured
electron–hole asymmetry therefore has to come from outside leading-order Born,
and two mechanisms compete in the literature. They are not reconciled, and
they do not even agree on which branch is worse.

| | channel origin | contact origin |
|---|---|---|
| source | Zion, Butenko, Sharoni, Kaveh, Shlimak, [arXiv:1801.03565](https://arxiv.org/abs/1801.03565) | Huard, Stander, Sulpizio, Goldhaber-Gordon, [arXiv:0804.2040](https://arxiv.org/abs/0804.2040) |
| geometry | Hall bar, non-invasive probes 30 µm from the contacts | invasive probes crossing the strip, vs. external side arms |
| finding | `μ_h < μ_e` at low density; asymmetry **shrinks** as `\|n\|` grows | `ρ_e/ρ_h = 1 ± 0.03` away from the metal, **≈ 1.20** through it |
| mechanism | beyond-Born (Novikov) relativistic scattering off **negatively charged acceptors**, ≈ 1.5 × 10¹² cm⁻² in the SiO₂, puddle density 0.25 × 10¹² cm⁻², potential amplitude ≈ 58 meV | the metal pins the graphene beneath it **p-type** (Pd ≈ 5.1 eV vs. graphene ≈ 4.5 eV), so the junction is p–n for an n-type channel and p–p′ for a p-type one, and the two junction resistances differ |
| scaling | with the channel | with the **metal/graphene interface width** `w`, not with channel length |
| which branch is worse | holes | electrons (for p-doping contacts) |

Two numbers from Huard et al. land directly on this repository's own
calibration. Their specific contact resistance is `λ = 135 Ω·µm`, which sits at
the bottom of the 110–500 Ω·µm band `graphene_fet_model.Rc_total` is calibrated
against — so the contact term this repository already ships is in the range
where their junction asymmetry is a measurable fraction of it. And their
mobility, 4600 cm²/V·s, is within 15 % of the shipped `μ = 4000 cm²/V·s`.

Zion et al.'s central observation is one this session's model **cannot**
represent: their asymmetry is density-dependent and dies off as the long-range
potential is screened (hole crossover at `p ≈ 3 × 10¹² cm⁻²`). A constant
`(μ_e, μ_h)` pair has no such behaviour. That is recorded as a follow-on item
rather than guessed at, and it is stated in the new module's own docstring
under "what this module does not claim" — because a constant pair that is
silently used where a density-dependent one belongs is exactly the kind of
scope error §3 below is about.

**Live search note.** Both papers were reached and read this session. PubMed
and PMC were not attempted, per the standing note (four logged reCAPTCHA
occurrences). `arxiv.org/pdf/0804.2040` fetched cleanly, and that is worth
recording precisely, because two *other* standing items — a lower bound on `ℓ`
and a `Δ_c` form for the chemisorbed metals — have been logged since 2026-10-18
as "blocked on Fig. 4, which is blocked on `arxiv.org`, refused by the egress
proxy".

What is now established is narrower than "arxiv.org is blocked" and narrower
than "arxiv.org is open":

- **`arxiv.org` PDF *text* is reachable** this session, through the web-fetch
  tool, which runs cloud-side. The 10-18 refusal was the **device VM's** egress
  proxy, which is a different network path with a possibly narrower allowlist.
  Those two items' wording conflates the two paths and should name which one.
- **Fig. 4 is still blocked, for a different reason.** A fetch returns text;
  Fig. 4 is a plotted curve, and no amount of text extraction reads a value off
  a curve. So the items do *not* become actionable today — but the reason they
  are blocked has changed from "the domain is refused" to "the datum is a
  picture", and the remedy named in the item (a copy of the PDF in the
  connected folder, where a page can be rendered and looked at) is the right
  remedy for the *second* reason and would not have helped with the first.

Corrected in today's log entry rather than left as it stands.

---

## 2. The derivation, which is three lines and gives a geometric mean

Mobility is a transport parameter; it does not enter the charge relation. So
keep the electrostatics exactly as shipped and put the asymmetry only in the
two-carrier Drude sum, on the signed decomposition the 2026-10-07 module
supplied:

```
  n_e, n_h = (S ± n_net)/2,          S = sqrt(n_net² + n_p²)
  σ(n_net) = e(μ_e n_e + μ_h n_h)
           = e(μ̄ S + (Δμ/2) n_net),  μ̄ = (μ_e+μ_h)/2,  Δμ = μ_e-μ_h   (A3)
```

`S` is **even** in `n_net` and `n_net` is **odd** in `dV`, so (A3) is a clean
split: the even part is the conductivity this repository already had, and the
odd part is new and linear in the mobility difference. That split is the whole
derivation. Minimising:

```
  dσ/dn_net = e(μ̄ n_net/S + Δμ/2) = 0
    ⇒  n_net/S = -r,         r ≡ Δμ/(2μ̄) = (μ_e-μ_h)/(μ_e+μ_h)
    ⇒  n_net* = -n_p r/sqrt(1-r²)                                      (A4)
    ⇒  σ_min  = e n_p μ̄ sqrt(1-r²) = e n_p sqrt(μ_e μ_h)               (A5)
```

The last step is an identity, not an approximation: `1-r² =
4μ_eμ_h/(μ_e+μ_h)²`. **The minimum sheet conductivity of a disordered graphene
channel is set by the geometric mean of the two branch mobilities.**

I did not expect the geometric mean and went looking for the arithmetic one; it
is worth saying that the reason it is geometric is visible in (A4) rather than
in (A5). The minimum does not sit at charge neutrality. It sits displaced onto
the **low-mobility** branch, at the density where the extra carriers of the
slow species exactly pay for the lost carriers of the fast one, and evaluating
(A3) at that displaced point is what converts an arithmetic mean into a
geometric one.

---

## 3. The correction of record, and why it moves no number

`graphene_fet_signed_carrier_model.sigma_min_ch4` returns `e·μ·n_puddle` and
its docstring calls that **EXACTLY**. On the `μ_e = μ_h` diagonal — where the
shipped model lives and where that function evaluates — it is exact, and (A5)
reproduces its value to **1 ulp**, checked. Off the diagonal the sentence is
wrong in two distinguishable ways:

1. `e·μ·n_puddle` **overestimates** the floor, by `μ̄/sqrt(μ_eμ_h)`.
2. The docstring's "at `dV = 0`" is no longer *where the minimum is*. At
   `dV = 0` the function returns a point **on** the curve which coincides with
   the curve's minimum only when `r = 0`.

So nothing computed changes and the Chapter 3 / Chapter 4 floor adjudication of
2026-10-09 stands as written. What changes is the **scope** of a word. Both
docstrings are annotated in place, with their numbers left where they are, per
the standing rule.

This is the fourth session in a row whose finding was a claim that was true of
the configuration it was evaluated in and false as stated. The pattern is now
specific enough to name: **a claim evaluated on a diagonal, asserted off it.**
`sigma_min_ch4` evaluates at `μ_e = μ_h`; `C7b`'s verdict (2026-10-26)
evaluated at one `gfet.C_ox`; the §7b census (2026-10-28) evaluated on a tree
with its own module skipped. In all three the *evaluation* was right and the
*quantifier* was missing.

---

## 4. The result that matters for a device, and the one that does not

The two orders are different, and that is the useful statement:

```
  apparent Dirac shift :  n_net* ∝ r          FIRST order
  floor change         :  sqrt(1-r²) = 1 - r²/2 + O(r⁴)   SECOND order
```

Inverting (A4) through this repository's own series-factor charge relation,
at the shipped 90 nm SiO₂ back gate (`C_ox = 3.8368 × 10⁻⁴ F/m²`):

| `μ_e/μ_h` | `r` | apparent `ΔV_dirac` | as fraction of the shipped `V_dirac = 0.80 V` | misattributed density |
|---|---|---|---|---|
| 1.05 | 0.0244 | −0.0524 V | 0.065 | 1.25 × 10¹⁰ cm⁻² |
| 1.10 | 0.0476 | −0.1011 V | 0.126 | 2.42 × 10¹⁰ cm⁻² |
| 1.20 | 0.0909 | −0.1922 V | 0.240 | 4.60 × 10¹⁰ cm⁻² |
| 1.30 | 0.1304 | −0.2763 V | **0.345** | 6.62 × 10¹⁰ cm⁻² |
| 1.50 | 0.2000 | −0.4278 V | 0.535 | 1.02 × 10¹¹ cm⁻² |

A 30 % mobility asymmetry — squarely inside the range Huard et al. report as
`ρ_e/ρ_h ≈ 1.20` and Zion et al. report qualitatively — walks the resistance
maximum **0.276 V** to lower `V_g`, which is **34.5 %** of the 0.80 V this
repository attributes *entirely* to built-in doping. The sign is informative
and falls out of (A4) without a sweep: `μ_e > μ_h` puts the minimum on the hole
side, i.e. the apparent Dirac point moves **down** in gate voltage.

The practical consequence is a measurement caveat, not a model change: **reading
a doping level off the position of a transfer-curve minimum assumes `r = 0`,
and nothing in a two-terminal transfer curve tests that assumption.** The same
split of orders explains why the floor literature looks so much cleaner than
the doping literature: at `μ_e/μ_h = 1.3` the floor moves −0.88 % and at 1.5 it
moves −2.02 %, both inside anybody's error bar, while the Dirac point moves by
a third of a volt.

One knock-on, quantified rather than hand-waved:
`n_puddle_matching_ch3()` divides Chapter 3's `4q_e²/h` floor by `e·μ` to get
the puddle density that would reproduce it. Under asymmetry the divisor is
`e·sqrt(μ_eμ_h)`, so the extracted puddle density is **larger** by
`μ̄/sqrt(μ_eμ_h)` — +0.88 % at `μ_e/μ_h = 1.3`. Second order, as everything
about the floor is.

---

## 5. Telling the two mechanisms apart, exactly

The literature conflict in §1 is not settled here, and it is not the kind of
thing this repository can settle. What *can* be done is to show the two
mechanisms are distinguishable in principle, and to show it exactly rather than
by eye. The odd part of the two-terminal resistance is

```
  R_odd(L) ≡ R(-dV) - R(+dV) = (L/W)(1/σ(-dV) - 1/σ(+dV)) + (Rc_p - Rc_n)
```

which is **exactly linear in L**. The slope carries the channel asymmetry; the
intercept carries the contact asymmetry. Neither can be absorbed into the
other at any length, because the contact term enters with no factor of `L` at
all. This is the same structure Huard et al. report from the measurement side:
their odd resistance collapses when scaled by interface width, and a bulk-only
reading of the same data gives a resistivity asymmetry that *drifts with probe
spacing* — which is what a nonzero intercept misread as a slope looks like.

Because the relation is exactly linear, a planted truth must come back out, and
that is the strongest check in the new module (`C7`). Planting
`(μ_e/μ_h = 1.3, ΔRc = 120 Ω)` and solving two points of a 25-length sweep
recovers the slope to **2.2 × 10⁻¹⁶ relative** and the intercept to the printed
digit, with residuals at **1.4 × 10⁻¹² Ω**. Each of the three scenarios —
channel only, contact only, both — is adjudicated correctly from the fitted
pair alone.

The honest cost is in the same section. At the shipped `L = 200 nm`:

```
  channel-only R_odd =  80.93 Ω
  contact-only R_odd = 120.00 Ω      ratio 1.48
```

The two mechanisms sit within a factor 1.48 of each other, so **a
single-length measurement cannot distinguish them at all**, and the shipped
device is a single length. The discriminator is a statement about an
experiment this repository does not model, which is why it is reported as a
method and not as a result.

---

## 6. The checks, and the one that was wrong first

43 of 43 green. Five are exact-value checks rather than plausible-range
checks, which is the standing requirement since the 2026-09-17 bug that a
range check passed:

| check | what it pins | measured |
|---|---|---|
| `C2` | `σ(+dV) - σ(-dV)` at `μ_e = μ_h` | **bitwise 0.0** |
| `C3` | `n_e·n_h = n_p²/4` for *every* mobility pair — catches a mobility leaking into the electrostatics | 5.40 × 10⁻¹⁶ relative |
| `C4` | analytic `dσ/dn` at the analytic `n*` | ≤ 4.4 × 10⁻³⁶ against a 3.2 × 10⁻⁴ scale |
| `C5` | exchanging `μ_e ↔ μ_h` flips `n*` and fixes `σ_min` | **both bitwise** |
| `C7` | planted `(μ_e, μ_h, Rc_n, Rc_p)` recovered from a length sweep | 2.2 × 10⁻¹⁶ relative |

`C1` — the symmetric limit reproducing the shipped sheet conductivity — was
written with a **1 ulp** threshold and **failed at 3**. The threshold was
wrong, not the code: 1 ulp was a number I picked because it looked strict, and
it is not a bound derived from anything. Counting instead: the shipped
`abs(n)*e*mu` rounds twice after `S`; the species-split form rounds five times
and the even/odd form four; seven roundings of `eps/2` bound the disagreement
at ~8 ulp. The budget is now counted and the measurement (3 ulp, 2.00 eps)
quoted beside it.

Loosening a failing threshold is the obvious way to fake this, so `C1` gained a
second check that a bound alone cannot provide. A rounding floor does not care
how finely the bias sweep is sampled; a small *model* difference, evaluated at
new bias points, generally would. Refining the sweep 10× (601 → 6001 points)
leaves the deviation at 6.661 × 10⁻¹⁶, inside the counted budget. That is the
check that distinguishes "my threshold was wrong" from "my model is wrong", and
it is the one that should have been there before the threshold was chosen.
