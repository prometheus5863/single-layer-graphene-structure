# Density-dependent electron/hole mobility, and why the gate stack decides whether the asymmetry matters

2026-10-30. Research notes behind `graphene_density_dependent_mobility.py`,
which closes the item ranked second on 2026-10-29.

## 1. What the literature actually supplies, and what it does not

Zion, Butenko, Sharoni, Kaveh and Shlimak, *Charge carrier transport asymmetry
in monolayer graphene*, [arXiv:1801.03565](https://arxiv.org/html/1801.03565v1)
— the paper the 10-29 module took mechanism (A) from — does **not** supply a
mobility function. It supplies three things, and the distinction matters because
the module built on it has to be honest about which of its content is measured
and which is interpolation:

1. **The two Boltzmann limits, with their density exponents.** Their eq. (4) is
   `mu = (e v_F^2 / E_F) tau(E_F)`, and since `E_F ~ sqrt(n)`:
   - long-range (charged impurities, screened Coulomb): `tau ~ n^(+1/2)`, so
     **`mu ~ const`**;
   - short-range (neutral defects — vacancies, dislocations):
     `tau ~ n^(-1/2)`, so **`mu ~ 1/n`**.
   The *opposite signs* are what make this a crossover rather than two parallel
   curves, and they are the only input the model in this repository needs.
2. **One crossover density.** The hole data "clearly demonstrate the transition
   from long-range to short-range scattering" at `p ≈ 3 × 10^12 cm^-2`. This is
   `N_CROSSOVER` in the module, converted once to `3e16 m^-2`.
3. **A direction and a decay, with no numbers.** `mu_h < mu_e` at low `|n|`, and
   "the asymmetry decreases with increasing `|n|`". **No numerical `mu_e/mu_h`
   is given in the text** — the mobility values appear only in their Fig. 5,
   which is a picture. This is the same obstruction already logged for the
   `ℓ` lower bound and the `Δ_c` form: the datum is a figure, not a number.

Supporting parameters they do give: acceptor density
`N_A^(S) ≈ 1.5 × 10^12 cm^-2`, from `N_A^(S) = C_ox V_NP / e` with `V_NP ≈ 20 V`
and `C_ox = 1.15 × 10^-4 F/m^2`; puddle density `⟨n_D⟩ = 0.25 × 10^12 cm^-2`;
potential amplitude `Δ ≈ 58 meV` from `Δ ≈ ħ v_F sqrt(π ⟨n_D⟩)` with
`v_F = 10^6 m/s`, against `Δ ≈ 54 meV` from their Ref. 16.

**Mechanism, and why the asymmetry belongs to one channel only.** They
attribute the asymmetry to the relativistic (beyond-Born) scattering
cross-section: massless Dirac fermions attracted to a charged impurity scatter
more strongly than those repelled, so negatively charged acceptors scatter holes
more strongly than electrons. That argument is *about a charged scatterer*.
Neutral defects have no sign for it to act on, so the short-range channel is
electron–hole symmetric **by construction, not by assumption** — this is the one
modelling decision in `(B1)` that the literature determines rather than leaves
open.

## 2. The closed form

Matthiessen on the two channels, asymmetry confined to the long-range one:

```
1/mu_b(n) = 1/mu_L,b + |n| / (n_c * mu_L_bar),    b in {e, h}        (B1)
```

Writing `a = 1/mu_L,e`, `b = 1/mu_L,h`, `s = |n|/(n_c mu_L_bar)`, the common
denominator `(a+s)(b+s)` cancels **exactly** between numerator and denominator of
`r = (mu_e-mu_h)/(mu_e+mu_h)`, leaving `r = (b-a)/(a+b+2s)`. Since `2/(a+b)` is
the *harmonic* mean `mu_H` and `mu_H/mu_bar = 1-r^2` identically for any pair —
the same identity that made the 10-29 conductivity floor a geometric mean —

```
r(n) = r0 / (1 + (1 - r0^2) * |n| / n_c)                             (B3)
```

The asymmetry decays **hyperbolically** in `|n|`, with scale `n_c/(1-r0^2)`, and
with no free parameter beyond the two the model already had. Checked against the
species computation in D1 to 0.75 eps absolute.

`(1-r0^2)` is 0.98 at the asymmetries of interest, which is exactly why D12
mutation-tests the closed form by dropping it: a factor of 0.98 is invisible to
any plausible-range check, and the mutant exceeds D1's budget by 3e13.

## 3. The device statement, and the independent corroboration I did not expect

`(B3)` turns `n_c` into a crossover gate overdrive, `dV_c = e n_c / C_ox`:

| gate stack | `C_ox` (F/m²) | `dV` at `n_c` | reached in operation? |
|---|---|---|---|
| 90 nm SiO₂ back gate (shipped) | 3.84e-4 | **12.53 V** | no — every sweep here is inside ±5 V |
| 10 nm HfO₂ top gate | 2.21e-2 | **0.217 V** | yes, essentially always |

> **Annotated 2026-11-02.** Both rows of this table use the *same* `n_c`, and
> §4 below already says they must not. `graphene_impurity_limited_crossover.py`
> now derives `n_c(N_imp)` and the top-gate row's honest overdrive is
> **0.478–0.912 V**, not 0.217 V — the uncorrected number is low by a factor of
> **2.2–4.2**. The back-gate row is unaffected: its dielectric is grown before
> the channel sits on it, so it adds no impurities the reference mobility does
> not already contain. **The numbers above are left as they were measured.**

So whether electron–hole mobility asymmetry is a first-order effect or a
negligible one is **a property of the gate stack, not of the graphene**. It is
large in the thick-oxide back-gated devices most transport papers measure and
small in the thin-oxide top-gated devices a semiconductor process would build.

> **ANNOTATION (2026-11-07): the second half of that sentence is withdrawn,
> and the open item at the foot of this note is what withdrew it.** The table
> above uses one `n_c` for both rows. Strand 4 of the Chapter 7 synthesis —
> `n_c = κ·N_imp` (11-02) with `N_imp ≥ N_A^(S) = C_ox·V_NP/e` (11-03) — makes
> `n_c` scale with the same `C_ox` that sets `|n|`, and the two cancel
> **exactly**: `ΔV_c = κ·f·V_NP`, independent of oxide thickness and
> permittivity. The two stacks then need the *same* overdrive, measured as a
> difference of exactly `0.0` V, and the 57.69× spread above is exactly the
> inverse `C_ox` ratio — a statement about the oxide, not about the graphene.
>
> With `κ ≤ 2.0889` and `V_NP = 20 V`, `κ·V_NP = 41.78 V`, which **neither**
> stack reaches. At 1 V of overdrive on the 10 nm HfO₂ top gate the surviving
> asymmetry is `r = 0.0977` — 97.7 % of the unscreened `r₀ = 0.10` — against
> the `0.0180` the fixed-`n_c` reading gives, a factor of 192 in `x`. So the
> asymmetry is **large in both** stacks, not large in one and small in the
> other. "A property of the gate stack" survives only in the sense that the
> stack enters through `V_NP`, how much net charge its deposition leaves
> behind, and not through `C_ox` at all.
>
> **What this does to the arXiv:1010.0913 corroboration below.** That paper's
> HfO₂-topped devices show electron and hole mobilities differing by less than
> 10 %, which was read here as independent support for screening in a
> thin-oxide stack. Under the corrected reading that support is **withdrawn as
> corroboration of this mechanism**: if `ΔV_c` is ~42 V regardless of stack,
> those devices are not past their crossover either, so their small asymmetry
> needs a different explanation — most plausibly a smaller `r₀`, i.e. a
> contact or interface difference, rather than screening. The observation
> stands; its attribution does not. This is logged as an open item rather than
> resolved here.
>
> Source: `graphene_screening_overdrive_invariance.py`, 13/13 green, with the
> cancellation checked bitwise over six decades of `C_ox` and four decades of
> `κ`, the predecessor reproduced exactly in its own regime (E4a), the
> overlapping limit bitwise at Zion's own `C_ox` (E5), and a mutation that
> breaks the cancellation and turns both invariance checks red (E7).
> Figure: `screening_overdrive_invariance.png`, rendered and inspected.
> Synthesis: Chapter 7 §7.3a.

The corroboration: looking for the mobility of HfO₂-topped graphene I found
*Dielectric Thickness Dependence of Carrier Mobility in Graphene with HfO₂ Top
Dielectric*, [arXiv:1010.0913](https://arxiv.org/pdf/1010.0913), and it reports,
in passing and for an entirely different purpose, that in its HfO₂-topped
devices **"the electron and hole mobilities differ by less than 10%"**. That is a
prediction of the paragraph above, made by a paper that was not asking the
question. It is weak evidence — one sentence, no density resolution, mobilities
extracted from a σ-vs-V_BG slope averaged over a 25 V window with the Dirac
region *excluded* — but it is evidence, and it points the right way.

## 4. The complication the same paper creates, which the module does not yet handle

That paper's actual finding cuts against the simple reading of the table above.
HfO₂ deposition **adds charged impurities**: it raises the charged-impurity
concentration by `1.5–4 × 10^12 cm^-2` and costs "about 50% of the initial
value" of the mobility, dropping steeply over the first 2–4 nm and then
plateauing. (Their own onset is stated inconsistently — "first 1-2 nm" in one
place, "first 2-4 nm" in another; recorded as stated rather than reconciled.)

`1.5–4 × 10^12 cm^-2` is **the same order as `n_c` itself**, so `n_c` is not a
constant of the material across the two rows of that table. In `(B1)`,
`A = n_c * mu_L_bar` depends only on the *neutral* defect density, so
`n_c = A / mu_L_bar`: halving `mu_L_bar` by adding charged impurities **doubles**
`n_c`. The honest version of the top-gate row is therefore `n_c ≈ 6e12 cm^-2`
and `dV_c ≈ 0.43 V`, not 0.217 V.

**The conclusion survives the correction, and the margin is worth writing down
rather than waving at.** A 10 nm HfO₂ top gate at 1 V reaches
`n = C_ox·1/e = 1.38 × 10^13 cm^-2`, which is 2.3× the corrected `n_c` and 4.6×
the uncorrected one. So the device is past the crossover either way. What the
correction removes is the *comfort* of the factor: the margin is 2.3, not 46, and
a thicker or dirtier dielectric could eat it.

> **Annotated 2026-11-02 — the hand calculation above is reproduced exactly and
> then narrowed.** `graphene_impurity_limited_crossover.py` D9 reproduces
> `n_c = 6.000e12 cm^-2`, `dV_c = 0.4343 V` and `M = 2.303` from this
> paragraph's own stated input (the mobility halved), so nothing here was
> arithmetically wrong. What was wrong was the *choice of input*: this
> paragraph used one of the two numbers arXiv:1010.0913 reports (the mobility
> loss) where the other one (the added impurity density, `1.5–4e12 cm^-2`) is
> the more direct route. Taking that route, **`n_c` = 6.6–12.6e12 cm^-2 and the
> margin is 1.10–2.09×, not 2.3×.** At the dirty end the device clears the
> crossover by **10 %**.
>
> So "a thicker or dirtier dielectric could eat it" was right, and it is nearer
> to eaten than this paragraph supposed. The falsifiable form is the design rule
> `N_imp* = C·C_ox·V/(e·A)`: this stack at 1 V tolerates `4.51e12 cm^-2` total,
> i.e. **12.7 % of headroom** above the paper's dirty end.
>
> The two routes disagree because **the paper's two numbers disagree with each
> other** once Chen's `C` is fixed — its ~50 % mobility loss implies
> `ΔN_imp = 1.25e12 cm^-2`, below even the bottom of its own `1.5–4e12` range
> (D10). That is a property of the source, not an error in this paragraph.

This is a new open item, not something the committed module handles: `n_c` is a
single module-level constant there, and it should be a function of the
charged-impurity density that the gate stack determines. Logged.

## 5. Fetch record, honestly

- `arxiv.org/html/1801.03565v1` — fetched cleanly, cloud-side. Consistent with
  10-29's finding that the `arxiv.org` block was the device VM's egress path,
  not the domain.
- `arxiv.org/pdf/1010.0913` — fetched cleanly, cloud-side.
- Zion *et al.* Fig. 5 (the `mu_e`, `mu_h` values) and Huard *et al.* Fig. 4
  (the `ℓ` bound, the `Δ_c` form) remain unavailable for the same reason as
  before: **the datum is a picture**, and no fetch of any kind supplies it. A
  copy of either PDF in the connected folder unblocks all three.
- PubMed/PMC were **not attempted**: four logged occurrences of reCAPTCHA pages
  make it a known dead end, and the 10-29 entry asks that the session not be
  spent re-confirming it.

## Sources

- [Zion, Butenko, Sharoni, Kaveh, Shlimak, *Charge carrier transport asymmetry in monolayer graphene*, arXiv:1801.03565](https://arxiv.org/html/1801.03565v1)
- [*Dielectric Thickness Dependence of Carrier Mobility in Graphene with HfO₂ Top Dielectric*, arXiv:1010.0913](https://arxiv.org/pdf/1010.0913)
- [Huard, Stander, Sulpizio, Goldhaber-Gordon, *Evidence of the role of contacts on the observed electron-hole asymmetry in graphene*, arXiv:0804.2040](https://arxiv.org/abs/0804.2040) (mechanism (B), carried from 10-29)
