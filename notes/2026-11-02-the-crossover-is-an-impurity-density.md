# The crossover density is a charged-impurity density in disguise

*2026-11-02. Closes the item ranked first on 2026-10-30. Sources: Chen et al.,
Nature Physics **4**, 377 (2008),
[arXiv:0708.2408](https://arxiv.org/pdf/0708.2408); Zion et al.,
[arXiv:1801.03565](https://arxiv.org/abs/1801.03565); Jang et al.,
[arXiv:1010.0913](https://arxiv.org/pdf/1010.0913).*

## 1. The thing that was wrong, which the 10-30 entry had already said in prose

`graphene_density_dependent_mobility.py` carried

```python
N_CROSSOVER = 3e16              # m^-2  (= 3e12 cm^-2)
```

as a module-level constant and used it for both rows of its gate-stack table.
The 10-30 entry's §4 identified the problem in words and computed the
consequence by hand in a notes file:

> In `(B1)`, `A = n_c * mu_L_bar` depends only on the *neutral* defect density,
> so halving `mu_L_bar` by adding charged impurities **doubles** `n_c`.

That is right, and it is the whole argument. The short-range (neutral-defect)
term of `(B1)` is written with two factors, `n_c` and `mu_L_bar`, and only one
of them is a neutral-defect property. The other is the long-range mobility,
which charged impurities certainly change. So the quantity that is actually the
material constant is the **product**

```
A ≡ n_c · mu_L_bar        [V^-1 s^-1]                                      (C1)
```

and `n_c` alone is not a constant of anything. What this session did was move
that from prose into the code, with checks, and then find that the hand
calculation's *input* was the wrong one of the two the source offers.

## 2. One piece of literature closes it

Chen, Jang, Adam, Fuhrer, Williams and Ishigami dose graphene with potassium
*in situ* and measure mobility against **deliberately added** charged-impurity
density — so `n_imp` is a controlled variable, not a fitted one, which is the
strongest form this input could take. Their scaling law:

```
mu · n_imp = C,     C = 5 × 10^15 V^-1 s^-1                                (C2)
```

with the dosing rate measured independently at `(2.6–3.2) × 10^15 m^-2 s^-1`
and `C` cross-checked against the shift in `V_g,min`. Eliminating `mu_L_bar`
between (C1) and (C2):

```
n_c(N_imp) = (A/C) · N_imp ≡ κ · N_imp,     κ = 2.40 (dimensionless)        (C3)
```

The crossover density is a **fixed multiple of the charged-impurity density**,
and the multiple is a ratio of two literature constants — one from Zion (the
crossover), one from Chen (the scaling law). No fitted parameter enters.

**The caveat that is load-bearing, stated here rather than buried.** Chen's
measured mobility is the *total* mobility, not the long-range channel in
isolation. Using (C2) for `mu_L` is legitimate exactly insofar as their dosed
samples are charged-impurity-limited, which is their central claim, but it means
`N_imp` inferred as `C/mu_total` is overstated by `(1 + |n|/n_c)` at whatever
density they measured, and **they do not state that density**. D7 measures the
size of the effect at this repository's own operating points: the long-range
channel carries **86 %** of `1/mu` at the puddle density and **71 %** at 5 V on
the back gate, so (C2) is being used inside its regime there. At the top-gate
operating point it carries only **32 %** — which is not a defect but a
restatement of `M > 1`, i.e. that device operates past the crossover, where
neutral defects dominate by construction.

## 3. Two consequences worth stating before any number

Both are checked rather than asserted, and neither is obvious from `(B1)` as it
was written.

1. **`r0` is invariant under `N_imp`.** Scaling the charged-impurity density
   scales `mu_L,e` and `mu_L,h` by the *same* factor, because the beyond-Born
   asymmetry is a ratio of per-scatterer cross-sections and carries no density.
   So the *amplitude* of the asymmetry does not move at all. D2: invariant over
   seven decades, to 1.4e-16.
2. **A dirtier device is MORE asymmetric.** More charged impurities push `n_c`
   up, so at fixed operating density the device sits further *below* the
   crossover and `r(n)` sits closer to `r0`. This is the right sign and it is a
   consistency check on the whole decomposition rather than a new result: the
   asymmetry lives in the charged-impurity channel, so more charged impurities
   must mean more asymmetry. D4-i: strictly monotonic over eight decades.

## 4. The reference point is a *pairing*, and that is the weakest link

```
A_ref     = n_c,Zion · mu_this-repo = 3e16 · 0.4   = 1.20e16 V^-1 s^-1
N_imp,ref = C / mu_this-repo        = 5e15 / 0.4   = 1.25e16 m^-2 (1.25e12 cm^-2)
κ         = A_ref / C                              = 2.40
```

`A` is **not measured anywhere.** It pairs Zion's crossover density, measured on
*their* CVD samples, with *this* repository's long-range mobility. If Zion's
samples were cleaner than `mu = 0.4 m^2/Vs`, `A_ref` is overestimated and so is
`n_c` at every impurity density. `N_imp,ref = 1.25e12 cm^-2` is an ordinary
value for graphene on SiO₂, and that is the only external check available on it.

**An attempt to resolve the pairing today reached only the abstract.** The
arXiv abstract page for 1801.03565 carries the asymmetry statements and the
T and B ranges but no mobility, substrate or impurity number. Recorded as a
fetch outcome, not worked around: the pairing stands unresolved and is the
first open item.

What survives the pairing and what does not, **measured as an invariance rather
than argued** (D5):

| quantity | carries the pairing? |
|---|---|
| the scaling `n_c ∝ N_imp` (C3) | **no** — follows from (C1)+(C2) alone |
| every ratio: `n_c(dirty)/n_c(clean)`, margin ratios, fractional change in `r` | **no** — `A_ref` cancels, bitwise under a ×4 perturbation |
| absolute `n_c`, absolute `dV_c`, absolute margin | **yes** — all move by exactly the perturbation |

A factor of **4** is used for the perturbation rather than a round 3, so that it
is exactly representable and the invariance can be tested bitwise instead of to
a tolerance.

## 5. The correction of record, and why the 10-30 arithmetic was not wrong

D9 reproduces the 10-30 hand calculation from **its own stated input** before
correcting it — `n_c = 6.000e12 cm^-2`, `dV_c = 0.4343 V`, `M = 2.303`. A
correction that cannot first reproduce the thing it corrects is not a
correction, and this one reproduces it to three decimals.

So the hand arithmetic was right. What was wrong was the **choice of input**:
arXiv:1010.0913 reports *two* numbers about HfO₂ deposition, and 10-30 used the
mobility loss where the added-impurity density is the more direct route.

| channel | `N_imp` (e12 cm⁻²) | `n_c` (e12 cm⁻²) | `dV_c` (V) | `M` at 1 V |
|---|---|---|---|---|
| reference, no deposition | 1.250 | 3.000 | 0.2171 | 4.605 |
| +1.5e12 (paper, clean end) | 2.750 | 6.600 | 0.4777 | **2.093** |
| mobility-halved (the 10-30 route) | 2.500 | 6.000 | 0.4343 | 2.303 |
| +4.0e12 (paper, dirty end) | 5.250 | 12.600 | 0.9120 | **1.096** |

**The conclusion survives and the margin does not.** The top-gated device is
still past the crossover across the whole range, but 10-30 reported 2.3× and the
honest range is **1.10–2.09×**. At the dirty end the device clears the crossover
by **10 %**, which is not a margin a process would design to.

## 6. The two numbers in the source disagree with each other

This is why the two routes differ, and it is a property of the source rather
than of either calculation. Given `C`, the paper's "~50 % of the mobility"
implies `ΔN_imp = 1.25e12 cm^-2` — **below even the bottom of its own stated
`1.5–4e12` range**, by a factor of 1.2, and below its top by a factor of 3.2.
So the paper's own mobility datum selects the *clean* end of the paper's own
impurity range (D10).

**The escape is real and is not assumed away.** If HfO₂ deposition also adds
*neutral* defects, then `A` grows too, the mobility loss is shared between the
two channels, and a larger `ΔN_imp` becomes consistent. That would change `A`,
which this module holds fixed, and it is logged as an open item rather than
modelled. The high end of the range is therefore *disfavoured*, not excluded.

## 7. The device statement, which is a design rule

Setting `M = 1` gives the dirtiest channel a gate stack can drive past the
crossover at a given overdrive:

```
N_imp* = C · C_ox · V / (e · A)
```

| gate stack and overdrive | `N_imp*` (e12 cm⁻²) | vs. the reference channel |
|---|---|---|
| 10 nm HfO₂ at 1 V | 5.757 | 4.61× dirtier allowed |
| 90 nm SiO₂ at 5 V | 0.499 | **0.40× — would need to be 2.5× cleaner** |

Two statements come out of this, and the second one is the counterintuitive one.

**(i) The HfO₂ deposition budget is real but thin.** That stack at 1 V tolerates
`ΔN_imp ≤ 4.51e12 cm^-2` of deposition damage against the paper's reported
maximum of `4.0e12` — **12.7 % of headroom**. This is the falsifiable form of
10-30's "a thicker or dirtier dielectric could eat it". Since
`N_imp* ∝ C_ox·V ∝ V/t_ox`, the budget tightens linearly with oxide thickness at
fixed drive voltage, and a process that thickens the dielectric to improve
leakage spends it directly.

**(ii) Cleaning the channel moves the crossover DOWN.** A cleaner channel has a
*lower* `n_c`, so it reaches the asymmetry-killing regime at *lower* gate
voltage. The back-gated device cannot reach `n_c` at 5 V not because it is too
dirty but because it is too *clean* for its own weak gate — it would need to be
2.5× cleaner still for `n_c` to fall to where 5 V can reach it. This inverts the
usual intuition that a cleaner device is the more ideal one in every respect:
for electron–hole symmetry, material quality and gate strength act in the same
direction, and the thick-oxide device is short of both.

## 8. What the checks caught, and one thing only looking caught

Twenty-one checks, 21 passing, three of them against values known exactly in
advance:

- **D1** requires the new model to *be* the 10-30 model at the reference point,
  not to approximate it: 0 ulp on `n_c`, 2 ulp on the mobility pair, 0 ulp on
  the conductivity over a 401-point ±5 V sweep. The asymmetry parameter gets an
  **absolute** eps budget rather than a relative one, for the reason 10-30's own
  D1 established — `r` is a difference of nearly equal mobilities over their
  sum, so cancellation amplifies relative error by `1/|r|` without bound.
- **D3** checks the high-density asymptote against `A·X/(1+X)`, an *exact*
  consequence of `(B1)` rather than a limit taken numerically, and the predicted
  value contains no `N_imp` at all — which is the sharpest available statement
  of (C1). Worst 1 ulp over seven decades.
- **D6** is a true bitwise symmetry: `r(2n, 2N_imp) == r(n, N_imp)` to the last
  bit, because scaling by a power of two is exact and commutes with rounding at
  every step. It says the model has exactly one dimensionless group.
- **D11** reinstates the pre-correction model — `n_c` held fixed while `mu_L`
  scales — as a *control*, and requires D3 to fail against it. It fails by
  3e16 ulp away from the reference point and passes at the reference point,
  where the two models coincide by construction.

**D4-ii was written as a single tolerance and it failed, by 83× its counted
budget.** The budget was not miscounted. At `N_imp = 1e12 × ref` the group
`n/n_c` is still `4.6e-12`, so the shortfall from `r0` is a real first-order
term and a tolerance loose enough to admit it would have been a number chosen to
pass. Split in two: the limit is now taken where it genuinely *is* the limit
(`1e20 × ref`, where the short-range term is `4.6e-20` of `1/mu_L` and cannot
move the float at all), and the finite-`N_imp` residual gets a **prediction** —
`r0·u/(1+u)` with `u = (1-r0²)n/n_c` — which accounts for it to `9.7e-5`
relative. Same lesson as 10-30's D4, reached the same way and on the same day
the check was written. Item 15 again.

**And one defect that no check in this module could have found.** The figure was
verified as a well-formed 1950×712 PNG and by its underlying numbers, and it
passed both — but rendering it and *looking* at it showed two colliding labels
in panel (a): the "after HfO₂ deposition" annotation sat underneath the legend
box, and the reference-pairing callout overlapped the green dashed line's label.
Both fixed and re-inspected. This is the first figure in the standing ~66-PNG
backlog to be actually looked at, and it took one render to find two real
defects, which is an argument for the item rather than against it.

## 9. Fetch record, honestly

- `arxiv.org/pdf/0708.2408` — **fetched cleanly**, and it is the paper this
  session needed: `C = 5e15 V^-1 s^-1` stated as the relation tested, with the
  dosing rate and the `V_g,min` cross-check.
- `arxiv.org/abs/1801.03565` — **abstract only.** No mobility, substrate or
  impurity number, so the §4 pairing is unresolved. The full text, or Fig. 5,
  would close it.
- No PubMed/PMC attempted: three logged reCAPTCHA occurrences make it a known
  dead end in this environment.
