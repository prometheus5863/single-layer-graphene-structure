# The diffusion term was a boundary term, and the potential was never named

**Date:** 2026-10-07
**Modules:** `graphene_diffusion_current_model.py`,
`graphene_diffusion_mutation.py`
**Transcripts:** `diffusion_current_output.txt`, `diffusion_mutation_output.txt`
**Figure:** `diffusion_current_model.png`
**Closes:** the top open item of 2026-10-06 (*a diffusion term in Eq. (4)*) —
and moves its premise

---

## 1. The item

2026-10-06 created this item and named it that session's top one, on the
strength of a literature result:

> **A DIFFUSION TERM IN EQ. (4).** Feijoo *et al.* 2020 measure diffusion as
> comparable to drift at the peak-`f_max` bias, in devices at ~45 % of
> `v_sat`. Eq. (4) is drift-only, so **every `V_ds` ladder in this
> repository's RF thread is a drift-only slice** and the peak-`f_max` bias is
> a bias this model cannot be asked about.

Three parts of that are now settled, and they do not settle the same way.

| the item said | the measurement says |
|---|---|
| Eq. (4) is drift-only | **it is not clear that it is** — §4 |
| diffusion is comparable to drift | **±0.50 % of `I_d` here** — §3 |
| the peak-`f_max` bias cannot be asked about | **true, for a different reason** — §5 |

## 2. The derivation, and the thing it turned out to be

Zebrev's generalised Einstein relation
([arXiv:1102.2348](https://arxiv.org/pdf/1102.2348), Eq. 20) is
`mu = e*D/eps_D` with `eps_D = n/(dn/dmu_c)`, and for graphene's linear
dispersion `eps_D = E_F/2` exactly in the degenerate limit, so
`D = mu*E_F/(2e)`. With `V_F(n) = E_F/e = A_F*sqrt(n)` built from this
repository's own `v_F`,

```
D * dn/dx  =  mu * n * dV_F/dx         exactly, no approximation
```

so the drift-diffusion current is `W*e*mu*n * d(V - V_F)/dx` — transport down
the gradient of the quasi-Fermi potential, which is the textbook statement and
is also §4's problem. The `beta = 1` separability that makes Eq. (4)
closed-form **survives**: with `Phi = V - lambda*V_F`,

```
I_d  =  mu*W*e*Q_D / (L + mu*S_D),      lambda = 0 recovers Eq. (4) bitwise
Q_D  =  Q + lambda*(A_F/3)*(n_s^{3/2} - n_d^{3/2})
S_D  =  S + int kappa(V)/v_sat(V) dV,     kappa = -lambda*dV_F/dV
```

**`int n dV_F` has an exact antiderivative.** `int n d(sqrt n) =
(1/3)[n^{3/2}]`, so the diffusion contribution to the channel integral is a
difference of **endpoint** values: a pure **boundary term**, with no
dependence whatsoever on the density profile between source and drain.

That is a stronger statement than "it is small", and it is worth more than the
number it produces. 2026-10-05's criterion C found that this model's
*profile-domain* handling carried a one-signed error 674× larger than the
Jensen gap it was looking for. The diffusion term **cannot inherit that
error**, because it does not read the profile. Check **X6a** asserts it
structurally rather than arguing it: a deliberately lopsided interior
quadrature grid with the same endpoints leaves the term **bitwise** unchanged,
and **X6c** is the contrast — the quadrature route to the same number moves by
0.4 % on the same regridding.

## 3. The size, and a sign nobody had asked about

At the two fixed biases 2026-10-06 reported at, `V_ds = 0.1 V`:

| bias | `I_d` (λ=0) | `I_d` (λ=1) | share | `Qd/Q` |
|---|---|---|---|---|
| saturated-peak, `V_g = +2.857143 V` | 9.007791e-05 | 9.052996e-05 | **+0.5018 %** | +1.1621 % |
| resistor-peak, `V_g = -1.112782 V` | 8.910103e-05 | 8.865209e-05 | **−0.5038 %** | −1.1546 % |

**The correction is signed, and its sign flips between branches** (G1b). It
follows `d|n|/dV`, and this model carries `n` as a magnitude
(`carrier_density()` returns `sqrt(n_eff^2 + n_puddle^2)`), so on the hole
branch `|n|` *rises* toward the drain and the term subtracts. Whether that is
the correct physical sign for holes cannot be decided by a model that has
discarded the sign of the carrier — this repository has a signed-carrier
treatment for photodetectors (`graphene_photodetector_signed_carrier_model.py`)
and nothing equivalent for the FET. **Recorded as a new item, not as a
result.**

The share of `I_d` is about half of `Qd/Q`, and X9c says why: `S_d > 0`, so
Eq. (12) **opposes** Eq. (11). The two halves of the diffusion term pull in
opposite directions, and `S_d/S = +1.1621 %` is the same magnitude as `Qd/Q`
with the opposite effect. That was found by this module's own mutation
harness — see §7.

### Zebrev's closed form, as an independent check

Zebrev Eq. 62 gives the local ratio in closed form, `kappa = C_ox/(C_Q +
C_it)`. With `C_it = 0` and the dispersion-consistent
`C_Q = (2e^3/(pi hbar^2 v_F^2)) V_F`, that is capacitor algebra against a
boundary term — two computations sharing no machinery:

| | source end | drain end |
|---|---|---|
| `kappa`, this model | 0.011699 | 0.011538 |
| `C_ox/C_Q`, Zebrev | 0.016677 | 0.016877 |
| `C_ox/(C_Q+C_ox)` | 0.016403 | 0.016597 |

Interior mean ratio **0.6927**. They agree to within 31 %, which for two
routes this different is agreement; the residual is this model's series factor
and puddle floor, and it is reported rather than tuned away.

## 4. The premise: which potential is `V_ch`?

The standard GFET compact model this repository's Eq. (4) descends from —
Pasadas and Jiménez, *IEEE TED* 2016,
[arXiv:1605.08235](https://arxiv.org/pdf/1605.08235) — writes the channel
current as `I = -W*Q_tot*v` with `v = mu*F`, `F = -dV/dx`, and states
explicitly that **"V(x) is the quasi-Fermi level along the graphene
channel"**.

A drift expression whose driving potential is the quasi-Fermi level **is
already the complete drift-diffusion current.** Under that reading, the term
added here is a *double count*.

`gfet.carrier_density(V_g, V_ch)` never says which potential `V_ch` is. It
enters an electrostatic charge relation `C_ox*(V_g - V_ch - V_dirac)`, which
reads electrostatic; and it is then corrected by a quantum-capacitance series
factor `C_q/(C_q+C_ox)`, which is exactly the correction one applies when
`V_ch` is the quasi-Fermi level and the graphene drop `E_F/e` is a separate
voltage. The function contains evidence for both readings.

**The two readings differ by exactly `lambda*V_F`.** So the ±0.50 % above is
simultaneously (a) the diffusion current under the electrostatic reading and
(b) **the size of an ambiguity that was already present in every committed
`I_d` in this repository.**

This is the 2026-10-04 fault — a quantity computed correctly under a name that
does not pin down what it is — located in a *variable name* rather than in a
figure of merit. And the repair is not the one the item asked for: **it is to
name the potential, not to add a term.** Recommended and recorded: declare
`V_ch` the quasi-Fermi potential, keep `lambda = 0` as the production path,
and keep this module as the *bound* on the choice. `lambda = 1` is
deliberately not wired into any other module.

What the item got right is narrower than it looks and still matters: it is no
longer true that this repository's RF numbers are *drift-only* slices. They are
slices of a model whose driving potential was never named, bounded at ±0.50 %.

## 5. A prediction, written down first, and refuted

Both literature sources point the same way. Zebrev's `kappa = C_ox/C_Q`
**diverges** as `C_Q -> 0` at charge neutrality; Feijoo *et al.* put the
peak-`f_max` bias near the onset of bipolar conduction. The prediction was
written into Section 3 of the module before the sweep was run:

> the share **peaks at the Dirac point**.

Measured: it **collapses** there. `−0.006925 %` at the nearest grid point to
the Dirac point, against `0.5053 %` at the maximum — a factor of **73** — and
the maximum sits at `V_g = −1.272 V`, two volts away.

The mechanism is this model's own regularisation, and G4b measures it rather
than arguing it. `carrier_density()` returns `n = sqrt(n_eff^2 + n_puddle^2)`,
so `d|n|/dV = (n_eff/n) * dn_eff/dV` is **exactly 0.0** at `n_eff = 0`, while
`C_Q` stays finite at `n_puddle`. Measured at the Dirac point: `d|n|/dV =
0.000000e+00` with the floor, against `−2.290589e+15 m^-2/V` unfloored.

**The puddle floor flattens precisely the region where the capacitor ratio
diverges.** The diffusion-dominated regime is not small in this model — it is
*absent*, and that is a statement about the regularisation rather than about
graphene. 2026-10-06 was right that the peak-`f_max` bias cannot be asked
about. The reason is the puddle floor, not the missing term, and that reason
does not go away by adding terms.

## 6. A second finding, annotated rather than absorbed

`gfet.quantum_capacitance()` evaluates `C_Q` at the **gate overdrive**
`(V_g - V_ch - V_dirac)` in the slot the dispersion wants `E_F/e = V_F`. At
the RF bias those are 2.057 V and 0.0977 V: a ratio of **21.05**. So the
repository's own `C_q` overstates `C_Q` by that factor and would understate
`kappa` by it, which is why §3's comparison uses the dispersion-consistent
form.

No committed number is withdrawn. `C_q` is used downstream only through the
series factor `C_q/(C_q+C_ox)`, where `C_q >> C_ox` makes the factor close to
1 either way. Chapter 4 §4.6.5 carries the annotation in place.

## 7. The harness, and the survivor that was not a hole

7 of 7 after a first run of 6 of 7; control A green, controls B and C
surviving. Battery 26/26 → 30/30. The mutants attack the *structure* of the
claim, not only its arithmetic: **M6** replaces the closed form with the
quadrature route to the same number — a <0.1 % numerical change that destroys
the only claim in the module title — and X6a is the only check that reads it.

Two checks the harness forced:

- **X9** — λ enters in **two** places, and M2 arrives at only the first. It was
  killed only by G5c, a check written for something else; an incidental
  detector is not a detector. X9 gave Eq. (12) its own check, and in doing so
  found the `S_d > 0` fact that explains §3's otherwise unexplained factor of
  two.
- **G5d** — the stencil is now pinned against the **committed** 0.683186
  (reproduced to 1.010e-07), not only against `vsat` on the same stencil. G5a
  compares λ=0 against `vsat` *on the same stencil*, so a stencil change moves
  both together and G5a is blind to it by construction. That is 2026-10-06's
  M5 fault exactly.

**And the part worth the most.** The first run's survivor was the stencil
*spacing*, `6/399 -> 6/400`. The reflex is to file a survivor as a hole and
add a check. Measured before filing: the change moves `f_max/f_T` by
**−6.560e-11** relative and `f_T` by 5.275e-08, against a value committed to
six figures. A central difference is `O(h^2)` accurate precisely so that its
answer does not depend on `h`, so the mutant is **semantically inert at the
precision the claim is made at**, and a check sensitive to it would be a check
on an artefact. It is now **control C** — it must survive — and M7 was replaced
by a mutation that is not inert (the stencil *centre* becoming its one-sided
*edge*).

That is the verification repository's 2026-10-01 rule arriving here for the
first time: *a constant that can be mutated with no observable effect is not a
suite weakness.* The two repositories' methodological series have run
separately since 08-21. **This is the first rule to cross between them, and it
crossed because the same reflex produced the same wrong filing in both.**

## 8. The RF consequence

| quantity | λ = 0 | λ = 1 | change |
|---|---|---|---|
| `f_T` | 2.073963e+10 | 2.085305e+10 | **+0.5469 %** |
| `f_max` | 1.416902e+10 | 1.421070e+10 | **+0.2941 %** |
| `f_max/f_T` | 0.683186 | 0.681468 | **−0.2514 %** |
| `g_ds` | 3.374491e-02 | 3.391519e-02 | +0.5046 % |

`f_max` rises and the ratio **falls**. That is 2026-10-06's finding — a
criterion written on a ratio is the wrong instrument when the mechanism scales
both arguments — in a **third** mechanism, and this time the two move in
*opposite* directions rather than merely by different amounts. §7.8.1a's
structural verdict survives: 0.681468 is 52.42 % of 1.3 either way.

## 9. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a way
that cannot drift. 10-02: when it agrees, that is a fact about two artefacts.
10-03: a PASS/FAIL at zero is silent about magnitude. 10-04: and a quantity can
be computed correctly under the wrong name. 10-05: and a criterion with a
number in it can be missed in a way that locates the real constraint. 10-06:
and a criterion can be decided correctly on the quantity it names and still
miss the finding, if that quantity is a ratio.

**10-07: AND AN ITEM CAN BE WELL-POSED, CORRECTLY MOTIVATED FROM THE
LITERATURE, AND BUILT EXACTLY AS WRITTEN, WHILE ITS PREMISE — THAT THE MODEL
LACKS THE TERM — IS NOT A FACT ABOUT THE MODEL BUT AN ARTEFACT OF A VARIABLE
WHOSE MEANING WAS NEVER FIXED.**

The checkable rule: before adding a term to a transport equation, ask what the
*potential* in it is. A drift expression in the quasi-Fermi potential and a
drift-diffusion expression in the electrostatic potential are the same
equation, and no amount of care about the term can distinguish them. The
mechanisable form, and the new standing item: **census every potential-like
variable in this repository and record which potential it is.** `V_ch` was the
first one asked and the answer was "both".

---

### Fetch failures, recorded honestly

- `link.springer.com` — **HTTP 429 (rate limited)**, with an explicit
  instruction not to retry, on *Drift-diffusion models for the simulation of a
  graphene field effect transistor* (Mathematics in Industry, 2022,
  DOI 10.1186/s13362-022-00120-3). Not retried. **Publishers are now the
  fourth distinct source of fetch refusals** in this repository and the first
  429 rather than a 403.
- PubMed/PMC — **not attempted**, per the standing reCAPTCHA note. One result
  (*Analytical Solution for the Potential Distribution in the Channel of a
  GFET*, PMC13173535) was visible in search results and skipped for that
  reason; it is plausibly relevant to §4 and is left for a session with
  another route to it.
- `researchgate.net` / `academia.edu` copies of Jiménez and Moldovan's 2011
  explicit drain-current model were not fetched; the 2016 Pasadas–Jiménez
  arXiv preprint carries the statement §4 needs, from the same group, and is a
  primary source rather than a mirror.

### A near-miss in this session's own method, recorded because it is the genre

An attempt to verify that `graphene_velocity_saturation_model.py` reproduces
its committed transcript on the **device VM** (Python 3.10.12, numpy 2.2.6)
was cut off by the device shell's 180 s ceiling at 175 s. `diff -q` against
the committed transcript then reported **no differences** — because the module
writes its transcript only at the end, so the file still held the committed
copy the run had never overwritten. **A diff against a file the run did not
write is not a reproduction check**, and it was one command away from being
logged as a PASS. The claim that is actually supported, and the one the
portability commit makes, is byte-identity on **Python 3.13.16** in the cloud
container, where the run completed in 2 m 11 s. The device VM's real
constraint is a *timeout*, not a `SyntaxError`, and it is narrower: this
session's new module runs there in 14.5 s.
