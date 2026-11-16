# The negatives were never sound, and v_F was never one name

*2026-10-25. Discharges the 2026-10-24 ranked-first item — audit the capture
surface of the other three replacements — and contradicts the sentence that
item was written under.*

---

## 1. The claim under test

2026-10-24 built `graphene_capture_surface_guard.py`, which adjudicates a
capture by **counterfactual** rather than by a read set, and which needs a
way to tell `PROTECTED` (the capture caused the bitwise zero) from
`UNREACHED` (the holder never looks the name up, so the zero is nobody's
doing). Its Section 0 grounds that split on one sentence:

> an over-approximate read set has reliable **negatives**, so *not in the
> static read surface* is trustworthy even though *in it* is not.

The sentence is plausible, it is the standard property of a sound
over-approximation, and in this repository it is **false**.

## 2. Why it is false

`rra.static_reads` walks bytecode and records a `(module, attr)` pair in
exactly two situations:

1. a module object is loaded and an attribute is then taken from it —
   `qca.A_FERMI`, `gfet.C_ox`;
2. a bare global is loaded and it resolves to a **repo function** —
   `cq_thermal`, `series_factor`.

A bare global that resolves to a module-level **scalar**, read from inside
the very module that defines it, matches neither branch and is recorded
**nowhere**. That is how every constant in this repository is read by its own
module's functions. So a negative does not mean *not read*; it means *not
read in a spelling the walker can see*.

The clean demonstration is one name read two ways:

| holder | on the old surface? | leverage of `qca.A_FERMI` |
|---|---|---|
| `qca._carrier_density_consistent` | **no** | +0.0121 % |
| `qca._quantum_capacitance_consistent` | **no** | **+29.8900 %** |
| `blast._cq_consistent_contact` | yes | +58.3417 % |
| `flip.make_contact_cq(live=False)` | yes | +58.3417 % |

One name, one value, one mechanism. `qca` reads it bare; `blast` spells it
`qca.A_FERMI`. **Its visibility to every audit in this repository is decided
by the spelling at the read site**, and nothing else. This is the sixteenth
running instance of the standing *resolve-on-spelling* item and the first
with the sign reversed: previously a test resolved on spelling, here the
**audited surface** does.

## 3. The standing question, answered — opposite ways

The item asked whether `gfet._thermal_rounding` and `qca.A_FERMI` carry
2026-10-24's ≈ +58 % onto the other three replacements, *"which would make
them repository-wide levers rather than a contact fact"*.

- **`gfet._thermal_rounding` is a contact fact.** Bitwise **exactly zero**
  on both channel replacements. Note the epistemics: that zero is a
  *measured counterfactual*, not the name's absence from a surface — after
  §2 an absence would prove nothing.
- **`qca.A_FERMI` is repository-wide**, and off-surface in both places where
  it matters.

So the hypothesis is half right, and the half that is right is right for a
reason the question did not anticipate.

## 4. What else was hiding there

Completing the surface on the four holders adds **18 (holder, name) pairs**,
8 distinct names, and **every one of them moves the value**:

| name | charge repl. | capacitance repl. | contact repl. |
|---|---|---|---|
| `qca.A_FERMI` | +0.0121 % | +29.8900 % | on old surface |
| `qca.CQ_PREFACTOR` | +0.0227 % | **+50.4145 %** | exact 0 |
| `qca.E_CHARGE` | −0.5361 % | −25.3101 % | −21.2518 % |
| `qca.KB` | +0.0173 % | +27.6751 % | — |
| `gfet.e` | — | — | **+278.4911 %** |
| `gfet.hbar` | — | — | −63.5273 % |
| `gfet.v_F` | exact 0 | exact 0 | −63.5273 % |
| `gfet.kB` | — | — | +0.0242 % |

Not one had ever been adjudicated, because the guard's missing-name check can
only demand a verdict for a name its surface *names*. The transcript records
the consequence directly: the 2026-10-24 guard **accepts** an adjudication of
all three replacements with every one of these omitted.

**Repository-wide, 75 of the 90 module-level scalars (84 %) are in the blind
class** — read bare inside their own module, and therefore invisible to
`rra.static_reads`, to the 2026-10-23 write/read guard, to
`bc.identity_class` and to the capture surface alike. The census's precision
is pinned from the other side too: `gfet.t_ox` and `gfet.epsilon_0` are
scalars that are *not* in the class, because they are read only at import
time to build `C_ox` and never inside a function body.

## 5. What this supersedes, and what it does not

2026-10-24's §4.6.20 claimed a **fully adjudicated** capture surface. That
claim is false. Its **verdicts are not**: all nine are re-measured here and
all nine still hold. What was missing was five further names on the same
holder, the largest at +278.4911 %. §4.6.20 (d) and (e) are annotated in
place; the numbers stand, the completeness claim and the sensitivity ranking
in (e) do not.

Concretely, (e) said *"the three names that do carry the contact C_q are
`gfet._thermal_rounding` and `qca.A_FERMI` at ≈ +58 %, and `cd.C_interface`
itself at +18 %, which is the sensitivity ranking a reader of §4.5 should
have."* The actual ranking is headed by `gfet.e` at +278.4911 %, with
`gfet.hbar` and `gfet.v_F` at −63.5273 % and `qca.E_CHARGE` at −21.2518 %
interleaved among the three names it lists.

## 6. The device-physics consequence: v_F is not one name

`C_q` carries the Fermi velocity twice, as
`CQ_PREFACTOR = 2e³/(π ħ² v_F²)` and as `A_FERMI = ħ v_F √π / e`, so a
consistent change of `v_F` by *k* must move one name as *k*⁻² and the other
as *k*⁺¹. It does not, because the two families reach it through two
different names:

| holder | `gfet.v_F` | `qca.CQ_PREFACTOR` | `qca.A_FERMI` |
|---|---|---|---|
| `qca._carrier_density_consistent` | **exact 0** | +0.000671 | +0.000247 |
| `qca._quantum_capacitance_consistent` | **exact 0** | +1.008232 | +0.519524 |
| `blast._cq_consistent_contact` | −2.323463 | **exact 0** | +1.161680 |

*(d ln value / d ln name; `exact 0` = a ×1.5 rebind leaves the value bitwise
unchanged.)*

The contact replacement reads `gfet`'s **live** prefactor and so follows
`gfet.v_F`. Both channel replacements read `qca.CQ_PREFACTOR` and
`qca.A_FERMI`, which `graphene_exact_charge_rederivation` computed from
`gfet.v_F` **at import time** and copied in. So `gfet.v_F` is bitwise
unreachable from every Chapter 4 channel number, and `qca.CQ_PREFACTOR` is
bitwise unreachable from the contact number.

**A v_F sensitivity study done by editing `gfet.v_F` would have reported
channel insensitivity to the Fermi velocity as a result.** This is the
2026-10-14 delivery fault in a new place. There the stale binding was a name
that had not been rebound; here it is a name that was **precomputed**, which
no audit in this repository looks for, and which `bc.identity_class` cannot
see because the two values are not the same object — they are a float and a
float derived from it.

### The perturbation turns out to be the physical case

`v_F` is not a constant of nature; substrate screening renormalises it.
Hwang *et al.* measure, by ARPES, 1.15 ± 0.02, 1.49 ± 0.08 and 2.49 ± 0.30
(×10⁶ m/s) for graphene on SiC(000-1) (ε = 7.26), BN (ε = 4.22) and quartz
(ε = 1.80), against a fully screened *ab initio* value of 0.85 × 10⁶ m/s,
and report a 1/ε dependence at fixed carrier concentration.

The committed `v_F = 1.00 × 10⁶` sits just above the **bare** value. A 90 nm
SiO₂ back gate has ε ≈ 3.9, between the BN and quartz cases, so the
substrate-appropriate value is near 1.5 × 10⁶ — which is **exactly the ×1.5
the probe applies by default**. The −63.5273 % it reports for the contact
`C_q` is therefore not a synthetic perturbation: it is the SiO₂-screened
Fermi velocity, and **Chapter 4 quotes the contact `C_q` at a v_F its own
substrate does not have.**

This is a real open question about what the thesis should claim, not a bug.
Using the bare value consistently is defensible (it is what most device
models do, and the renormalisation is carrier-density dependent, so a single
number is already an approximation); quoting it as *the* contact `C_q` while
the same chapter models a specific SiO₂ stack is less so. What is not
defensible is the present state, where the question cannot even be *asked*
by changing one name.

## 7. Validated against an exactly known symmetry

The whole of §6 rests on a probe, so the probe is anchored to a value known
in advance: `gfet.hbar` and `gfet.v_F` must have **bitwise identical**
leverage, because the prefactor depends on them only through `(ħ v_F)²`.
They do — both move the contact `C_q` to `0x1.87bc5bc0c26d0p-6`. A probe
reaching the value by any other route would not reproduce that coincidence.

This is the fourth session running in which an exact check, rather than a
plausible-range check, is what made a result trustworthy.

## 8. A defect in this session's own checks, caught by this session's own mutation test

The first draft of check V10c read
`(blind - b_ex) == set(_PLANTED_APPARATUS)`. Since `b_ex` is
`blind - _PLANTED_APPARATUS` by construction, that is **true for any value
of the exclusion** — a tautology checking the exclusion against itself.
Mutation M7 widens the exclusion to swallow `qca.A_FERMI` and
`qca.CQ_PREFACTOR`, i.e. hides the day's two findings behind the apparatus
carve-out the 2026-10-17 standing item warned about, and the first draft
scored **31/31** against it.

V10c now compares the exclusion with the census's own result for this
module, and V10f requires every excluded name to belong to this module *as
named by `__file__`*. Neither can be satisfied by editing the exclusion. The
unmutated suite is 36/36 and all 9 injected defects are detected.

One refusal branch has **no** control and the transcript says so: `EXPOSED`
on an added name whose rebind leaves the value bitwise unchanged. No added
pair has exact-zero leverage, and the only candidates are the planted
controls, whose sites `csg.leverage` cannot resolve. Sixth instance of the
standing *every verdict needs a planted control* item, in this session's own
code.

## 9. The lesson, in the form the others take

- §4.6.13 — a repair to one binding is not a repair to the class.
- §4.6.14 — a class can span two operating points, so it must be partitioned.
- §4.6.17 — a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- §4.6.18 — a capture that fixes one name does not fix the names the captured
  object reads.
- §4.6.19 — a refusal that nothing reaches is not a refusal.
- §4.6.20 — a name outside every rebind set is not thereby harmless, only
  unwatched.
- **§4.6.21 — a name's absence from a read surface is a fact about the
  surface's vocabulary, not about the name. Two names for one parameter are
  not one knob, and a precomputed copy is unreachable from the name it was
  computed from.**

## Sources

- C. Hwang, D. A. Siegel, S.-K. Mo, W. Regan, A. Ismach, Y. Zhang, A. Zettl,
  A. Lanzara, *Fermi velocity engineering in graphene by substrate
  modification*, Sci. Rep. **2**, 590 (2012),
  <https://doi.org/10.1038/srep00590> —
  <https://www.nature.com/articles/srep00590>

### Fetch record

Live search was available. One search and one fetch, both successful; no
retries and no blocked hosts this session. PubMed/PMC were not attempted, in
line with the standing note that they return reCAPTCHA in this environment
(three logged occurrences).
