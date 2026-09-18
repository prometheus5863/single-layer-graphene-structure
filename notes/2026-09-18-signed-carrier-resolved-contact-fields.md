# Signed, carrier-resolved contact fields: the crossover is not at graphene's work function

**Date:** 2026-09-18
**Closes:** the top "not yet covered" item of 2026-09-17 — *"signed,
carrier-resolved two-contact treatment ... the most consequential open
item in Chapter 6"* — and Simplification 1 of
`graphene_photodetector_two_contact_model.py`.
**Code:** `graphene_photodetector_signed_carrier_model.py`
**Figure:** `photodetector_signed_carrier_response.png`

---

## 1. What the magnitude convention was hiding

The 2026-09-17 two-contact model wrote both contact fields with
`|W_metal - W_graphene|`. It flagged that as a simplification. Working
through it properly shows it was concealing **two** separate physical
facts, not one.

### 1.1 The sign

A contact that *n*-dopes graphene and one that *p*-dopes it build
electrostatic fields that point the **same way** along the channel, not
opposite ways. Forcing both to a magnitude makes every contact behave
like a p-type one, so in the old model every pair partially cancelled
and every pair developed an interior stagnation point. That is an
artefact of the convention, not physics.

### 1.2 The crossover — the part that was not anticipated

The sign of the contact doping is **not** `sign(W_metal - W_graphene)`.
Graphene's own work function is 4.5 eV, but the n/p crossover for a
metal *on* graphene sits near **5.4 eV**:

> Giovannetti, Khomyakov, Brocks, Karpan, van den Brink & Kelly,
> "Doping Graphene with Metal Contacts", *Phys. Rev. Lett.* **101**,
> 026803 (2008). <https://arxiv.org/abs/0802.2267>
> Crossover at "a metal work function of ~5.4 eV"; the ~0.9 eV offset
> from graphene's own 4.5 eV is the short-range metal–graphene chemical
> interaction adding an interface dipole on top of vacuum-level
> alignment.

`graphene_contact_doping_model.py`'s `METAL_WORK_FUNCTIONS` table
**already knew this** — its inline comments call Cu an "n-type dopant"
even though 4.65 > 4.5, and Pt "clearly above the ~5.4 eV crossover ->
p-type". The magnitude convention simply could not express it, because
`|W - 4.5|` discards the crossover entirely.

Under the physical crossover, **six of the seven metals in the table
n-dope graphene and only Pt p-dopes it**:

| metal | W (eV) | ΔW vs 4.5 | type | ΔW vs 5.4 | type |
|-------|--------|-----------|------|-----------|------|
| Ti | 4.33 | −0.17 | n | −1.07 | n |
| Cr | 4.50 | +0.00 | n | −0.90 | n |
| Cu | 4.65 | +0.15 | **p** | −0.75 | n |
| Ni | 5.04 | +0.54 | **p** | −0.36 | n |
| Au | 5.10 | +0.60 | **p** | −0.30 | n |
| Pd | 5.12 | +0.62 | **p** | −0.28 | n |
| Pt | 5.65 | +1.15 | p | +0.25 | p |

Pt is the only available partner for a genuine n/p pair. Because this
changes which metals are which *type*, the model carries the crossover
as an explicit parameter (`w_cross`) and every result below is reported
under both values.

---

## 2. The model

Contact A at *x*=0, contact B at *x*=L, signed offset
ΔW = W_metal − w_cross. The Dirac point in the channel is rigidly
shifted near each contact with the same saturating profile the rest of
the repo uses (λ = 250 nm):

```
E_D(x) = ΔW_A/(1 + x/λ) + ΔW_B/(1 + (L−x)/λ)
```

A rigid band shift is an electron potential energy, so
E(x) = (1/e) dE_D/dx, and both carriers drift in that single field in
opposite directions: v_h = +μE, v_e = −μE. Each is transported by the
same stagnation-aware logic as 2026-09-17. The figure of merit is the
net charge delivered to contact A per absorbed photon:

```
N = (1/L) ∫ [ (+1)·s_h(x)·p_h(x) + (−1)·s_e(x)·p_e(x) ] dx
```

with s = +1 collected at A, −1 at B, 0 stalled. **N now runs over
[−2, +2], not [−1, +1]** — both carriers can be collected, which the
single-species magnitude convention structurally could not represent.

---

## 3. Validation: three exactly-known values

The 2026-09-17 session recorded the lesson that an exactly-known
validation beats a plausible-range one; it was an exact zero that
exposed an `inf − inf` bug a range check had waved through. All three
checks here are against exact values.

| # | Check | Worst deviation |
|---|-------|-----------------|
| 1 | Hole-only + \|ΔW\| + w_cross = 4.5 eV must reproduce the 2026-09-17 model on all 49 ordered pairs at both biases (98 comparisons) | **0.000e+00** (bitwise) |
| 2 | Identical contacts at zero bias → N = 0 exactly (antisymmetry; Weiss & Duan's "net zero photocurrent") | 6.6e−17 |
| 3 | Charge conjugation: at zero bias, flipping the sign of both offsets must give N(−ΔW) = −N(ΔW) exactly | **0.000e+00** |

Check 1 matters because it establishes the two modules are the *same
expression* under the old restrictions — algebraically E(x) = −F_old,
so "hole moves toward A" is exactly "F_old > 0". It is a real
reduction, not an approximate agreement. Check 3 is the one that only
becomes statable once the model is signed: the magnitude convention
could not even express it, and a model that mixed carriers wrongly or
leaked a magnitude anywhere would fail it while still passing check 2.

---

## 4. Result 2 (asymmetric pairs) — confirmed, and robust

**This confirms the specific prediction 2026-09-17 made about its own
open item.** That entry said the signed treatment "could invert Result
2's ordering by letting an n/p pair such as Ti/Pt add rather than
partially cancel". It does:

| | best zero-bias pair | \|N\| |
|---|---|---|
| magnitude convention (2026-09-17) | **Cr/Pt** | 0.917 |
| signed model, w_cross = 5.4 eV | **Ti/Pt** | 1.832 |
| signed model, w_cross = 4.5 eV | **Ti/Pt** | 1.830 |

The ordering inverts (Ti/Pt overtakes Cr/Pt) and the magnitude roughly
**doubles**. Crucially the result is **insensitive to the crossover
choice** — 1.832 vs 1.830 — because Ti and Pt sit on opposite sides of
*both* candidate crossovers. This is the strongest result of the run.

### Mechanism (stagnation audit, zero bias)

| pair | N | E(x) nulls | holes stalled | electrons stalled | h→A | h→B | e→A | e→B |
|------|------|-----|------|------|------|------|------|------|
| Pt/Pt (p/p) | +0.000 | 2 | 0.00 | **1.00** | 0.357 | 0.357 | 0 | 0 |
| Pd/Pd (n/n) | +0.000 | 2 | **1.00** | 0.00 | 0 | 0 | 0.368 | 0.368 |
| Ti/Pd (n/n) | −1.599 | **0** | 0.00 | 0.00 | 0 | 0.720 | 0.879 | 0 |
| Ti/Pt (n/p) | −1.832 | **0** | 0.00 | 0.00 | 0 | 0.907 | 0.925 | 0 |
| Cr/Pt (n/p) | −1.810 | **0** | 0.00 | 0.00 | 0 | 0.896 | 0.914 | 0 |

The interior null **disappears entirely** for an unequal pair, and both
carrier species are collected, at opposite ends. In a same-type
symmetric pair one whole carrier species is 100% stalled at the null
while the other splits evenly and cancels — which is why N = 0 there
for two independent reasons at once.

Note Ti/Pd reaches 1.599 despite both metals being n-type under the
physical crossover: what matters is the *difference* in ΔW, not the
sign pair. A large-ΔW / small-ΔW same-type pair is already most of the
way to an n/p pair. That is the practically useful version of the
statement, since Pt is the only p-type metal in the table.

### Experimental anchor

This is not hypothetical. Mueller, Xia & Avouris, *Nature Photonics*
**4**, 297 (2010) — already cited in this repo for τ = 1 ps — built
exactly this device: "One set of fingers was made of palladium/gold
(20/25 nm in thickness), and the other of titanium/gold (20/25 nm)",
precisely because "if both electrodes ... consist of the same metal,
the built-in electric field profile in the channel between two
neighbouring fingers is symmetric, and the total photocurrent is zero."
Asymmetric metallisation gave **6.1 mA/W at 1.55 μm, a 15-fold
improvement**. The Ti/Pd pair they chose is the second-best pair in the
table above (1.599), and the best that avoids Pt.

Supporting sign/band-bending picture, for the contact-doping profile
itself: Mueller, Xia, Freitag, Tsang & Avouris,
*"The role of contacts in graphene transistors: a scanning photocurrent
study"* (arXiv:0902.1479) — measured potential step ≈ **0.12 eV** at
the graphene/electrode interface, with charge-transfer doping extending
**0.2–0.3 μm** into the channel, consistent with the λ = 250 nm used
here; and explicitly "a p–n junction forms close to the
electrode/graphene interface" depending on gate bias.

---

## 5. Result 1 (symmetric ranking) — NOT robust; 2026-09-17's Result 1 must be downgraded

This is the uncomfortable half, and it is recorded rather than
smoothed over.

Net response of a **symmetric** device at V_bias = 0.1 V:

| metal | magnitude (2026-09-17) | signed, w_cross = 5.4 | signed, w_cross = 4.5 |
|-------|------------------------|------------------------|------------------------|
| Ti | 0.672 | 0.180 | 0.987 |
| Cr | 0.632 | 0.210 | **1.264** |
| Cu | 0.673 | 0.245 | 1.101 |
| Ni | 0.321 | 0.436 | 0.321 |
| Au | 0.295 | 0.495 | 0.295 |
| Pd | 0.287 | 0.518 | 0.287 |
| Pt | **0.169 (worst)** | **0.556 (best)** | **0.169 (worst)** |

Three different answers for the same question:

* 2026-09-07 (single contact): **Pt best**, 1.47× enhancement.
* 2026-09-17 (two contact, magnitude): **Pt worst**, and the whole
  ranking reversed.
* 2026-09-18 (signed, physical crossover): **Pt best again** — but
  **Pt worst** if the crossover is placed at 4.5 eV instead.

The honest conclusion is not "Pt is best after all". It is that
**the symmetric-device metal ranking is not an established result of
this repo**: it flips with a modelling choice (where the n/p crossover
sits) that the 2026-09-07 and 2026-09-17 models never made explicit
because they used a magnitude. 2026-09-17's Result 1 should be read as
*crossover-dependent*, not as a reversal that supersedes 2026-09-07.

The structure behind the 5.4 eV column is clean and worth keeping
independently of the ranking: under the physical crossover the signed
response is **monotone in |ΔW|** — Ti (|ΔW| = 1.07) worst, Pt
(|ΔW| = 0.25) best. For a *symmetric* device the contact doping field
is purely parasitic: it is the same at both ends, so it contributes
nothing to the net current and only creates a null that strands
carriers. The winner is whichever metal is **closest to the
crossover**, i.e. dopes graphene least. That statement is
crossover-independent in form, and it is also why Cr wins the 4.5 eV
column (ΔW = 0 exactly there).

**Design rule that survives both columns:** for a symmetric
two-terminal device, pick the metal that perturbs graphene least; for
an asymmetric device, maximise |ΔW_A − ΔW_B|. These are opposite
instructions, and the asymmetric device wins outright — 1.83 vs at best
0.56.

---

## 6. What this does and does not claim

* It does **not** claim the ΔW → doping-profile relation is linear. It
  is only roughly linear in Giovannetti et al., and not linear at all
  for the chemisorbed metals (Ti, Ni, Pd). Only the sign structure and
  the crossover were changed here; the profile magnitude convention is
  inherited unchanged from 2026-08-26.
* Uniform illumination, drift only, no diffusion, μ_e = μ_h, no
  photogain / photo-thermoelectric / bolometric contribution.
* The 0.12 eV measured potential step (arXiv:0902.1479) is well below
  the 0.25–1.07 eV offsets the table assumes, because the measurement
  is of the *residual* step in a gated device rather than the flat-band
  charge transfer. Reconciling the two is not attempted here and is a
  genuine open item.

## 7. Fetch log (honest record)

* `arxiv.org/abs/0802.2267` (Giovannetti et al.) — **OK**. Note: the
  fetch summary returned the n/p labels *inverted* (it said metals
  above 5.4 eV n-dope). The correct physics is the opposite — a
  high-work-function metal withdraws electrons and p-dopes graphene —
  and the crossover value 5.4 eV was the number actually used. Recorded
  because an uncorrected paraphrase would have flipped every sign in
  this module.
* `arxiv.org/pdf/0902.1479` (Mueller et al., scanning photocurrent) — **OK**.
* Mueller, Xia & Avouris 2010 review PDF (bilkent mirror) — **OK**.
* PubMed/PMC — **not attempted**. Three consecutive reCAPTCHA blocks
  (2026-09-05, 09-07, 09-17) established this as a standing environment
  limitation, so the session did not spend time on it.
