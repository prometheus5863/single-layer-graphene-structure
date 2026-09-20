# The dW -> doping relation is not linear, and it is least linear exactly where this thesis leans on it

**Date:** 2026-09-20
**Open item closed (partially):** "A non-linear dW -> doping-profile relation",
the top open item in Chapter 6 since 2026-09-19.
**Status:** literature settled; model implemented in
`graphene_contact_doping_nonlinear_model.py`; consequences for Chapter 6
computed and reported in the same file.

---

## 1. What every model in this repo currently assumes

Three modules take the contact-induced doping profile to be **linear in the
work-function offset**, with a single decay length for every metal:

- `graphene_contact_doping_model.py` -- `lambda_decay = 250 nm` for all metals
- `graphene_photodetector_signed_carrier_model.py` -- `total_field()` builds
  the field as `dW / lam / (1 + x/lam)**2`, i.e. the profile *magnitude* is
  strictly proportional to `dW = W_metal - 5.4 eV`
- `graphene_photodetector_nonuniform_illumination_model.py` -- inherits the
  above through `total_field`, so the 2026-09-19 `max|k|` ceiling inherits it too

Everything downstream -- the Ti/Pt figure of 1.832, the `max|k|/|N_uniform| <
1.02` bound, the per-metal collection-efficiency ranking -- is a statement
about `dW`, not about the physical Fermi-level shift, *if and only if the two
are proportional*. They are not.

## 2. The actual relation (Khomyakov et al. 2009)

Khomyakov, Giovannetti, Rusu, Brocks, van den Brink and Kelly,
*First-principles study of the interaction and charge transfer between
graphene and metals*, **Phys. Rev. B 79, 195425 (2009)**,
[arXiv:0902.1203](https://arxiv.org/abs/0902.1203),
[full text PDF](https://www.ifw-dresden.de/uploads/users/59/uploads/publications/PhysRevB_79_195425.pdf).
This is the companion paper to the more-cited PRL,
Giovannetti *et al.*, *Doping Graphene with Metal Contacts*,
**Phys. Rev. Lett. 101, 026803 (2008)**,
[arXiv:0802.2267](https://arxiv.org/abs/0802.2267), which states the headline
(|dE_F| <~ 0.5 eV for Al, Ag, Cu, Au, Pt; p/n crossover at ~5.4 eV, not at
graphene's own 4.5 eV) but leaves the derivation to the PRB.

The PRB's analytic model (its Eq. 7) is:

```
                sqrt( 1 + 2 D0 (d - d0) dW' ) - 1
  dE_F(d)  =   -----------------------------------      ,   dW' = W_M - W_G - D_c(d)
                        D0 (d - d0)
```

with the sign carried through `dW'`, and where

- `W_G = 4.5 eV` is free-standing graphene's work function,
- `D_c(d)` is a short-range **chemical** interface term, fitted as
  `D_c(d) = exp(-d/a) (a0 + a1 d + a2 d^2)`, parametrised on Cu,
- `d0 = 2.4 A` is the separation below which the parallel-plate charge-
  transfer term switches off,
- `D0` is the slope of graphene's linear DOS, `D(E) = D0 |E|`.

Two things matter for this thesis.

**(a) The p/n crossover is a consequence, not an input.** Neutrality needs
`dW' = 0`, i.e. `W_0(d) = W_G + D_c(d)`. At the physisorbed equilibrium
separation `d_eq ~ 3.3 A`, `D_c ~ 0.9 eV`, so `W_0 ~ 5.4 eV`. This repo has
been *hard-coding the answer* (`W_CROSS_CHEM = 5.4`) since 2026-09-18 and
citing the PRL for it; the PRB shows it is the value of a function at one
separation. A metal sitting at a different separation has a different
crossover, which is the first hint that a single `w_cross` for all seven
metals is an approximation, not a constant of nature.

**(b) The square root is the whole point.** Expanding Eq. 7 for small `dW'`:

```
  dE_F  =  dW'  -  (D0 (d-d0) / 2) dW'^2  +  O(dW'^3)
```

so the **linear model is the leading term and is exact only at dW' -> 0**.
For large `|dW'|` the relation goes as `sqrt(2 dW' / (D0 (d-d0)))` --
*sublinear*, because pushing graphene's Fermi level further costs quadratically
more transferred charge (`n ~ E_F^2` from the linear DOS) and therefore
quadratically more electrostatic potential drop across the gap. Graphene
resists being doped, increasingly, the harder you push.

The physical origin is exactly the quantum-capacitance argument
`graphene_contact_doping_model.py` already makes in its docstring for the
carrier density at the contact edge -- but the *profile magnitude* used by the
photodetector chain never inherited it. That is the inconsistency this note
closes.

## 3. The constant, derived rather than fitted

The repo does not need the PRB's fitted `D0`; the combination can be built
from first principles. Writing `dW' = phi + (alpha/2) phi^2` with `phi` the
Fermi-level shift in eV:

```
  n(phi)   = (e phi)^2 / (pi hbar^2 v_F^2)          graphene's linear DOS, g=4
  D_tr     = e n (d - d0) / eps_0                   parallel-plate step, in eV
  =>  alpha = 2 e^3 (d - d0) / (eps_0 pi hbar^2 v_F^2)
```

which for the physisorbed `d - d0 = 3.3 - 2.4 = 0.9 A` gives

```
  alpha = 2.39 eV^-1        (v_F = 1e6 m/s)
```

and Eq. 7 in the form actually implemented:

```
  dE_F(dW') = sgn(dW') * ( sqrt(1 + 2 alpha |dW'|) - 1 ) / alpha
```

Cross-check against the PRB: Pt, `W_M = 6.13 eV` (their calculated clean-slab
value, not the 5.65 eV in this repo's `METAL_WORK_FUNCTIONS`, which is an
experimental polycrystalline number), gives `dW' = 0.73 eV` and
`dE_F = 0.47 eV`, against their reported `~0.32-0.33 eV`. **Agreement is to
about 1.4x, i.e. the right size and the right sign, not quantitative.**

**Honest reporting of a data-extraction failure.** Two separate WebFetch
passes over the *same* PRB PDF returned **mutually contradictory** Table I
contents: the first gave `Al -0.57, Ag -0.32, Cu -0.17, Au +0.19, Pt +0.33`,
the second gave `Al +0.19, Ag +0.33, Au -0.17, Pt -0.32` and a different set
of work functions. At least one is a transcription artefact of PDF-to-text
conversion (the sign convention and the metal column appear to have been
shuffled). **No per-metal `dE_F` from that table is used anywhere in the
model or in Chapter 6.** Only quantities that agreed across both passes, and
are independently well known, are used: `d0 = 2.4 A`, `D_c(3.3 A) ~ 0.9 eV`,
`W_0 = 5.4 eV`, and the physisorbed/chemisorbed separations below. The
calibration constant `alpha` is derived, not fitted to those numbers --
which is why the derivation in Section 3 is worth having rather than a
convenience.

## 4. The part that should worry Chapter 6 most: chemisorption

Both fetches agreed on the equilibrium separations, and these are the
well-established qualitative result of the paper:

| bonding | metals | `d_eq` (A) | `d_eq - d0` (A) |
|---|---|---|---|
| physisorbed | Al, Cu, Ag, Au, Pt | ~3.3 | **+0.9** |
| chemisorbed | Pd | ~2.3 | **-0.1** |
| chemisorbed | Ti | ~2.1 | **-0.3** |
| chemisorbed | Ni, Co | ~2.05 | **-0.35** |

For Ni, Co, Ti and Pd, **`d_eq < d0`**. Eq. 7's parallel-plate term is
therefore not merely small, it is *outside the regime the equation was
written for*: the charge-transfer dipole cannot be described as a gap
capacitance when the metal d-states and graphene's pi band have hybridised
and the Dirac cone is destroyed under the contact. This is precisely
Giovannetti *et al.*'s stated conclusion -- that the chemisorbed metals are
**not** characterised by their work function alone.

The uncomfortable consequence for this repo: **the four metals whose
`METAL_WORK_FUNCTIONS` entries produce the largest `|dW|` under the 5.4 eV
crossover -- Ti (-1.07), Cr (-0.90), Cu (-0.75), Ni (-0.36) -- include the two
worst-behaved chemisorbed cases, and the Ti/Pt and Ti/Pd pairs that carry
every headline result of Chapter 6 are Ti-based.** Ti's `dW` is the single
largest in the table and Ti is chemisorbed. The honest position is not that
the Chapter 6 numbers are wrong; it is that Ti's offset is the least
defensible number in the table on grounds the repo has been citing all along
without acting on.

## 5. What was computed (results live in the model file, not here)

`graphene_contact_doping_nonlinear_model.py` does three things:

1. implements `fermi_shift(dW, alpha)` with three **exactly-known**
   validations: `fermi_shift(0) = 0` bitwise; odd symmetry
   `fermi_shift(-x) = -fermi_shift(x)` bitwise; and `alpha -> 0` reproducing
   the linear model to within the analytic `O(alpha dW^2)` residual;
2. re-runs the 2026-09-19 collection kernel with the profile magnitude
   replaced by `fermi_shift(dW)` and reports, side by side with the linear
   result, whether **Ti/Pt's 1.832** and the **`max|k| < 1.02` ceiling**
   survive;
3. sweeps `alpha` from 0 (linear) to 5 eV^-1 so that the conclusion does not
   rest on the 2.39 derived above, which is good to about 1.4x at best.

The prediction, written before running (per the 2026-09-19 practice, which
caught two errors that day): **compression, not reversal.** The square root
is monotone, so no p/n assignment can flip and the symmetric-pair zero is
untouched; but it compresses the ratio of the two contacts' offsets, so the
Ti/Pt asymmetry should shrink by a few tens of percent, and `max|k|` -- a
*ratio* of kernel to its own mean -- should move very little. If the model
instead reverses a ranking, this prediction is wrong and will be recorded as
wrong in the same file.

## 6. Sources

- Giovannetti, Khomyakov, Brocks, Karpan, van den Brink, Kelly,
  "Doping Graphene with Metal Contacts", Phys. Rev. Lett. **101**, 026803
  (2008). https://arxiv.org/abs/0802.2267 --
  https://link.aps.org/doi/10.1103/PhysRevLett.101.026803
- Khomyakov, Giovannetti, Rusu, Brocks, van den Brink, Kelly,
  "First-principles study of the interaction and charge transfer between
  graphene and metals", Phys. Rev. B **79**, 195425 (2009).
  https://arxiv.org/abs/0902.1203 --
  https://www.ifw-dresden.de/uploads/users/59/uploads/publications/PhysRevB_79_195425.pdf
- Khomyakov, Starikov, Brocks, Kelly, "Nonlinear screening of charges induced
  in graphene by metal contacts", Phys. Rev. B **82**, 115437 (2010),
  https://arxiv.org/abs/0911.2027 -- already cited in this repo for the `x^-1`
  spatial falloff; note its title also says *nonlinear*, about the
  **spatial** screening, a separate nonlinearity from the one treated here.

**Fetch log (honest).** 1 WebSearch + 3 WebFetch for this note. `arxiv.org/abs/0802.2267`
returned the PRL abstract (usable). `arxiv.org/abs/0902.1590` was fetched on a
mis-guessed identifier and returned an unrelated paper -- recorded rather than
silently discarded. The IFW Dresden PRB mirror returned usable text for the
equations and separations, but **contradicted itself between two passes on
Table I** (Section 3). PubMed carried the PRL record but was not fetched:
five consecutive reCAPTCHA blocks (2026-09-05, 09-07, 09-17, 09-18) make it a
standing environment limitation, and the arXiv copy is open access anyway.
