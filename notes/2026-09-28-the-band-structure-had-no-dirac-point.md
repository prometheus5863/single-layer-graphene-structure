# The band structure had no Dirac point, and three figures beside it could not have noticed

**2026-09-28.** Code, output and figure:
`graphene_band_structure_audit.py`, `band_structure_audit_output.txt`,
`band_structure_audit.png`. Fix and annotations in
`graphene_band_structure.py`. Chapter: `thesis_draft/02-electronic-structure.md`
§2.10.

Not pre-registered. This was not the session's planned work — it was found
while drafting Chapter 2, which is itself the finding (§5 below).

---

## 1. What was wrong

`k_path_graphene()` contained, three lines apart:

```python
k_point = np.array([4.0/(3*np.sqrt(3)), 0.0])   # K point   <-- no pi
m_point = np.array([np.pi/(3*np.sqrt(3)), np.pi/3])  # M point  <-- pi
```

`graphene_hamiltonian()` builds its phases from nearest-neighbour vectors of
**unit length** — the bond-length convention — in which the zone corner is at
`|K| = 4π/(3√3) = 2.41840`. The shipped `K` was `0.76980`: exactly a factor of
π too small (measured ratio **3.141592654**). The shipped `M` was wrong in both
components.

The file mixes three conventions across 294 lines: lattice vectors in Ångström,
Hamiltonian phases in bond-length units, high-symmetry points in neither. The
inconsistency that mattered was **internal to one function**, which is why "it
was a units problem between modules" is not the right description.

## 2. What it cost

| quantity | as shipped | corrected | exact |
|---|---|---|---|
| `\|φ(K)\|` | 2.571775 | 4.4e-16 | 0 |
| gap at the point labelled K | **14.4019 eV** | 2.5e-15 eV | 0 |
| minimum gap over the whole Γ–K–M–Γ path | **11.2000 eV** | 2.5e-15 eV | 0 |
| gap at the point labelled M | 11.2000 eV | 5.6000 eV | `2t` = 5.6 |

Graphene's one defining electronic property was absent from
`band_structure.png` — the repository's oldest figure — for five weeks of daily
automated sessions. The superseded values are kept as `HIGH_SYMMETRY_LEGACY`
and the broken spectrum is still reproducible via
`calculate_band_structure(legacy=True)`, per the 09-24 rule.

## 3. Why five weeks of daily audits never touched it

**Three band/DOS figures exist in this repository and not one of them was a
calculation that could have failed.**

| figure | how it was produced | could it have failed? |
|---|---|---|
| `band_structure.png` | computed from `H(k)` | yes — and it did |
| `band_structure_simple.png` | draws `E = ±v_F\|k\|` | **no** — it plots the known answer |
| `density_of_states.png` | returns `\|E\|/(πt²)` analytically | **no** — never calls `H(k)` at all |

The third is demonstrated rather than asserted, by mutation: replace
`graphene_hamiltonian` with the k-independent gapped matrix `diag(+7, −7)` —
no Dirac point, no dispersion, nothing resembling graphene — and the shipped
DOS output is **bitwise identical**. It would draw the correct V shape for any
band structure whatsoever, including the broken one.

A DOS *computed from the bands* would have caught this on the very first run:
no states near zero energy, van Hove peaks in the wrong place. The check that
was available and would have worked was replaced by a picture of the expected
answer.

Both non-discriminating functions are annotated in place and **retained**, per
09-25 item 11 and 09-26 item 10, with pointers to `dos_from_bands()`.

## 4. The validations, including the two that failed first

Five exact checks, all now passing. The two that failed in their first form are
more informative than the three that did not, and both first forms are kept in
the source.

**V1 — `|φ|` = 3, 1, 0 at Γ, M, K.** Three independent integers, so one check
has three chances to fail. It had three chances for five weeks and took none,
because nothing ever evaluated `|φ|` at a *correct* K.

**V2 — `v_F` from the slope at K. FAILED FIRST, and the failure was physics.**
The first form assumed second-order convergence of the one-sided slope and
failed by a factor of 2 × 10⁴. The measured errors were 2.512e-3, 2.501e-4,
2.500e-5 at `dk` = 1e-2, 1e-3, 1e-4 — clean **first** order. Derived in closed
form along Γ–K, with `u = √3q/2`:

```
|φ| = |1 − cos u − √3 sin u| = √3u − u²/2 + O(u³) = (3/2)q(1 − q/4 + O(q²))
```

so the leading correction is linear in `q` with coefficient exactly **1/4**.
That is trigonal warping. Rewritten to check the derived constant: measured
**0.250013**; correcting the slope by `(1 + q/4)` recovers `v_F =
9.062708e5 m/s` to 1.9e-9.

*An instrument reporting failure can mislead exactly as one reporting success
does.* This is the second instance in one day — the morning's session had a
Richardson oracle withdraw a true claim — and in both cases what settled it was
deriving the expected behaviour in closed form and testing **that** rather than
a decay rate.

**V3 — exact particle–hole symmetry.** `E₊ = −E₋` bitwise at **4000 of 4000**
random k. The strongest check in the chapter and the cheapest.

**V4 — DOS slope from the bands. FAILED FIRST, TWICE, AND ONE OF THE FAILURES
WAS THE AUDITED FAULT COMMITTED BY THE AUDITOR.**

- The analytic constant was first written `2q_e/(πħ²v_F²)`, missing one factor
  of `q_e`.
- More seriously: the reciprocal vectors used to sample the zone were
  `(4π/3)(1,0)` and `(4π/3)(½,√3/2)`. They have the **right magnitude** and
  span the **right area**, but they are **not reciprocal lattice vectors** of
  this lattice — `|φ|` is not periodic under them, because `b·δᵢ` are not all
  equal mod 2π. The parallelogram they span is not a fundamental domain: it
  holds **one** Dirac cone instead of two, which halved the low-energy DOS.
- **And the obvious check did not notice.** "The DOS integrates to 4 states per
  unit-cell area" passed to **sixteen digits** with the wrong vectors in place,
  because the normalisation divides by the sample count and therefore returns
  4 states per cell whatever the Hamiltonian or the sampling region is. It is
  the same fault as the shipped DOS figure, committed inside the module written
  to audit that fault, within the hour.
- It is retained and **printed as vacuous**, beside the two checks that replaced
  it and do have power: periodicity of `|φ|` under `b₁, b₂` (1.6e-15 against a
  derived floor of 1.4e-14), and **a count of the Dirac cones inside the
  sampled cell**, which must be 2 and was 1.

Corrected, the computed DOS matches the exact analytic slope
`2q_e²/(πħ²v_F²) = 1.7891e18 eV⁻¹m⁻²` to **3.6 %** — the residual being the
same `q/4` warping term, against a derived window bound of 19.6 % — and places
the van Hove peak at **2.771 eV** against the exact `t = 2.800 eV`.

**V5 — the diagnosis must reproduce the shipped number.** The 14.401937320702
eV gap equals `2t|φ(K_shipped)|` **bitwise**. Without this the diagnosis would
be a story about a different program.

## 5. The new physics, and it is small but was missing

`result_3_validity_of_the_linear_cone` gives the **first quantitative validity
window for the linear dispersion this thesis has ever had**. Carrier density
from the full nearest-neighbour bands against `n = E_F²/(πħ²v_F²)`:

| `E_F` | full bands | linear cone | ratio |
|---|---|---|---|
| 0.05 eV | 2.2614e11 | 2.2364e11 | 1.0112 |
| 0.10 eV | 8.9044e11 | 8.9455e11 | 0.9954 |
| 0.30 eV | 8.0682e12 | 8.0509e12 | 1.0021 |
| 0.50 eV | 2.2478e13 | 2.2364e13 | 1.0051 |
| 1.00 eV | 9.1423e13 | 8.9455e13 | 1.0220 |
| 1.50 eV | 2.1213e14 | 2.0127e14 | 1.0539 |

**1 % over 0.1–0.5 eV, 5 % to 1.0 eV.** Stated as a *window* and not an upper
limit because the deviation is **not monotone** — 1.0112 at 0.05 eV is worse
than 1.0051 at 0.5 eV — and that is k-sampling, not physics: at 0.05 eV the
occupied cones cover 5.9e-5 of the cell, so an 1800 × 1800 grid places ~190
samples inside them. The audit prints that count so the artefact is legible
rather than smoothed over.

The **sign** is what the device chapters should carry: the full bands hold
*more* states than the cone above ~0.25 eV, because warping flattens the band.
A linear-DOS device model therefore **understates** carrier density at high
bias — conservative for drive current, and against itself for quantum
capacitance.

## 6. Nothing in Chapters 4–7 moves, and the reason is the finding

Those chapters import `v_F`, the linear DOS, and the absence of a gap. They
import the k-path and the 2 × 2 `H(k)` **nowhere** (§2.11 tabulates this). That
is simultaneously:

- why a 14 eV error in the band structure was survivable for five weeks, and
- why no device number changes now that it is fixed.

`v_F` itself is confirmed to 1e-6 relative by the corrected computation.

## 7. What this session concludes, stated as narrowly as the evidence allows

Two things, and the first is the one worth keeping.

**Prose is a detector, and it is the one this repository did not have.** The bug
was not found by an audit, a validation, a convergence study, or a
pre-registration. It was found by **writing the chapter**, because a chapter
cannot contain the sentence "the two bands touch at K" without someone checking
whether they do. Chapter 2 had been marked "computational results complete"
since 2026-08-23 and had never been drafted; **the undrafted chapter is exactly
where the undetected error was.** Five weeks of automated sessions had
systematically preferred new analysis to writing up old analysis — Chapter 7 was
displaced seven consecutive times before 09-27 — and this is the bill for that
preference, stated in eV.

**A check that cannot fail is worse than no check, because it certifies
whatever sits beside it.** Three figures, three unfalsifiable instruments, one
of them created this morning inside the audit module. The defence is not care:
it is to *make the check fail on purpose* — mutate the Hamiltonian and require
the DOS to change; break the reciprocal vectors and require the cone count to
drop; feed the census a snippet containing the pattern it reports absent (which
is what this morning's other session did). Every non-discriminating check found
today was found by that move and by no other.

---

## Open, created or sharpened by today

- **Whether any DEVICE figure in this repository is drawn from a known answer
  rather than computed.** Three of three band-structure figures were. The
  question has not been asked of Chapters 4–6's figures, and it is now the most
  important open item in the repository.
- **Angular trigonal warping.** §2.5 derives the coefficient along Γ–K only;
  the `cos 3θ` structure is what matters for anisotropic transport.
- **Finite-temperature carrier density** `n(E_F, T)`; §2.7 is the `T = 0`
  integral and Chapter 6 runs at room temperature.
- **Whether `t′` can be excluded quantitatively** as the source of Chapter 4's
  electron–hole asymmetry, rather than argued away on magnitude.
- **Chapter 3 remains undrafted**, and after today that is no longer a writing
  backlog item but a **risk** item: its computational results carry the same
  "complete" label Chapter 2's did this morning.
