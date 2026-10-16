# Chapter 2: The Electronic Structure of Single-Layer Graphene

*Draft section — first drafted 2026-09-28. Chapter 2 had been listed as
"computational results complete, chapter not started" in the Chapter 1 status
table since 2026-08-23, which made it, together with Chapter 3, the largest
remaining block of pure writing in this thesis. It is drafted here after the
band-structure audit of the same day (`graphene_band_structure_audit.py`,
`band_structure_audit_output.txt`,
`notes/2026-09-28-the-band-structure-had-no-dirac-point.md`), which found that
the computational results were **not** in fact complete: the shipped band
structure had no Dirac point. Section 2.10 is the correction of record, and
Section 2.11 states exactly what the device chapters import from here.*

---

## 2.1 What a device thesis needs from a band structure, and what it does not

This chapter has a narrower job than a condensed-matter thesis's Chapter 2
would. Chapters 4 through 7 are about graphene as a **device material** —
field-effect transistors, interconnects, contacts, photodetectors — and every
one of them reaches back into this chapter for the same short list of
quantities:

1. the **Fermi velocity** `v_F`, which sets the dispersion and hence the
   density of states;
2. the **linear density of states** `g(E) ∝ |E|`, which sets the quantum
   capacitance and the carrier density at a given Fermi level;
3. the **absence of a gap**, which is why every one of those chapters has to
   talk about off-state leakage, ambipolar transfer characteristics, and the
   lack of current saturation;
4. the statement that graphene is a **two-dimensional** conductor with one
   atom of thickness, which is what makes contact resistance and edge
   scattering first-order effects rather than corrections.

What those chapters do **not** use is the full Brillouin-zone band structure,
the k-path, or the 2 × 2 Hamiltonian. This division of labour is worth stating
at the outset for two reasons. The first is honest scoping: a device thesis
that computed a full band structure and then never used it would be padding.
The second is a matter of record, and is the subject of Section 2.10 — the
part of this chapter that nothing downstream depended on is precisely the part
that turned out to be wrong, and it stayed wrong for five weeks of daily
automated work.

## 2.2 Lattice, and the two constants everything else is built from

Graphene is a single sheet of carbon atoms on a honeycomb lattice: a
triangular Bravais lattice with a two-atom basis, conventionally labelled the
A and B sublattices. Two lengths fix everything:

| quantity | symbol | value used throughout this thesis |
|---|---|---|
| lattice constant | `a` | 2.4600 Å |
| C–C bond length | `a_cc = a/√3` | 1.420282 Å |

The Bravais lattice vectors, and the nearest-neighbour vectors from an A site
to its three B neighbours, are

```
a₁ = a(√3/2,  1/2)        δ₁ = a_cc( 0,     1  )
a₂ = a(√3/2, -1/2)        δ₂ = a_cc( √3/2, -1/2)
                          δ₃ = a_cc(-√3/2, -1/2)
```

The three `δ` are of equal length and 120° apart, which is the entire source
of graphene's electronic peculiarity: the sublattice sum below can vanish
identically, and it does so at isolated points rather than along lines.

The reciprocal lattice is triangular, rotated 30° from the real-space lattice,
with primitive vectors of length `4π/(3a_cc) = 4π/3` in units of `1/a_cc`. The
first Brillouin zone is a hexagon whose **six corners fall into two
inequivalent classes**, K and K′, which is the origin of the twofold valley
degeneracy that appears in every carrier-density formula in this thesis. This
is a genuine degeneracy and not a bookkeeping convention: K and K′ are not
connected by any reciprocal lattice vector.

Two points on this lattice are used often enough below to be worth writing
down in the bond-length convention the computation actually uses:

```
K = (4π/(3√3), 0)     |K| = 2.418399
M = (2π/3)(cos30°, sin30°)   |M| = 2.094395
```

Section 2.10 is about what happened because the repository's code had these
written down in a third, inconsistent convention.

## 2.3 The nearest-neighbour tight-binding model

Each carbon contributes one half-filled `p_z` orbital. Keeping only
nearest-neighbour hopping between sublattices, the Bloch Hamiltonian is the
2 × 2 matrix

```
        ⎡   0      t φ(k) ⎤                    3
H(k) =  ⎢                 ⎥ ,      φ(k)  =  Σ  exp(i k·δⱼ)
        ⎣ t φ*(k)    0    ⎦                  j=1
```

with hopping `t = 2.8 eV` throughout this thesis. The eigenvalues are

```
E±(k) = ± t |φ(k)| .
```

Three features follow immediately, and all three are used later as exact
checks on the computation rather than as decoration:

- **Exact particle–hole symmetry.** `E₊(k) = −E₋(k)` for every `k`, with no
  approximation, because the diagonal of `H` is zero. In
  `graphene_band_structure_audit.v3_particle_hole_symmetry` this holds
  **bitwise at 4000 of 4000 random k-points** — not to a tolerance, but to the
  last bit. It is the cleanest validation available in this chapter, and it is
  also the first thing that fails when next-nearest-neighbour hopping is added
  (Section 2.9).
- **A total bandwidth of `6t = 16.8 eV`**, since `|φ|` runs from 0 to 3.
- **The zeros of `φ` are the Dirac points**, and there is nothing else in the
  spectrum at the Fermi level. Graphene is a zero-gap semiconductor with a
  Fermi *surface* that has shrunk to two *points*.

## 2.4 Three exactly known values, and why this chapter leans on them

`|φ|` takes exactly integer values at the three high-symmetry points:

| point | `|φ|` (exact) | `E±` at `t = 2.8 eV` | measured in the audit |
|---|---|---|---|
| Γ | 3 | ± 8.4 eV | 3.000000000000 |
| M | 1 | ± 2.8 eV | 1.000000000000 |
| K | 0 | 0 (degenerate) | 4.4 × 10⁻¹⁶ |

These are not approximations to be checked against a plausible range; they are
integers. That matters methodologically as much as physically. The audit
practice this thesis has converged on over the last two weeks is to test every
computation against a value known in closed form and to report the measured
value beside a **derived** floor rather than against a round tolerance, and
this chapter is unusually well supplied: three independent exact values at
three different points, so a single check has three chances to fail. It had
three chances and, until 2026-09-28, took none of them, because nothing in the
repository ever evaluated `|φ|` at a *correct* K.

The value at M deserves a note of its own. `E(M) = ±t` is where the two bands
are separated by exactly `2t = 5.6 eV`, and it is a **saddle point** of the
dispersion. A saddle point in two dimensions gives a logarithmically divergent
density of states — the van Hove singularity of Section 2.6 — so the exact
value at M and the shape of the DOS are the same statement seen twice.

## 2.5 The Dirac cone, the Fermi velocity, and the exact size of the first
correction to it

Expanding `φ` about K for small `q = k − K` gives the linear dispersion that
the rest of this thesis uses:

```
E±(q) ≈ ± ħ v_F |q| ,        v_F = 3 t a_cc / (2ħ) .
```

With `t = 2.8 eV` and `a_cc = 1.420282 Å` this is

```
v_F = 9.062708 × 10⁵ m/s ,       ħ v_F = 0.5965 eV·nm .
```

— the number quoted loosely as "about 10⁶ m/s, roughly c/300" and used
quantitatively in Chapters 4 and 6. The audit confirms it from the computed
slope to **1 × 10⁻⁶ relative**.

That confirmation required getting one thing right that the audit initially
got wrong, and the correction is a result in its own right. Along the Γ–K
direction the expansion can be done in closed form. With `k = K + q x̂` and
`u = √3 q/2`,

```
|φ| = |1 − cos u − √3 sin u| = √3 u − u²/2 + O(u³)
    = (3/2) q ( 1 − q/4 + O(q²) ) .
```

So the leading correction to the linear cone along this direction is **first
order in `q`, with the exact coefficient 1/4**, not second order. This is
trigonal warping — the hexagonal anisotropy of the bands away from K — seen in
its simplest form. The audit's first validation of `v_F` assumed quadratic
convergence, failed by a factor of 2 × 10⁴, and was rewritten to check the
derived constant instead: the measured coefficient is **0.250013**, and
correcting the finite-difference slope by `(1 + q/4)` recovers the closed-form
`v_F` to 1.9 × 10⁻⁹.

The episode is recorded here rather than tidied away because the failure mode
is instructive and recurs in Section 2.10: a convergence rate that looked
wrong was **physics the check had not accounted for**. An instrument that
reports a failure can be as misleading as one that reports success, and the
only way the distinction was settled was by deriving the coefficient in closed
form and testing *that*.

## 2.6 Density of states, computed rather than assumed

Near the Dirac points, a linear dispersion in two dimensions gives a density
of states linear in energy. Including both spins and both valleys, per unit
area,

```
g(E) = 2|E| / (π ħ² v_F²) ,
```

whose slope in the units this thesis uses is

```
dg/d|E| = 2 q_e² / (π ħ² v_F²) = 1.7891 × 10¹⁸ states eV⁻¹ m⁻²
                               = 1.7891 × 10¹⁴ states eV⁻¹ cm⁻² .
```

Two features of `g(E)` matter downstream. It **vanishes at the Dirac point**,
which is why an ungated graphene sheet has a small but non-zero conductivity
set by disorder rather than by thermal activation; and it is **linear rather
than constant**, unlike a conventional two-dimensional electron gas, which is
what makes the quantum capacitance of Section 2.8 bias-dependent and therefore
a real term in a graphene FET's gate stack rather than a constant absorbed
into the oxide capacitance.

Computed over the full Brillouin zone from `E = ±t|φ|`
(`graphene_band_structure_audit.dos_from_bands`), the density of states
reproduces both the low-energy form and the high-energy structure the
analytic expression cannot describe:

| feature | computed | exact | agreement |
|---|---|---|---|
| low-energy slope | 1.8538 × 10¹⁸ eV⁻¹ m⁻² | 1.7891 × 10¹⁸ | 3.6 % |
| van Hove peak position | 2.771 eV | `t` = 2.800 eV | 1.0 % |

The 3.6 % excess in the slope is not numerical error: it is the same trigonal
warping term of Section 2.5, averaged over the fit window, and the derived
bound for that window is 19.6 %. The van Hove peak is the saddle point at M
(Section 2.4), and its position is a second, independent confirmation that the
high-symmetry points are now consistent with the Hamiltonian.

Section 2.10 explains why this paragraph could not have been written before
2026-09-28: the density of states this repository plotted for five weeks was
not computed from the bands at all.

## 2.7 Carrier density, and the validity range of the linear model

Integrating `g(E)` to a Fermi level `E_F` gives the carrier density the device
chapters use:

```
n = E_F² / (π ħ² v_F²) ,          E_F = ħ v_F √(π n) .
```

| `E_F` | `n` (cm⁻²) | `g(E_F)` (eV⁻¹cm⁻²) |
|---|---|---|
| 0.05 eV | 2.24 × 10¹¹ | 8.95 × 10¹² |
| 0.10 eV | 8.95 × 10¹¹ | 1.79 × 10¹³ |
| 0.20 eV | 3.58 × 10¹² | 3.58 × 10¹³ |
| 0.30 eV | 8.05 × 10¹² | 5.37 × 10¹³ |
| 0.50 eV | 2.24 × 10¹³ | 8.95 × 10¹³ |

The square-root relation `E_F ∝ √n` is worth flagging: to move graphene's
Fermi level by 0.3 eV requires ~8 × 10¹² cm⁻², which is a perfectly ordinary
electrostatic doping level, and this is the reason a graphene FET's channel can
be swung between electron and hole conduction by a gate at all.

**The validity range of the linear form had never been quantified in this
repository, and is quantified here for the first time**
(`result_3_validity_of_the_linear_cone`). Comparing the carrier density from
the full nearest-neighbour bands against the Dirac-cone formula:

| `E_F` | `n` from full bands | `n` from linear cone | ratio |
|---|---|---|---|
| 0.10 eV | 8.9044 × 10¹¹ | 8.9455 × 10¹¹ | 0.9954 |
| 0.30 eV | 8.0682 × 10¹² | 8.0509 × 10¹² | 1.0021 |
| 0.50 eV | 2.2478 × 10¹³ | 2.2364 × 10¹³ | 1.0051 |
| 1.00 eV | 9.1423 × 10¹³ | 8.9455 × 10¹³ | 1.0220 |
| 1.50 eV | 2.1213 × 10¹⁴ | 2.0127 × 10¹⁴ | 1.0539 |

**The linear form is good to 1 % over `0.1 ≤ E_F ≤ 0.5 eV` and to 5 % up to
1.0 eV.** This is stated as a *window* rather than an upper limit because the
deviation is not monotone: at `E_F = 0.05 eV` the ratio is 1.0112, worse than
at 0.5 eV. That is a k-sampling artefact — at 0.05 eV the occupied cones
cover 5.9 × 10⁻⁵ of the reciprocal cell, so an 1800 × 1800 grid places only
about 190 samples inside them — and the audit prints that sample count so the
artefact is legible rather than smoothed over. The clean physical trend is the
monotone rise above 0.2 eV.

The **sign** of the deviation is the part the device chapters should carry
forward. The full band structure has *more* states than the linear cone above
about 0.25 eV, because trigonal warping flattens the band. A device model
built on the linear density of states therefore **understates** carrier density
at high bias: it errs conservatively for drive current, and against itself for
quantum capacitance. Every gate bias in Chapter 4 and every photon energy in
Chapter 6 sits inside the 1 % window, so no number in those chapters is
affected — but the direction of the error is now on record rather than assumed.

## 2.8 Quantum capacitance: where this chapter enters the device chapters most directly

Because `g(E)` is finite and bias-dependent, charging the graphene channel
costs electrostatic energy *and* chemical-potential energy. The second appears
in series with the oxide as a quantum capacitance

```
C_Q = q_e² g(E_F) = 2 q_e² |E_F| / (π ħ² v_F²) .
```

| `E_F` | `C_Q` |
|---|---|
| 0.05 eV | 1.43 µF/cm² |
| 0.10 eV | 2.87 µF/cm² |
| 0.20 eV | 5.73 µF/cm² |
| 0.30 eV | 8.60 µF/cm² |
| 0.50 eV | 14.33 µF/cm² |

These numbers are the reason Chapter 4 cannot treat the gate capacitance as
the oxide capacitance alone. A 1 nm-equivalent high-κ oxide gives a
`C_ox` of order 3–4 µF/cm², which is **comparable to `C_Q` at the Fermi levels
a real device operates at**, not negligible beside it. The series combination
`C_g⁻¹ = C_ox⁻¹ + C_Q⁻¹` therefore degrades gate control worst exactly where
`E_F` is smallest — near the Dirac point, which is where a graphene FET is
supposed to switch. `C_Q ∝ |E_F| ∝ √n` vanishing at the Dirac point is a
statement about the band structure that turns directly into a
transconductance limit, and it is the single place where Chapter 2's physics is
most visible in Chapter 4's numbers.

## 2.9 What this model leaves out, and which omissions matter for later chapters

The nearest-neighbour model above is the right level of description for this
thesis, but the omissions should be named, since three of them are raised again
in Chapters 6 and 7.

- **Next-nearest-neighbour hopping (`t′ ≈ 0.1–0.3 eV`)** adds a diagonal term
  and **breaks the exact particle–hole symmetry** of Section 2.3. Its main
  observable effect is an asymmetry between the electron and hole branches of a
  transfer characteristic. Chapter 4 observes such asymmetry and attributes it
  to contact doping rather than to `t′`; that attribution is a *choice* made on
  the grounds that the measured asymmetry is far larger than `t′/t` would give,
  and it is flagged in Chapter 7's open items rather than claimed as settled.
- **Substrate interaction.** Graphene on SiO₂ is doped, strained and
  disordered by its substrate; graphene on hexagonal boron nitride much less
  so. Nothing in this chapter is a statement about graphene on a substrate,
  which is why Chapters 4 and 5 carry substrate-dependent parameters
  (mobility, residual carrier density, edge roughness) as inputs rather than
  deriving them.
- **Electron–electron interactions.** The Dirac-cone picture is a
  single-particle one. Interactions renormalise `v_F` logarithmically near the
  Dirac point; the value used throughout is the effective one appropriate at
  the carrier densities of Section 2.7, and no chapter in this thesis resolves
  the renormalisation.
- **Spin–orbit coupling**, which is very small in carbon and is neglected
  everywhere here.
- **Trigonal warping**, which is *not* neglected: Section 2.5 gives its exact
  leading coefficient and Section 2.7 gives the resulting error in carrier
  density. It is the one omission of the linear model that this chapter
  quantifies rather than merely names.

## 2.10 Correction of record: this repository's band structure had no Dirac point

**What was wrong.** From the creation of `graphene_band_structure.py` until
2026-09-28, the module's high-symmetry points were written as

```
k_point = [4.0/(3·√3), 0]          ← no π
m_point = [π/(3·√3), π/3]          ← π
```

three lines apart in the same function, while `graphene_hamiltonian()` builds
its phases from nearest-neighbour vectors of unit length — the bond-length
convention of Section 2.2, in which the zone corner is at `4π/(3√3)`. The
shipped K was therefore **exactly a factor of π too small** (the measured
ratio `|K_correct|/|K_shipped|` is 3.141592654), and the shipped M was wrong in
both components.

**What it cost.** `|φ(K_shipped)| = 2.5718` instead of 0. The consequences,
measured:

| quantity | as shipped | corrected | exact value |
|---|---|---|---|
| gap at the point labelled K | **14.4019 eV** | 2.5 × 10⁻¹⁵ eV | 0 |
| minimum gap over the whole Γ–K–M–Γ path | **11.2000 eV** | 2.5 × 10⁻¹⁵ eV | 0 |
| gap at the point labelled M | 11.2000 eV | 5.6000 eV | `2t` = 5.6 eV |

Graphene's single defining electronic property — the gapless linear crossing —
was **absent from `band_structure.png`**, the oldest figure in the repository
and the one this chapter is built on, through five weeks of daily automated
sessions. The superseded values are retained in the module as
`HIGH_SYMMETRY_LEGACY` and the broken spectrum is still reproducible via
`calculate_band_structure(legacy=True)`, per the standing rule that deleting a
superseded implementation destroys the only oracle available for judging its
replacement.

**Why nothing caught it, which is the more useful half.** The repository
contained three band-structure or density-of-states figures and **not one of
them was a calculation that could have failed**:

- `band_structure.png` was computed, and was wrong;
- `band_structure_simple.png` *draws* `E = ±v_F|k|` — the known answer,
  plotted rather than derived;
- `density_of_states.png` came from a function that returns `|E|/(πt²)`
  analytically and **never calls the Hamiltonian at all**. This is demonstrated
  in the audit by mutation: replacing the Hamiltonian with the k-independent
  gapped matrix `diag(+7, −7)` — which has no Dirac point and no dispersion
  whatsoever — leaves the shipped DOS output **bitwise identical**.

A density of states computed from the bands would have caught this on the
first run, because it would have shown no states near zero energy and its van
Hove peaks in the wrong place. The figure that could have falsified the
calculation was instead the figure that certified it. That function is now
annotated in place and retained, labelled as non-discriminating, with a pointer
to `dos_from_bands()`.

**What does and does not move.** Nothing in Chapters 4 through 7 changes. Those
chapters import `v_F` and the linear density of states (Section 2.11) and never
the k-path, which is simultaneously the reason the error was survivable for
five weeks and the reason no device number depends on it. `v_F` itself is
confirmed to 1 × 10⁻⁶ relative by the corrected computation.

**The methodological point, stated narrowly.** The audit that found this also
committed the same fault inside itself, twice, within an hour: its first
normalisation check on the computed DOS — that the density of states integrates
to four states per unit cell — passed to sixteen digits *while the reciprocal
lattice vectors were wrong*, because the normalisation divides by the sample
count and returns four states per cell whatever the Hamiltonian is. It is kept
in the code and printed as **vacuous**, beside the two checks that replaced it
and do have power: periodicity of `|φ|` under the reciprocal vectors, and a
count of the Dirac cones inside the sampled cell (which must be two, and was
one). The lesson this chapter carries forward is therefore not "check your
conventions" but something more specific and more actionable: **a check that
cannot fail is worse than no check, because it certifies whatever sits beside
it**, and the way to tell the two apart is to make the check fail on purpose.

## 2.11 What the later chapters import from this chapter

Stated explicitly, so that the dependency is auditable rather than implied:

| imported quantity | value | used in |
|---|---|---|
| Fermi velocity `v_F` | 9.0627 × 10⁵ m/s | Ch. 4 (transconductance, `f_T`), Ch. 6 (photoresponse) |
| `ħ v_F` | 0.5965 eV·nm | Ch. 5 (edge-scattering length scales) |
| DOS slope `2q_e²/(πħ²v_F²)` | 1.7891 × 10¹⁸ eV⁻¹m⁻² | Ch. 4 (quantum capacitance) |
| `n(E_F) = E_F²/(πħ²v_F²)` | table in §2.7 | Ch. 4, Ch. 6 |
| zero gap | — | Ch. 4 (off-state, ambipolarity), Ch. 6 (broadband absorption) |
| validity window of the above | 1 % for 0.1–0.5 eV | **new, §2.7** |
| the k-path and 2 × 2 `H(k)` | — | **nothing** |

The last row is the one this chapter's correction turns on, and the one before
it is its only genuinely new physical content. Together they describe a useful
and slightly uncomfortable fact about how this thesis is built: the chapter
everything depends on was depended on through **two numbers and one absence**,
and a chapter can be wrong in every part that nothing reads.

---

### Follow-on items generated by this chapter

- The next-nearest-neighbour term `t′` and Chapter 4's electron–hole asymmetry:
  whether `t′` can be *excluded* quantitatively as the source, rather than
  argued away on magnitude (Chapter 7 open item).
- Trigonal warping in the **angular** direction, not just along Γ–K: §2.5
  derives the coefficient for one direction only, and the `cos 3θ` structure is
  the thing that actually matters for anisotropic transport.
- A finite-temperature carrier density, `n(E_F, T)`, since §2.7 is the `T = 0`
  integral and Chapter 6 operates at room temperature.
- Whether any other figure in this repository is drawn from a known answer
  rather than computed — §2.10 found three in the band-structure family alone,
  and the question has not been asked of the device figures.
