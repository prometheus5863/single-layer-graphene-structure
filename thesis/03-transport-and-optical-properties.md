# Chapter 3: Transport and Optical Properties of Single-Layer Graphene

*Draft section — first drafted 2026-09-29. Chapter 3 had been listed as
"Computational results complete, chapter not started" in the Chapter 1 status
table since 2026-08-23 — the same label, entered on the same day, that Chapter
2 carried until drafting it on 2026-09-28 revealed that the band structure it
referred to had no Dirac point. The 2026-09-28 log entry recorded that label
as meaning nothing until tested, and made Chapter 3 a risk item rather than a
writing-backlog item on that basis.*

*It was right to. Drafting this chapter found **four** defects in
`graphene_transport_properties.py`, the module the label refers to. They are
set out in Section 3.9, the correction of record, and the audit that found them
is `graphene_transport_optical_audit.py` /
`transport_optical_audit_output.txt`. As in Chapter 2, every one was exposed by
comparison against an **exactly known** value and none by a plausible range,
and — also as in Chapter 2 — every one was in a part of the repository that
nothing downstream imports.*

---

## 3.1 What this chapter is for

Chapter 2 established what graphene's electronic structure gives the device
chapters: a Fermi velocity, a linear density of states, and the absence of a
gap. This chapter takes those three things and asks what they imply for the two
classes of measurement a device engineer actually makes on a material —
**how it carries current**, and **how it interacts with light**.

The framing is deliberately narrower than a survey. Chapters 4 through 7 need
five specific things from this chapter:

1. the **room-temperature mobility** and what limits it, because Chapter 4's
   transfer characteristics are linear in it and Chapter 6's transit time is
   inversely proportional to it;
2. the **minimum conductivity** at the Dirac point, because it is the physical
   origin of the off-state floor that makes a graphene FET a poor switch;
3. the **universal optical absorption**, because Chapter 6's entire
   responsivity budget starts from it;
4. the **Pauli-blocking edge**, because it is the one gate-tunable knob
   graphene's optics has and the reason graphene modulators exist at all;
5. the **half-integer quantum Hall effect**, because it is the cleanest
   experimental confirmation that the Dirac-fermion picture of Chapter 2 is
   physically real and not merely a convenient model.

Everything else — Klein tunneling, weak antilocalisation, hydrodynamic flow —
is noted where relevant and not developed, on the same scoping grounds Section
2.1 gives.

## 3.2 Diffusive transport and the mobility budget

At the carrier densities a device operates at (`10¹² – 10¹³ cm⁻²`, Chapter 4
§4.3) and at room temperature, transport in supported graphene is diffusive,
and the sheet conductivity is

    σ_s = n q_e μ

with `n` from §2.7's `n(E_F) = E_F²/(πħ²v_F²)`. Everything device-relevant
about graphene transport is therefore a statement about `μ`.

Three scattering mechanisms set it, and they combine by Matthiessen's rule:

| mechanism | temperature dependence | typical magnitude on SiO₂ |
|---|---|---|
| charged-impurity / Coulomb | approximately T-independent | 10 000 – 20 000 cm²/V·s |
| longitudinal acoustic phonon (intrinsic) | `μ ∝ 1/T` | ~10⁵ cm²/V·s at 300 K, intrinsic |
| remote (surface) polar phonons of the substrate | strongly activated above ~200 K | caps SiO₂-supported devices near 40 000 cm²/V·s |

The second row is the one this repository had wrong, and it is worth stating
why the correct exponent is `1` and not the `3/2` that a three-dimensional
semiconductor would give. For longitudinal acoustic phonons in graphene the
resistivity is

    ρ_LA = π D_A² k_B T / (4 q_e² ħ ρ_s v_s² v_F²)

which is **linear in T and independent of carrier density** (Hwang and Das
Sarma, *Phys. Rev. B* **77**, 115449 (2008)). Linear-in-T resistivity at fixed
`n` means `μ ∝ 1/T`. The `T^{-3/2}` law is the bulk deformation-potential
result, and it arises from the three-dimensional density of states; graphene's
density of states is the linear one of §2.6, and the exponent changes with it.
Using the 3D exponent understates the 500 K mobility by a factor 1.291 —
small enough to look plausible, which is exactly why it survived (§3.9, T4).

The device-relevant consequence is the ordering. Because the intrinsic phonon
limit is around 10⁵ cm²/V·s at 300 K, **graphene's room-temperature mobility in
a real device is set by its environment, not by graphene.** Charged impurities
in the oxide and remote polar phonons of the substrate are what Chapter 4's
`μ = 4000 cm²/V·s` reflects; the number is a property of the *stack*, and the
factor-of-25 gap between it and the intrinsic limit is the single largest
engineering lever in this thesis. It is also why hBN encapsulation appears in
every high-performance graphene device result in the literature.

## 3.3 The minimum conductivity, and why graphene FETs cannot switch off

The linear density of states vanishes at the Dirac point, so a naive reading of
`σ_s = n q_e μ` predicts zero conductivity there. Measurement does not find
zero. It finds a floor, and the floor is of order the conductance quantum.

Two numbers matter and they are not the same number:

| quantity | value | resistance |
|---|---|---|
| conductance quantum `G₀ = 2q_e²/h` | 7.7481 × 10⁻⁵ S | 12.906 kΩ |
| theoretical `σ_min = 4q_e²/(πh)` | 4.9326 × 10⁻⁵ S | 20.27 kΩ/□ |
| measured `σ_min ≈ 4q_e²/h` | 1.5496 × 10⁻⁴ S | 6.45 kΩ/□ |

The measured floor sits a factor of exactly **π** above the ballistic
self-consistent theoretical value — exactly, because the two expressions differ
only by that factor, so the ratio is a statement about which formula is being
quoted and not an experimental finding. The discrepancy itself is real and
long-standing, and is generally attributed to the charge puddles induced by the
same charged impurities that cap the mobility in §3.2: the sample is never
uniformly at the Dirac point, so the measurement averages over regions of
finite local density. This repository encodes that picture as
`n_puddle = 5 × 10¹⁵ m⁻²` in `graphene_fet_model.py`, and §4.4's on/off ratio
is a direct consequence of it.

Both values are now returned side by side by
`calculate_quantum_conductance()` rather than conflated, because a caller who
needs the off-state floor wants the measured one and a caller checking the
theory wants the other.

**For the device chapters this is the central negative result of the thesis.**
A switch is characterised by its on/off ratio, and graphene's off-state is not
off: it is a 6.45 kΩ/□ resistor. Chapter 4 finds an on/off ratio in the single
digits for exactly this reason, and Chapter 7 §7.2's conclusion — that
graphene's device future is in interconnects, RF and photodetection rather than
in logic — descends from this section more directly than from any other.


Two later results bear on this section's floor: Chapter 4's disorder floor is 2.07× this one, and, because the floor is a geometric mean of the two branch mobilities, a puddle density matched to it is an underestimate (§4.11). The comparisons are in [Appendix A, §3.3.1](appendix-A-verification-record.md#s3-3-1) and [§3.3.2](appendix-A-verification-record.md#s3-3-2).

## 3.4 The half-integer quantum Hall effect

In a perpendicular field the Dirac spectrum gives Landau levels

    E_N = ± v_F √(2 q_e ħ B |N|),  N = 0, ±1, ±2, …

with two features no conventional two-dimensional electron gas has: a level
pinned at exactly **zero** energy, shared between electrons and holes, and a
`√(NB)` spacing rather than the equally spaced `(N + ½)ħω_c`. The zero-energy
level is what produces the half-integer sequence

    σ_xy = 4 q_e²/h × (N + ½)

where the `4` is the combined spin and valley degeneracy. The first plateau sits
at `2q_e²/h`, not at `0` or `4q_e²/h`.

The device-relevant number here is the level spacing, because it determines
whether any of this survives to room temperature:

| B | `E₁` | vs. `k_BT` at 300 K (25.9 meV) |
|---|---|---|
| 1 T | 32.9 meV | 1.27 × |
| 10 T | 104.0 meV | 4.02 × |

Even at 1 T the first Landau gap exceeds the room-temperature thermal energy,
which is why the quantum Hall effect in graphene is observable at 300 K where in
GaAs it demands millikelvins. That is a direct consequence of the large `v_F`
and the `√B` spacing, and it is the most unambiguous experimental confirmation
available that the carriers really are massless Dirac fermions. For this thesis
its role is evidential rather than applied: it is the reason Chapter 2's model
is trusted well enough for Chapters 4–6 to be built on it.

## 3.5 Universal optical absorption

The interband optical conductivity of graphene in the collisionless,
zero-temperature limit is frequency-independent:

    σ₀ = π q_e² / (2h) = q_e² / (4ħ) = 6.0853 × 10⁻⁵ S

Normal-incidence absorption of a free-standing conducting sheet is
`A = σ₀/(ε₀c)` to first order, and substituting gives

    A = π α = 2.29253 %

with `α` the fine-structure constant. This is the celebrated result that a
single atomic layer absorbs a fixed 2.3% of visible light, determined by a pure
number and containing no material parameter at all — not the Fermi velocity, not
the lattice constant, not the hopping integral. The two routes to it,
`σ₀/(ε₀c)` and `πα`, agree to 3.0 × 10⁻¹⁶: they are the same identity written
twice, which is what makes this the strongest anchor in the chapter.

*(A methodological aside worth recording. Checking these two routes against each
other at a tolerance of 10⁻¹² makes the check FAIL, at a residual of 4.5 ×
10⁻¹². The residual is not an error in either route: it is the slack between
CODATA's independently measured `α` and the `α` implied by its own `q_e`, `h`,
`ε₀` and `c`, and it sits well inside the ~1.6 × 10⁻¹⁰ relative uncertainty
CODATA quotes. The correct response was to check the algebraic identity and
name the CODATA slack separately, not to loosen a tolerance until the check
passed.)*

For Chapter 6 the important consequence is a budget rather than a number.
2.3% is simultaneously remarkable — for one atom of material — and ruinous:
97.7% of the incident light passes straight through. Every architecture in
Chapter 6 (plasmonic enhancement §6.5, waveguide integration, microcavities) is
an attempt to buy back that 97.7%, and the fact that the bare figure is a
*constant* is what makes the enhancement factor, rather than the absorption
itself, the design variable.

The corresponding sheet resistance at the universal conductivity,
`1/σ₀ = 16.4 kΩ/□`, is also the reason graphene is not a transparent electrode
for displays without doping: the undoped material trades its excellent
transparency for a sheet resistance two orders of magnitude worse than ITO.

## 3.6 Pauli blocking: the one tunable knob in graphene's optics

The universality of §3.5 holds only while the interband transition is available.
An interband transition at photon energy `ħω` connects states at `−ħω/2` and
`+ħω/2`; if the Fermi level is raised to `E_F`, every final state below `E_F` is
occupied and the transition is blocked. The condition for absorption is
therefore

    ħω > 2 E_F

and gating graphene switches its interband absorption off below that edge. The
edge is where the device interest lies, because `E_F` is electrostatically
controllable through exactly the machinery of Chapter 4:

| `E_F` | onset `2E_F` | cutoff wavelength |
|---|---|---|
| 0.1 eV | 0.20 eV | 6199 nm |
| 0.2 eV | 0.40 eV | 3100 nm |
| 0.3 eV | 0.60 eV | 2066 nm |
| 0.4 eV | 0.80 eV | 1550 nm |

The last row is the one that matters commercially: at `E_F = 0.4 eV` — a gate
swing well within Chapter 4's range — the Pauli edge lands exactly on the
telecom C-band. That coincidence, and not the 2.3% of §3.5, is what makes
graphene electro-absorption modulators a real technology. A material whose
absorption can be switched at a wavelength set by a gate voltage is doing
something no bulk semiconductor with a fixed band gap can do.

For Chapter 6 the same physics is a constraint rather than an opportunity: a
photodetector operated at high `E_F` for low contact resistance is a
photodetector that has begun to switch off its own absorption. That tension is
not currently modelled in Chapter 6, and is recorded below as a follow-on item.

## 3.7 What this chapter does not develop

Three effects are real, frequently cited, and deliberately left out:

- **Klein tunneling.** Normally-incident carriers cross an electrostatic
  barrier with unit transmission, because pseudospin conservation forbids
  backscattering. Device-relevant chiefly as an explanation of *why* §3.3's
  off-state floor cannot be fixed by putting a barrier in the channel, which is
  the form in which Chapter 4 uses it.
- **Weak antilocalisation.** A low-temperature quantum-correction signature of
  the Berry phase; no room-temperature device consequence.
- **Hydrodynamic electron flow.** Requires electron–electron scattering to
  dominate both impurity and phonon scattering, which §3.2 shows does not happen
  in a supported device at 300 K.

## 3.8 What the later chapters import from this chapter

| imported quantity | value | used in |
|---|---|---|
| room-temperature mobility (SiO₂-supported) | 4000 cm²/V·s = 0.4 m²/V·s | Ch. 4 (σ_s, g_m), Ch. 6 (τ_transit) |
| mobility is substrate-limited, not intrinsic | — | Ch. 7 (the hBN lever) |
| measured `σ_min ≈ 4q_e²/h` | 1.5496 × 10⁻⁴ S | Ch. 4 (off-state floor, on/off ratio) |
| puddle density `n_puddle` | 5 × 10¹⁵ m⁻² | Ch. 4 (§4.4) |
| universal absorption `πα` | 2.29253 % | Ch. 6 (responsivity budget, §6.1) |
| Pauli edge `ħω = 2E_F` | table in §3.6 | **nothing yet — see follow-on items** |
| half-integer QHE | — | Ch. 2 (evidential), nothing downstream |
| Landau spacing `v_F√(2q_eħB)` | 32.9 meV at 1 T | nothing downstream |

As in §2.11, the last rows are the honest ones: two of this chapter's results
are cited by nothing. Unlike Chapter 2, where the unused part was also the wrong
part, here the unused rows are correct and the **used** rows are where three of
the four defects lived — `σ₀`, `σ_min` and the mobility exponent are all in the
first five rows. The protection was weaker than Chapter 2's, and it held only
because the device modules re-typed the right answers instead of importing them
(§3.9).

## 3.9 Correction of record (moved)

Four defects in `graphene_transport_properties.py` (T1–T4; e.g. the universal optical conductivity, shipped a factor 2π too large, which implied 14.4044% single-layer absorption instead of 2.2925%) were found and corrected on drafting. Every number in this chapter is the corrected one; see [Appendix A, §3.9](appendix-A-verification-record.md#s3-9).
