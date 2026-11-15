# Low-frequency (1/f) noise in graphene FETs: mechanisms, gate dependence, and a sign convention that does not survive ambipolarity

**Date:** 2026-11-15
**Code:** `graphene_low_frequency_noise_model.py`
**Transcript:** `low_frequency_noise_output.txt`
**Figure:** `low_frequency_noise_gate_dependence.png`

---

## 0. Why this topic, and why it was missing

Before today this repository had no treatment of low-frequency noise at all:
`grep -ril flicker` returned nothing, and so did `noise figure` and `NF_min`.
That is a real gap for a thesis aimed at device/semiconductor roles, because
1/f noise is one of the few figures of merit on which a transistor technology
is accepted or rejected for analog and RF front ends:

- it up-converts to oscillator **phase noise**, contributing a `(1/f)^3` region
  near the carrier (Balandin review, Introduction, citing its Ref. 19);
- the **1/f corner** is quoted as a figure of merit, and for graphene devices
  the review gives it as roughly **1 kHz to 100 kHz** (Sec. II; the
  Introduction gives the wider "a few Hz to tens of kHz");
- unlike thermal noise it cannot be averaged away: the review states that for
  a spectral exponent `gamma >= 1` the accuracy of a 1/f-limited measurement
  "cannot be improved by extending the measuring time" (p. 4);
- the review is explicit that this is a commercial gate, not an academic one:
  "The level of 1/f noise is important for determining the competitiveness of
  graphene technology for cell phones, radars or other communication
  applications" (p. 4).

Graphene sharpens the question relative to silicon because it has **no gap**.
The carrier number `N` passes through a minimum, not through zero, at the
charge-neutrality point (CNP), and both candidate mechanisms are functions of
`N`. So the **gate dependence** of the noise is itself the experiment that
discriminates mechanism, which is why the literature classifies graphene noise
by the SHAPE of `S_I/I^2(V_g)` rather than by one number.

## 1. Sources, with their strength stated

Two sources were fetched and read today. Provenance is recorded at two
strengths throughout, because one of the two equations this module needs was
**not** in the part of its paper that was read verbatim.

### 1.1 Balandin, "Review of the Low-Frequency 1/f Noise in Graphene Devices"
arXiv:1307.4797 — <https://arxiv.org/pdf/1307.4797>
Read: Introduction through Ref. [46] verbatim (about 40k of 50k characters);
reference tail and figure captions summary-only.

VERBATIM from it:

- Normalised spectral density is written `S_I/I^2`; current fluctuations from
  `I = q N mu` (Introduction).
- Hooge, Box I: `S_R/R^2 = S_I/I^2 = S_V/V^2`, with the Hooge parameter `H`
  called "the most commonly used figure of merit", read as
  `S_I/I^2 = alpha_H/(f N)`. **The review warns that applying Hooge's
  volume-noise figure of merit to 2D graphene "presents conceptual
  difficulties"** (Box I) — worth carrying, since this module uses exactly
  that form.
- Magnitude: "Different groups reported consistent values of `S_I/I^2` in the
  range from `1e-9` to `1e-7` Hz^-1 at `f = 10 Hz`" (Sec. II). Area-normalised,
  `(S_I/I^2)(L x W)` is about `1e-8` to `1e-7` um^2/Hz.
- Gate dependence: near the CNP the noise follows a **V-shape with its
  minimum at the Dirac point where the resistance is at its maximum**
  (Sec. II, Fig. 1(c)), reported independently by several groups. "In some
  graphene devices, V-shape becomes **M-shape** dependence over the extended
  bias range."
- Attribution of the M-shape: **spatial charge inhomogeneity** — electron and
  hole puddles (its Ref. 28, Xu 2010). An M-to-V change after annealing is
  attributed to "the interplay between the long- and short-range scattering
  mechanisms" (its Ref. 37, Kaverzin 2012).
- Graphene "does not follow the McWhorter model conventionally used for Si
  CMOS devices", which "predicts that `S_I/I^2` decreases in the inversion
  regime as `~(1/nC)^2`"; deviations are read as contacts, inhomogeneous trap
  distribution, or mobility fluctuations (p. 6).
- Contacts: the noise "mostly originates from graphene itself and is not
  dominated by metal contact contributions", from an area-variation study
  (Sec. II, its Ref. 36).

### 1.2 Mavredakis et al., "Understanding the Bias Dependence of Low Frequency Noise in Single Layer Graphene FETs"
arXiv:1806.06589 — <https://arxiv.org/pdf/1806.06589>
Read: main text through the Conclusions verbatim (about 40k of 68k
characters); **Supplementary Information A–F summary-only**.

VERBATIM from it:

- The CNF local PSD, its Eq. (3), in terms of `(e/Q_gr) * C_q/(C_q + C_top +
  C_back)` and `k T lambda N_T/(W dx f)`, with `N_T` a *volumetric* trap
  density per unit energy in `1/(eV cm^3)`, used as a fitting parameter, and
  `lambda ~ 0.1 nm` the tunnelling attenuation distance.
- A Hooge-type `dmu` term, its Eqs. (6)–(7), with dimensionless `alpha_H`.
- Total is the sum, its Eq. (8).
- Table 1 (electrolyte top-gated SLG, `W = 40 um`): `alpha_H` from `1.1e-3`
  (L = 5.5 um) to `1.5e-3` (L = 43 um); `N_T` from `5.5e19` to `1.3e20`
  `1/(eV cm^3)`; `mu` 3250–3600 cm^2/Vs; other-source rows give `alpha_H` of
  `3e-4` (its Ref. 33) and `7e-3` (its Ref. 32).
- The paper states its fitted `N_T` is **higher than the Si range of
  `1e16`–`1e18` 1/(eV cm^3)**.
- Shapes: `dN` with substantial residual charge `rho_0` gives the **M-shape**;
  `dmu` alone **always** gives a Lambda-shape with a maximum at the CNP;
  classical Hooge alone cannot produce an M-shape. Contact noise is negligible
  near the CNP and matters at high `V_GS`.
- Normalised noise amplitude `B = f * Area * S_I/I^2` is about `1e-7` to
  `1e-6` um^2.

SUMMARY-ONLY from it, and therefore **used but flagged everywhere**:

- Supp. Eq. A24, the simplified CNF form `(g_m/I_D)^2 q^2 kT lambda N_T/(W L
  C_t^2)` — i.e. the `(g_m/I_D)^2 S_Vfb` shape this module implements.
- Supp. Eq. A25, the CNF+CMF form `(g_m/I_D + alpha_c mu C_t)^2 ...` with
  `alpha_c` the "Coulomb scattering coefficient in V.s/C" and
  `alpha_c = 450 V s/C`.

The sub-agent that read the paper flagged, unprompted, that the familiar
`(1 +- alpha_c mu C_ox I_D/g_m)^2` rearrangement is **not** in the paper's
text and should not be attributed to it. It is not used here.

### 1.3 Fetch record

Both fetches returned HTTP 200 on the first attempt and neither needed a
retry. No PubMed/PMC URL was attempted, per the standing note that they
return reCAPTCHA in this environment (now four occurrences logged).
`scipy` is still absent from the device VM and was installed to a `/tmp`
target this session, because the session filesystem was full (Section 5).

## 2. What was built

`graphene_low_frequency_noise_model.py` evaluates three models **on this
repository's own `I_d(V_g)`** — `graphene_fet_model.transfer_characteristic`
and `carrier_density` are imported unchanged, so the device is the same
0.2 um^2, 90 nm-SiO2, 4000 cm^2/Vs GFET the FET and RF chapters already use:

1. CNF: `S_I/I^2 = (g_m/I_d)^2 S_Vfb`, `S_Vfb = q^2 kT lambda N_T/(W L C_ox^2 f)`
2. CNF+CMF: `(g_m/I_d + alpha_c mu C_ox)^2 S_Vfb`
3. Hooge: `alpha_H/(f N)`, `N = n_sheet W L`

with `S_Vfb(N_T = 1.1e20) = 2.479e-5 V^2/Hz` and
`S_Vfb(N_T = 1e17) = 2.254e-8 V^2/Hz` at 10 Hz.

## 3. Finding 1 — the imported sign convention notches the hole branch

The CNF+CMF form is written for a **unipolar** MOSFET, where `g_m > 0` at
every usable bias; the sign of the correlated term is then a convention, and
the literature writes it `+-`. Graphene is **ambipolar**: `g_m` changes sign
at the CNP, so on the hole branch `g_m/I_d` is negative and the sum
`g_m/I_d + alpha_c mu C_ox` passes through **zero**.

Transcribed verbatim, the form therefore predicts a gate voltage at which a
trap-limited graphene FET has **no 1/f noise at all from its dominant
mechanism**. In this device:

- `alpha_c mu C_ox = 6.90627e-2 1/V`, and it is **bias-independent** — it is
  a property of the oxide and the mobility only;
- `g_m/I_d` crosses `-6.906e-2 1/V` at **`V_g = -2.90 V`**, on the hole side;
- there the verbatim form gives `6.17e-13 1/Hz` against `4.72e-7 1/Hz` for the
  branch-corrected form — a factor **7.6e5**;
- its mirror bias `V_g = +4.50 V` (same `|V_g - V_dirac|`, electron side) gives
  `4.79e-7 1/Hz`. So the verbatim form makes two mirror-image bias points
  differ by **7.8e5x**, with nothing physical distinguishing them;
- in the TOTAL (CMF + Hooge), which is what a measurement sees, the Hooge
  floor hides most of it: `6.39e-8` against `5.36e-7 1/Hz`, a factor **8.4**.
  Visible to a fit, easy to miss by eye on a log plot.

The minimal repair is `(|g_m/I_d| + alpha_c mu C_ox)^2`: Coulomb scattering
from a trapped charge raises the scattering rate on **both** branches, so the
correlated term adds in magnitude. Both forms are kept in the module —
`s_cnf_cmf_signed()` so the artefact stays reproducible and checkable, and
`s_cnf_cmf()` for every result.

This is a *modelling* finding, not a claim about the paper: the paper's
devices are electrolyte-gated and it fits them over a gate range where the
cancellation need never be reached. The failure is in transcribing a unipolar
form to an ambipolar device without re-deriving the sign.

## 4. Finding 2 — the M-shape threshold is INSIDE the ordinary SiO2 trap range, which is the opposite of what this module was drafted expecting

Define the central **dip depth** as `log10( min(max_left, max_right) /
S(CNP) )`: positive means an M (noise dips at the CNP), zero means a Lambda
(the maximum IS the CNP). Then:

| `N_T` [1/(eV cm^3)] | dip depth | shape |
|---|---|---|
| 1.1e20 (fitted, electrolyte-gated) | **+0.517 decades** | M |
| 1.0e17 (mid Si/SiO2) | 0 | Lambda |

and the threshold, by bisection, is

> **`N_T* = 7.474e17 1/(eV cm^3)`** — i.e. **inside** the conventional Si/SiO2
> interval of `1e16`–`1e18`, at **75 %** of its top end.

The draft prose of this module asserted the opposite — that `N_T*` would lie
far above the conventional range, and hence that an M-shape measured on a
SiO2-gated GFET could not be homogeneous CNF. **The computation refuted it,
and the prose was corrected rather than the number reinterpreted.** The
corrected statement:

> On this device, homogeneous charge-number fluctuation at a perfectly
> ordinary oxide trap density already produces the M-shape. **The M-shape is
> therefore not by itself evidence of the electron/hole-puddle inhomogeneity
> that the Balandin review attributes it to.** The shape discriminates
> mechanism only in combination with an independent trap-density measurement.

Robustness: `N_T*` is where CNF overtakes the Hooge floor on the flanks, so it
moves with `alpha_H`. Over the full factor-23 span of `alpha_H` that the same
Table 1 reports for other devices,

| `alpha_H` | `N_T*` |
|---|---|
| 3.0e-4 (its Ref. 33) | 1.725e17 |
| 1.3e-3 (used here) | 7.474e17 |
| 7.0e-3 (its Ref. 32) | 4.024e18 |

— the threshold straddles the top of the conventional interval and stays two
orders below the fitted `1.1e20` throughout.

### The Section-3 artefact would have flipped this conclusion

Recomputing `N_T*` on the **verbatim signed** total instead of the corrected
one:

| total used | `N_T*` [1/(eV cm^3)] | ratio |
|---|---|---|
| branch-corrected | **7.474e17** | 1 |
| verbatim signed | **5.466e21** | **7314** |

The dip depth takes the MINIMUM of the two flank maxima, so the notched flank
is the one that decides — and the notch suppresses it. With the imported sign
convention one would therefore have concluded that the M-shape needs a trap
density about **three orders of magnitude above** the conventional SiO2 range,
i.e. **exactly the claim this module was drafted asserting and Section 4
refutes**. A sign convention transcribed from a unipolar device would have
produced a plausible-looking, confidently-stated, wrong physical conclusion,
and nothing in the shape of the curves would have looked odd.

### Bisection on a predicate, not a root

`dip_depth` is **identically zero** over the whole Lambda region (the maximum
IS the CNP there, so the ratio is exactly 1). So the threshold is the boundary
of a region, not a sign change, and a sign-bracketing bisection cannot start —
the first draft tried and raised `ValueError`. The committed version bisects
on the predicate `dip_depth > 0`, and prints the dip depth on both sides of
the final bracket (`+4.06e-8` above, exactly `0` below) so the degeneracy is
visible rather than hidden.

## 5. Validation: ten checks, eight against exactly known values

| # | kind | check | result |
|---|---|---|---|
| V1 | EXACT | `g_m/I_d` of an analytic linear `I(V_g)` equals `1/(V_g-V_0)` (central differences are exact on a line) | 4.09e-14 rel |
| V2 | EXACT | CNF noise is **even about the CNP in the `V_ds -> 0` limit** — `carrier_density` depends on `(V_g - V_dirac)` through its magnitude, so `I_d` is exactly even and `g_m` exactly odd | 2.49e-14 |
| V2b | range | finite-`V_ds` asymmetry halves when `V_ds` halves | ratio 2.0000 |
| V3 | EXACT | both CMF variants at `alpha_c = 0` are **bitwise** CNF | residue 0.0 |
| V4 | EXACT | every term is 1/f: `S(1 Hz)/S(100 Hz) = 100` | 4.26e-16 rel |
| V5 | EXACT | Hooge scales as `1/(W L)`: doubling both divides by exactly 4 | 4.000000000000000 |
| V6 | range | area-normalised `S_I/I^2` overlaps the published `1e-8`–`1e-7` um^2/Hz band | spans 1.78e-9 … 2.60e-8 |
| V7 | EXACT | the signed form is **exactly 0** at `g_m/I_d = -alpha_c mu C_ox` — the artefact pinned as an identity, not described in prose | 0.0 |
| V7b | EXACT | the branch-corrected form there is exactly `4(alpha_c mu C_ox)^2 S_Vfb` | 0 rel err |
| V8 | EXACT | the branch correction is a **bitwise no-op** wherever `g_m >= 0` (117 points) | 0.0 |

**V2 is the check that earned its keep.** The first draft asserted an exact
even-ness on `transfer_characteristic()` itself and failed at `3.4e-6`. The
failure was not a bug: `transfer_characteristic` averages the local channel
resistance over a `0..V_ds` quasi-Fermi profile, so its output is only
*asymptotically* even. The repair was to test the limit itself
(`id_zero_vds_limit`, `V_ch = 0` exactly), which is exactly even and passes at
2.5e-14 — and to keep the finite-`V_ds` behaviour as the separate, honestly
weaker check V2b, which confirms the residue is linear in `V_ds`. A
plausible-range check would have accepted `3.4e-6` without ever asking which
of the two it was.

V6 is deliberately labelled a RANGE check: the model spans `1.78e-9` to
`2.6e-8` um^2/Hz against a published band of `1e-8`–`1e-7`, which **overlaps**
rather than lands inside, and the transcript says "overlaps".

## 6. Environment note

The device VM's session filesystem was **100 % full** at the start of this
session (`/sessions` 9.3G used, 0 available, with only 29 MB of it this
session's), so the standard `$HOME/work` clone location failed with "No space
left on device". Work moved to a per-session directory under `/tmp` (on a
different filesystem, 3.1 GB free), and `scipy` — still not preinstalled —
was installed with `pip --target /tmp/pylibs`, with `MPLCONFIGDIR` and
`TMPDIR` redirected likewise. Both repositories cloned and ran normally after
that. The sibling repository logged a related `TMPDIR` repair on 2026-11-14;
this is the same condition seen from the other side.

## 7. Not yet covered, for a future session

- **Measure the 1/f corner frequency**, which is the review's own stated
  figure of merit and the quantity an RF designer actually asks for. It needs
  the thermal + shot noise floor of the same device, which the repository does
  not have, and crosses the existing `rf_small_signal_model.py`.
- **A puddle-FLUCTUATION model.** `n_puddle` here regularises the density; it
  does not fluctuate. Until it does, this module cannot weigh the two
  mechanisms the review names against each other — it can only show that one
  of them suffices.
- **Phase-noise up-conversion.** `(1/f)^3` near the carrier is the reason an
  RF chapter should care; nothing here converts `S_I/I^2` into `L(f)`.
- **Area/oxide scaling of `N_T*`.** `S_Vfb` carries `1/(W L C_ox^2)`, so the
  threshold moves with geometry; today's number is one device.
- **The contact-noise term** that arXiv:1806.06589 adds at high `|V_g|`;
  without it the flanks here are a lower bound.
- **Where else in this repository is a unipolar form transcribed to an
  ambipolar device?** Section 3 is one instance of a class. Any imported
  expression carrying a `g_m` (not `|g_m|`), a `V_gs - V_th` overdrive, or a
  "+-" sign the source left as a convention is a candidate. `rf_small_signal_model.py`
  and the Chapter 4 FET material are where to look first. Not audited today.
- **A puddle-fluctuation model is now the ranked item**, because Section 4
  removed the only argument this repository had for *needing* one.
