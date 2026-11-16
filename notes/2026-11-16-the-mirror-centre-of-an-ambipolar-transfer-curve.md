# The mirror centre of an ambipolar transfer curve: V_Dirac + V_ds/2

*Research notes, 2026-11-16. Sources fetched today; HTTP outcomes recorded in
Section 5 rather than assumed.*

Companion to `graphene_ambipolar_mirror_audit.py` and its output
`ambipolar_mirror_audit_output.txt`. Written for the 2026-11-15 entry's
ranked-second item ("audit the repository for other unipolar forms
transcribed to an ambipolar device") and to score that entry's prediction 1.

---

## 1. Why a unipolar device never raises this question

A silicon n-MOSFET has one branch. `V_gs - V_th` is positive wherever the
device is on, `g_m > 0` at every usable bias, and there is no gate voltage
about which the transfer curve is expected to reflect. Every reflection-
sensitive quantity in a unipolar small-signal model is therefore unconstrained
by symmetry, and a model can be internally inconsistent about *where* the
origin of the gate axis is without anything looking wrong.

Graphene removes that freedom. The transfer curve is V-shaped and ambipolar,
so for a device with equal electron and hole mobility there **is** a
reflection symmetry, and any expression built out of symmetric ingredients
must respect it. That makes the symmetry an instrument: it is a property that
has to hold exactly, so a violation is a defect rather than a modelling
choice. This is the same leverage the 2026-11-15 noise module used in its
V2, and the 2026-09-30 branch-asymmetry module used in its C2.

The question this note answers is narrower than "is the model symmetric":
**about what?**

## 2. The answer in the literature

D. Jiménez, *Explicit drain current, charge and capacitance model of graphene
field-effect transistors*, arXiv:1105.2734 (DOI 10.1109/TED.2011.2168960),
gives the relation directly. Writing `V_gs0` for the gate voltage at which
the channel is neutral at the source end, the paper names

> `V_gs,D = V_gs0 + V_ds/2`

the **Dirac voltage**, and justifies it not by a current minimum but by a
charge balance over the whole channel:

> "For this particular gate voltage, named Dirac voltage, electron and hole
> charges are exactly balanced over the channel length."

Three things in that sentence matter for this repository.

1. **The centre is a channel-integral property, not a terminal property.**
   It is `V_ds/2` because the channel potential runs linearly from `0` at the
   source to `V_ds` at the drain, so the channel's *average* local
   neutrality condition sits at half the drain bias. Anything compared
   against it has to be a channel-averaged quantity too.
2. **It is a statement about charge, not current.** Jiménez does *not* say
   the current minimum is at half the drain bias; he adds that `Q_d` "goes
   across zero at `V_gs` slightly smaller than `V_gs,D`". In the compact
   model in this repository the current minimum and the charge-balance
   point coincide exactly, because the drain current is built as the
   harmonic-style average of a local resistance that depends only on the
   magnitude of the local overdrive — the audit module measures this
   (its V9) rather than importing it.
3. **`V_ds` there is the internal drain-source voltage**, after the contact
   drops `V_ds = V_ds,ext - I_ds(R_s + R_d)` have been subtracted. This
   repository's `transfer_characteristic` integrates over a profile built
   from the *terminal* `Vds` and then adds `Rc_total` in series outside the
   integral, so the two definitions differ by the contact drop. The audit
   module's law is exact for the model as written; it would acquire a
   correction in a model that re-solved the profile self-consistently with
   the contacts. Recorded as an open item, not as a defect.

The same paper is also the reference for the capacitance half of the finding.
Its gate charge is

> `Q_g = -W * integral_0^L Q_c(y) dy`

and the small-signal capacitances are obtained by differentiating these
*channel-integrated* charges with respect to the channel-end potentials,
with Ward–Dutton partitioning weighting charge by position along the channel.
So in the literature model, `C_gs` at finite drain bias is a derived,
non-reciprocal, channel-integrated quantity. It is **not** the source-end
series capacitance evaluated at `V_ch = 0`, which is what
`rf_small_signal_model.gate_capacitance` returns.

## 3. The contrast: a treatment where the question is not raised

K. M. Lee, *Design Considerations of Graphene FETs for RF Applications* (PhD
dissertation, UC San Diego, 2016) states the gate-referenced version
throughout:

> "Figure 1.6 displays Id-Vg characteristic of graphene FET is a V-shaped
> ambipolar curve, with minimum current conduction point at Vg=VDirac."

and

> "When Vtg reaches VDirac (charge neutral point, or Dirac point), the channel
> region becomes intrinsic and the drain current reaches a minimum."

Both sentences identify the minimum-current gate voltage with a single
parameter `V_Dirac` and carry no `V_ds` term. The dissertation uses the
gradual-channel approximation ("A model with the Gradual Channel Approximation
describes the I-V characteristics") with an "effective gate capacitance"
`C_eff`, and does not discuss whether `C_eff` is a source-end or a
channel-averaged object in the portion of the text read.

This is not an error in that dissertation — at the small `V_ds` of a transfer
sweep the `V_ds/2` shift is tens of millivolts and the distinction is
immaterial to its argument. It is, however, **exactly the wording a compact
model inherits**, and it is the wording this repository had: a parameter named
`V_dirac`, described as the gate voltage where the current is minimal, used in
one place as the curve's centre.

## 4. What the repository had, and what is now measured

`gfet.V_dirac = 0.8 V` is a **flat-band parameter**: the built-in doping
offset. It is exactly what its docstring says and nothing about it is wrong.
What was missing is that the *curve's* centre is a different number,
`V_dirac + V_ds/2 = 0.85 V` at the RF bias, and the audit module measures:

- `I_d` is even about `V_dirac + V_ds/2` to **4.5e-16** relative, at six drain
  biases spanning a factor 500 — exact, not asymptotic;
- about `V_dirac` the same curve is asymmetric by **1.05 %** at `V_ds = 0.1 V`;
- `rf.gate_capacitance` is **bitwise** even about `V_dirac`, i.e. it carries
  the other centre as an exact algebraic property;
- so `f_T = |g_m| / (2 pi C_gs)` is even about neither, and the electron- and
  hole-branch peak `f_T` of a device whose `I_d` is exactly mirror-symmetric
  come out unequal by **4.24e-5** (20.27790 vs 20.27876 GHz).

The peak-`f_T` consequence is small and is reported as small: at the peak,
`C_gs` is `C_ox`-dominated and flat, so the centre mismatch has little to act
on. **No committed Chapter 4 number is withdrawn.** Where it is not small is
near neutrality, where `C_q` carries the bias dependence: the committed
`C_gs` differs between the two branches by **2.85 %** at 50 mV from the
centre and **1.74 %** at 100 mV.

## 5. Fetch record

| Source | URL | Outcome |
|---|---|---|
| Jiménez, arXiv:1105.2734 | `arxiv.org/abs/1105.2734` | HTTP 200, abstract too thin for the equations |
| same, full text | `arxiv.org/pdf/1105.2734` | HTTP 200, 29,491 chars, read in full |
| Lee dissertation | `escholarship.org/content/qt2rk2j7rb/...pdf` | HTTP 200, 101,694 chars; chars 0–79,770 read verbatim, 79,770–101,694 **not read** |

Honest limitation on the second source: the fetched copy appears to stop
around Chapter 3 §3.3, so Chapter 5 (compact GFET models, and the figure
captioned "parasitic capacitance (Cgs, Cgd) versus effective device width")
was **not** read. Section 3's claim is therefore about the text that was read
and is not an absence claim about the dissertation as a whole.

Arithmetic caveat carried from the Jiménez paper rather than silently fixed:
the paper's worked example states "At some specific `V_gs,D = V_gs0 + V_ds/2`
(=0.65 V)", while its own Fig. 4 values are `V_gs0 = 0.85 V` and
`V_ds = -0.5 V`, which give 0.60 V. Its companion relation
`V_gs,p-off-d = V_ds + V_gs0 = 0.35 V` is consistent with `V_ds = -0.5 V`. The
*relation* is used here; the paper's 0.65 V is not, and the discrepancy is
recorded because a number taken on trust from a worked example is how a
wrong constant enters a repository.

PubMed/PMC was not attempted. Three prior sessions logged reCAPTCHA pages
from it and the 2026-11-15 entry records that as settled.

## 6. Open after this note

- **Does the centre survive velocity saturation and diffusion?** Both
  `graphene_velocity_saturation_model.py` and
  `graphene_diffusion_current_model.py` integrate over the same symmetric
  profile but with a velocity field that is not even in the local overdrive.
  The audit does not test them.
- **The contact-drop correction.** Jiménez's `V_ds` is internal; this model's
  is terminal. The size of `I_ds(R_s + R_d)` at the RF bias is already
  computed elsewhere in the repository and the two could be reconciled.
- **A Ward–Dutton `C_gs`.** The audit's channel-averaged capacitance is a
  midpoint average of the source-end formula, which is enough to restore the
  shared centre and is *not* the non-reciprocal `C_gs` of the literature
  model. Building the latter is a separate piece of work.
- **Branch-asymmetric mobility.** `mu_e != mu_h` destroys the premise on
  purpose (`graphene_branch_asymmetric_mobility.py`). What replaces the
  centre there — if anything — is not known.
