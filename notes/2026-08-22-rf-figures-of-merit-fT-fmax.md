# Device Physics Notes: RF Figures of Merit (f_T, f_max) for Graphene FETs

**Date:** 2026-08-22
**Focus:** Cutoff frequency (f_T) and maximum oscillation frequency (f_max) —
the two standard RF figures of merit used to benchmark graphene FETs against
III-V and Si RF technologies, and the ones flagged as follow-on work in
`notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`. This note
builds directly on that day's contact-resistance and quantum-capacitance
material and sets up the small-signal model added in
`rf_small_signal_model.py`.

## 1. Why graphene is interesting for RF, and why f_max (not just f_T) is the honest metric

Graphene's zero bandgap makes it unattractive for digital logic (no good
on/off ratio — this repo's `graphene_fet_model.py` transfer curves show the
V-shaped, non-saturating ambipolar characteristic that is the direct
consequence of this), but the same gapless, extremely high-mobility band
structure is attractive for analog/RF applications, where a large on/off
ratio is not required and raw transconductance-to-capacitance ratio matters
more. This has made RF/microwave/mm-wave circuits (frequency multipliers,
mixers, phase shifters) one of the most credible near-term application
spaces for graphene FETs, alongside interconnects and photodetectors.

Two figures of merit are reported, and it is important to distinguish them:

- **f_T (current-gain cutoff frequency):** the frequency at which the
  magnitude of the small-signal current gain |h21| extrapolated with a
  -20 dB/decade slope crosses unity gain.
- **f_max (maximum oscillation / power-gain cutoff frequency):** the
  frequency at which the unilateral power gain U extrapolates to unity —
  the frequency above which the transistor can no longer provide power
  gain into a matched load, i.e. the true ceiling for oscillator/amplifier
  use.

State-of-the-art graphene devices have historically shown a **large gap
between f_T and f_max**, with f_max typically about an order of magnitude
lower than f_T (Wu et al., *Nano Letters* "State-of-the-Art Graphene
High-Frequency Electronics," https://pubs.acs.org/doi/10.1021/nl300904k).
This is the central practical problem in graphene RF electronics: f_T
alone flatters the technology, while f_max is set by the *output
conductance* and *gate resistance*, which the gapless band structure and
contact/interconnect parasitics make hard to control.

## 2. Representative literature numbers

- Intrinsic f_T as high as ~400 GHz (with a widely cited 427 GHz intrinsic
  cutoff frequency figure) has been reported for graphene FETs, and >300 GHz
  intrinsic f_T was demonstrated on wafer-scale CVD-grown graphene and on
  epitaxial graphene on SiC (Wu et al., *Nano Letters*, as above).
- **f_max lags badly behind f_T.** A 60 nm gate-length graphene transistor
  achieved a record (at time of publication) de-embedded **f_max = 200 GHz**
  using a T-shaped/mushroom gate to simultaneously shrink gate length and
  gate resistance (Feng et al., "200 GHz Maximum Oscillation Frequency in
  CVD Graphene Radio Frequency Transistors," *ACS Appl. Mater. Interfaces*,
  https://pubs.acs.org/doi/10.1021/acsami.6b05791). Related device families
  have reported f_max in the ~105 GHz range
  (https://www.nature.com/articles/srep35717, "Deep-submicron Graphene
  Field-Effect Transistors with State-of-Art fmax").
- On flexible substrates (a use case distinct from raw performance
  records, but relevant if graphene's mechanical flexibility is part of
  its value proposition), f_T = 39 GHz was demonstrated while maintaining
  mechanical robustness under bending
  (https://pubmed.ncbi.nlm.nih.gov/27396243/).
- For more conservative, foundry-style multi-finger devices with gate
  lengths from 0.5-2 um (i.e. not exotic record-chasing geometries),
  extrinsic f_T and f_max of ~34 GHz and ~37 GHz were measured at the
  shortest gate length studied, with extrapolation predicting **extrinsic
  f_T, f_max ~ 100 GHz at L_g = 50 nm** — a more realistic "if you actually
  built this at a foundry" number than the exotic record devices above
  (ResearchGate 329285885, "Graphene field-effect transistors with high
  extrinsic f_T and f_max").
- Contact resistance has been directly identified (via benchmarking against
  MoS2 devices) as one of the dominant parasitics separating *intrinsic*
  from *extrinsic* f_T/f_max in graphene and other 2D-material FETs
  (ResearchGate 310661245, "Impact of Contact Resistance on the f_T and
  f_max of Graphene vs. MoS2 Transistors") — this directly connects to the
  R_c ~ 110-500 Ohm.um literature range already gathered in the
  2026-08-21 notes and used in `graphene_fet_model.py`.

## 3. Small-signal model used for the estimate in this repo

The standard hybrid-pi small-signal FET model, applied to a GFET, gives:

```
f_T   = g_m / (2*pi*C_gs)                         [intrinsic]

f_max = f_T / (2 * sqrt( g_ds*(R_g + R_s + R_i) + g_m*R_g*C_gd/C_gs ))
        (Kim et al./standard FET RF approximation; simplified further
         below to the common textbook form used for a first estimate)
```

with the commonly used simplified/practical form (valid when the
gate-resistance-dominated term dominates over drain-conductance terms,
appropriate for a first estimate rather than a full two-port S-parameter
fit):

```
f_max ~= f_T / (2 * sqrt( g_ds * R_g + 2*pi*f_T*C_gd*R_g ))
```

Ingredients, and where they come from in *this* repo's existing model:

- **g_m (transconductance):** dI_d/dV_g extracted numerically from the
  existing `transfer_characteristic()` function in `graphene_fet_model.py`
  — no new physics needed, this is a direct reuse of Monday's GFET model.
- **C_gs (gate-to-channel capacitance):** the *same* series
  C_ox/C_q combination (`quantum_capacitance()`, already implemented)
  that sets the channel charge — i.e. f_T is directly linked to the
  quantum-capacitance physics from the 2026-08-21 notes, not an
  independent parameter. This is the connection flagged as "planned
  follow-on work" in that note.
- **R_g, R_c (gate and contact/access resistance):** literature contact
  resistance range (110-500 Ohm.um) already used for R_c; gate resistance
  is modeled with a simple sheet-resistance-times-geometry estimate for a
  metal gate finger, following the qualitative point made in Feng et al.
  that gate resistance (not just contact resistance) is what limits f_max.
- **g_ds (output conductance):** approximated from the channel resistance
  model's sensitivity to V_ds, since graphene's lack of current saturation
  (visible in the non-flattening `gfet_transfer_characteristics.png`
  curves) is exactly what keeps g_ds large and f_max low relative to f_T —
  this is the physical origin of the f_T/f_max gap documented in the
  literature above.

## 4. What this note sets up

`rf_small_signal_model.py` (companion code commit) implements the
g_m and C_gs extraction from the existing GFET model, computes f_T(V_g)
and a simplified f_max(V_g) across the gate-voltage sweep, and plots both
against the literature benchmarks summarized in Section 2 so the model's
output can be sanity-checked against real reported numbers rather than
just "looking physically reasonable."

## References

1. Wu, Y. et al. "State-of-the-Art Graphene High-Frequency Electronics."
   *Nano Letters* 12, 3062-3067 (2012).
   https://pubs.acs.org/doi/10.1021/nl300904k
2. Feng, Z. et al. "200 GHz Maximum Oscillation Frequency in CVD Graphene
   Radio Frequency Transistors." *ACS Appl. Mater. Interfaces* (2016).
   https://pubs.acs.org/doi/10.1021/acsami.6b05791
3. "Deep-submicron Graphene Field-Effect Transistors with State-of-Art
   fmax." *Scientific Reports* (2016). https://www.nature.com/articles/srep35717
4. "Mechanically robust 39 GHz cut-off frequency graphene field effect
   transistors on flexible substrates." PubMed 27396243.
   https://pubmed.ncbi.nlm.nih.gov/27396243/
5. "Graphene field-effect transistors with high extrinsic f_T and f_max."
   ResearchGate 329285885.
   https://www.researchgate.net/publication/329285885
6. "Impact of Contact Resistance on the f_T and f_max of Graphene vs. MoS2
   Transistors." ResearchGate 310661245.
   https://www.researchgate.net/publication/310661245
7. "A physics-based, small-signal model for graphene field effect
   transistors." *Solid-State Electronics*, ScienceDirect
   S003811011100298X.
   https://www.sciencedirect.com/science/article/abs/pii/S003811011100298X
8. "Small-Signal Capacitance and Current Parameter Modeling in Large-Scale
   High-Frequency Graphene Field-Effect Transistors." ResearchGate
   236341003.
   https://www.researchgate.net/publication/236341003
