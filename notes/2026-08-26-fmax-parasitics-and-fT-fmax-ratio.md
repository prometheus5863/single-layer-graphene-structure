# Device Physics Notes: Resolving the f_max > f_T Artifact with Realistic Parasitics

**Date:** 2026-08-26 (second session)
**Focus:** Closing the open caveat flagged in `rf_small_signal_model.py` since
2026-08-22 — the simplified analytic f_max estimate exceeding f_T at some
bias points, which that day's notes claimed "does NOT match the literature
trend." This note re-examines that claim against more recent, specific RF
graphene-FET literature, and grounds a revised f_max model that adds source/
drain access resistance and GSG pad capacitance, both currently missing from
`rf_small_signal_model.compute_fT_fmax()`.

## 1. The 2026-08-22 claim was too broad, and the literature is more specific

The 2026-08-22 note asserted (citing Wu et al., *Nano Letters* 2012,
"State-of-the-Art Graphene High-Frequency Electronics",
https://pubs.acs.org/doi/10.1021/nl300904k) that f_max is "typically about
an order of magnitude lower than f_T" for graphene FETs, and treated any
model result with f_max > f_T as unphysical on that basis.

Re-reading the deep-submicron f_max literature already cited that same day
(Feijoo, Pasadas et al., "Deep-submicron Graphene Field-Effect Transistors
with State-of-Art fmax," *Scientific Reports* 6, 35717 (2016),
https://www.nature.com/articles/srep35717) shows this is not a fixed rule —
it is a consequence of *how large gds*(Rg+Rs) is*, which is a design choice,
not an intrinsic graphene property. That paper's own de-embedded numbers:

| Gate length | f_T (de-embedded) | f_max (de-embedded) | f_max / f_T |
|---|---|---|---|
| 200 nm | 35.4 GHz | 50 GHz | 1.41 |
| 300 nm | 34.2 GHz | 45 GHz | 1.32 |
| 400 nm | 25.5 GHz | 35.5 GHz | 1.39 |

f_max *exceeds* f_T at every gate length in this device family, the
opposite of the 2012-era Wu et al. trend. The paper reports a fitted gate
resistance of only ~15 Ω (a low value, deliberately engineered) and states
explicitly that "f_max is inversely proportional to the square root of
r_g" while f_T is independent of r_g — i.e. f_max/f_T is controlled almost
entirely by how well the gate resistance (and, by the same
gds*(R_g+R_s) term, the source/drain access resistance) is engineered down,
not by any graphene-specific ceiling. A companion device family in the same
literature space, "High f_max/f_T ratio in multi-finger embedded T-shaped
gate graphene transistors" (attempted to fetch for this note but the source
returned a rate-limit error; referenced by title only, not used as a
citation for numbers below), corroborates the same qualitative point from
its title alone: multi-finger/T-gate layouts specifically target a high,
not low, f_max/f_T ratio by minimizing R_g.

Before de-embedding, the same 200 nm device in Feijoo et al. measured
f_T/f_max = 21.3/42 GHz — f_max is still above f_T, but both numbers are
substantially suppressed relative to the de-embedded values, which is
direct, paper-reported evidence that un-de-embedded (i.e. "extrinsic,
including pad parasitics") measurements pull both figures down, f_T more
than f_max in this case.

**Revised understanding:** the *sign* of f_max − f_T is not itself a
correctness check for a compact model (2026-08-22's framing was wrong to
treat it as one). What the model needs to get right is the physical
mechanism — gds*(R_g + R_s) and the R_g*C_gd feedback term — with
realistic parameter values, and then whichever ratio falls out is
whichever ratio falls out. This note therefore does not "fix" the
f_max > f_T *result* — it fixes the *model*, by adding the R_s term that
was structurally missing from the f_max formula (only R_g appeared) and a
pad-capacitance term that separates "intrinsic" (de-embedded-equivalent)
from "extrinsic" (raw, pad-parasitic-included) estimates, consistent with
how the literature itself reports both numbers.

## 2. Source/drain access resistance in the f_max formula

The standard hybrid-pi f_max formula used throughout the RF-FET literature
(see e.g. the general delay/parasitic treatment in Wang, Hsu, Wu, "Delay
Analysis of Graphene Field-Effect Transistors," arXiv:1112.4831,
https://arxiv.org/pdf/1112.4831, which identifies source and drain access
resistances R_S, R_D as first-order parasitics that "play a key role in
limiting f_T of short channel GFETs" through the parasitic delay term
τ_par = C_gd(R_S+R_D)[1 + (1+C_gs/C_gd)·g_ds/g_m]) includes an access
resistance contribution to the gds*R term that the f_max formula in
`rf_small_signal_model.py` (added 2026-08-22) omitted — that module's
`compute_fT_fmax()` only includes `gds*Rg`, not `gds*(Rg+Rs)`.

This repo already computes a literature-calibrated total contact
resistance, `graphene_fet_model.Rc_total` (2026-08-21 notes), as a lumped
two-contact series resistance. Splitting it evenly between source and
drain gives R_s = Rc_total / 2 — the natural access-resistance term to add,
reusing an already-calibrated number rather than introducing a new free
parameter.

## 3. GSG pad capacitance for back-gated devices on a conductive substrate

`graphene_fet_model.py`'s device is explicitly back-gated (global Si
back-gate through a 90 nm oxide, i.e. a conductive substrate under the
whole chip). The arXiv:1112.4831 delay-analysis paper notes this exact
substrate choice matters for RF parasitics: devices "on conductive Si
substrates exhibit large parasitics of its GSG probe pads," while
sapphire-substrate devices (highly resistive, >10^16 Ohm.cm) show minimal
pad-to-substrate coupling. Since this repo's device sits on a conductive
(back-gate) substrate by construction, GSG pad capacitance is a relevant,
not a negligible, parasitic for an "extrinsic" (as-measured, pre-de-
embedding) estimate — directly analogous to the raw-vs-de-embedded gap
in Feijoo et al.'s 21.3/42 GHz vs. 35.4/50 GHz numbers above.

A precise pad-capacitance number for this exact geometry was not found in
the sources checked this session (the de-embedding-methodology paper
that would likely report one, "Open-Thru de-embedding for Graphene RF
devices," https://www.researchgate.net/publication/260851587, was not
fetchable in this session — ResearchGate blocks the fetch tool used
here). In the absence of a specific number, this model uses a
representative GSG pad capacitance of C_pad = 15 fF per pad (gate pad and
drain pad each), the low end of the range typically quoted for on-wafer
GSG pads on lossy/conductive substrates in III-V and Si RF-CMOS
de-embedding literature generally; this is flagged in the code as an
assumed, not literature-pinned, value, and is applied only to the
"extrinsic" case for contrast, not to the "intrinsic" curve.

## 4. What changes in the model

`rf_small_signal_model.py` is updated to:

1. Add `Rs = gfet.Rc_total / 2.0` (source access resistance) to the f_max
   denominator: `gds*(Rg + Rs) + 2*pi*fT*Cgd*Rg`, per the standard formula
   and the arXiv:1112.4831 access-resistance discussion above.
2. Compute two variants: an **intrinsic** estimate (no pad capacitance,
   R_s included) and an **extrinsic** estimate (R_s and C_pad = 15 fF/pad
   both included, added to C_gs and C_gd respectively), matching the
   de-embedded-vs-raw structure reported in Feijoo et al.
3. Replace the old blanket "f_max > f_T is unphysical" caveat with the
   corrected framing from Section 1: report both f_T and f_max without
   asserting a required sign, and separately report the R_g and R_s
   contributions so the *mechanism* is checkable even when the resulting
   ratio isn't compared against a fixed rule of thumb.

## References

- Feijoo, Pasadas et al., "Deep-submicron Graphene Field-Effect Transistors
  with State-of-Art fmax," *Scientific Reports* 6, 35717 (2016).
  https://www.nature.com/articles/srep35717
- Wang, Hsu, Wu, "Delay Analysis of Graphene Field-Effect Transistors,"
  arXiv:1112.4831. https://arxiv.org/pdf/1112.4831
- Wu et al., "State-of-the-Art Graphene High-Frequency Electronics,"
  *Nano Letters* 12, 3062 (2012). https://pubs.acs.org/doi/10.1021/nl300904k
  (context/contrast only — this is the source of the "order of magnitude"
  claim revised above)
- "High fMAX/fT ratio in multi-finger embedded T-shaped gate graphene
  transistors" — referenced by title only (fetch attempt returned a
  rate-limit error this session); a follow-up session with working access
  should confirm its reported R_g values and fold them in.

**WebSearch/WebFetch availability:** both were available and used this
session. Two fetch attempts failed (an IEEE Xplore PDF returned HTTP 418;
a ResearchGate page returned a rate limit) and are noted above rather than
silently omitted or worked around with an unsourced guess.
