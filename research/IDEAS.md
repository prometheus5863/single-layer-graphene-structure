# Research ideas

One row per idea. A daily research task reads and updates this table; keep it to the columns below.
`status` is one of: new, checking, running, parked, killed, done.

| id | idea | status | novelty check | next experiment | kill condition | last touched |
|---|---|---|---|---|---|---|
| I1 | Simulate the TLM extraction on top of the doping profile, so `R_extra` and `R_c,measured` are defined on the same region | new | none yet | fit total R vs contact spacing generated from the §4.5 profile; compare the extracted R_c with R_extra | the extracted R_c absorbs < 10 % of R_extra for every physisorbed metal | 2026-11-17 |
| I2 | Au's doping sign: polycrystalline W = 5.10 eV vs DFT Au(111) at 5.54 eV | running (1 day) | Giovannetti 2008 and Khomyakov 2009 give the (111) value; result alone not novel (consistency fix) | 2026-11-17: Au flips n→p, R_extra −161→+947 Ω·µm, 5/5 signs vs DFT, RMS ΔE_F 0.074 eV (`research/I2-au-facet-work-function/`). Next: second reference (Khomyakov 2009 tabulated ΔE_F) + one experimental Au-contact sign, then close | no sign changes for any physisorbed metal (not hit: Au flips) | 2026-11-17 |
| I3 | A doping description for chemisorbed Ni, Pd, Ti that does not go through work function (§7.9 item 2) | parked | not done | literature search for a d_eq- or hybridisation-based doping model | no tabulated model with a stated error bar exists | 2026-11-17 |
| I4 | Single-layer contact model as the normal-metal side of a TBG–superconductor junction (Ch 8) | new | not done | list what Ch 4's contact model must supply to a BTK-type interface (barrier Z, doping, mismatch) | Ch 8's outline needs none of Ch 4's quantities | 2026-11-17 |
