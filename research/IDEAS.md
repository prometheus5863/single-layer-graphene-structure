# Research ideas

One row per idea. A daily research task reads and updates this table; keep it to the columns below.
`status` is one of: new, checking, running, parked, killed, done.

| id | idea | status | novelty check | next experiment | kill condition | last touched |
|---|---|---|---|---|---|---|
| I1 | Simulate the TLM extraction on top of the doping profile, so `R_extra` and `R_c,measured` are defined on the same region | new | none yet | fit total R vs contact spacing generated from the §4.5 profile with I2's calc (111) W (Au p-type); compare the extracted R_c with R_extra | the extracted R_c absorbs < 10 % of R_extra for every physisorbed metal | 2026-11-17 |
| I2 | Au's doping sign: polycrystalline W = 5.10 eV vs DFT Au(111) at 5.54 eV | parked (closed 2026-11-18 after 2 days; consistency fix, no standalone novelty; output feeds I1/I4) | Giovannetti 2008 + Khomyakov 2009 give the (111) values; not novel alone | done: 5/5 signs vs Khomyakov 2009 printed ΔE_F (RMS 0.053 eV); Au p-type confirmed by STS (Sławińska 2012, +0.25–0.55 eV; model +0.13); with experimental W, Au wrongly goes n, so W and W0 must share a source (`research/I2-au-facet-work-function/RESULT-2026-11-18.md`). No further experiment | no sign changes for any physisorbed metal (not hit) | 2026-11-18 |
| I3 | A doping description for chemisorbed Ni, Pd, Ti that does not go through work function (§7.9 item 2) | parked | not done | literature search for a d_eq- or hybridisation-based doping model | no tabulated model with a stated error bar exists | 2026-11-17 |
| I4 | Single-layer contact model as the normal-metal side of a TBG–superconductor junction (Ch 8) | new | not done | list what Ch 4's contact model must supply to a BTK-type interface (barrier Z, doping, mismatch) | Ch 8's outline needs none of Ch 4's quantities | 2026-11-17 |
