# Path map: 2026-10-16 reorganisation

Every file that moved, old → new. All moves were made with `git mv`, so `git log --follow <new path>` shows the full history. Outputs listed under *new outputs* did not exist in the repository before; they are what the canonical modules produce and were previously either uncommitted or shadowed by `generate_plots.py` (legacy).

| old | new |
|---|---|
| `band_structure_simple.png` | `legacy/figures/band_structure_simple.png` |
| `density_of_states_simple.png` | `legacy/figures/density_of_states_simple.png` |
| `honeycomb_lattice.png` | `legacy/figures/honeycomb_lattice.png` |
| `honeycomb_lattice_simple.png` | `legacy/figures/honeycomb_lattice_simple.png` |
| `optical_absorption.png` | `legacy/figures/optical_absorption.png` |
| `optical_absorption_simple.png` | `legacy/figures/optical_absorption_simple.png` |
| `quantum_hall_effect.png` | `legacy/figures/quantum_hall_effect.png` |
| `quantum_hall_simple.png` | `legacy/figures/quantum_hall_simple.png` |
| `test_plot.png` | `legacy/figures/test_plot.png` |
| `generate_plots.py` | `src/generate_plots.py` (legacy) |
| `initialize_repo.py` | `src/initialize_repo.py` (legacy) |
| `main.py` | `src/main.py` (legacy) |
| `simple_graphene_plots.py` | `src/simple_graphene_plots.py` (legacy) |
| `summary.py` | `src/summary.py` (legacy) |
| `default_scale_audit.png` | `results/figures/audits/default_scale_audit.png` |
| `default_scale_probe_audit.png` | `results/figures/audits/default_scale_probe_audit.png` |
| `figure_provenance_audit.png` | `results/figures/audits/figure_provenance_audit.png` |
| `import_binding_sweep.png` | `results/figures/audits/import_binding_sweep.png` |
| `log_sensitivity_step_audit.png` | `results/figures/audits/log_sensitivity_step_audit.png` |
| `band_structure.png` | `results/figures/ch2/band_structure.png` |
| `band_structure_audit.png` | `results/figures/ch2/band_structure_audit.png` |
| `density_of_states.png` | `results/figures/ch2/density_of_states.png` |
| `contact_argument_blast_radius.png` | `results/figures/ch4/contact_argument_blast_radius.png` |
| `contact_consistent_argument.png` | `results/figures/ch4/contact_consistent_argument.png` |
| `contact_doping_profile.png` | `results/figures/ch4/contact_doping_profile.png` |
| `contact_resistance_crossover.png` | `results/figures/ch4/contact_resistance_crossover.png` |
| `cross_module_delivery_audit.png` | `results/figures/ch4/cross_module_delivery_audit.png` |
| `diffusion_current_model.png` | `results/figures/ch4/diffusion_current_model.png` |
| `double_regularisation.png` | `results/figures/ch4/double_regularisation.png` |
| `edge_vs_top_contact.png` | `results/figures/ch4/edge_vs_top_contact.png` |
| `exact_charge_rederivation.png` | `results/figures/ch4/exact_charge_rederivation.png` |
| `fet_signed_carrier_model.png` | `results/figures/ch4/fet_signed_carrier_model.png` |
| `fmax_shortfall_decomposition.png` | `results/figures/ch4/fmax_shortfall_decomposition.png` |
| `gds_quadrature_audit.png` | `results/figures/ch4/gds_quadrature_audit.png` |
| `gfet_transfer_characteristics.png` | `results/figures/ch4/gfet_transfer_characteristics.png` |
| `perfect_contact_counterfactual.png` | `results/figures/ch4/perfect_contact_counterfactual.png` |
| `quantum_capacitance.png` | `results/figures/ch4/quantum_capacitance.png` |
| `quantum_capacitance_argument.png` | `results/figures/ch4/quantum_capacitance_argument.png` |
| `rc_recalibration.png` | `results/figures/ch4/rc_recalibration.png` |
| `rf_figures_of_merit.png` | `results/figures/ch4/rf_figures_of_merit.png` |
| `sensitivity_audit.png` | `results/figures/ch4/sensitivity_audit.png` |
| `velocity_saturation_model.png` | `results/figures/ch4/velocity_saturation_model.png` |
| `interconnect_resistivity_vs_linewidth.png` | `results/figures/ch5/interconnect_resistivity_vs_linewidth.png` |
| `contact_doping_nonlinear_relation.png` | `results/figures/ch6/contact_doping_nonlinear_relation.png` |
| `crossover_sensitivity.png` | `results/figures/ch6/crossover_sensitivity.png` |
| `differential_crossover.png` | `results/figures/ch6/differential_crossover.png` |
| `per_metal_crossover.png` | `results/figures/ch6/per_metal_crossover.png` |
| `photodetector_collection_efficiency_by_metal.png` | `results/figures/ch6/photodetector_collection_efficiency_by_metal.png` |
| `photodetector_nonuniform_illumination.png` | `results/figures/ch6/photodetector_nonuniform_illumination.png` |
| `photodetector_responsivity_gain_tradeoff.png` | `results/figures/ch6/photodetector_responsivity_gain_tradeoff.png` |
| `photodetector_signed_carrier_response.png` | `results/figures/ch6/photodetector_signed_carrier_response.png` |
| `photodetector_two_contact_net_response.png` | `results/figures/ch6/photodetector_two_contact_net_response.png` |
| `plasmonic_photodetector_enhancement.png` | `results/figures/ch6/plasmonic_photodetector_enhancement.png` |
| `second_anchor_bound.png` | `results/figures/ch7/second_anchor_bound.png` |
| `default_scale_audit_output.txt` | `results/transcripts/audits/default_scale_audit_output.txt` |
| `default_scale_probe_audit_output.txt` | `results/transcripts/audits/default_scale_probe_audit_output.txt` |
| `figure_provenance_audit_output.txt` | `results/transcripts/audits/figure_provenance_audit_output.txt` |
| `import_binding_sweep_output.txt` | `results/transcripts/audits/import_binding_sweep_output.txt` |
| `log_sensitivity_step_audit_output.txt` | `results/transcripts/audits/log_sensitivity_step_audit_output.txt` |
| `mutation_arrival_probe_output.txt` | `results/transcripts/audits/mutation_arrival_probe_output.txt` |
| `pristine_transcript_audit_output.txt` | `results/transcripts/audits/pristine_transcript_audit_output.txt` |
| `band_structure_audit_output.txt` | `results/transcripts/ch2/band_structure_audit_output.txt` |
| `transport_optical_audit_output.txt` | `results/transcripts/ch3/transport_optical_audit_output.txt` |
| `audit_output.txt` | `results/transcripts/ch4/audit_output.txt` |
| `contact_argument_blast_radius_output.txt` | `results/transcripts/ch4/contact_argument_blast_radius_output.txt` |
| `contact_consistent_argument_output.txt` | `results/transcripts/ch4/contact_consistent_argument_output.txt` |
| `cross_module_delivery_audit_output.txt` | `results/transcripts/ch4/cross_module_delivery_audit_output.txt` |
| `dead_name_sweep_output.txt` | `results/transcripts/ch4/dead_name_sweep_output.txt` |
| `diffusion_current_output.txt` | `results/transcripts/ch4/diffusion_current_output.txt` |
| `diffusion_mutation_output.txt` | `results/transcripts/ch4/diffusion_mutation_output.txt` |
| `double_regularisation_output.txt` | `results/transcripts/ch4/double_regularisation_output.txt` |
| `exact_charge_rederivation_output.txt` | `results/transcripts/ch4/exact_charge_rederivation_output.txt` |
| `fet_signed_carrier_output.txt` | `results/transcripts/ch4/fet_signed_carrier_output.txt` |
| `fmax_shortfall_decomposition_output.txt` | `results/transcripts/ch4/fmax_shortfall_decomposition_output.txt` |
| `gds_quadrature_audit_output.txt` | `results/transcripts/ch4/gds_quadrature_audit_output.txt` |
| `perfect_contact_counterfactual_output.txt` | `results/transcripts/ch4/perfect_contact_counterfactual_output.txt` |
| `perfect_contact_mutation_output.txt` | `results/transcripts/ch4/perfect_contact_mutation_output.txt` |
| `potential_census_output.txt` | `results/transcripts/ch4/potential_census_output.txt` |
| `potential_declaration_output.txt` | `results/transcripts/ch4/potential_declaration_output.txt` |
| `quantum_capacitance_argument_output.txt` | `results/transcripts/ch4/quantum_capacitance_argument_output.txt` |
| `quantum_capacitance_overflow_output.txt` | `results/transcripts/ch4/quantum_capacitance_overflow_output.txt` |
| `rf_small_signal_output.txt` | `results/transcripts/ch4/rf_small_signal_output.txt` |
| `sensitivity_converged_output.txt` | `results/transcripts/ch4/sensitivity_converged_output.txt` |
| `velocity_saturation_mutation_output.txt` | `results/transcripts/ch4/velocity_saturation_mutation_output.txt` |
| `velocity_saturation_output.txt` | `results/transcripts/ch4/velocity_saturation_output.txt` |
| `differential_crossover_output.txt` | `results/transcripts/ch6/differential_crossover_output.txt` |
| `rootfinder_audit_output.txt` | `results/transcripts/ch6/rootfinder_audit_output.txt` |
| `second_anchor_bound_output.txt` | `results/transcripts/ch7/second_anchor_bound_output.txt` |
| `bracketed_root.py` | `src/bracketed_root.py` |
| `contact_resistance_crossover.py` | `src/contact_resistance_crossover.py` |
| `graphene_band_structure.py` | `src/graphene_band_structure.py` |
| `graphene_band_structure_audit.py` | `src/graphene_band_structure_audit.py` |
| `graphene_binding_class.py` | `src/graphene_binding_class.py` |
| `graphene_contact_argument_blast_radius.py` | `src/graphene_contact_argument_blast_radius.py` |
| `graphene_contact_consistent_argument.py` | `src/graphene_contact_consistent_argument.py` |
| `graphene_contact_doping_model.py` | `src/graphene_contact_doping_model.py` |
| `graphene_contact_doping_nonlinear_model.py` | `src/graphene_contact_doping_nonlinear_model.py` |
| `graphene_covariance_probe_audit.py` | `src/graphene_covariance_probe_audit.py` |
| `graphene_cross_module_delivery_audit.py` | `src/graphene_cross_module_delivery_audit.py` |
| `graphene_crossover_sensitivity_model.py` | `src/graphene_crossover_sensitivity_model.py` |
| `graphene_dead_name_sweep.py` | `src/graphene_dead_name_sweep.py` |
| `graphene_default_scale_audit.py` | `src/graphene_default_scale_audit.py` |
| `graphene_differential_crossover_model.py` | `src/graphene_differential_crossover_model.py` |
| `graphene_diffusion_current_model.py` | `src/graphene_diffusion_current_model.py` |
| `graphene_diffusion_mutation.py` | `src/graphene_diffusion_mutation.py` |
| `graphene_double_regularisation.py` | `src/graphene_double_regularisation.py` |
| `graphene_edge_contact_model.py` | `src/graphene_edge_contact_model.py` |
| `graphene_exact_charge_rederivation.py` | `src/graphene_exact_charge_rederivation.py` |
| `graphene_fet_model.py` | `src/graphene_fet_model.py` |
| `graphene_fet_signed_carrier_model.py` | `src/graphene_fet_signed_carrier_model.py` |
| `graphene_figure_provenance_audit.py` | `src/graphene_figure_provenance_audit.py` |
| `graphene_fmax_shortfall_decomposition.py` | `src/graphene_fmax_shortfall_decomposition.py` |
| `graphene_gds_quadrature_audit.py` | `src/graphene_gds_quadrature_audit.py` |
| `graphene_import_binding_sweep.py` | `src/graphene_import_binding_sweep.py` |
| `graphene_interconnect_model.py` | `src/graphene_interconnect_model.py` |
| `graphene_log_sensitivity_step_audit.py` | `src/graphene_log_sensitivity_step_audit.py` |
| `graphene_mutation_arrival_probe.py` | `src/graphene_mutation_arrival_probe.py` |
| `graphene_per_metal_crossover_model.py` | `src/graphene_per_metal_crossover_model.py` |
| `graphene_perfect_contact_counterfactual.py` | `src/graphene_perfect_contact_counterfactual.py` |
| `graphene_perfect_contact_mutation.py` | `src/graphene_perfect_contact_mutation.py` |
| `graphene_photodetector_collection_model.py` | `src/graphene_photodetector_collection_model.py` |
| `graphene_photodetector_model.py` | `src/graphene_photodetector_model.py` |
| `graphene_photodetector_nonuniform_illumination_model.py` | `src/graphene_photodetector_nonuniform_illumination_model.py` |
| `graphene_photodetector_signed_carrier_model.py` | `src/graphene_photodetector_signed_carrier_model.py` |
| `graphene_photodetector_two_contact_model.py` | `src/graphene_photodetector_two_contact_model.py` |
| `graphene_plasmonic_photodetector_model.py` | `src/graphene_plasmonic_photodetector_model.py` |
| `graphene_potential_census_audit.py` | `src/graphene_potential_census_audit.py` |
| `graphene_potential_declaration.py` | `src/graphene_potential_declaration.py` |
| `graphene_pristine_transcript_audit.py` | `src/graphene_pristine_transcript_audit.py` |
| `graphene_quantum_capacitance_argument.py` | `src/graphene_quantum_capacitance_argument.py` |
| `graphene_quantum_capacitance_overflow.py` | `src/graphene_quantum_capacitance_overflow.py` |
| `graphene_rootfinder_audit.py` | `src/graphene_rootfinder_audit.py` |
| `graphene_second_anchor_bound.py` | `src/graphene_second_anchor_bound.py` |
| `graphene_sensitivity_audit.py` | `src/graphene_sensitivity_audit.py` |
| `graphene_transport_optical_audit.py` | `src/graphene_transport_optical_audit.py` |
| `graphene_transport_properties.py` | `src/graphene_transport_properties.py` |
| `graphene_velocity_saturation_model.py` | `src/graphene_velocity_saturation_model.py` |
| `graphene_velocity_saturation_mutation.py` | `src/graphene_velocity_saturation_mutation.py` |
| `rf_small_signal_model.py` | `src/rf_small_signal_model.py` |
| `thesis_draft` | `thesis/` |
| `graphene_frozen_default_semantics.py` (added on main 2026-10-16, moved at merge) | `src/graphene_frozen_default_semantics.py` |
| `frozen_default_semantics_output.txt` (same) | `results/transcripts/ch4/frozen_default_semantics_output.txt` |

## New outputs (added 2026-10-16)

Produced by the canonical modules at HEAD and committed for the first time. Before the move the first five
were written to the root but never committed; the last three names were taken by `generate_plots.py`'s
versions, which are now in `legacy/figures/`.

| file | producer |
|---|---|
| `results/figures/ch2/fermi_surface.png` | `src/graphene_band_structure.py` |
| `results/figures/ch3/conductance_map.png` | `src/graphene_transport_properties.py` |
| `results/figures/ch3/klein_tunneling.png` | `src/graphene_transport_properties.py` |
| `results/figures/ch3/linear_dispersion.png` | `src/graphene_transport_properties.py` |
| `results/figures/ch3/mobility_vs_temperature.png` | `src/graphene_transport_properties.py` |
| `results/figures/ch2/honeycomb_lattice.png` | `src/graphene_band_structure.py` |
| `results/figures/ch3/quantum_hall_effect.png` | `src/graphene_transport_properties.py` |
| `results/figures/ch3/optical_properties.png` | `src/graphene_transport_properties.py` |

## Merge of `main` @ `e4fbc26` (2026-11-16)

Files that main added after `77a2f3a`, old → new. Modules go to `src/` (flat, as in §6 of
`docs/INVENTORY.md`); outputs go to `results/<figures|transcripts>/<chapter>/`, the chapter being the first
thesis chapter that cites the producer (audits of the verification layer go to `audits/`). Two committed
records that no module regenerates (`census_environment_verdict_in_anger.txt`,
`standing_list_normalisation_1108.txt`) go to `results/transcripts/audits/` and are listed in
`scripts/reproduce_all.py`'s `NOT_REGENERATED`. Notes stay in `notes/`.

| old | new |
|---|---|
| `ambipolar_mirror_audit.png` | `results/figures/ch4/ambipolar_mirror_audit.png` |
| `ambipolar_mirror_audit_output.txt` | `results/transcripts/ch4/ambipolar_mirror_audit_output.txt` |
| `areal_thermal_resistance.png` | `results/figures/ch4/areal_thermal_resistance.png` |
| `areal_thermal_resistance_output.txt` | `results/transcripts/ch4/areal_thermal_resistance_output.txt` |
| `attribute_default_capture_output.txt` | `results/transcripts/audits/attribute_default_capture_output.txt` |
| `branch_asymmetric_mobility.png` | `results/figures/ch4/branch_asymmetric_mobility.png` |
| `branch_asymmetric_mobility_output.txt` | `results/transcripts/ch4/branch_asymmetric_mobility_output.txt` |
| `capture_surface_guard_output.txt` | `results/transcripts/audits/capture_surface_guard_output.txt` |
| `capture_surface_leverage.png` | `results/figures/audits/capture_surface_leverage.png` |
| `capture_surface_leverage_output.txt` | `results/transcripts/audits/capture_surface_leverage_output.txt` |
| `capture_surface_mutation_output.txt` | `results/transcripts/audits/capture_surface_mutation_output.txt` |
| `census_environment_verdict_in_anger.txt` | `results/transcripts/audits/census_environment_verdict_in_anger.txt` |
| `census_split_equivalence_output.txt` | `results/transcripts/audits/census_split_equivalence_output.txt` |
| `chapter_status_audit_output.txt` | `results/transcripts/audits/chapter_status_audit_output.txt` |
| `chord_decay_length_output.txt` | `results/transcripts/ch4/chord_decay_length_output.txt` |
| `chord_vs_local_decay_length.png` | `results/figures/ch4/chord_vs_local_decay_length.png` |
| `citation_census_output.txt` | `results/transcripts/audits/citation_census_output.txt` |
| `default_flip.png` | `results/figures/ch4/default_flip.png` |
| `default_flip_output.txt` | `results/transcripts/ch4/default_flip_output.txt` |
| `density_dependent_mobility.png` | `results/figures/ch4/density_dependent_mobility.png` |
| `density_dependent_mobility_output.txt` | `results/transcripts/ch4/density_dependent_mobility_output.txt` |
| `graphene_ambipolar_mirror_audit.py` | `src/graphene_ambipolar_mirror_audit.py` |
| `graphene_areal_thermal_resistance.py` | `src/graphene_areal_thermal_resistance.py` |
| `graphene_attribute_default_capture.py` | `src/graphene_attribute_default_capture.py` |
| `graphene_branch_asymmetric_mobility.py` | `src/graphene_branch_asymmetric_mobility.py` |
| `graphene_capture_leverage_figure.py` | `src/graphene_capture_leverage_figure.py` |
| `graphene_capture_surface_guard.py` | `src/graphene_capture_surface_guard.py` |
| `graphene_capture_surface_mutation.py` | `src/graphene_capture_surface_mutation.py` |
| `graphene_census_split_equivalence.py` | `src/graphene_census_split_equivalence.py` |
| `graphene_chapter_status_audit.py` | `src/graphene_chapter_status_audit.py` |
| `graphene_chord_decay_length.py` | `src/graphene_chord_decay_length.py` |
| `graphene_citation_census_audit.py` | `src/graphene_citation_census_audit.py` |
| `graphene_default_flip.py` | `src/graphene_default_flip.py` |
| `graphene_density_dependent_mobility.py` | `src/graphene_density_dependent_mobility.py` |
| `graphene_impurity_limited_crossover.py` | `src/graphene_impurity_limited_crossover.py` |
| `graphene_installed_surface_audit.py` | `src/graphene_installed_surface_audit.py` |
| `graphene_installed_surface_mutation.py` | `src/graphene_installed_surface_mutation.py` |
| `graphene_kappa_impurity_bound.py` | `src/graphene_kappa_impurity_bound.py` |
| `graphene_low_frequency_noise_model.py` | `src/graphene_low_frequency_noise_model.py` |
| `graphene_mixed_capture_repair_audit.py` | `src/graphene_mixed_capture_repair_audit.py` |
| `graphene_replacement_read_audit.py` | `src/graphene_replacement_read_audit.py` |
| `graphene_screening_overdrive_invariance.py` | `src/graphene_screening_overdrive_invariance.py` |
| `graphene_self_heating_model.py` | `src/graphene_self_heating_model.py` |
| `graphene_self_heating_mutation.py` | `src/graphene_self_heating_mutation.py` |
| `graphene_standing_list_audit.py` | `src/graphene_standing_list_audit.py` |
| `graphene_tolerance_census_audit.py` | `src/graphene_tolerance_census_audit.py` |
| `graphene_tolerance_census_figure.py` | `src/graphene_tolerance_census_figure.py` |
| `graphene_value_read_mutation.py` | `src/graphene_value_read_mutation.py` |
| `graphene_value_read_surface.py` | `src/graphene_value_read_surface.py` |
| `graphene_write_read_guard.py` | `src/graphene_write_read_guard.py` |
| `impurity_limited_crossover.png` | `results/figures/ch4/impurity_limited_crossover.png` |
| `impurity_limited_crossover_output.txt` | `results/transcripts/ch4/impurity_limited_crossover_output.txt` |
| `installed_surface_audit_output.txt` | `results/transcripts/audits/installed_surface_audit_output.txt` |
| `installed_surface_mutation_output.txt` | `results/transcripts/audits/installed_surface_mutation_output.txt` |
| `kappa_impurity_bound.png` | `results/figures/ch4/kappa_impurity_bound.png` |
| `kappa_impurity_bound_output.txt` | `results/transcripts/ch4/kappa_impurity_bound_output.txt` |
| `low_frequency_noise_gate_dependence.png` | `results/figures/ch4/low_frequency_noise_gate_dependence.png` |
| `low_frequency_noise_output.txt` | `results/transcripts/ch4/low_frequency_noise_output.txt` |
| `mixed_capture_repair_output.txt` | `results/transcripts/audits/mixed_capture_repair_output.txt` |
| `replacement_read_audit.png` | `results/figures/audits/replacement_read_audit.png` |
| `replacement_read_audit_output.txt` | `results/transcripts/audits/replacement_read_audit_output.txt` |
| `screening_overdrive_invariance.png` | `results/figures/ch7/screening_overdrive_invariance.png` |
| `screening_overdrive_invariance_output.txt` | `results/transcripts/ch7/screening_overdrive_invariance_output.txt` |
| `self_heating_electrothermal.png` | `results/figures/ch4/self_heating_electrothermal.png` |
| `self_heating_mutation_output.txt` | `results/transcripts/audits/self_heating_mutation_output.txt` |
| `self_heating_output.txt` | `results/transcripts/ch4/self_heating_output.txt` |
| `standing_list_audit_output.txt` | `results/transcripts/audits/standing_list_audit_output.txt` |
| `standing_list_normalisation_1108.txt` | `results/transcripts/audits/standing_list_normalisation_1108.txt` |
| `tolerance_census_licensed_vs_tested.png` | `results/figures/audits/tolerance_census_licensed_vs_tested.png` |
| `tolerance_census_output.txt` | `results/transcripts/audits/tolerance_census_output.txt` |
| `value_read_mutation_output.txt` | `results/transcripts/audits/value_read_mutation_output.txt` |
| `value_read_surface_output.txt` | `results/transcripts/audits/value_read_surface_output.txt` |
| `write_read_guard_output.txt` | `results/transcripts/audits/write_read_guard_output.txt` |
