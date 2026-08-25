# Biological control-augmented dynamics

This analysis evaluates candidate module-level biological control features in the form:

```text
dz/dtau = f_meth(z,tau) + sum_m beta_m u_m(tau) b_m
```

Current feature modalities and leakage status:
- ATAC: n=5, statuses=methylation_non_leaking_motif_TF_ATAC_gated_no_morula_ATAC

Top model results:
- measured_missing_correction_upper_bound: occupancy=1.000, cosine=1.000, PC3_recovery=1.000, status=measured_residual_upper_bound
- ATAC_ridge_beta_diagnostic: occupancy=1.000, cosine=1.000, PC3_recovery=0.996, status=uses_morula_methylation_residual_for_beta
- combined_external_ridge_beta_diagnostic: occupancy=1.000, cosine=1.000, PC3_recovery=0.996, status=uses_morula_methylation_residual_for_beta
- all_available_ridge_beta_diagnostic: occupancy=1.000, cosine=1.000, PC3_recovery=0.996, status=uses_morula_methylation_residual_for_beta
- ATAC_unit_beta: occupancy=0.222, cosine=0.452, PC3_recovery=0.094, status=external_feature_defined
- combined_external_unit_beta: occupancy=0.222, cosine=0.452, PC3_recovery=0.094, status=external_feature_defined
- all_available_unit_beta: occupancy=0.222, cosine=0.452, PC3_recovery=0.094, status=includes_internal_proxy_features
- methylation_only_strict_baseline: occupancy=0.044, cosine=nan, PC3_recovery=-0.000, status=baseline
- ATAC_unit_beta_sign_flip: occupancy=0.000, cosine=-0.452, PC3_recovery=-0.094, status=sign_flip_control
- combined_external_unit_beta_sign_flip: occupancy=0.000, cosine=-0.452, PC3_recovery=-0.094, status=sign_flip_control

Interpretation boundary: models whose beta was fit to the measured correction are diagnostic upper bounds. A real u_bio result requires external feature-defined controls that improve occupancy without using morula methylation to define beta or direction.
