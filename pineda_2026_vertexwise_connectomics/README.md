# Vertex-wise structural connectome
Code and derived data to reproduce the figures of

> **Diffusion-based vertex-wise connectomics encodes individual differences in cognition, affect and pathology**
> Pineda J, Cabello-Toscano M, Sarri E, Pascual-Diaz S, for the Alzheimer's Disease Neuroimaging Initiative.

The pipeline that builds the vertex-wise connectome is documented in the interactive toolbox:
<https://saulpascualdiaz.github.io/understanding_vertexwise_connectome>

## Layout

```
scripts/   one script per figure (s01–s12) + two shared helpers (s00_*)
figures/   output (PNG 300 dpi + PDF)
```

| Script | Figure | Data folder |
|---|---|---|
| `s01_cohort_qc.R` | Fig. 1 — cohort characteristics and quality control | `cohort/` |
| `s02_acquisition_ladder.R` | Fig. 2 — identifiability vs. streamline count and shell sampling | `acquisition_ladder/` |
| `s03_distance_bands.R` | Fig. 3 — connection-length bands | `distance_bands/` |
| `s04_community_hierarchy.R` | Fig. 4 — nested community hierarchy | `hierarchy/` |
| `s05_fingerprint_scale.R` | Fig. 5 — fingerprint across spatial scales and atlases | `fingerprint/` |
| `s06_clinical_balance.R` | Fig. 6 — short-to-long balance and clinical measures | `clinical/` |
| `s07_pathology_burden.R` | Fig. 7 — short-to-long balance and pathology burden | `pathology/` |
| `s08_supp_sample_sizes.R` | Supp. Fig. 1 — sample sizes | `cohort/` |
| `s09_supp_short_band_cut.R` | Supp. Fig. 2 — short-band definition | `acquisition_ladder/` |
| `s10_supp_band_cuts.R` | Supp. Fig. 3 — tertile boundaries by cohort | `distance_bands/` |
| `s11_supp_scale_coverage.R` | Supp. Fig. 4 — identifiability, scale and coverage | `fingerprint/` |
| `s12_supp_clinical_models.R` | Supp. Fig. 5 — clinical models and site adjustment | `clinical/` |

`s00_setup.R` holds the palette, theme and input/output helpers; `s00_surface.R` draws vertex maps on the cortical surface (Figs. 6 and 7).
