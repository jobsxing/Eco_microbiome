# OH-MicroDynamics analysis code
Analysis code accompanying the OH-MicroDynamics (One Health Microbiome
Dynamics) project, a year-long longitudinal metagenomic study of microbiome,
resistome, strain-sharing, candidate horizontal gene transfer (HGT), and
population-genomic dynamics across a human-animal-environment farming
ecosystem.
## Repository status
- **Study:** OH-MicroDynamics
- **Release type:** manuscript-associated analysis code
- **Primary languages:** R and Python
Replace the bracketed records above when the GitHub repository, Zenodo archive, and manuscript record are public.
## Analyses included
The deposited workflow covers:
1. incidence-based sampling completeness and gamma-diversity estimation;
2. sensitivity to the combined sampling-round/extraction/sequencing
   `batch/time` factor;
3. repeated-measures-corrected richness and Shannon-diversity comparisons at
   genus, family, and species levels;
4. taxonomic and HUMAnN pathway association models, including transformation
   sensitivity analyses;
5. ARG-class association models, abundance-scale sensitivity analyses, and
   microbiome-resistome Procrustes analysis;
6. canonical same-round, intra-habitat, inter-habitat, denominator-aware, and
   temporally lagged strain-sharing analyses;
7. ribosomal-feature and taxonomic-relatedness sensitivity analyses for
   candidate HGT fragments; and
8. multiple-testing-corrected McDonald-Kreitman (MK), asymptotic MK, and
   stratified gene-recurrence analyses for the focal lineage.
## Repository structure
|-- code/
|   |-- analysis/
|   |   |-- 00_validate_inputs.py
|   |   |-- 01_sampling_sufficiency/
|   |   |-- 02_batch_effects/
|   |   |-- 03_alpha_diversity/
|   |   |-- 04_feature_associations/
|   |   |-- 05_arg_associations/
|   |   |-- 06_strain_sharing/
|   |   |-- 07_hgt_sensitivity/
|   |   `-- 08_population_genomics/
|   |-- check_release.py
|   `-- run_manifest.tsv
|-- data/
|   `-- metadata_template.csv
