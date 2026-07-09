# Supplementary table descriptions

## Table S19. Domain-architecture reassessment of *Stentor coeruleus* CaMKII homolog candidates.

Per-protein InterPro / member-database domain evidence for the 14 *S. coeruleus* proteins
previously annotated as CaMKII homologs (Rajan et al., 2026). Columns: **Protein_ID**,
the source (StentorDB) gene model (Slabodnick et al., 2017); **new_gene_id_if_mapped**, the corresponding
SteCoeF_v1 gene model (this study), assigned by BLASTp of the 2017 sequences against
the SteCoeF_v1 predicted proteome; **Length_aa**; **StentorDB_annotation**, the
description in the source annotation; **Kinase_domain_evidence**,
**EF_hand_CaM_like_evidence**, and **InterPro_PANTHER_family**, the supporting
member-database signatures with match coordinates; **CaMKII_association_domain**, the
result of the Pfam PF08332 search (see Methods); **Revised_annotation**, the reassessed
description; and **Notes**, manual-curation remarks. All 14 candidates carry a protein
kinase domain, a C-terminal EF-hand/calmodulin-like region, and PANTHER/InterPro assignment
to the calcium-dependent protein kinase (CDPK) family (PTHR24349 / IPR050205), and none
carries a detectable CaMKII association domain.


## Table S20. CDPK-versus-CaMKII kinase-domain HMM score comparison for *Stentor coeruleus* CaMKII homolog candidates.

PF00069-aligned kinase-domain slices from the 14 *S. coeruleus* candidates were scored against
custom HMMs built from reference CDPK/CDPK-like and CaMKII-family kinase domains. Columns:
candidate protein ID; CDPK/CDPK-like model bit score and E-value; CaMKII-family model bit
score and E-value; Δ bit score (CDPK/CDPK-like model bit score minus CaMKII-family model
bit score); and preferred model. Positive Δ bit scores indicate stronger similarity to the
CDPK/CDPK-like model. Reference training sets, custom-HMM construction, and the full
workflow are described in Methods and in `README_cdpk_vs_camkii_reassessment.md`. The *Stentor*
candidates were used only as test sequences and were not included in either training set.
