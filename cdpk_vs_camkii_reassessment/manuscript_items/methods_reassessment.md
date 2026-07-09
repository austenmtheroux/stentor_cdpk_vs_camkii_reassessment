# Methods text for reassessment of *Stentor* CaMKII homolog candidates as CDPK-like kinases

### Reassessment of *Stentor coeruleus* CaMKII homolog candidates

The 14 *S. coeruleus* proteins previously described as CaMKII homologs (named in text and Table S2
in Rajan et al., 2026) were reassessed in three complementary ways: by domain architecture, 
by a targeted search for the CaMKII association domain, and by kinase-domain HMM scoring.

**Domain-architecture annotation.** Domain, family, and functional-site annotations for the
14 candidate proteins were obtained from InterPro (release 5.70-109.0; member databases Pfam 38.2,
PANTHER 19.0, PROSITE Patterns and Profiles 2026_01, SMART 9.0, CDD 3.21, CATH-Gene3D 4.3.0, and
SUPERFAMILY 1.75).
For each protein, kinase-domain, EF-hand/calmodulin-like, and CDPK-family evidence was
tabulated from the member-database matches with their sequence coordinates, and each protein
was checked for the CaMKII association domain (Pfam PF08332 / InterPro IPR013543). Proteins
in which the Ser/Thr protein-kinase active-site signature (PROSITE PS00108) was not detected
were flagged as predicted pseudokinases. The resulting per-protein summary was curated
manually to add SteCoeF_v1 gene-model identifiers, assigned by BLASTp (v2.16.0+; e-value 1e-100,
identity 90%, coverage 99%) of the prior genome (Slabodnick et al., 2017) candidate sequences against 
the SteCoeF_v1 predicted proteome (this study). Every candidate comprised a protein kinase domain followed 
by a C-terminal EF-hand/calmodulin-like region, with PANTHER/InterPro assignment to the calcium-dependent
protein kinase (CDPK) family (PTHR24349 / IPR050205), and none carried a detectable CaMKII
association domain (Table S19; Figure S7).

**CaMKII association-domain search.** To test whether any *S. coeruleus* protein carries the
CaMKII association domain, the Pfam CaMKII association-domain HMM (PF08332, CaMKII_AD) was
searched against both the previously published 2017 annotation (Slabodnick et al., 2017;
34,506 proteins) and the SteCoeF_v1 annotation (this study; 32,306 proteins) using HMMER
`hmmsearch` (v3.4) at the Pfam gathering threshold (`--cut_ga`). No protein in either
annotation returned a hit (best per-proteome domain E-value 2 and 60, respectively, both far above
significance). At profile-HMM sensitivity, the CaMKII association domain is therefore
undetectable in the *S. coeruleus* proteome under either annotation.

**Kinase-domain HMM score comparison.** Full-length sequences of the 14 candidates and of
reference CDPK/CDPK-like and CaMKII-family proteins from UniProt were searched against the
Pfam protein kinase domain HMM (PF00069) with HMMER `hmmsearch --cut_ga` (v3.4); for each
protein the highest-scoring PF00069 match was retained and the kinase-domain region was
extracted using the target alignment coordinates (`ali_from`–`ali_to`). Reference kinase
domains were split into CDPK/CDPK-like and CaMKII-family training sets by predefined UniProt
accession lists, aligned separately with MAFFT L-INS-i (v7.526; `--localpair --maxiterate 1000`),
and used to build custom kinase-domain HMMs with `hmmbuild`. The 14 candidate kinase domains
were then scored against both custom HMMs with `hmmsearch --max -E 1000 --domE 1000`, and model
preference was summarized as Δ bit score (CDPK/CDPK-like minus CaMKII-family bit score); the
candidates were used only as test sequences and were not included in either training set. All
14 candidates scored higher against the CDPK/CDPK-like model (positive Δ bit score; Table S20).
Full workflow and software versions are provided with the deposited code 
(`README_cdpk_vs_camkii_reassessment.md`; `software_versions.txt`).

---

### Data and code availability

Code and inputs for the CaMKII/CDPK reassessment (the domain-architecture summary script,
the CaMKII association-domain (PF08332) search, and the kinase-domain HMM score comparison,
together with the candidate FASTA, InterPro annotation files, custom HMMs, and reference
accession lists) are available at [repository URL]. Software versions used are recorded
in `software_versions.txt`.
