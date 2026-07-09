# Reassessment of *Stentor coeruleus* CaMKII homolog candidates as CDPK-like kinases

## Analyses performed

1. **Domain-architecture reassessment**: (InterPro member database evidence) the
   reannotation summary table + the domain-architecture schematic.
2. **CaMKII association-domain search**: the Pfam PF08332 HMM searched proteome-wide
   against the 2017 and present-study annotations (negative result; no hits).
3. **Kinase-domain HMM score comparison**: PF00069-aligned kinase-domain slices scored against
   custom CDPK vs CaMKII HMMs.

## Directory layout

```
cdpk_vs_camkii_reassessment/
├── README_cdpk_vs_camkii_reassessment.md         # this file
├── software_versions.txt              # software versions used in analyses
│
├── 01_domain_architecture/
│   ├── data/
│   │   ├── stentor_CamKII_homolog_candidates_full_length.fasta   # 14 candidates
│   │   └── interpro_files/            # one InterPro TSV per candidate (14 files)
│   ├── scripts/
│   │   └── make_reannotation_summary.py     # make reannotation_summary_AUTO.tsv
│   └── results/
│       └── reannotation_summary_AUTO.tsv    # auto output (pre-curation)
│
├── 02_association_domain_search/       # PF08332 (CaMKII association domain) proteome-wide search
│   ├── data/                           # PF08332.hmm; 2017 + present-study (SteCoeF_v1) proteomes
│   ├── scripts/                        
│       └── run_PF08332_association_domain_search.sh      # hmmsearch PF08332 against the proteomes
│   └── results/                        # hmmsearch output hit tables
│
├── 03_hmm_score_comparison/           
│   ├── data/                                        
│       ├── stentor_CaMKII_homolog_candidates_full_length.fasta       # 14 candidates
│       ├── camk_reference_uniprot_full_length.fasta           # UniProt reference proteins used to build the custom CDPK/CDPK-like and CaMKII-family HMMs
│       ├── PF00069.hmm                              # Pfam protein kinase domain HMM used to identify/extract kinase-domain regions
│       ├── cdpk_training_accessions.txt             # UniProt accessions for CDPK/CDPK-like reference training set
│       └── camkii_training_accessions.txt           # UniProt accessions for CaMKII-family reference training set
│   ├── scripts/                                     
│       ├── run_hmm_score_comparison.sh              # main script; runs the full HMM comparison workflow
│       ├── extract_best_hmmer_domain.py             # extracts the highest-scoring PF00069 kinase-domain slice from each protein
│       ├── split_kinase_domain_sets.py              # splits extracted kinase domains into CDPK training, CaMKII training, and Stentor test sets
│       └── make_cdpk_vs_camkii_score_table.py       # parses HMMER outputs and creates the final CDPK-vs-CaMKII score table
│   ├── generated/                                   # intermediate files generated during workflow - concatenated input FASTAs
│   ├── results/                                     # main analysis outputs; extracted kinase domains, custom HMMs, HMMER tables, final score table
│   └── logs/                                        # runtime logs and diagnostic messages
│
└── manuscript_items/                   # what goes into the paper
    ├── figures/
    │   └── Figure_S7_domain_architecture.svg    # made by hand
    ├── tables/
    │   ├── Table_S19_reannotation_summary_CURATED.tsv    # main output, curated, from analysis 1
    │   └── Table_S20_CDPK_vs_CaMKII_KD_scores.tsv      # main output from analysis 3; stentor_CaMKII_homolog_candidates_CDPK_vs_CaMKII_KD_scores.tsv
    ├── figure_descriptions.md          # Figure S7 legend
    ├── table_descriptions.md           # Table S19 / S20 legends
    └── methods_reassessment.md         # Methods paragraphs + version reporting
```

## mapping to the manuscript (placeholder numbers)

1. Domain architecture, **Figure S7**
1. Domain architecture, curated summary table, **Table S19**
2. Association-domain (PF08332) search, negative result, Methods
3. HMM score comparison, score table, **Table S20** 
all three, Methods paragraphs, `methods_reassessment.md`
all three, code + inputs, paper's *Data and code availability* section

## Analysis 1: Domain-architecture reassessment

In `01_domain_architecture`, run: `python scripts/make_reannotation_summary.py`

AUTO vs. CURATED reannotation table:

`01_domain_architecture/results/reannotation_summary_AUTO.tsv` is produced by
`make_reannotation_summary.py` and features 8 columns: Protein_ID, Length_aa,
StentorDB_annotation, Kinase_domain_evidence, EF_hand_CaM_like_evidence,
InterPro_PANTHER_family, CaMKII_association_hub_domain, Revised_annotation

`manuscript_items/tables/Table_S19_..._CURATED.tsv` is the version in the paper. It
adds (by hand) `new_gene_id_if_mapped` (present-study gene models, by BLASTp of the 
Slabodnick et al., 2017 sequences against the present-study proteome) and a `Notes` column.

## Analysis 2: CaMKII association-domain search

In `02_association_domain_search`, run: `bash scripts/run_PF08332_association_domain_search.sh`

## Analysis 3: Kinase-domain HMM score comparison

In `03_hmm_score_comparison`, run: `bash scripts/run_hmm_score_comparison.sh`

# Overview

The workflow performs the following steps:

1. Concatenates the full-length reference proteins and the full-length *Stentor* candidate proteins
2. Searches all full-length proteins against the Pfam PF00069 protein kinase domain HMM using HMMER `hmmsearch --cut_ga`
3. Extracts the highest-scoring PF00069 kinase-domain slice from each protein using the HMMER target alignment coordinates, `ali_from`–`ali_to`
4. Splits the extracted reference kinase domains into CDPK/CDPK-like and CaMKII-family training sets using predefined UniProt accession lists
5. Aligns each reference training set separately with MAFFT L-INS-i
6. Builds custom CDPK/CDPK-like and CaMKII-family kinase-domain HMMs with HMMER `hmmbuild`
7. Scores the 14 *S. coeruleus* kinase-domain sequences against both custom HMMs using HMMER `hmmsearch`
8. Produces a final table comparing CDPK/CDPK-like and CaMKII-family HMM bit scores for each candidate

Main output: `results/stentor_CaMKII_homolog_candidates_CDPK_vs_CaMKII_KD_scores.tsv`
