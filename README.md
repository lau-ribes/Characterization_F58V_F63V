# Genomic and in silico analyses — F58V and F63V bacteriophages

## Overview

This repository contains all bioinformatic input and output files supporting
the manuscript:

> **A conserved surface polysaccharide degradation module in Phapecoctavirus
> is associated with differential anti-biofilm activity against
> ESBL-producing *Escherichia coli***
> AUTORES — to be submitted to *Communications Biology*, 2026

Two lytic bacteriophages isolated from hospital wastewater (Hospital
Universitario Infanta Elena, Valdemoro, Madrid) were characterized: F58V
(Phapecoctavirus, Straboviridae, ~149 kb myovirus) and F63V (novel genus,
Drexlerviridae, ~46 kb siphovirus), both targeting ESBL-producing *E. coli*
clinical isolates.

The central argument is that a predicted tripartite colanic acid depolymerase
module in F58V (CDS_0012, absent in F63V) is associated with superior
activity against preformed biofilms, supported by phenotypic, structural,
and evolutionary evidence.

---

## Repository structure

### `/genomes`
Assembled genome sequences in FASTA format.
- `F58V.fasta` — F58V complete genome (149,321 bp, GC 39.0%)
- `F63V.fasta` — F63V complete genome (45,748 bp, GC 44.1%)
- `all_phapecoctavirus_accessions.txt` — GenBank accession numbers for the
  24 *Phapecoctavirus* reference genomes used in comparative analyses

### `/annotation`
Pharokka v1.7.1 outputs (PHANOTATE gene prediction, PHROGs functional
annotation, tRNAscan-SE 2.0, MASH/INPHARED taxonomy).
- `F58V_pharokka.gbk` / `F63V_pharokka.gbk` — annotated genomes (GenBank format)
- `F58V_proteins.faa` / `F63V_proteins.faa` — predicted protein sequences
- `F58V_CDS_nt.ffn` / `F63V_CDS_nt.ffn` — CDS nucleotide sequences

#### `/annotation/InterProScan`
InterProScan 5.77 domain architecture analysis of all predicted proteins.
Applications run: Pfam, Gene3D, SUPERFAMILY, PANTHER, CDD, PRINTS, NCBIfam, 
ProSiteProfiles, SMART, PIRSF, MobiDBLite, Coils, TMHMM, Phobius, SignalP. 
Reclassification threshold: E-value < 1×10⁻⁵ consistent across ≥2 member 
databases.
- `F58V_interproscan.tsv` — 313 CDS (F58V)
- `F63V_interproscan.tsv` — 83 CDS (F63V)

#### `/annotation/HHpred`
Remote homology validation (HHpred, PDB_mmCIF70 + Pfam-A, March 2026) for
three proteins insufficiently characterized by InterProScan alone. HHpred
profile-vs-profile HMM alignment provides superior sensitivity for remote
homology and direct structural comparison against PDB.

- `HHpred_F63V_CDS0004_portal.hhr` — F63V CDS_0004 (425 aa). InterProScan
  assigned PF06381/IPR024459 (anti-CBASS Acb1-like). HHpred top hit:
  9KMH (portal protein; Prob=100%, E=1.3×10⁻⁴²). InterProScan assignment
  incorrect: PF06381 corresponds to a phage portal domain, not an anti-CBASS
  effector. CDS_0004 is annotated as portal protein in the final manuscript.

- `HHpred_F63V_CDS0056_tailfiber.hhr` — F63V CDS_0056 (844 aa, tail fiber).
  InterProScan detected pectin lyase fold (SSF51126) but not WcaM; catalytic
  vs. structural role ambiguous. HHpred top hit: 7VYV_B (depolymerase;
  Prob=99.9%, E=1.7×10⁻²¹), with multiple additional glycoside hydrolase and
  pectin lyase hits. Confirms the SSF51126 fold is catalytic.

- `HHpred_F58V_CDS0015_tailfiber.hhr` — F58V CDS_0015 (353 aa, secondary
  tail fiber). No domains detected by InterProScan. HHpred top hit: 5YVQ_A
  (tail fiber protein S, phage Mu; Prob=97.7%, E=1.9×10⁻⁴). Confirms tail
  fiber identity; receptor target remains undetermined.

### `/crispr`
CRISPR spacer screening with MinCED (via Pharokka pipeline).
- `pharokka_minced_spacersF58V.txt`
- `pharokka_minced_spacersF63V.txt`

### `/blast_results`
tBLASTn results of F58V CDS against 24 *Phapecoctavirus* reference genomes.
All TSV files follow the format: `qseqid sseqid pident length qlen slen
evalue bitscore`.
- `blast_all_cds.tsv` — full CDS conservation screen (313 CDS × 24 genomes)
- `blast_depol_all.tsv` — depolymerase module conservation specifically
- `blast_depolymerase.tsv` — CDS_0012 (predicted colanic acid depolymerase)
- Remaining TSV files — individual gene screens used for annotation
  validation and comparative analysis

### `/phylogenetics`
One subdirectory per gene. Each contains:
- `[gene].fna` — protein sequence alignment (MAFFT L-INS-i v7.505)
- `[gene]_codon.fna` — back-translated codon alignment
- `[gene].phy` — PHYLIP-format alignment for PAML input
- `[gene].tree` — guide tree (Newick) for PAML input
- `[gene].iqtree` — IQ-TREE log (model selection, likelihood scores)
- `[gene].treefile` — maximum-likelihood tree with UFBoot2 support values

Genes analyzed: depolymerase (CDS_0012), main_tail_fiber, TerL, baseplate_hub,
endolysin, portal_protein, tail_sheath, major_head_protein, methyltransferase,
DNA_helicase, DNA_primase.

Tree inference: IQ-TREE 2.3.6, GTR+G4 model, UFBoot2 (1000 replicates).
11 taxa: F58V + 10 representative *Phapecoctavirus* genomes spanning genus
diversity (MZ726793, OP172792, MN850598, JX561091, MH252123, MT944117,
OR352955, OM386656, KX664695, OR437326).

### `/selection_analysis`
PAML 4.10.7 codon-level selection analysis. One subdirectory per gene.
Each contains the control file (`.ctl`) and output (`.out`) for each model run.

- **All 11 genes**: M0 (one-ratio) model only
- **depolymerase and main_tail_fiber**: full site-model comparisons —
  M0, M1a (nearly neutral), M2a (positive selection), M7 (beta), M8 (beta+ω)

Key results:
- Depolymerase CDS_0012: ω = 0.020 (M0), M2a/M8 not significant → strong
  purifying selection, no evidence of positive selection
- Main tail fiber: ω = 0.093 (M0), M2a significant → evidence of positive
  selection at receptor-binding sites

#### `/structural`
AlphaFold3 structural modeling of the CDS_0012 C-terminal effector domain
(residues 344–1017), modeled as a C3-symmetric homotrimer (5 seeds).

Input:
- `F58V_depolymerase_AF3_input.fasta` — submitted sequence (674 aa × 3 chains)

Outputs (top-ranked model by ipTM):
- `fold_2026_03_18_19_32_model_1.cif` — predicted structure (CIF format)
- `fold_2026_03_18_19_32_full_data_1.json` — per-residue pLDDT and PAE matrix
- `fold_2026_03_18_19_32_summary_confidences_1.json` — ipTM, pTM, summary scores
- `fold_2026_03_18_19_32_template_hit_1_chains_c_a_b.cif` — template used
- `fold_2026_03_18_19_32_paired_msa_chains_c_a_b.a3m` — paired MSA
- `fold_2026_03_18_19_32_unpaired_msa_chains_c_a_b.a3m` — unpaired MSA

Key quality metrics: ipTM = 0.90. Structural superposition against PDB 6E0V
(colanidase gp150, phage Phi92) in UCSF ChimeraX 1.11.1: RMSD = 0.905 Å.
Candidate catalytic residues: D554, Y562, W578, D656 (pLDDT > 95).

AlphaFold3 modeling was performed on the public server (alphafoldserver.com).
Top-ranked model: ipTM = 0.90. Structural superposition against PDB 6E0V
(colanidase gp150, phage Phi92) in UCSF ChimeraX 1.11.1: RMSD = 0.905 Å.
Four candidate catalytic residues identified: D554, Y562, W578, D656
(pLDDT > 95). AlphaFold3 model output files (.cif, confidence JSON) are
PENDIENTE deposited at Zenodo, DOI:.

---

## Software versions

| Tool | Version | Reference |
|---|---|---|
| Pharokka | 1.7.1 | Bouras et al. 2023 |
| MAFFT | 7.505 | Katoh & Standley 2013 |
| BLAST+ | 2.12.0 | Camacho et al. 2009 |
| IQ-TREE | 2.3.6 | Nguyen et al. 2015 |
| trimAl | 1.5.rev1 | Capella-Gutierrez et al. 2009 |
| UFBoot2 | — | Hoang et al. 2018 |
| PAML (codeml) | 4.10.7 | Yang 2007 |
| VIRIDIC | 1.1 | Moraru et al. 2020 |
| AlphaFold3 | server | Abramson et al. 2024 |
| ChimeraX | 1.11.1 | Meng et al. 2023 |
| InterProScan | 5.77-108.0 | Jones et al. 2014 |
| FoldSeek | server | van Kempen et al. 2023 |
| HHpred | PDB_mmCIF70_20_Feb Pfam-A_v38.2 | Zimmermann et al. 2018 |
| clinker | 0.0.32 | Gilchrist & Chooi 2021 |
| CARD RGI | RGI 6.0.5, CARD 4.0.1 | Alcock et al. 2023 |
| VFDB | server | Liu et al. 2022 |
| MinCED | — | via Pharokka |

---

## Contact

Laura Ribes-Martinez — lauraribes@outlook.com
