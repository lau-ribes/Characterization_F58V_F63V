# Genomic and in silico analyses — F58V and F63V bacteriophages

## Overview

This repository contains the bioinformatic input and output files supporting
the manuscript:

> **Phage depolymerase specificity shapes treatment outcomes for preformed
> biofilms in ESBL-producing *Escherichia coli***
> Ribes-Martinez et al. — submitted to *Scientific Reports*, 2026

Two lytic bacteriophages isolated from hospital wastewater (Hospital
Universitario Infanta Elena, Valdemoro, Madrid) were characterized: F58V
(*Phapecoctavirus*, Straboviridae, ~149 kb myovirus) and F63V (proposed new
genus *Aunosvirus* gen. nov., Drexlerviridae/Braunvirinae, ~46 kb siphovirus),
both targeting ESBL-producing *E. coli* clinical isolates.

The central argument is that the distinct substrate specificities of the
tail-associated polysaccharide depolymerases of F58V (a predicted colanic acid
depolymerase, CDS_0012) and F63V (a predicted K2-type capsule depolymerase,
CDS_0056) are associated with their differential activity against preformed
biofilms, supported by phenotypic, structural, evolutionary and comparative
genomic evidence.

> **Note on reference-genome sets.** Three different genome sets are used in
> this work and should not be confused:
> - **PAML codon selection** (F58V genes): F58V + 10 representative
>   *Phapecoctavirus* genomes.
> - **tBLASTn conservation screen**: F58V CDS vs. 24 *Phapecoctavirus* genomes.
> - **VIRIDIC intergenomic similarity**: run separately for F58V
>   (*Phapecoctavirus* set) and for F63V (Braunvirinae set).

---

## Repository structure

### `/genomes`
Assembled genome sequences in FASTA format.
- `F58V.fasta` — F58V complete genome (149,321 bp, GC 39.0%)
- `F63V.fasta` — F63V complete genome (45,748 bp, GC 44.1%)
- `all_phapecoctavirus_accessions.txt` — GenBank accessions for the
  *Phapecoctavirus* reference genomes used in comparative analyses.

### `/annotation`
Pharokka v1.7.1 outputs (PHANOTATE, PHROGs, tRNAscan-SE 2.0, MASH/INPHARED).
- `F58V_pharokka.gbk` / `F63V_pharokka.gbk` — annotated genomes (GenBank).
- `F58V_proteins.faa` / `F63V_proteins.faa` — predicted proteins.
- `F58V_CDS_nt.ffn` / `F63V_CDS_nt.ffn` — CDS nucleotide sequences.

#### `/annotation/phold`
Structure-aware functional reannotation with phold v1.2.0 (Foldseek + ProstT5)
on top of the Pharokka output. Increases the proportion of CDSs with a
specific functional assignment from 78 (24.9%) to 102 (32.6%) in F58V, and
from 35 (42.2%) to 42 (50.6%) in F63V. Zero CDSs were assigned to the PHROG
category "integration and excision" in either genome, supporting a predicted virulent lifestyle.
- `F58V-phold_per_cds_predictions.tsv`
- `F63V-phold_per_cds_predictions.tsv`

#### `/annotation/interproScan`
InterProScan 5.77-108.0 domain analysis of all predicted proteins.

#### `/annotation/HHpred`
Remote-homology validation (HHpred, PDB_mmCIF70 + Pfam-A) of three proteins:
- `HHpred_F63V_CDS0004_portal.hhr` — CDS_0004 reclassified as portal protein
  (9KMH; Prob 100%, E = 1.3e-42); supersedes the InterProScan anti-CBASS call.
- `HHpred_F63V_CDS0056_tailfiber.hhr` — CDS_0056 K2-type capsule depolymerase
  (7VYV_B; Prob 99.9%, E = 1.7e-21).
- `HHpred_F58V_CDS0015_tailfiber.hhr` — CDS_0015 secondary tail fiber
  (5YVQ_A; Prob 97.7%).

### `/crispr`
CRISPR spacer screening with MinCED (via Pharokka).

### `/blast_results`
tBLASTn screens of F58V CDS against *Phapecoctavirus* reference genomes.
Format: `qseqid sseqid pident length qlen slen evalue bitscore`.

### `/phylogenetics`
F58V gene trees (one subdirectory per gene). Each contains the protein
alignment (`.fna`, MAFFT L-INS-i), codon alignment (`_codon.fna`), PHYLIP
alignment (`.phy`), guide tree (`.tree`), IQ-TREE log (`.iqtree`) and ML tree
(`.treefile`). Genes: depolymerase (CDS_0012), main_tail_fiber, TerL,
baseplate_hub, endolysin, portal_protein, tail_sheath, major_head_protein,
methyltransferase, DNA_helicase, DNA_primase. F58V + 10 *Phapecoctavirus*
genomes (MZ726793, OP172792, MN850598, JX561091, MH252123, MT944117, OR352955,
OM386656, KX664695, OR437326). IQ-TREE2 2.3.6, GTR+F+G4.

`/phylogenetics/F63V_single_gene` — TerL and MCP single-gene ML trees
(Supplementary Figure 2). IQ-TREE2 2.2.0, LG+F+G4, 1000 UFBoot + 1000 SH-aLRT.
TerL groups F63V with ES10/AN_ECEAS/IMM-001; MCP groups F63V with Loudonvirus
DTL (94.2% MCP amino-acid identity).

`/phylogenetics/F63V_whole_genome` — whole-genome tree of F63V and 22
Braunvirinae genomes (`F63V_whole_genome_input_23genomes.fasta`,
`F63V_whole_genome_tree.png`).

### `/viridic`
VIRIDIC v1.1 intergenomic-similarity analyses (BLASTn; R v4.5.2, BLAST+ 2.17).
- `VIRIDIC_F63V_sim-dist_table.tsv` + `VIRIDIC_F63V_input.fasta` +
  `VIRIDIC_F63V_matrix.pdf` — F63V vs. Braunvirinae set. F63V shares
  78.7% / 78.8% with ES10 / AN_ECEAS (99.8% between them), and <= 62% with any
  classified Braunvirinae genus exemplar (nearest: *Guelphvirus*), a 16-point
  gap supporting genus-level demarcation under ICTV criteria.
- `VIRIDIC_F58V_matrix.pdf` (+ `VIRIDIC_F58V_input_31genomes.fasta`).

### `/selection_analysis`
PAML 4.10.7 codon-level selection analysis (one subdirectory per gene).
- All 11 F58V genes: M0 (one-ratio) model.
- depolymerase (CDS_0012) and main_tail_fiber: full site models M0/M1a/M2a/M7/M8.

Key results: depolymerase omega = 0.020 (M0; 0.011-0.020 across reference-set
iterations), M1a/M2a not significant. Main tail fiber omega = 0.093.

`/selection_analysis/M3_per_site` — M3 (K = 3) NEB per-site dN/dS profiles
(Supplementary Figure 1). `*_M3.out`, `*_M3_rst_persite.txt`,
`Suppl_Fig_dNdS_persite.{pdf,svg}`. Catalytic residues (D554, Y562, W578, D656)
fall in the strongly constrained class; the two NEB-positive sites (Y383,
D434) lie outside the catalytic pocket.

### `/structural`
AlphaFold3 model of CDS_0012 C-terminal effector domain (residues 344-1017),
C3-symmetric homotrimer. ipTM = 0.90, pTM = 0.91. Superposition vs. PDB 6E0V
(ChimeraX 1.11.1): RMSD 0.905 A full chain, 0.732 A catalytic domain, crystal
control 0.248 A.

---

## Software versions

| Tool | Version | Reference |
|---|---|---|
| Pharokka | 1.7.1 | Bouras et al. 2023 |
| phold | 1.2.0 | Bouras et al. (in prep) |
| BACPHLIP | 0.9.6 | Hockenberry & Wilke 2021 |
| MAFFT | 7.520 | Katoh & Standley 2013 |
| BLAST+ | 2.17 | Camacho et al. 2009 |
| IQ-TREE | 2.3.6 | Minh et al. 2020 |
| trimAl | 1.5.rev1 | Capella-Gutierrez et al. 2009 |
| UFBoot2 / SH-aLRT | - | Hoang et al. 2018 |
| PAML (codeml) | 4.10.7 | Yang 2007 |
| VIRIDIC | 1.1 | Moraru et al. 2020 |
| AlphaFold3 | server | Abramson et al. 2024 |
| ChimeraX | 1.11.1 | Meng et al. 2023 |
| InterProScan | 5.77-108.0 | Jones et al. 2014 |
| HHpred | PDB_mmCIF70 + Pfam-A v38.2 | Zimmermann et al. 2018 |
| CARD RGI | RGI 6.0.5, CARD 4.0.1 | Alcock et al. 2023 |
| VFDB | server | Liu et al. 2022 |
| MinCED | - | via Pharokka |

---

## Contact

Laura Ribes-Martinez — lauraribes@outlook.com
