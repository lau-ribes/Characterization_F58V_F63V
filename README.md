# Genomic and in silico analyses: bacteriophages F58V and F63V

## Overview

This repository holds the input files and primary outputs of the bioinformatic
analyses supporting the manuscript:

> **Depolymerase architecture and activity against preformed biofilms in two phages targeting ESBL-producing *Escherichia coli***
> Ribes-Martinez et al., *Scientific Reports* (in revision, 2026)

Two lytic bacteriophages isolated from hospital wastewater (Hospital
Universitario Infanta Elena, Valdemoro, Madrid) were characterized against
ESBL-producing *E. coli* clinical isolates: **F58V** (149,321 bp, GC 39.0%,
313 CDS; *Phapecoctavirus*, Straboviridae; contractile tail) and **F63V**
(45,748 bp, GC 44.1%, 83 CDS; Braunvirinae, Drexlerviridae; proposed new genus
*Aunosvirus* gen. nov.).

The central argument is that the two phages carry tail-associated
polysaccharide depolymerases with different predicted substrates, F58V CDS_0012
a colanic acid depolymerase and F63V CDS_0056 a structurally distinct enzyme
whose closest structural match is a *Klebsiella* K2 capsular depolymerase, and
that this difference parallels their divergent activity against preformed
biofilms. The genomic, structural, phylogenetic and evolutionary evidence
deposited here supports the architectural characterization; it does not
establish causality, and the manuscript treats depolymerase architecture as a
correlate rather than a demonstrated cause.

Genome sequences are in NCBI GenBank under **PZ445148** (F58V) and **PZ445149**
(F63V). Library preparation, sequencing and assembly were performed externally
and are described in the manuscript Methods; this repository covers the in
silico analyses only.

---

## Note on reference-genome sets

Several different reference sets are used in this work and should not be
conflated. Each analysis states its own set:

| Analysis | Set | Files |
|---|---|---|
| VIRIDIC intergenomic similarity, F58V (Figure 4c) | F58V + 24 *Phapecoctavirus* genomes (25 sequences), the set identified in Supplementary Table 5 | `/viridic` |
| VIRIDIC intergenomic similarity, F58V (extended run) | F58V + 30 *Phapecoctavirus* genomes (31 sequences), which covers every ratified genus exemplar | `/viridic` |
| VIRIDIC intergenomic similarity, F63V (Figure 4d) | F63V + 22 Braunvirinae and outgroup genomes (23 sequences) | `/viridic` |
| Whole-genome ML tree, F58V (Figure 4a) | F58V + the same 24 *Phapecoctavirus* genomes used in Figure 4c (25 sequences) | `/phylogenetics/F58V_whole_genome` |
| Whole-genome ML tree, F63V (Figure 4b) | F63V + 24 Braunvirinae and outgroup genomes (25 sequences), the Figure 4d set plus MK907243 and MK907257 | `/phylogenetics/F63V_whole_genome` |
| tBLASTn conservation screen and genomic null distribution | 46 *Phapecoctavirus* genomes | `/blast_results`, `genomes/all_phapecoctavirus_accessions.txt` |
| Codon-level selection and single-gene trees, F58V | F58V + 10 *Phapecoctavirus* genomes (nine for TerL) | `/phylogenetics`, `/selection_analysis` |

One further figure in the manuscript has no input set of its own: a tBLASTn
search of CDS_0012 returned 65 *Phapecoctavirus* subject genomes. Thirteen of those are RefSeq copies of
genomes already present under their GenBank accession, leaving **52
non-redundant genomes**, of which 46 carry the protein at or above 99.0%
identity and 23 at or above 99.7%. Six fall below 99.0%: five (LC739533,
LC739534, LC739535, LC739537, LC739539) carry a more divergent copy at 84.8 to
85.1%, and OR067838 falls just short of the cut at 98.7%. The raw output is
`blast_results/conservation/tblastn_CDS_0012_65subjects.tsv` (identical to
`blast_results/blast_depolymerase.tsv`; the same hits in 8-column format are in
`blast_results/blast_depol_all.tsv`), and the duplicate-detection table is
`blast_results/conservation/duplicados_en_tblastn_65subjects.tsv`. The BLASTp
search of CDS_0012 against the NCBI non-redundant protein database is a separate
output, `blast_results/conservation/blastp_CDS_0012_vs_nr.tsv`.

---

## Repository structure

### `/genomes`

- `F58V.fasta`: F58V complete genome, 149,321 bp, GC 39.0%.
- `F63V.fasta`: F63V complete genome, 45,748 bp, GC 44.1%.
- `all_phapecoctavirus_accessions.txt`: the 46 GenBank accessions of the
  tBLASTn reference set. These are exactly the 46 subject genomes present in
  `blast_results/blast_all_cds.tsv`.

### `/annotation`

Pharokka v1.3.2 outputs (PHANOTATE gene calling, PHROGs functional assignment,
tRNAscan-SE, mash-based placement against INPHARED).

- `F58V_pharokka.gbk` / `F63V_pharokka.gbk`: annotated genomes (GenBank).
- `F58V_proteins.faa` / `F63V_proteins.faa`: predicted proteins.
- `F58V_CDS_nt.ffn` / `F63V_CDS_nt.ffn`: CDS nucleotide sequences.

Of the 313 F58V CDS, 276 (88.2%) matched a PHROG entry and 78 (24.9%) carried a
specific functional product at the Pharokka stage. For F63V the corresponding
figures are 83 CDS and 35 (42.2%).

#### `/annotation/phold`

Structure-aware reannotation with phold v1.2.0 (Foldseek against the PHROG
structural database plus ProstT5 embeddings) on top of the Pharokka output.
Specific functional assignments rise from 78 (24.9%) to **102 (32.6%)** in F58V
and from 35 (42.2%) to **42 (50.6%)** in F63V; 211 F58V CDS retain an unknown
function. No CDS in either genome was assigned to the PHROG category
"integration and excision", consistent with the predicted virulent lifestyle.

- `F58V-phold_per_cds_predictions.tsv`
- `F63V-phold_per_cds_predictions.tsv`

#### `/annotation/hmmer_pfam`

HMMER v3.4 searches of both proteomes against Pfam release 37.2, used to
validate the PHROG and phold assignments. `*_pfam.tbl` holds the per-sequence
hits and `*_pfam.domtbl` the per-domain hits, in the standard HMMER tabular
formats.

#### `/annotation/interproScan`

InterProScan 5.77-108.0 domain analysis of all predicted proteins
(`F58V_interproscan.tsv`, `F63V_interproscan.tsv`). These files are the source
of the domain boundaries drawn in Figure 2.

#### `/annotation/HHpred`

Profile-profile searches (HHpred, Tübingen toolkit; PDB_mmCIF70 and Pfam-A) for
the proteins whose assignment required remote-homology evidence. Top hit of
each report:

| File | Protein | Top hit | Probability | E-value | Coverage |
|---|---|---|---|---|---|
| `HHpred_F63V_CDS0004_portal.hhr` | F63V CDS_0004 | 9KMH, portal protein | 100.0 | 1.3e-42 | 396 cols, query 1-400 |
| `HHpred_F63V_CDS0056_tailfiber.hhr` | F63V CDS_0056 | 7VYV_B, K2 capsular depolymerase | 99.9 | 1.7e-21 | 411 cols, query 1-485 |
| `HHpred_F58V_CDS0015_tailfiber.hhr` | F58V CDS_0015 | 5YVQ_A, tail fiber protein S of phage Mu | 97.7 | 1.9e-04 | 68 cols, query 278-353 |
| `HHpred_F58V_CDS0019_endosialidase.hhr` | F58V CDS_0019 | 4HIZ_A, endosialidase of *Escherichia* phage phi92 | 100.0 | 4.2e-97 | 682 cols, query 75-759 |
| `HHpred_F58V_CDS0018_DUF2644.hhr` | F58V CDS_0018 | PF10841, DUF2644 | 99.9 | 2.0e-25 | 59 cols, query 14-73 |
| `HHpred_F58V_CDS0218.hhr` | F58V CDS_0218 | PF09964, DUF2198 | 89.8 | 3.2 | 42 cols, query 13-54 |

The CDS_0004 report supersedes the InterProScan anti-CBASS prediction
(Supplementary Note 1). The CDS_0015 report places the Mu tail fiber match on
the C-terminal 76 residues only, with no enzymatic fold anywhere in the
protein. The CDS_0018 and CDS_0218 reports are the two holin candidates of
Supplementary Note 2; the CDS_0218 hit is weak and neither prediction has been
tested experimentally.

### `/crispr`

CRISPR spacer screening with MinCED, run through Pharokka. MinCED detected no
CRISPR array and no spacer in either genome, so both output files are empty.

### `/blast_results`

tBLASTn searches of F58V proteins against *Phapecoctavirus* genomes.
`blast_all_cds.tsv` covers all F58V CDS against the 46-genome reference set,
with columns `qseqid sseqid pident length qlen slen evalue bitscore`. The
per-gene files (`blast_<gene>.tsv`) are searches of single proteins that
returned between 52 and 65 subject genomes. Most use the 14-column format
`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue
bitscore qlen slen`; `blast_depol_all.tsv` and `blast_methyltransferase.tsv`
use the 8-column format above.

`blast_all_cds.tsv` is the genome-wide screen from which the null distribution
of protein identity was built, over the 165 core CDS with an ortholog in every
one of the 46 reference genomes. The colanic acid depolymerase sits at the 50th
percentile of that distribution.

#### `/blast_results/conservation`

The computation behind the conservation claim for CDS_0012, deposited so that
the figures in the Results can be reproduced step by step.

- `blastp_CDS_0012_vs_nr.tsv`: BLASTp of CDS_0012 against the NCBI
  non-redundant protein database.
- `tblastn_CDS_0012_65subjects.tsv`: the tBLASTn search of CDS_0012 that
  returned the 65 *Phapecoctavirus* subject genomes.
- `duplicados_en_tblastn_65subjects.tsv`: the genome-length matching that
  identifies RefSeq copies of genomes already present under their GenBank
  accession, which reduces those 65 to the 52 non-redundant genomes reported.
- `conservacion_CDS_0012_46genomas.tsv`: per-genome identity and summed-HSP
  query coverage across the 46-genome reference set.
- `calcular_conservacion_CDS_0012.py`, `null_distribution_metodo.py`: the two
  scripts that produce the tables above and the null distribution. Usage is in
  each script's docstring; both use only the Python standard library.
- `null_distribution_summary.txt`, `null_distribution_per_cds.tsv`: the output
  of `null_distribution_metodo.py` run from this directory on
  `../blast_all_cds.tsv` and `../../annotation/F58V_pharokka.gbk`.

### `/phylogenetics`

Single-gene maximum-likelihood trees, one subdirectory per gene. Each holds the
ortholog CDS sequences (`.fna`), the codon alignment (`_codon.fna`), the same
alignment in PHYLIP format with PAML-compatible labels (`.phy`), the IQ-TREE
report (`.iqtree`), the ML tree (`.treefile`) and the same tree with
PAML-compatible labels, used as the codeml guide tree (`.tree`). The protein
alignments from which the codon alignments were back-translated are not
deposited. Eleven genes: depolymerase
(CDS_0012), main_tail_fiber, TerL, baseplate_hub, endolysin, portal_protein,
tail_sheath, major_head_protein, methyltransferase, DNA_helicase, DNA_primase.
Alignments contain F58V plus ten *Phapecoctavirus* genomes, except TerL with
nine. Inferred with IQ-TREE v2.3.6 under GTR+F+G4, selected by ModelFinder for
all eleven genes.

`/phylogenetics/F63V_single_gene`: TerL and major capsid protein trees for
F63V and its candidate relatives (Supplementary Figure 2). IQ-TREE v2.3.6,
1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates. Seven taxa: F63V, ES10,
AN_ECEAS, IMM-001, Loudonvirus DTL and the two Braunvirinae outgroups NC_047810
and NC_047812. The files in the top level are the fixed LG+F+G4 runs used for
Supplementary Figure 2.

`modelfinder_2026-09-23/` holds the model-selection reanalysis, described in the
`README.md` of that directory. `*_mfp.*` are the ModelFinder runs reported in the
Methods, which select LG+R2 for TerL and Q.pfam+I for MCP; `*_rep.*` reproduce
the original fixed-model analysis, matching its log-likelihood to the fourth
decimal. Both recover the same topology as the fixed-model trees. TerL groups F63V with ES10, AN_ECEAS and IMM-001 (MF630922); MCP
instead groups F63V with *Loudonvirus* DTL, the single-gene discordance
discussed in the manuscript.

`/phylogenetics/F58V_whole_genome` and `/phylogenetics/F63V_whole_genome`: the
two whole-genome maximum-likelihood analyses of Figure 4a and Figure 4b. Each
directory holds the IQ-TREE report (`.iqtree`), the run log, the consensus tree
(`.contree`), the ML tree (`.tree`) and the rendered figure as PNG, PDF and SVG.

Both were run with IQ-TREE v1.6.12 as
`iqtree -s <alignment> -m TESTNEW -bb 1000 -alrt 1000`, that is ModelFinder
selection with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates. The
*Phapecoctavirus* alignment is 25 sequences and 286,370 nucleotide sites, with
TIM+F+I+G4 selected by BIC; the Braunvirinae alignment is 25 sequences and
107,300 sites, with GTR+F+R3. Node labels in the figures are SH-aLRT support /
ultrafast bootstrap support.

The query genomes appear in these files under their internal working labels,
`58V-B-fago` for F58V and `63V-K-fago-r` for F63V. The input alignments
themselves are not deposited; their dimensions, composition and model selection
are recorded in the `.iqtree` reports.

### `/viridic`

VIRIDIC v1.1 intergenomic similarity (BLASTn; word size 7, reward 2, penalty
-3, filter on; R v4.5.2, BLAST+ v2.17). Species and genus demarcation follow
the ICTV thresholds for *Caudoviricetes*, 95% and 70%.

**F58V, Figure 4c.** `VIRIDIC_F58V_Fig4c_input_25genomes.fasta`,
`VIRIDIC_F58V_Fig4c_sim-dist_table.tsv`, `VIRIDIC_F58V_Fig4c_cluster_table.tsv`
and `VIRIDIC_F58V_Fig4c_heatmap.pdf` are the run shown in the manuscript: F58V
and the 24 reference genomes flagged in Supplementary Table 5. The lowest
pairwise value in the matrix is 83.997% and all 25 genomes fall in a single
genus cluster. F58V reaches 94.0% against the genus exemplar
phAPEC8 (JX561091 / NC_020079), below the 95% species threshold, and 95.4%
against the unclassified genome PNJ1809-36 (MT944117), marginally above it.

**F58V, extended run.** `VIRIDIC_F58V_extended31_input.fasta`,
`VIRIDIC_F58V_extended31_sim-dist_table.tsv` and
`VIRIDIC_F58V_extended31_matrix.pdf` are a wider comparison over 31 sequences,
which includes every ratified *Phapecoctavirus* species exemplar. It is
deposited because it is what supports the statement in the Results that F58V
was compared against all genus exemplars: the highest similarity to F58V across
them is 94.0%, and the exemplars absent from Figure 4c all fall lower, between
89.6% and 92.6%.

**F63V, Figure 4d.** `VIRIDIC_F63V_Fig4d_input_23genomes.fasta`,
`VIRIDIC_F63V_Fig4d_sim-dist_table.tsv` and `VIRIDIC_F63V_Fig4d_matrix.pdf` are
the run shown in the manuscript: F63V and 22 Braunvirinae and outgroup genomes.
F63V shares 78.8% and 78.7% with AN_ECEAS (PP944860) and ES10 (OR466728), which
share 99.8% with each other. Every other genome in the matrix falls at or below
67.0%, well under the 70% genus demarcation threshold, the most distant being
MF630922 at 54.3%. The gap between the intra-clade similarities and the next
nearest genome is the basis of the genus-level proposal.

### `/selection_analysis`

Codon-level selection analysis with PAML codeml v4.9j, one subdirectory per
gene, each with the control file and the codeml output. The control files read
`paml_<gene>.phy` and `paml_<gene>.tree`, which are
`phylogenetics/<gene>/<gene>.phy` and `<gene>.tree` under the names codeml
expects.

- All eleven F58V genes: one-ratio model M0.
- Depolymerase (CDS_0012) and main tail fiber (CDS_0011): the site models
  M1a/M2a and M7/M8 in addition.

Gene-wide ω under M0 spans two orders of magnitude, from DNA primase at 0.002
to endolysin at 0.182. The depolymerase is at **ω = 0.020**, sixth of eleven
and within the range of essential structural and replicative genes; the main
tail fiber is at 0.093. In the depolymerase, M1a/M2a does not support a
positively selected site class (2ΔlnL = 2.03, p = 0.36). M7/M8 is nominally
significant (2ΔlnL = 15.89, p < 0.001), with 0.4% of sites at ω = 2.99, but no
site reaches a Bayes empirical Bayes posterior above 0.95; the manuscript
discusses this result. Absolute ω values depend on the
composition of the reference set, so comparisons should be read as
within-analysis contrasts between genes.

`/selection_analysis/M3_per_site`: M3 (K = 3) NEB per-site dN/dS profiles
(Supplementary Figure 1): `*_M3.out`, `*_M3_rst_persite.txt`,
`Suppl_Fig_dNdS_persite.{pdf,svg}`. In the depolymerase, two of the three site
classes are under strong purifying selection (ω = 0.011 and 0.012) and the
third, 0.36% of sites, has ω = 3.03. The four predicted catalytic residues
(D554, Y562, W578, D656) have a posterior probability of at most 0.002 for that
third class; the only two sites assigned to it (Y383, D434; posterior > 0.99)
lie outside the catalytic pocket. The M3 control file is not deposited; the M0
control file with `NSsites = 3` and `ncatG = 3` reproduces the run (lnL
-5742.779 against -5742.780 in `depolymerase_M3.out`).

### `/structural`

AlphaFold3 model of the CDS_0012 C-terminal effector domain (residues
344-1017) as a C3-symmetric homotrimer, ipTM 0.90, pTM 0.91. Full-length
monomer prediction failed to converge, which is why the trimeric effector
domain was modeled instead. Superposition against PDB 6E0V, the apo gp150
colanidase tail spike of *Escherichia* phage phi92, with UCSF ChimeraX v1.11.1
(matchmaker, default settings): Cα RMSD 0.905 Å over the 354 of 595 pairs
retained after iterative pruning, 0.732 Å over 222 of 264 pairs for the
catalytic domain alone, against a crystal-to-crystal control (6E0V chain A vs.
6E0W chain A) of 0.248 Å over all 602 pairs.

- `F58V_depolymerase_AlphaFold.fasta`: input sequence.
- `fold_2026_03_18_19_32_model_1.cif`: predicted structure.
- `fold_2026_03_18_19_32_summary_confidences_1.json`,
  `fold_2026_03_18_19_32_full_data_1.json`: confidence metrics, including the
  per-residue pLDDT values quoted in the manuscript.
- `fold_2026_03_18_19_32_paired_msa_chains_c_a_b.a3m`,
  `..._unpaired_msa_chains_c_a_b.a3m`, `..._template_hit_1_chains_c_a_b.cif`:
  MSAs and template hit returned by the server.

`/structural/foldseek` holds the structural similarity search that identified
6E0V as the closest structural homolog: `foldseek_CDS_0012_vs_PDB100_*.tsv` is
the result table and `foldseek_run_info.txt` records the query chain, the
database build date, the Foldseek commit and the exact command.

### `/figures`

- `figura2.R`: the script that draws the Figure 2 domain architecture diagram,
  ggplot2 v3.4.4 under R v4.3.3.
- `dominios.tsv`: the domain and fold boundaries it draws, taken from the
  InterProScan and HHpred outputs in `/annotation`.

---

## Software versions

| Tool | Version | Reference |
|---|---|---|
| Pharokka | 1.3.2 | Bouras et al. 2023 |
| PHANOTATE | via Pharokka | McNair et al. 2019 |
| phold | 1.2.0 | Bouras et al. |
| Foldseek | via phold | van Kempen et al. 2024 |
| ProstT5 | via phold | Heinzinger et al. 2023 |
| HMMER | 3.4 (Pfam 37.2) | Eddy 2011 |
| InterProScan | 5.77-108.0 | Jones et al. 2014 |
| HHpred | PDB_mmCIF70 + Pfam-A | Zimmermann et al. 2018 |
| BLAST+ | 2.17 | Camacho et al. 2009 |
| MAFFT | 7.520 | Katoh & Standley 2013 |
| PAL2NAL | 14 | Suyama et al. 2006 |
| trimAl | - | Capella-Gutierrez et al. 2009 |
| IQ-TREE (single-gene) | 2.3.6 | Minh et al. 2020 |
| IQ-TREE (whole-genome) | 1.6.12 | Nguyen et al. 2015 |
| ModelFinder | via IQ-TREE | Kalyaanamoorthy et al. 2017 |
| UFBoot2 / SH-aLRT | via IQ-TREE | Hoang et al. 2018 |
| PAML (codeml) | 4.9j | Yang 2007 |
| VIRIDIC | 1.1 | Moraru et al. 2020 |
| AlphaFold3 | AlphaFold Server | Abramson et al. 2024 |
| ChimeraX | 1.11.1 | Meng et al. 2023 |
| BACPHLIP | 0.9.6 | Hockenberry & Wilke 2021 |
| CARD RGI | RGI 6.0.5, CARD 4.0.1 | Alcock et al. 2023 |
| VFDB | core dataset | Liu et al. 2022 |
| MinCED | via Pharokka | Bland et al. 2007 |
| ARAGORN | via Pharokka | Laslett & Canback 2004 |
| R | 4.5.2 (VIRIDIC), 4.3.3 (figures) | R Core Team |
| ggplot2 | 3.4.4 | Wickham 2016 |

---

## Data availability

Genome sequences: NCBI GenBank **PZ445148.1** (F58V), **PZ445149.1** (F63V). All
other data supporting the findings of the study are in the manuscript, its
supplementary information, and this repository. Phenotypic raw data, including
the biofilm, efficiency of plating and lysis-curve measurements, are tabulated
in the supplementary information.

## Contact

Laura Ribes-Martinez, lauraribes@outlook.com
