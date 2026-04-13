# Characterization\_F58V\_F63V



\# Genomic and in silico analyses — F58V and F63V bacteriophages



\## Overview



This repository contains all bioinformatic input and output files supporting

the manuscript:



> \\\*\\\*A conserved surface polysaccharide degradation module in Phapecoctavirus

> is associated with differential anti-biofilm activity against

> ESBL-producing \\\*Escherichia coli\\\*\\\*\\\*

> AUTORES. — intended to be submitted to \\\*Communications Biology\\\*, 2026



Two lytic bacteriophages isolated from hospital wastewater (Hospital

Universitario Infanta Elena, Valdemoro, Madrid) were characterized: F58V

(Phapecoctavirus, Straboviridae, \~149 kb myovirus) and F63V (novel genus,

Drexlerviridae, \~46 kb siphovirus), both targeting ESBL-producing *E. coli*

clinical isolates.



The central argument is that a predicted tripartite colanic acid depolymerase

module in F58V (CDS\_0012, absent in F63V) is associated with superior

activity against preformed biofilms, supported by phenotypic, structural,

and evolutionary evidence.



\---



\## Repository structure



\### `/genomes`

Assembled genome sequences in FASTA format.

\- `F58V.fasta` — F58V complete genome (149,321 bp, GC 39.0%)

\- `F63V.fasta` — F63V complete genome (45,748 bp, GC 44.1%)

\- `all\\\_phapecoctavirus\\\_accessions.txt` — GenBank accession numbers for the

&#x20; 24 \*Phapecoctavirus\* reference genomes used in comparative analyses



\### `/annotation`

Pharokka v1.3.2 outputs (PHANOTATE gene prediction, PHROGs functional

annotation, tRNAscan-SE 2.0, MASH/INPHARED taxonomy).

\- `F58V\\\_pharokka.gbk` / `F63V\\\_pharokka.gbk` — annotated genomes (GenBank format)

\- `F58V\\\_proteins.faa` / `F63V\\\_proteins.faa` — predicted protein sequences

\- `F58V\\\_CDS\\\_nt.ffn` / `F63V\\\_CDS\\\_nt.ffn` — CDS nucleotide sequences



\#### `/annotation/InterProScan`

InterProScan 5.64 domain architecture analysis of all predicted proteins.

Applications run: Pfam, Gene3D, SUPERFAMILY, PANTHER, CDD, PRINTS, NCBIfam, 

ProSiteProfiles, SMART, PIRSF, MobiDBLite, Coils, TMHMM, Phobius, SignalP. 

Reclassification threshold: E-value < 1×10⁻⁵ consistent across ≥2 member databases.

\- `F58V\\\_interproscan.tsv` — 313 CDS (F58V)

\- `F63V\\\_interproscan.tsv` — 83 CDS (F63V)



\#### `/annotation/HHpred`

Remote homology validation (HHpred, PDB\_mmCIF70 + Pfam-A, March 2026) for

three proteins insufficiently characterized by InterProScan alone. HHpred

profile-vs-profile HMM alignment provides superior sensitivity for remote

homology and direct structural comparison against PDB.



\- `HHpred\\\_F63V\\\_CDS0004\\\_portal.hhr` — F63V CDS\_0004 (425 aa). InterProScan

&#x20; assigned PF06381/IPR024459 (anti-CBASS Acb1-like). HHpred top hit:

&#x20; 9KMH (portal protein; Prob=100%, E=1.3×10⁻⁴²). InterProScan assignment

&#x20; incorrect: PF06381 corresponds to a phage portal domain, not an anti-CBASS

&#x20; effector. CDS\_0004 is annotated as portal protein in the final manuscript.



\- `HHpred\\\_F63V\\\_CDS0056\\\_tailfiber.hhr` — F63V CDS\_0056 (844 aa, tail fiber).

&#x20; InterProScan detected pectin lyase fold (SSF51126) but not WcaM; catalytic

&#x20; vs. structural role ambiguous. HHpred top hit: 7VYV\_B (depolymerase;

&#x20; Prob=99.9%, E=1.7×10⁻²¹), with multiple additional glycoside hydrolase and

&#x20; pectin lyase hits. Confirms the SSF51126 fold is catalytic.



\- `HHpred\\\_F58V\\\_CDS0015\\\_tailfiber.hhr` — F58V CDS\_0015 (353 aa, secondary

&#x20; tail fiber). No domains detected by InterProScan. HHpred top hit: 5YVQ\_A

&#x20; (tail fiber protein S, phage Mu; Prob=97.7%, E=1.9×10⁻⁴). Confirms tail

&#x20; fiber identity; receptor target remains undetermined.



\### `/crispr`

CRISPR spacer screening with MinCED (via Pharokka pipeline).

\- `pharokka\\\_minced\\\_spacersF58V.txt`

\- `pharokka\\\_minced\\\_spacersF63V.txt`



\### `/blast\\\_results`

tBLASTn results of F58V CDS against 24 \*Phapecoctavirus\* reference genomes.

All TSV files follow the format: `qseqid sseqid pident length qlen slen

evalue bitscore`.

\- `blast\\\_all\\\_cds.tsv` — full CDS conservation screen (313 CDS × 24 genomes)

\- `blast\\\_depol\\\_all.tsv` — depolymerase module conservation specifically

\- `blast\\\_depolymerase.tsv` — CDS\_0012 (predicted colanic acid depolymerase)

\- Remaining TSV files — individual gene screens used for annotation

&#x20; validation and comparative analysis



\### `/phylogenetics`

One subdirectory per gene. Each contains:

\- `\\\[gene].fna` — protein sequence alignment (MAFFT L-INS-i v7.505)

\- `\\\[gene]\\\_codon.fna` — back-translated codon alignment

\- `\\\[gene].phy` — PHYLIP-format alignment for PAML input

\- `\\\[gene].tree` — guide tree (Newick) for PAML input

\- `\\\[gene].iqtree` — IQ-TREE log (model selection, likelihood scores)

\- `\\\[gene].treefile` — maximum-likelihood tree with UFBoot2 support values



Genes analyzed: depolymerase (CDS\_0012), main\_tail\_fiber, TerL, baseplate\_hub,

endolysin, portal\_protein, tail\_sheath, major\_head\_protein, methyltransferase,

DNA\_helicase, DNA\_primase.



Tree inference: IQ-TREE 2.3.6, GTR+G4 model, UFBoot2 (1000 replicates).

11 taxa: F58V + 10 representative \*Phapecoctavirus\* genomes spanning genus

diversity (MZ726793, OP172792, MN850598, JX561091, MH252123, MT944117,

OR352955, OM386656, KX664695, OR437326).



\### `/selection\\\_analysis`

PAML 4.10.7 codon-level selection analysis. One subdirectory per gene.

Each contains the control file (`.ctl`) and output (`.out`) for each model run.



\- \*\*All 11 genes\*\*: M0 (one-ratio) model only

\- \*\*depolymerase and main\_tail\_fiber\*\*: full site-model comparisons —

&#x20; M0, M1a (nearly neutral), M2a (positive selection), M7 (beta), M8 (beta+ω)



Key results:

\- Depolymerase CDS\_0012: ω = 0.020 (M0), M2a/M8 not significant → strong

&#x20; purifying selection, no evidence of positive selection

\- Main tail fiber: ω = 0.093 (M0), M2a significant → evidence of positive

&#x20; selection at receptor-binding sites



\### `/structural`

AlphaFold3 structural modeling of the CDS\_0012 C-terminal effector domain

(residues 344–1017), modeled as a C3-symmetric homotrimer (5 seeds).



Input:

\- `F58V\\\_depolymerase\\\_AF3\\\_input.fasta` — submitted sequence (674 aa × 3 chains)



Outputs (top-ranked model by ipTM):

\- `fold\\\_2026\\\_03\\\_18\\\_19\\\_32\\\_model\\\_1.cif` — predicted structure (CIF format)

\- `fold\\\_2026\\\_03\\\_18\\\_19\\\_32\\\_full\\\_data\\\_1.json` — per-residue pLDDT and PAE matrix

\- `fold\\\_2026\\\_03\\\_18\\\_19\\\_32\\\_summary\\\_confidences\\\_1.json` — ipTM, pTM, summary scores

\- `fold\\\_2026\\\_03\\\_18\\\_19\\\_32\\\_template\\\_hit\\\_1\\\_chains\\\_c\\\_a\\\_b.cif` — template used

\- `fold\\\_2026\\\_03\\\_18\\\_19\\\_32\\\_paired\\\_msa\\\_chains\\\_c\\\_a\\\_b.a3m` — paired MSA

\- `fold\\\_2026\\\_03\\\_18\\\_19\\\_32\\\_unpaired\\\_msa\\\_chains\\\_c\\\_a\\\_b.a3m` — unpaired MSA



Key quality metrics: ipTM = 0.90. Structural superposition against PDB 6E0V

(colanidase gp150, phage Phi92) in UCSF ChimeraX 1.8: RMSD = 0.905 Å.

Candidate catalytic residues: D554, Y562, W578, D656 (pLDDT > 95).

AlphaFold3 modeling was performed on the public server (alphafoldserver.com).

Top-ranked model: ipTM = 0.90. Structural superposition against PDB 6E0V

(colanidase gp150, phage Phi92) in UCSF ChimeraX 1.8: RMSD = 0.905 Å.

Four candidate catalytic residues identified: D554, Y562, W578, D656

(pLDDT > 95). AlphaFold3 model output files (.cif, confidence JSON) are 

**PENDIENTE deposited at Zenodo, DOI:** .



\---



\## Analyses not represented by output files in this repository



The following analyses were performed via web servers; results are described

in the manuscript and supplementary methods but raw output files were not

systematically archived:



\- \*\*VIRIDIC v1.1\*\* — intergenomic similarity (F58V vs. 24 \*Phapecoctavirus\*;

&#x20; F63V vs. 24 Drexlerviridae)

\- \*\*InterProScan 5.64\*\* — domain architecture of all CDS (Pfam, TIGRFAM,

&#x20; CDD, SUPERFAMILY, Gene3D, SMART, PANTHER)

\- \*\*HHpred\*\* — remote homology validation for three proteins of uncertain

&#x20; function (PDB\_mmCIF70 + Pfam-A)

\- \*\*FoldSeek\*\* — structural similarity search (PDB100 + SwissProt)

\- \*\*clinker v0.0.28\*\* — synteny analysis across 5 representative genomes



\---



\## Software versions



| Tool | Version | Reference |

|---|---|---|

| Pharokka | 1.3.2 | Bouras et al. 2023 |

| MAFFT | 7.505 | Katoh \& Standley 2013 |

| IQ-TREE | 2.3.6 | Nguyen et al. 2015 |

| UFBoot2 | — | Hoang et al. 2018 |

| PAML (codeml) | 4.10.7 | Yang 2007 |

| VIRIDIC | 1.1 | Moraru et al. 2020 |

| AlphaFold3 | server | Abramson et al. 2024 |

| UCSF ChimeraX | 1.8 | Meng et al. 2023 |

| clinker | 0.0.28 | Gilchrist \& Chooi 2021 |

| CARD RGI | — | Alcock et al. 2023 |

| VFDB | — | Liu et al. 2022 |

| MinCED | — | via Pharokka |



\---



\## Contact



Laura Ribes-Martinez — lauraribes@outlook.com

