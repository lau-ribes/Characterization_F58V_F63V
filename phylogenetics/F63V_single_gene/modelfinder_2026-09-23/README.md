# F63V single-gene trees: model-selection reanalysis

Rerun with IQ-TREE 2.3.6 on the TerL and MCP protein alignments of the parent
directory, `TerL_F63V.aln` and `MCP_F63V.aln` (named `terL.aln` and `mcp.aln`
in the commands recorded in the logs). Same seven taxa as the parent directory.

## `*_rep.*`: reproduction of the fixed-model trees

    iqtree2 -s <gene>.aln -m LG+F+G4 -B 1000 -alrt 1000

The log-likelihoods match those of the trees in the parent directory to the
fourth decimal (TerL -2545.4229, MCP -1650.3290), and the topology is the same.

## `*_mfp.*`: model selection with ModelFinder

    iqtree2 -s <gene>.aln -m MFP -B 1000 -alrt 1000

BIC selects LG+R2 for TerL and Q.pfam+I for MCP. The topology is the same as
under LG+F+G4 in both genes; only the support values change (SH-aLRT / UFBoot):

| Node | LG+F+G4 | ModelFinder |
|---|---|---|
| TerL: F63V, ES10, AN_ECEAS and IMM-001 | 99.9 / 100 | 100 / 100 |
| TerL: F63V with IMM-001 | 65.2 / 71 | 81.6 / 79 |
| MCP: F63V with DTL and IMM-001 | 92.5 / 85 | 92.7 / 87 |

Each `.log` holds the full command line, and each `.iqtree` the version, input,
model ranking and support values.
