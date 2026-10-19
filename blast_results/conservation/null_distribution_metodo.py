#!/usr/bin/env python3
"""
Genomic null distribution of protein identity for F58V CDS across Phapecoctavirus.

Recomputes the null distribution against which the conservation of the
tail-associated colanic acid depolymerase (CDS_0012) is assessed, directly from
the deposited tBLASTn output. Standard library only; no dependencies.

INPUTS
  blast_all_cds.tsv   tBLASTn of all F58V CDS against the Phapecoctavirus
                      reference genomes, outfmt 6 with the columns:
                      qseqid sseqid pident length qlen slen evalue bitscore
  F58V_pharokka.gbk   Pharokka annotation of F58V, used only to attach a
                      product name to each locus_tag.
  --refset FILE       optional; one subject accession per line. If omitted,
                      every subject present in the tBLASTn output is used.

METHOD (every choice below is a parameter, not a hidden default)
  1. Ortholog identity. For each (query CDS, subject genome) pair, the identity
     of the LONGEST HSP is taken as the identity of the ortholog. The longest
     HSP, not the one with the highest percent identity: where tBLASTn splits a
     gene into a long alignment and a short one, the short fragment often scores
     higher and is not representative. For CDS_0012 in OR067838.1, for example,
     the HSPs are 924 aa at 98.701% and 103 aa at 99.029%.
  2. Ortholog coverage. tBLASTn fragments long alignments into several HSPs, so
     coverage is the SUM of HSP lengths for that pair divided by the query
     length. Using a single HSP instead would spuriously drop long proteins
     such as CDS_0012 (1017 aa).
  3. An ortholog is called when coverage >= MIN_COV and identity >= MIN_ID.
  4. Core CDS are those with a called ortholog in EVERY genome of the
     reference set. Only core CDS enter the null distribution.
  5. The statistic per core CDS is the MEDIAN ortholog identity across the
     reference set.
  6. The percentile of a focal gene is the proportion of core CDS whose median
     identity is <= that of the focal gene, expressed as a percentage. This
     convention is stated because the alternative (strictly <) gives a
     materially different number when many CDS are at 100% identity.

OUTPUTS
  null_distribution_per_cds.tsv   one row per CDS with a called ortholog in at
                                  least one reference genome
  null_distribution_summary.txt   parameters, counts and the focal-gene and
                                  housekeeping comparisons

USAGE
  python3 null_distribution.py blast_all_cds.tsv F58V_pharokka.gbk
  python3 null_distribution.py blast_all_cds.tsv F58V_pharokka.gbk --refset refs.txt

Ribes-Martinez et al. Script deposited alongside the data it reads.
"""

import argparse
import re
import statistics as st
import sys
from collections import defaultdict

# --- Parameters -------------------------------------------------------------
MIN_COV = 0.80          # minimum summed-HSP coverage of the query
MIN_ID = 70.0           # minimum longest-HSP amino acid identity (%)
FOCAL = "CDS_0012"      # the gene whose percentile is reported
# The ten housekeeping genes curated for the comparison, by locus_tag. Listed
# explicitly rather than matched on product name, because several further
# baseplate and polymerase CDS also match by name and were not part of the set.
HOUSEKEEPING = [
    "CDS_0001",  # DNA polymerase I
    "CDS_0025",  # baseplate hub
    "CDS_0029",  # tail length tape measure protein
    "CDS_0033",  # tail sheath
    "CDS_0039",  # major head protein
    "CDS_0043",  # portal protein
    "CDS_0044",  # terminase large subunit
    "CDS_0275",  # DNA helicase
    "CDS_0302",  # endolysin
    "CDS_0311",  # DNA primase/helicase
]


def read_products(gbk_path):
    """locus_tag -> product, from a Pharokka GenBank file."""
    products, current = {}, None
    with open(gbk_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.search(r'/locus_tag="([^"]+)"', line)
            if m:
                current = m.group(1)
            m = re.search(r'/product="([^"]+)"', line)
            if m and current:
                products[current] = m.group(1)
    return products


def read_blast(tsv_path):
    """Return longest-HSP identity and summed coverage per (query, subject)."""
    best_id = defaultdict(dict)
    best_len = defaultdict(dict)
    summed = defaultdict(lambda: defaultdict(float))
    qlen = {}
    with open(tsv_path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            q, s = f[0], f[1]
            pident, length, q_len = float(f[2]), float(f[3]), float(f[4])
            summed[q][s] += length
            qlen[q] = q_len
            if length > best_len[q].get(s, -1.0):
                best_len[q][s] = length
                best_id[q][s] = pident
    return best_id, summed, qlen


def short(tag):
    """BJJLIZOP_CDS_0012 -> CDS_0012"""
    i = tag.find("CDS_")
    return tag[i:] if i >= 0 else tag


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("blast_tsv")
    ap.add_argument("pharokka_gbk")
    ap.add_argument("--refset", help="file with one subject accession per line")
    ap.add_argument("--min-cov", type=float, default=MIN_COV)
    ap.add_argument("--min-id", type=float, default=MIN_ID)
    args = ap.parse_args()

    products = read_products(args.pharokka_gbk)
    best_id, summed, qlen = read_blast(args.blast_tsv)

    all_subjects = sorted({s for d in best_id.values() for s in d})
    if args.refset:
        with open(args.refset, encoding="utf-8") as fh:
            refset = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        missing = [s for s in refset if s not in all_subjects]
        if missing:
            sys.exit("Reference genomes absent from the tBLASTn output: "
                     + ", ".join(missing))
    else:
        refset = all_subjects
    ref = set(refset)

    def is_ortholog(q, s):
        if s not in best_id[q] or best_id[q][s] < args.min_id:
            return False
        ql = qlen.get(q, 0)
        return bool(ql) and (summed[q][s] / ql) >= args.min_cov

    rows, core = [], {}
    for q in sorted(best_id):
        present = [s for s in ref if is_ortholog(q, s)]
        median = st.median([best_id[q][s] for s in present]) if present else None
        if len(present) == len(ref):
            core[q] = median
        if present:
            rows.append((q, present, median))

    if not core:
        sys.exit("No core CDS under the current thresholds.")

    def percentile(value):
        return 100.0 * sum(1 for v in core.values() if v <= value) / len(core)

    focal = next((q for q in core if short(q) == FOCAL), None)

    with open("null_distribution_per_cds.tsv", "w", encoding="utf-8") as out:
        out.write("locus_tag\tproduct\tn_orthologs\tn_reference_genomes"
                  "\tmedian_identity\tin_core\tpercentile\n")
        for q, present, median in rows:
            in_core = len(present) == len(ref)
            out.write("%s\t%s\t%d\t%d\t%.3f\t%s\t%s\n" % (
                short(q), products.get(q, "NA"), len(present), len(ref),
                median, "yes" if in_core else "no",
                "%.1f" % percentile(median) if in_core else "NA"))

    lines = []
    lines.append("Genomic null distribution of protein identity")
    lines.append("=" * 60)
    lines.append("tBLASTn input      : %s" % args.blast_tsv)
    lines.append("Annotation         : %s" % args.pharokka_gbk)
    lines.append("Reference genomes  : %d%s" % (
        len(ref), "" if args.refset else " (all subjects in the tBLASTn output)"))
    lines.append("Ortholog criteria  : summed-HSP coverage >= %.0f%%, "
                 "longest-HSP identity >= %.0f%%" % (args.min_cov * 100, args.min_id))
    lines.append("Percentile         : proportion of core CDS with median "
                 "identity <= focal")
    lines.append("")
    lines.append("CDS with >= 1 ortholog : %d" % len(rows))
    lines.append("Core CDS (all genomes) : %d" % len(core))
    lines.append("Median of distribution : %.3f%%" % st.median(list(core.values())))
    lines.append("")
    if focal is None:
        lines.append("%s is NOT in the core set under these thresholds." % FOCAL)
        tag = next((q for q in best_id if short(q) == FOCAL), None)
        if tag is not None:
            n = sum(1 for s in ref if is_ortholog(tag, s))
            lines.append("  orthologs called in %d of %d reference genomes"
                         % (n, len(ref)))
        else:
            lines.append("  %s is absent from the tBLASTn output." % FOCAL)
    else:
        lines.append("%s (%s)" % (FOCAL, products.get(focal, "NA")))
        lines.append("  median identity : %.3f%%" % core[focal])
        lines.append("  percentile      : %.1f" % percentile(core[focal]))
    lines.append("")
    lines.append("Curated housekeeping genes in the core set")
    lines.append("-" * 60)
    hk = []
    for tag in HOUSEKEEPING:
        q = next((k for k in core if short(k) == tag), None)
        if q is None:
            lines.append("  %-10s  not in core set" % tag)
        else:
            hk.append(core[q])
            lines.append("  %-10s  %-42s %7.3f%%" % (
                tag, products.get(q, "NA")[:42], core[q]))
    if hk:
        lines.append("")
        lines.append("  median, housekeeping : %.3f%%" % st.median(hk))
        if focal is not None:
            lines.append("  median, %s      : %.3f%%" % (FOCAL, core[focal]))
        lines.append("")
        lines.append("  No hypothesis test is reported. A single gene cannot be")
        lines.append("  tested against a reference group of ten: the smallest")
        lines.append("  attainable two-sided Mann-Whitney p for n1 = 1, n2 = 10")
        lines.append("  is 2/11 = 0.182, so any such p-value reports the sample")
        lines.append("  sizes rather than the data.")

    text = "\n".join(lines) + "\n"
    with open("null_distribution_summary.txt", "w", encoding="utf-8") as out:
        out.write(text)
    print(text)


if __name__ == "__main__":
    main()
