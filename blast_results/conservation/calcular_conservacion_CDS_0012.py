#!/usr/bin/env python3
"""
Intra-genus conservation of F58V CDS_0012 (colanic acid depolymerase),
recomputed from the deposited tBLASTn output. Standard library only.

Writes the two tables behind the intra-genus block of Supplementary Table 6
and shows why the second tBLASTn output in this folder is not used for it.

Inputs, paths relative to this script
  ../blast_all_cds.tsv
      tBLASTn of the 313 F58V CDS against the 46 Phapecoctavirus reference
      genomes, outfmt 6: qseqid sseqid pident length qlen slen evalue bitscore.
  tblastn_CDS_0012_65subjects.tsv
      tBLASTn of CDS_0012 alone against a different, broader database,
      outfmt 6 with 14 columns: qseqid sseqid pident length mismatch gapopen
      qstart qend sstart send evalue bitscore qlen slen.
      Not used for the table; processed only to document why.

Method, the same as in null_distribution_metodo.py so that both sets of
figures are comparable
  identity  longest HSP. tBLASTn splits long genes into several HSPs, and the
            short fragment often scores higher without being representative.
  coverage  summed HSP length / query length. It exceeds 100% when HSPs
            overlap and is reported uncapped.
  E-value   best HSP.

Outputs
  conservacion_CDS_0012_46genomas.tsv
  duplicados_en_tblastn_65subjects.tsv
"""
import collections
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
QUERY = "CDS_0012"
QLEN_ESPERADA = 1017


def leer(nombre, col_subject, col_pident, col_alen, col_evalue, col_slen, col_qlen=None):
    ruta = os.path.join(AQUI, nombre)
    if not os.path.exists(ruta):
        sys.exit("No encuentro %s en %s" % (nombre, AQUI))
    hsps = collections.defaultdict(list)
    slen = {}
    qlen = set()
    for linea in open(ruta):
        linea = linea.rstrip("\n")
        if not linea or linea.startswith("#"):
            continue
        c = linea.split("\t")
        if QUERY not in c[0]:
            continue
        s = c[col_subject]
        hsps[s].append((int(c[col_alen]), float(c[col_pident]), c[col_evalue]))
        slen[s] = int(c[col_slen])
        if col_qlen is not None:
            qlen.add(int(c[col_qlen]))
    return hsps, slen, qlen


def resumir(hsps, qlen_consulta):
    filas = []
    for s, h in hsps.items():
        h.sort(key=lambda x: -x[0])
        ident_mas_largo = h[0][1]
        cobertura = sum(x[0] for x in h) / qlen_consulta * 100.0
        mejor_e = min(h, key=lambda x: float(x[2]))[2]
        filas.append((s, ident_mas_largo, cobertura, len(h), mejor_e))
    filas.sort(key=lambda r: (-r[1], -r[2]))
    return filas


def main():
    hsps46, slen46, qlens = leer(
        os.path.join("..", "blast_all_cds.tsv"),
        col_subject=1, col_pident=2, col_alen=3, col_evalue=6, col_slen=5, col_qlen=4)
    if not hsps46:
        sys.exit("No hay filas de %s en el fichero de 46 genomas" % QUERY)
    qlen = qlens.pop() if len(qlens) == 1 else QLEN_ESPERADA
    if qlen != QLEN_ESPERADA:
        print("AVISO: la longitud de consulta en el fichero es %d, no %d" % (qlen, QLEN_ESPERADA))

    filas = resumir(hsps46, qlen)
    salida = os.path.join(AQUI, "conservacion_CDS_0012_46genomas.tsv")
    with open(salida, "w") as f:
        f.write("# Conservacion de F58V CDS_0012 (%d aa) en los 46 genomas de referencia\n" % qlen)
        f.write("# identidad: del HSP mas largo. cobertura: suma de HSP / longitud de consulta.\n")
        f.write("subject\tgenome_length_bp\tidentity_pct\tcoverage_pct\tn_HSP\tbest_evalue\n")
        for s, ident, cob, n, ev in filas:
            f.write("%s\t%d\t%.3f\t%.1f\t%d\t%s\n" % (s, slen46[s], ident, cob, n, ev))

    idents = [r[1] for r in filas]
    cobs = [r[2] for r in filas]
    orden = sorted(idents)
    mediana = orden[len(orden) // 2] if len(orden) % 2 else (orden[len(orden) // 2 - 1] + orden[len(orden) // 2]) / 2
    print("Genomas con ortologo        : %d" % len(filas))
    print("Identidad minima y maxima   : %.3f  %.3f" % (min(idents), max(idents)))
    print("Identidad mediana           : %.3f" % mediana)
    print("Cobertura >= 99.5%%          : %d de %d" % (sum(1 for c in cobs if c >= 99.5), len(cobs)))
    print("Escrito                     : %s" % os.path.basename(salida))

    # 65-subject file: an extra subject with the same genome length as one of
    # the 46 is counted as a RefSeq duplicate of it.
    hsps65, slen65, _ = leer(
        "tblastn_CDS_0012_65subjects.tsv",
        col_subject=1, col_pident=2, col_alen=3, col_evalue=10, col_slen=13)
    if not hsps65:
        print("\n(No encuentro el fichero de 65 subjects; omito el bloque de duplicados.)")
        return
    extra = sorted(set(hsps65) - set(hsps46))
    por_longitud = collections.defaultdict(list)
    for s, L in slen46.items():
        por_longitud[L].append(s)
    salida2 = os.path.join(AQUI, "duplicados_en_tblastn_65subjects.tsv")
    n_dup = 0
    with open(salida2, "w") as f:
        f.write("# Subjects presentes en tblastn_CDS_0012_65subjects.tsv y ausentes del conjunto de 46.\n")
        f.write("# Se empareja por longitud de genoma: una coincidencia exacta indica que el subject\n")
        f.write("# es la copia RefSeq de un genoma ya presente bajo su accesion GenBank.\n")
        f.write("subject_extra\tgenome_length_bp\tduplicado_de\tveredicto\n")
        for s in extra:
            L = slen65[s]
            iguales = [x for x in por_longitud.get(L, []) if x != s]
            if iguales:
                n_dup += 1
                f.write("%s\t%d\t%s\tduplicado\n" % (s, L, ";".join(sorted(iguales))))
            else:
                f.write("%s\t%d\t-\tgenoma nuevo\n" % (s, L))
    print("\nSubjects del fichero de 65    : %d" % len(hsps65))
    print("Ausentes del conjunto de 46   : %d" % len(extra))
    print("  de ellos, duplicados RefSeq : %d" % n_dup)
    print("  genomas realmente nuevos    : %d" % (len(extra) - n_dup))
    print("Escrito                       : %s" % os.path.basename(salida2))


if __name__ == "__main__":
    main()
