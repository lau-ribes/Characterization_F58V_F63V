# Figure 2: domain architecture of CDS_0012, CDS_0019 and CDS_0015 (F58V) and
# CDS_0056 (F63V), with every feature drawn at its reported coordinates.
# HHpred structural homology goes in a dashed lower lane because it overlaps
# the domains; that lane is drawn for CDS_0056 and CDS_0015 only.
#
# Input: dominios.tsv, compiled from
#   InterProScan v5.77  ../annotation/interproScan/F58V_interproscan.tsv
#                       ../annotation/interproScan/F63V_interproscan.tsv
#   HHpred PDB_mmCIF70  ../annotation/HHpred/HHpred_F58V_CDS0015_tailfiber.hhr
#                       ../annotation/HHpred/HHpred_F63V_CDS0056_tailfiber.hhr
# CDS_0015 has no InterProScan hit.
#
# Usage, from this folder: Rscript figura2.R

suppressPackageStartupMessages(library(ggplot2))

TIPO <- "Liberation Sans"   # metric-compatible with Arial

# FALSE drops CDS_0019 and writes the files without the _con_CDS_0019 suffix.
INCLUIR_0019 <- TRUE

d <- read.delim("dominios.tsv", stringsAsFactors = FALSE)
if (!INCLUIR_0019) d <- subset(d, proteina != "CDS_0019")

prot <- data.frame(
  proteina = c("CDS_0012", "CDS_0019", "CDS_0056", "CDS_0015"),
  fago     = c("F58V", "F58V", "F63V", "F58V"),
  longitud = c(1017, 921, 844, 353),
  y        = c(0, -1.30, -2.55, -4.20),
  stringsAsFactors = FALSE)
if (!INCLUIR_0019) {
  prot <- subset(prot, proteina != "CDS_0019")
  prot$y <- c(0, -1.00, -2.65)
}
prot$etiqueta <- sprintf("%s\n%s (%s aa)", prot$proteina, prot$fago,
                         formatC(prot$longitud, big.mark = ","))

relleno <- c(fibra = "#B5D4F4", tsp1 = "#9FE1CB", pectina = "#C0DD97",
             wcam  = "#FAC775", endocat = "#F0B3A0", endocat2 = "#F0B3A0",
             endoflank = "#F8E0D7", endoflank2 = "#F8E0D7",
             hh5yvq = "#CECBF6", hh9cuy = "#E4E2FA", hh7vyv = "#D7EAB8")
borde   <- c(fibra = "#378ADD", tsp1 = "#1D9E75", pectina = "#639922",
             wcam  = "#BA7517", endocat = "#C2674C", endocat2 = "#C2674C",
             endoflank = "#C2674C", endoflank2 = "#C2674C",
             hh5yvq = "#534AB7", hh9cuy = "#8E88D6", hh7vyv = "#639922")
nom <- c(fibra   = "Tail fibre N-terminal (PF26208)",
         tsp1    = "TSP1/Gp66 (PF18668)",
         pectina = "Pectin lyase-like fold (SSF51126)",
         wcam    = "WcaM colanic acid (PF27800)",
         endocat   = "Endosialidase catalytic (PF12217, PF12219)",
         endoflank = "Endosialidase flanking (PF12218, PF13884)",
         hh7vyv  = "HHpred: depolymerase 7VYV (chain B)",
         hh9cuy  = "HHpred: short tail fibre 9CUY",
         hh5yvq  = "HHpred: tail fibre protein S 5YVQ (chain A)")

# Heights in y-axis data units.
ALTO   <- c(dominio = 0.70, plegamiento = 0.70, estructura = 0.50)
ESQ    <- 0.105            # half-height of the grey backbone
DESPL  <- -0.62            # HHpred lane offset below the row centre

d <- merge(d, prot[, c("proteina", "y")], by = "proteina")
d$y <- d$y + ifelse(d$pista == "estructura", DESPL, 0)
d$h <- ALTO[d$pista]
d$clave[d$clave == "endocat2"] <- "endocat"
d$clave[d$clave == "endoflank2"] <- "endoflank"
d$clave  <- factor(d$clave, levels = names(nom))
d$rotulo <- sprintf("%s-%s", d$inicio, d$fin)
d$ancho  <- d$fin - d$inicio

# WcaM lies inside the pectin lyase-like fold of CDS_0012, so the fold label
# moves to the free stretch.
d$x_rot <- (d$inicio + d$fin)/2
d$x_rot[d$proteina == "CDS_0012" & d$clave == "pectina"] <- 445
dentro <- subset(d, ancho >= 105)

# Single-line callouts, so they do not take height from the bars.
y19  <- if (INCLUIR_0019) prot$y[prot$proteina == "CDS_0019"] else NA
yF63 <- prot$y[prot$proteina == "CDS_0056"]
yF15 <- prot$y[prot$proteina == "CDS_0015"]
llam <- data.frame(
  clave  = c("tsp1", "wcam", "endoflank", "endoflank", "hh7vyv", "hh9cuy", "hh5yvq"),
  x      = c(331, 582, 109, 790, 243, 35, 315),
  texto  = c("TSP1/Gp66  302-361", "WcaM (PF27800)  502-662",
             "PF12218  76-142", "PF13884  760-819",
             "7VYV (K2-type)", "9CUY  1-70", "5YVQ chain A  278-353"),
  arriba = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  yb     = c(0.35, 0.35, y19 - 0.35, y19 - 0.35,
             yF63 + DESPL - 0.25, yF15 + DESPL - 0.25, yF15 + DESPL - 0.25),
  yt     = c(0.62, 0.62, y19 - 0.52, y19 - 0.52,
             yF63 + DESPL - 0.40, yF15 + DESPL - 0.40, yF15 + DESPL - 0.40),
  stringsAsFactors = FALSE)
if (!INCLUIR_0019) llam <- subset(llam, clave != "endoflank")
llam$clave <- factor(llam$clave, levels = names(nom))

notas <- data.frame(
  x = c(700, 372), y = c(yF63, yF15),
  texto = c("WcaM absent", "no InterProScan assignment; 71-277 unassigned"),
  stringsAsFactors = FALSE)

p <- ggplot() +
  geom_rect(data = prot,
            aes(xmin = 1, xmax = longitud, ymin = y - ESQ, ymax = y + ESQ),
            fill = "#E8E8E8", colour = "#C8C8C8", linewidth = 0.3) +
  geom_rect(data = subset(d, pista != "estructura"),
            aes(xmin = inicio, xmax = fin, ymin = y - h/2, ymax = y + h/2,
                fill = clave, colour = clave), linewidth = 0.7) +
  geom_rect(data = subset(d, pista == "estructura"),
            aes(xmin = inicio, xmax = fin, ymin = y - h/2, ymax = y + h/2,
                fill = clave, colour = clave),
            linewidth = 0.6, linetype = "22") +
  geom_text(data = dentro, aes(x = x_rot, y = y, label = rotulo),
            size = 3.3, colour = "black", family = TIPO) +
  geom_segment(data = llam, aes(x = x, xend = x, y = yb, yend = yt,
                                colour = clave),
               linewidth = 0.35, linetype = "22", show.legend = FALSE) +
  geom_label(data = llam, aes(x = x, y = yt, label = texto, colour = clave),
             vjust = ifelse(llam$arriba, 0, 1), size = 3.1, label.size = 0.35,
             label.r = unit(0.12, "lines"), label.padding = unit(0.25, "lines"),
             fill = "white", show.legend = FALSE, family = TIPO) +
  geom_text(data = notas, aes(x = x, y = y, label = texto),
            size = 3.0, colour = "grey45", fontface = "italic", hjust = 0,
            family = TIPO) +
  scale_fill_manual(values = relleno, labels = nom, breaks = names(nom),
                    name = NULL, drop = FALSE) +
  scale_colour_manual(values = borde, labels = nom, breaks = names(nom),
                      name = NULL, drop = FALSE) +
  scale_y_continuous(breaks = prot$y, labels = prot$etiqueta,
                     expand = expansion(add = c(0, 0))) +
  scale_x_continuous(name = "Amino acid position", limits = c(1, 1035),
                     breaks = seq(0, 1000, 200),
                     expand = expansion(mult = c(0.015, 0.015))) +
  # Explicit limits: automatic expansion clips the top and bottom callouts.
  coord_cartesian(ylim = c(min(prot$y) + DESPL - 0.75, 1.00), clip = "off") +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE,
                             override.aes = list(linetype = "solid")),
         colour = guide_legend(ncol = 2, byrow = TRUE)) +
  theme_classic(base_size = 12, base_family = TIPO) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(face = "bold", size = 11, colour = "black",
                                    lineheight = 1.15, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        axis.line.x  = element_line(colour = "black", linewidth = 0.6),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.6),
        axis.text.x  = element_text(size = 10.5, colour = "black"),
        axis.title.x = element_text(size = 11.5, colour = "black",
                                    margin = margin(t = 6)),
        legend.position = "bottom",
        legend.key.size = unit(0.48, "cm"),
        legend.text = element_text(size = 9.5, colour = "black"),
        legend.margin = margin(t = 4),
        legend.spacing.x = unit(0.3, "cm"),
        plot.margin = margin(16, 14, 4, 4))

SUF <- if (INCLUIR_0019) "_con_CDS_0019" else ""
ALT <- if (INCLUIR_0019) 5.5 else 4.7
ggsave(paste0("Figure_2_Domain_Architecture", SUF, ".png"), p,
       width = 7.6, height = ALT, dpi = 600, bg = "white")
ggsave(paste0("Figure_2_Domain_Architecture", SUF, ".tiff"), p,
       width = 7.6, height = ALT, dpi = 600, bg = "white", compression = "lzw")
cat("hecho:", R.version.string, "/ ggplot2",
    as.character(packageVersion("ggplot2")), "\n")
