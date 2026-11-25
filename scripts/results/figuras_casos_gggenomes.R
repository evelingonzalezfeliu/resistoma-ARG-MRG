#!/usr/bin/env Rscript

invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
.libPaths(c(
  "/Users/efeliu/Documents/Mlatorre/analisis/.Rlib",
  .libPaths()
))

suppressPackageStartupMessages({
  library(gggenomes)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
analysis_dir <- file.path(base_dir, "analisis_arg_mrg_20261022")
figure_dir <- file.path(analysis_dir, "figuras")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, comment.char = "") {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = comment.char,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

attr_value <- function(x, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]*)")
  hit <- regexec(pattern, x, perl = TRUE)
  value <- regmatches(x, hit)
  vapply(value, function(z) {
    if (length(z) >= 2) URLdecode(z[[2]]) else NA_character_
  }, character(1))
}

read_gff_region <- function(sample, contig, region_start, region_end) {
  path <- file.path(base_dir, "gff", paste0(sample, ".gff"))
  gff <- read.delim(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  names(gff) <- c(
    "seq_id", "source", "type", "start", "end", "score",
    "strand", "phase", "attributes"
  )
  gff %>%
    filter(
      seq_id == contig,
      type == "CDS",
      end >= region_start,
      start <= region_end
    ) %>%
    mutate(
      feat_id = attr_value(attributes, "ID"),
      gene = attr_value(attributes, "gene"),
      product = attr_value(attributes, "product")
    )
}

format_bp <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
}

base_theme <- theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold", size = 16, hjust = 0.5, colour = "#24323B"
    ),
    plot.subtitle = element_text(
      size = 11, hjust = 0.5, colour = "#4C5961", margin = margin(b = 8)
    ),
    plot.caption = element_text(size = 8.5, colour = "#59666E"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E1E6E9", linewidth = 0.4),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(18, 24, 12, 24)
  )

save_plot <- function(plot, stem) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot,
    width = 14,
    height = 6.6,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot,
    width = 14,
    height = 6.6,
    units = "in",
    dpi = 300,
    device = ragg::agg_png,
    background = "white"
  )
}

make_c04 <- function() {
  sample <- "B28_C04"
  contig <- "B28_C04|contig_1"
  region_start <- 370000L
  region_end <- 440000L
  offset <- region_start - 1L
  seq_id_plot <- paste0(sample, "_mobiloma")

  genes <- read_gff_region(sample, contig, region_start, region_end) %>%
    mutate(
      categoria = case_when(
        gene %in% c("tbtM", "tbtB", "tbtA") ~ "ARG",
        gene == "chrA" ~ "MRG",
        grepl(
          "transposase|integrase|recombinase|TniQ",
          product,
          ignore.case = TRUE
        ) ~ "Movilización",
        TRUE ~ "Otros genes"
      ),
      etiqueta = case_when(
        gene %in% c("tbtM", "tbtB", "tbtA", "chrA") ~ gene,
        TRUE ~ ""
      ),
      file_id = sample,
      seq_id = seq_id_plot,
      type = "CDS",
      start = pmax(start, region_start) - offset,
      end = pmin(end, region_end) - offset
    ) %>%
    filter(categoria != "Otros genes") %>%
    select(file_id, seq_id, start, end, strand, type, categoria, etiqueta)

  mef_path <- file.path(
    base_dir, "resultados_mobileelementfinder", sample, paste0(sample, ".csv")
  )
  mef <- read.csv(
    mef_path,
    comment.char = "#",
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    mutate(contig_clean = sub(" .*", "", contig)) %>%
    filter(
      contig_clean == .env$contig,
      end >= region_start,
      start <= region_end,
      name %in% c("Tn5563", "cn_14636_ISPsme1", "cn_31491_ISPsme1")
    ) %>%
    mutate(
      file_id = sample,
      seq_id = seq_id_plot,
      start = pmax(start, region_start) - offset,
      end = pmin(end, region_end) - offset,
      strand = NA_character_,
      feature_type = ifelse(
        type == "unit transposon", "Transposón unitario", "Transposón compuesto"
      ),
      label = case_when(
        name == "cn_14636_ISPsme1" ~ "Tn compuesto 1 (ISPsme1)",
        name == "cn_31491_ISPsme1" ~ "Tn compuesto 2 (ISPsme1)",
        TRUE ~ name
      )
    ) %>%
    select(file_id, seq_id, start, end, strand, feature_type, label)

  seqs <- tibble(
    file_id = sample,
    seq_id = seq_id_plot,
    seq_desc = " ",
    length = region_end - region_start + 1L
  )

  chr_a <- genes %>% filter(etiqueta == "chrA")
  tbt_m <- genes %>% filter(etiqueta == "tbtM")
  d1 <- chr_a$end[[1]]
  d2 <- tbt_m$start[[1]]
  dmid <- (d1 + d2) / 2

  gggenomes(genes = genes, seqs = seqs) +
    geom_seq(colour = "#53636D", linewidth = 1.1) +
    geom_gene(
      aes(fill = categoria),
      position = "strand",
      colour = "#43515A",
      size = 5.6
    ) +
    geom_gene_tag(
      aes(label = etiqueta),
      position = "strand",
      angle = 25,
      nudge_y = 0.10,
      size = 4,
      fontface = "bold",
      check_overlap = FALSE
    ) +
    annotate(
      "segment",
      x = mef$start,
      xend = mef$end,
      y = c(0.82, 0.68, 0.54),
      yend = c(0.82, 0.68, 0.54),
      linewidth = 5,
      colour = c("#008837", "#7B3294", "#C51B7D"),
      lineend = "butt"
    ) +
    annotate(
      "text",
      x = (mef$start + mef$end) / 2,
      y = c(0.70, 0.56, 0.42),
      label = mef$label,
      size = 3.1,
      fontface = "bold",
      colour = c("#006D2C", "#6A1B75", "#A50F62")
    ) +
    annotate(
      "segment",
      x = d1, xend = d2, y = 1.48, yend = 1.48,
      linewidth = 0.55, colour = "#343E45"
    ) +
    annotate(
      "segment",
      x = c(d1, d2), xend = c(d1, d2),
      y = 1.44, yend = 1.52,
      linewidth = 0.55, colour = "#343E45"
    ) +
    annotate(
      "text",
      x = dmid, y = 1.58,
      label = "11.206 pb",
      size = 3.8, fontface = "bold", colour = "#343E45"
    ) +
    scale_fill_manual(
      values = c(
        "ARG" = "#B2182B",
        "MRG" = "#2166AC",
        "Movilización" = "#E69F00"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      labels = function(x) label_number(
        big.mark = ".", decimal.mark = ",", accuracy = 1
      )(x + offset),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(ylim = c(0.38, 1.82), clip = "off") +
    labs(
      title = "B28_C04: ARG-MRG en una región rica en elementos móviles",
      subtitle = paste0(
        "chrA y el operón tbtM-tbtB-tbtA; ",
        "dos transposones compuestos superpuestos"
      ),
      x = NULL,
      y = NULL,
      caption = NULL
    ) +
    base_theme
}

make_d07 <- function() {
  sample <- "B28_D07"
  contig <- "B28_D07|contig_1"
  region_start <- 3790000L
  region_end <- 3867000L
  offset <- region_start - 1L
  seq_id_plot <- paste0(sample, "_integron")
  integron_start <- 3854473L
  integron_end <- 3860361L
  arr_end_abs <- 3800100L
  compressed_gap <- 12000L

  compress_x <- function(x) {
    ifelse(
      x <= arr_end_abs,
      x - region_start + 1L,
      ifelse(
        x >= integron_start,
        arr_end_abs - region_start + 1L + compressed_gap +
          (x - integron_start),
        arr_end_abs - region_start + 1L +
          (x - arr_end_abs) / (integron_start - arr_end_abs) * compressed_gap
      )
    )
  }

  genes <- read_gff_region(sample, contig, region_start, region_end) %>%
    mutate(
      categoria = case_when(
        gene == "arr" | feat_id == "B28_D07|IFLEDJ_17655" ~ "ARG",
        grepl("Integron integrase", product, ignore.case = TRUE) ~ "Integrasa intI",
        grepl("transposase", product, ignore.case = TRUE) ~ "Transposasa",
        start >= integron_start & end <= integron_end ~ "Cassette del integrón",
        TRUE ~ "Otros genes"
      ),
      etiqueta = case_when(
        feat_id == "B28_D07|IFLEDJ_17655" ~ "arr",
        grepl("Integron integrase", product, ignore.case = TRUE) ~ "intI",
        grepl("IS3 family transposase", product, ignore.case = TRUE) ~ "IS3",
        TRUE ~ ""
      ),
      file_id = sample,
      seq_id = seq_id_plot,
      type = "CDS",
      start = compress_x(pmax(start, region_start)),
      end = compress_x(pmin(end, region_end))
    ) %>%
    filter(categoria != "Otros genes") %>%
    select(file_id, seq_id, start, end, strand, type, categoria, etiqueta)

  integron_path <- file.path(
    base_dir,
    "resultados_integronfinder",
    sample,
    paste0("Results_Integron_Finder_", sample),
    paste0(sample, ".integrons")
  )
  int_rows <- read_tsv(integron_path, comment.char = "#") %>%
    filter(type == "complete")
  attc <- int_rows %>%
    filter(type_elt == "attC") %>%
    transmute(
      file_id = sample,
      seq_id = seq_id_plot,
      start = compress_x(pos_beg),
      end = compress_x(pos_end),
      strand = ifelse(strand == "1", "+", "-"),
      feature_type = "Sitio attC",
      label = ""
    )
  integron_feature <- tibble(
    file_id = sample,
    seq_id = seq_id_plot,
    start = compress_x(integron_start),
    end = compress_x(integron_end),
    strand = NA_character_,
    feature_type = "Integrón completo",
    label = "Integrón completo (7 sitios attC)"
  )
  seqs <- tibble(
    file_id = sample,
    seq_id = seq_id_plot,
    seq_desc = " ",
    length = compress_x(region_end)
  )

  arr_gene <- genes %>% filter(etiqueta == "arr")
  dist_start <- arr_gene$end[[1]]
  dist_end <- compress_x(integron_start)
  dist_mid <- (dist_start + dist_end) / 2

  gggenomes(genes = genes, seqs = seqs) +
    geom_seq(colour = "#53636D", linewidth = 1.1) +
    geom_gene(
      aes(fill = categoria),
      position = "strand",
      colour = "#43515A",
      size = 5.6
    ) +
    geom_gene_tag(
      aes(label = etiqueta),
      position = "strand",
      angle = 25,
      nudge_y = 0.10,
      size = 4,
      fontface = "bold",
      check_overlap = FALSE
    ) +
    annotate(
      "segment",
      x = integron_feature$start,
      xend = integron_feature$end,
      y = 0.72,
      yend = 0.72,
      linewidth = 7,
      colour = "#1B7837",
      lineend = "butt"
    ) +
    annotate(
      "segment",
      x = (attc$start + attc$end) / 2,
      xend = (attc$start + attc$end) / 2,
      y = 0.64,
      yend = 0.80,
      linewidth = 1,
      colour = "#2166AC"
    ) +
    annotate(
      "text",
      x = (integron_feature$start + integron_feature$end) / 2,
      y = 0.54,
      label = "Integrón completo (7 sitios attC; marcas azules)",
      size = 3.2,
      fontface = "bold",
      colour = "#146B2D"
    ) +
    annotate(
      "segment",
      x = dist_start, xend = dist_end,
      y = 1.48, yend = 1.48,
      linewidth = 0.55, colour = "#343E45"
    ) +
    annotate(
      "segment",
      x = c(dist_start, dist_end), xend = c(dist_start, dist_end),
      y = 1.44, yend = 1.52,
      linewidth = 0.55, colour = "#343E45"
    ) +
    annotate(
      "text",
      x = dist_mid, y = 1.58,
      label = "54.373 pb",
      size = 3.8, fontface = "bold", colour = "#343E45"
    ) +
    annotate(
      "label",
      x = dist_mid, y = 1.02,
      label = "//",
      size = 6,
      fontface = "bold",
      colour = "#53636D",
      fill = "white",
      linewidth = 0
    ) +
    scale_fill_manual(
      values = c(
        "ARG" = "#B2182B",
        "Integrasa intI" = "#542788",
        "Cassette del integrón" = "#C2A5CF",
        "Transposasa" = "#E69F00"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = NULL,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(ylim = c(0.40, 1.82), clip = "off") +
    labs(
      title = "B28_D07: ARG arr próximo a un integrón completo",
      subtitle = paste0(
        "Integrón de 5,9 kb con siete sitios attC, transposasas e integrasa intI ",
        "(intervalo comprimido)"
      ),
      x = NULL,
      y = NULL,
      caption = NULL
    ) +
    base_theme
}

c04_plot <- make_c04()
d07_plot <- make_d07()

save_plot(c04_plot, "figura_C_B28_C04_ARG_MRG_mobiloma")
save_plot(d07_plot, "figura_D_B28_D07_ARG_integron")

cat("Figuras creadas en:", figure_dir, "\n")
