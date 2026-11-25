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

base_dir <- "/Users/efeliu/Documents/Mlatorre/analisis"
analysis_dir <- file.path(base_dir, "analisis_arg_mrg_20261022")
pairs_path <- file.path(
  analysis_dir,
  "colocalizacion/pares_arg_mrg_mismo_contig_alta_confianza.tsv"
)
mob_root <- file.path(base_dir, "resultados_mobsuite")
is_root <- file.path(base_dir, "resultados_isescan")
output_dir <- file.path(
  analysis_dir,
  "colocalizacion/plots/atlas_ARG_MRG_gggenomes"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

interval_distance <- function(a_start, a_end, b_start, b_end) {
  ifelse(
    b_end < a_start,
    a_start - b_end - 1L,
    ifelse(b_start > a_end, b_start - a_end - 1L, 0L)
  )
}

format_bp <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
}

pairs <- read_tsv(pairs_path) %>%
  arrange(sample_id, distancia_bordes_bp, arg_symbol, mrg_feature_label) %>%
  group_by(sample_id) %>%
  slice(1L) %>%
  ungroup()

make_plot <- function(pair) {
  sample <- pair$sample_id[[1]]
  target_contig <- pair$arg_contig[[1]]
  arg_start_abs <- pair$arg_start[[1]]
  arg_end_abs <- pair$arg_stop[[1]]
  mrg_start_abs <- pair$mrg_start[[1]]
  mrg_end_abs <- pair$mrg_stop[[1]]
  pair_start <- min(arg_start_abs, mrg_start_abs)
  pair_end <- max(arg_end_abs, mrg_end_abs)
  padding <- 50000L
  region_start <- max(1L, pair_start - padding)
  region_end <- pair_end + padding
  coordinate_offset <- region_start - 1L
  seq_id_plot <- paste0(sample, "_region")

  mob_path <- file.path(mob_root, sample, "contig_report.txt")
  mob <- read_tsv(mob_path) %>%
    mutate(contig_clean = sub(" .*", "", contig_id)) %>%
    filter(contig_clean == target_contig)
  if (nrow(mob) != 1L) {
    stop("No se pudo resolver MOB-suite para ", sample, " / ", target_contig)
  }
  molecule <- mob$molecule_type[[1]]
  molecule_es <- recode(
    molecule,
    chromosome = "cromosoma",
    plasmid = "plasmido",
    .default = molecule
  )
  plasmid_mobility <- ""
  if (molecule == "plasmid") {
    typer_path <- file.path(mob_root, sample, "mobtyper_results.txt")
    if (file.exists(typer_path)) {
      typer <- read_tsv(typer_path) %>%
        filter(primary_cluster_id == mob$primary_cluster_id[[1]])
      if (nrow(typer) >= 1L) {
        plasmid_mobility <- typer$predicted_mobility[[1]]
      }
    }
  }

  genes <- tibble(
    file_id = sample,
    seq_id = seq_id_plot,
    start = c(arg_start_abs, mrg_start_abs) - coordinate_offset,
    end = c(arg_end_abs, mrg_end_abs) - coordinate_offset,
    strand = c(pair$arg_strand[[1]], pair$mrg_strand[[1]]),
    type = "CDS",
    categoria = c("ARG", "MRG"),
    etiqueta = c(pair$arg_symbol[[1]], pair$mrg_feature_label[[1]])
  )

  is_path <- file.path(is_root, sample, "fna", paste0(sample, ".fna.tsv"))
  is_features <- tibble()
  if (file.exists(is_path)) {
    is_features <- read_tsv(is_path) %>%
      filter(
        seqID == target_contig,
        isEnd >= region_start,
        isBegin <= region_end
      ) %>%
      mutate(
        distance_arg = interval_distance(
          arg_start_abs, arg_end_abs, isBegin, isEnd
        ),
        distance_mrg = interval_distance(
          mrg_start_abs, mrg_end_abs, isBegin, isEnd
        ),
        distance_target = pmin(distance_arg, distance_mrg)
      ) %>%
      transmute(
        file_id = sample,
        seq_id = seq_id_plot,
        start = pmax(isBegin, region_start) - coordinate_offset,
        end = pmin(isEnd, region_end) - coordinate_offset,
        strand = as.character(
          ifelse(strand %in% c("+", "-"), strand, NA_character_)
        ),
        feature_type = "IS (ISEScan)",
        label = as.character(
          ifelse(
            distance_target <= 10000,
            paste0(family, ifelse(type == "p", " (parcial)", "")),
            ""
          )
        ),
        distance_target
      )
  }

  mge_path <- file.path(mob_root, sample, "mge.report.txt")
  mge_features <- tibble()
  if (file.exists(mge_path) && file.info(mge_path)$size > 0L) {
    mge_features <- read_tsv(mge_path)
    if (nrow(mge_features) > 0L) {
      mge_features <- mge_features %>%
        mutate(contig_clean = sub(" .*", "", contig_id)) %>%
        filter(
          contig_clean == target_contig,
          contig_end >= region_start,
          contig_start <= region_end
        ) %>%
        transmute(
          file_id = sample,
          seq_id = seq_id_plot,
          start = pmax(contig_start, region_start) - coordinate_offset,
          end = pmin(contig_end, region_end) - coordinate_offset,
          strand = as.character(
            ifelse(sstrand %in% c("+", "-"), sstrand, NA_character_)
          ),
          feature_type = "Hit MGE (MOB-suite)",
          label = as.character(paste0(mge_type, " / ", mge_subtype)),
          distance_target = pmin(
            interval_distance(
              arg_start_abs, arg_end_abs, contig_start, contig_end
            ),
            interval_distance(
              mrg_start_abs, mrg_end_abs, contig_start, contig_end
            )
          )
        )
    }
  }

  features <- bind_rows(is_features, mge_features)
  if (nrow(features) == 0L) {
    features <- tibble(
      file_id = character(), seq_id = character(), start = integer(),
      end = integer(), strand = character(), feature_type = character(),
      label = character(), distance_target = integer()
    )
  }

  seqs <- tibble(
    file_id = sample,
    seq_id = seq_id_plot,
    seq_desc = paste(sample, target_contig, molecule_es, sep = " | "),
    length = region_end - region_start + 1L
  )

  arg_gene <- genes %>% filter(categoria == "ARG")
  mrg_gene <- genes %>% filter(categoria == "MRG")
  distance_x1 <- min(arg_gene$end, mrg_gene$end)
  distance_x2 <- max(arg_gene$start, mrg_gene$start)
  if (arg_gene$start < mrg_gene$start) {
    distance_x1 <- arg_gene$end
    distance_x2 <- mrg_gene$start
  } else {
    distance_x1 <- mrg_gene$end
    distance_x2 <- arg_gene$start
  }
  distance_mid <- (distance_x1 + distance_x2) / 2

  context_label <- paste0("MOB-suite: ", molecule_es)
  if (nzchar(plasmid_mobility)) {
    context_label <- paste0(context_label, " (", plasmid_mobility, ")")
  }
  evidence_label <- paste0(
    "ISEScan: ", nrow(is_features), " IS en la region | ",
    "MOB-suite: ", nrow(mge_features), " hit(s) MGE"
  )

  p <- gggenomes(genes = genes, seqs = seqs, feats = features) +
    geom_seq(colour = "#53636D", linewidth = 1.1) +
    geom_gene(
      aes(fill = categoria),
      position = "strand",
      colour = "#43515A",
      size = 5.8
    ) +
    geom_gene_tag(
      aes(label = etiqueta),
      position = "strand",
      angle = 30,
      nudge_y = 0.05,
      size = 4,
      fontface = "bold",
      check_overlap = FALSE
    ) +
    geom_feat(
      aes(colour = feature_type),
      linewidth = 3.5,
      position = "pile"
    ) +
    geom_feat_tag(
      aes(label = label, colour = feature_type),
      angle = 0,
      nudge_y = -0.14,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    annotate(
      "segment",
      x = distance_x1, xend = distance_x2,
      y = 1.47, yend = 1.47,
      linewidth = 0.55, colour = "#343E45"
    ) +
    annotate(
      "segment",
      x = c(distance_x1, distance_x2),
      xend = c(distance_x1, distance_x2),
      y = 1.43, yend = 1.51,
      linewidth = 0.55, colour = "#343E45"
    ) +
    annotate(
      "text",
      x = distance_mid, y = 1.57,
      label = paste0(format_bp(pair$distancia_bordes_bp[[1]]), " pb"),
      size = 3.8, fontface = "bold", colour = "#343E45"
    ) +
    annotate(
      "label",
      x = (region_end - region_start + 1L) * 0.02,
      y = 1.70,
      label = context_label,
      hjust = 0,
      size = 3.4,
      fontface = "bold",
      colour = "#34434C",
      fill = "#F3F6F7",
      linewidth = 0.25
    ) +
    annotate(
      "text",
      x = (region_end - region_start + 1L) * 0.98,
      y = 1.70,
      label = evidence_label,
      hjust = 1,
      size = 3.2,
      colour = "#4C5961"
    ) +
    scale_fill_manual(
      values = c("ARG" = "#D95F3B", "MRG" = "#16877A"),
      labels = c("ARG", "MRG"),
      name = NULL
    ) +
    scale_colour_manual(
      values = c(
        "IS (ISEScan)" = "#8C4A91",
        "Hit MGE (MOB-suite)" = "#E69F00"
      ),
      name = NULL,
      drop = FALSE
    ) +
    scale_x_continuous(
      labels = function(x) label_number(
        big.mark = ".", decimal.mark = ",", accuracy = 1
      )(x + coordinate_offset),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(ylim = c(0.62, 1.80), clip = "off") +
    labs(
      title = paste0(
        sample, ": ", pair$arg_symbol[[1]], " - ",
        pair$mrg_feature_label[[1]]
      ),
      x = paste0("Posicion en ", target_contig, " (pb)"),
      y = NULL,
      caption = paste0(
        "IS cercana = proximidad; no confirma transposon. ",
        "Hit MGE = coincidencia reportada por MOB-suite."
      )
    ) +
    theme_minimal(base_size = 12, base_family = "sans") +
    theme(
      plot.title = element_text(
        face = "bold", size = 16, hjust = 0.5, colour = "#24323B"
      ),
      plot.caption = element_text(size = 8.5, colour = "#59666E"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E1E6E9", linewidth = 0.4),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = margin(20, 26, 12, 26)
    )

  summary <- tibble(
    sample_id = sample,
    arg = pair$arg_symbol[[1]],
    mrg = pair$mrg_feature_label[[1]],
    contig = target_contig,
    distance_bp = pair$distancia_bordes_bp[[1]],
    molecule_type = molecule,
    plasmid_mobility = ifelse(nzchar(plasmid_mobility), plasmid_mobility, "NA"),
    is_in_region = nrow(is_features),
    is_within_10kb_of_gene = sum(is_features$distance_target <= 10000),
    mge_hits_in_region = nrow(mge_features),
    transposon_status = "No evaluable: sin MobileElementFinder"
  )

  list(plot = p, summary = summary)
}

atlas_pdf <- file.path(output_dir, "atlas_49_parejas_ARG_MRG_gggenomes.pdf")
grDevices::pdf(atlas_pdf, width = 14, height = 6.4, useDingbats = FALSE)

summaries <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  result <- make_plot(pairs[i, , drop = FALSE])
  print(result$plot)
  summaries[[i]] <- result$summary
  png_path <- file.path(
    output_dir,
    paste0(pairs$sample_id[[i]], "_ARG_MRG_gggenomes.png")
  )
  ggsave(
    png_path,
    result$plot,
    width = 14,
    height = 6.4,
    units = "in",
    dpi = 180,
    device = ragg::agg_png,
    background = "white"
  )
}
grDevices::dev.off()

summary_table <- bind_rows(summaries)
summary_path <- file.path(output_dir, "resumen_atlas_ARG_MRG.tsv")
write.table(
  summary_table,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Parejas representadas:", nrow(summary_table), "\n")
cat("PDF multipagina:", atlas_pdf, "\n")
cat("Resumen:", summary_path, "\n")
