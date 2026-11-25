#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(file.path(getwd(), "analisis_arg_mrg_20261022"), mustWork = TRUE)
mode <- if (length(args) >= 2) args[2] else "all"
tabdir <- file.path(base, "tablas")
figdir <- file.path(base, "figuras")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

freq <- read.delim(file.path(tabdir, "frecuencias_arg_mrg.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
mat <- read.delim(file.path(tabdir, "matriz_presencia_ausencia.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

meta_cols <- c("sample_id", "analysis_amrfinder", "analysis_diamond", "taxonomia", "taxonomia_grupo", "procedencia", "fenotipo_cobre")
feature_cols <- setdiff(names(mat), meta_cols)

freq$key <- paste(freq$block, freq$feature, sep = ":")
freq <- freq[order(freq$block, -freq$positive_isolates_n, freq$feature), ]
arg_order <- freq$key[freq$block == "ARG"]
mrg_part <- freq[freq$block == "MRG", , drop = FALSE]
mrg_part <- mrg_part[order(mrg_part$category, -mrg_part$positive_isolates_n, mrg_part$feature), ]
mrg_order <- mrg_part$key
feature_cols <- c(arg_order, mrg_order)
feature_cols <- feature_cols[feature_cols %in% names(mat)]
arg_cols <- feature_cols[startsWith(feature_cols, "ARG:")]
mrg_cols <- feature_cols[startsWith(feature_cols, "MRG:")]

short_tax <- function(x) {
  x <- gsub("^1a: ", "DISCORD: ", x)
  x <- gsub("; 2a: ", " / ", x)
  ifelse(nchar(x) > 42, paste0(substr(x, 1, 39), "..."), x)
}

tax_order <- order(mat$taxonomia_grupo, mat$taxonomia, mat$sample_id)
mat <- mat[tax_order, , drop = FALSE]
rownames(mat) <- paste(mat$sample_id, short_tax(mat$taxonomia), sep = " | ")

ann <- data.frame(
  Procedencia = mat$procedencia,
  Cobre = mat$fenotipo_cobre,
  stringsAsFactors = FALSE,
  row.names = rownames(mat)
)

origin_levels <- sort(unique(ann$Procedencia))
origin_palette <- setNames(grDevices::hcl.colors(length(origin_levels), "Dark 3"), origin_levels)
copper_palette <- c("Alta" = "#7F0000", "Media" = "#D95F0E", "Baja" = "#FEC44F", "No evaluado/no informado" = "#D9D9D9")
ann_colors <- list(
  Procedencia = origin_palette,
  Cobre = copper_palette,
  Bloque = c("ARG" = "#B2182B", "MRG" = "#2166AC")
)

make_heatmap <- function(
  cols, filename_stem, title, width, height, fontsize_col,
  fontsize_row = 7.5, cellwidth = NA, cellheight = 8.5, angle_col = 45,
  show_annotation_names = FALSE,
  row_annotation = ann, annotation_palette = ann_colors
) {
  if (length(cols) == 0) return(invisible(NULL))
  x <- as.matrix(mat[, cols, drop = FALSE])
  storage.mode(x) <- "numeric"
  display_colnames <- sub("^[^:]+:", "", colnames(x))
  block <- sub(":.*$", "", cols)
  gaps <- if (length(unique(block)) > 1) max(which(block == "ARG")) else NULL
  ann_col <- data.frame(Bloque = block, row.names = cols, stringsAsFactors = FALSE)
  palette <- colorRampPalette(c("#F2F2F2", "#2166AC"))(2)
  draw <- function(path, device) {
    if (device == "pdf") pdf(path, width = width, height = height, useDingbats = FALSE)
    else ragg::agg_png(path, width = width, height = height, units = "in", res = 300)
    pheatmap(
      x, cluster_rows = FALSE, cluster_cols = FALSE,
      color = palette, breaks = c(-0.01, 0.5, 1.01),
      annotation_row = row_annotation,
      annotation_col = ann_col,
      annotation_colors = annotation_palette,
      annotation_names_row = show_annotation_names,
      gaps_col = gaps, border_color = "#DDDDDD",
      cellwidth = cellwidth, cellheight = cellheight,
      fontsize = 9.5, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
      angle_col = angle_col, labels_col = display_colnames, main = title,
      legend_breaks = c(0, 1), legend_labels = c("No detectado", "Detectado"),
      na_col = "#8C8C8C"
    )
    dev.off()
  }
  draw(file.path(figdir, paste0(filename_stem, ".pdf")), "pdf")
  draw(file.path(figdir, paste0(filename_stem, ".png")), "png")
}

# Resumen: >=5% de prevalencia; si un bloque excede 35, conserva los 35 más frecuentes.
select_summary <- function(block) {
  z <- freq[freq$block == block & freq$positive_isolates_n >= 5, , drop = FALSE]
  if (nrow(z) > 35) z <- z[seq_len(35), , drop = FALSE]
  z$key
}
summary_cols <- c(select_summary("ARG"), select_summary("MRG"))
summary_cols <- summary_cols[summary_cols %in% feature_cols]
if (length(summary_cols) == 0) summary_cols <- feature_cols[seq_len(min(50, length(feature_cols)))]
make_heatmap(summary_cols, "figura_A_heatmap_resumen", "ARG y MRG por aislado - resumen (prevalencia >=5%; max. 35 por bloque)", 14.5, 15.5, 9, 7.5, 11, 8.5, 45)

# Heatmap con contexto genómico y señal de movilidad potencial. La pareja
# representativa es la pareja ARG-MRG confiable más cercana de cada aislado.
context_path <- file.path(base, "colocalizacion", "contexto_movil_por_muestra.tsv")
if (file.exists(context_path)) {
  context <- read.delim(
    context_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  context <- context[match(mat$sample_id, context$sample_id), , drop = FALSE]
  if (any(is.na(context$sample_id))) {
    stop("Faltan muestras en contexto_movil_por_muestra.tsv")
  }

  ann_context <- data.frame(
    Procedencia = mat$procedencia,
    Cobre = mat$fenotipo_cobre,
    `Dist. ARG-MRG` = factor(
      context$arg_mrg_distance_class,
      levels = c("Sin pareja confiable", ">50 kb", "10-50 kb", "<=10 kb")
    ),
    Molecula = factor(
      context$genomic_context,
      levels = c("Sin pareja confiable", "Cromosoma", "Plasmido", "No clasificado")
    ),
    `IS cercana` = factor(
      context$is_proximity,
      levels = c(
        "Sin pareja confiable", "Sin IS <=50 kb", "IS 10-50 kb", "IS <=10 kb"
      )
    ),
    Integron = factor(
      context$integron_status,
      levels = c(
        "Sin integron", "In0", "CALIN", "Integron completo", "No analizado"
      )
    ),
    MGE = factor(
      context$mge_status,
      levels = c("Sin MGE", "MGE detectado", "No analizado")
    ),
    `MGE-ARG/MRG` = factor(
      context$mef_context,
      levels = c(
        "Sin pareja confiable", "No analizado", "Sin MGE <=50 kb",
        "MGE 10-50 kb", "MGE <=10 kb",
        "ARG o MRG en MGE", "ARG y MRG en MGE"
      )
    ),
    check.names = FALSE,
    row.names = rownames(mat)
  )

  ann_colors_context <- c(
    ann_colors,
    list(
      `Dist. ARG-MRG` = c(
        "Sin pareja confiable" = "#E6E6E6",
        ">50 kb" = "#BDBDBD",
        "10-50 kb" = "#FDB863",
        "<=10 kb" = "#D73027"
      ),
      Molecula = c(
        "Sin pareja confiable" = "#E6E6E6",
        "Cromosoma" = "#2C7FB8",
        "Plasmido" = "#7A0177",
        "No clasificado" = "#969696"
      ),
      `IS cercana` = c(
        "Sin pareja confiable" = "#E6E6E6",
        "Sin IS <=50 kb" = "#F7F7F7",
        "IS 10-50 kb" = "#FDAE61",
        "IS <=10 kb" = "#D01C8B"
      ),
      Integron = c(
        "Sin integron" = "#F7F7F7",
        "In0" = "#D8B365",
        "CALIN" = "#C2A5CF",
        "Integron completo" = "#762A83",
        "No analizado" = "#969696"
      ),
      MGE = c(
        "Sin MGE" = "#F7F7F7",
        "MGE detectado" = "#1B7837",
        "No analizado" = "#969696"
      ),
      `MGE-ARG/MRG` = c(
        "Sin pareja confiable" = "#E6E6E6",
        "No analizado" = "#969696",
        "Sin MGE <=50 kb" = "#F7F7F7",
        "MGE 10-50 kb" = "#A6DBA0",
        "MGE <=10 kb" = "#5AAE61",
        "ARG o MRG en MGE" = "#1B7837",
        "ARG y MRG en MGE" = "#00441B"
      )
    )
  )

  make_heatmap(
    summary_cols,
    "figura_A_heatmap_resumen_contexto_movil",
    paste0(
      "ARG y MRG por aislado: contexto movil e integrones"
    ),
    15.2, 14.7, 9, 7.5, 11, 8.5, 45, TRUE,
    row_annotation = ann_context,
    annotation_palette = ann_colors_context
  )
}

# Versión completa suplementaria; ancho adaptado al número de columnas.
full_width <- 17.5
make_heatmap(feature_cols, "figura_A_heatmap_completo_suplementario", "ARG y MRG por aislado - todas las caracteristicas", full_width, 15.5, 7.5, 7.5, 6.5, 8.5, 90)

# Figura de frecuencia: máximo 25 características por bloque; tabla TSV conserva todas.
if (mode != "heatmap-only") {
topfreq <- do.call(rbind, lapply(c("ARG", "MRG"), function(b) {
  z <- freq[freq$block == b, , drop = FALSE]
  z[seq_len(min(25, nrow(z))), , drop = FALSE]
}))
topfreq$panel <- factor(topfreq$block, levels = c("ARG", "MRG"), labels = c("ARG (AMRFinderPlus)", "MRG (BacMet/DIAMOND)"))
topfreq$label <- sprintf("%d/%d (%.1f%%)", topfreq$positive_isolates_n, topfreq$analyzed_isolates_N, topfreq$frequency_pct)
topfreq$fill_group <- ifelse(topfreq$block == "ARG", "ARG", ifelse(topfreq$category == "", "Metal no resuelto", topfreq$category))
topfreq$feature_panel <- paste(topfreq$feature, topfreq$block, sep = "___")
topfreq$feature_panel <- factor(topfreq$feature_panel, levels = rev(topfreq$feature_panel))

fill_levels <- unique(topfreq$fill_group)
fill_colors <- setNames(grDevices::hcl.colors(length(fill_levels), "Set 2"), fill_levels)
fill_colors["ARG"] <- "#B2182B"

p <- ggplot(topfreq, aes(x = positive_isolates_n, y = feature_panel, fill = fill_group)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = label), hjust = -0.05, size = 3) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(values = fill_colors, name = "Categoria") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(
    title = "Frecuencia descriptiva de ARG y MRG en la coleccion",
    subtitle = "Hasta 25 caracteristicas por panel; cada aislado cuenta una vez por caracteristica",
    x = "Aislados positivos (n)", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 8), legend.position = "bottom"
  )

ggsave(file.path(figdir, "figura_B_frecuencias_ARG_MRG.pdf"), p, width = 12, height = 15, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(figdir, "figura_B_frecuencias_ARG_MRG.png"), p, width = 12, height = 15, dpi = 300, device = ragg::agg_png)
}

writeLines(capture.output(sessionInfo()), file.path(base, "qc", "sessionInfo_R.txt"))
