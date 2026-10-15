# Human adaptation of revised 08.1 Mff.expression.youngold.R.
# Paste ALL code in SERVER R, or source a server copy.
# Main outputs: 3 columns x 5 rows of violins + cell-mean log2FC dots.
# Scopes: scRNA + snRNA, scRNA, snRNA. No counts/pseudobulk or cell exclusion.
local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Cairo)
  })

  # 1. Settings: panel sizes/fonts match revised mouse 08.1.
  gene <- "MFF"
  assay <- "RNA"
  layer <- "data"
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  rds.file <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  source.h5ad <- "/BiO/Store/UNIST-Aging-MFF-2026-09/HLMA_scRNA/OMIX004308-03.h5ad"
  output.root <- file.path(project, "analysis/3.Mffexpression")
  analysis.folder <- "young_old"     # Fixed folder; reruns overwrite matching output files.
  group.order <- c("Young", "Old")
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  min.cells.per.group <- 2L
  fc.pseudocount <- 1                 # Same formula as revised mouse 08.1.
  p.cutoff <- 0.05                    # Figures use raw cell-level p.
  expression.scale <- "auto"          # auto, ln1p, log2p1, or linear.
  allow.assumed.ln1p <- TRUE          # If log metadata unavailable: FLAG ln1p assumption.
  save.individual.violins <- TRUE
  panel.width <- 900L
  panel.height <- 850L
  plot.res <- 160L
  n.columns <- 5L
  n.rows <- 3L
  direction.colors <- c("Old lower" = "#D95F02", "Old higher" = "#000000", "No change" = "#777777")

  # 2. Input: preserve stored expression and author Young/Old labels.
  x <- if (exists("hlma.umap", envir = .GlobalEnv, inherits = FALSE)) {
    get("hlma.umap", envir = .GlobalEnv)
  } else {
    if (!file.exists(rds.file)) stop("RDS not found: ", rds.file)
    readRDS(rds.file)
  }
  stopifnot(inherits(x, "Seurat"), assay %in% Assays(x), layer %in% Layers(x[[assay]]))
  meta <- x[[]]
  required <- c("Annotation", "sample", "tech", "age_pop", "age", "Country")
  stopifnot(all(required %in% names(meta)), identical(rownames(meta), colnames(x)),
            !anyDuplicated(colnames(x)))
  for (nm in required) {
    z <- as.character(meta[[nm]])
    if (anyNA(z) || any(trimws(z) == "")) stop("Missing metadata: ", nm)
  }
  stopifnot(all(as.character(meta$tech) %in% c("scRNA", "snRNA")),
            all(as.character(meta$age_pop) %in% c("young_pop", "old_pop")))
  rna <- GetAssayData(x, assay = assay, layer = layer)
  stopifnot(sum(rownames(rna) == gene) == 1L, identical(colnames(rna), colnames(x)))
  cells <- data.frame(
    donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
    tech = as.character(meta$tech), age = suppressWarnings(as.numeric(as.character(meta$age))),
    Country = as.character(meta$Country),
    group = factor(ifelse(meta$age_pop == "young_pop", "Young", "Old"), levels = group.order),
    expression = as.numeric(rna[gene, ]), stringsAsFactors = FALSE
  )
  rm(rna)
  stopifnot(all(is.finite(cells$age)), all(is.finite(cells$expression)), all(cells$expression >= 0))
  donor.info <- unique(cells[, c("donor", "age", "Country", "group")])
  if (anyDuplicated(donor.info$donor)) stop("Inconsistent age/age_pop/Country within donor.")
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE))
  if (length(celltypes) != n.columns * n.rows) stop("Expected 15 cell types for the configured overview layout.")

  # Read log metadata only; this does not perform raw-count matching.
  import <- x@misc$HLMA_import
  if (!is.null(import$source_h5ad)) source.h5ad <- import$source_h5ad
  read.log.base <- function(path) {
    unavailable <- list(base = NA_real_, evidence = "no_usable_h5ad_log1p_metadata")
    if (!file.exists(path) || !requireNamespace("hdf5r", quietly = TRUE)) return(unavailable)
    h <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h$close_all(), add = TRUE)
    if (!h$exists("uns/log1p")) return(unavailable)
    g <- h[["uns/log1p"]]
    if (!inherits(g, "H5Group") || !g$exists("base")) {
      return(list(base = NA_real_, evidence = "log1p_record_base_omitted"))
    }
    b <- g[["base"]]
    encoding <- tryCatch(hdf5r::h5attr(b, "encoding-type"), error = function(e) "")
    if (identical(as.character(encoding), "null")) {
      return(list(base = exp(1), evidence = "h5ad_log1p_base_null_natural_log"))
    }
    value <- suppressWarnings(as.numeric(b$read()))
    if (length(value) == 1L && is.finite(value) && value > 1) {
      return(list(base = value, evidence = "h5ad_log1p_numeric_base"))
    }
    unavailable
  }
  stopifnot(expression.scale %in% c("auto", "ln1p", "log2p1", "linear"),
            is.finite(fc.pseudocount), fc.pseudocount > 0)
  log.info <- if (expression.scale == "auto") {
    tryCatch(read.log.base(source.h5ad), error = function(e) {
      list(base = NA_real_, evidence = paste0("metadata_read_error: ", conditionMessage(e)))
    })
  } else list(base = switch(expression.scale, ln1p = exp(1), log2p1 = 2, linear = NA_real_),
              evidence = paste0("manual_setting_", expression.scale))
  scale.assumed <- expression.scale != "auto"
  if (expression.scale == "auto" && !is.finite(log.info$base)) {
    if (!allow.assumed.ln1p) stop("Log base unconfirmed. Set a documented expression.scale before log2FC.")
    log.info$base <- exp(1)
    scale.assumed <- TRUE
    log.info$evidence <- paste0(log.info$evidence, "; assumed_ln1p_from_paper_and_Scanpy_default")
  }
  cells$normalized <- if (expression.scale == "linear") cells$expression else
    expm1(cells$expression * log(log.info$base))
  if (any(!is.finite(cells$normalized)) || any(cells$normalized < 0)) stop("Invalid inverse-transform for log2FC.")
  scale.note <- if (scale.assumed) paste0("FC scale: ",
    if (expression.scale == "auto") "ln1p assumed" else paste(expression.scale, "manual setting")) else
    paste0("FC log base from h5ad: ", signif(log.info$base, 5))
  cat("LOG2FC_SCALE: ", scale.note, "\nEvidence: ", log.info$evidence, "\n", sep = "")
  if (scale.assumed) warning("LOG2FC_SCALE_ASSUMED: ", scale.note,
    ". Verify before final interpretation; assumption is marked on FC plots.", call. = FALSE)
  warning("Cell-level Wilcoxon p-values do not adjust donor dependence, technology, country or batch; these are exploratory, not donor-level ageing tests.", call. = FALSE)

  # 3. Statistics and donor summaries, separately for each scope.
  selected <- list(combined = cells, scRNA = cells[cells$tech == "scRNA", ],
                   snRNA = cells[cells$tech == "snRNA", ])
  stopifnot(nrow(selected$scRNA) + nrow(selected$snRNA) == nrow(cells))
  test.cells <- function(d) {
    young <- d$expression[d$group == "Young"]
    old <- d$expression[d$group == "Old"]
    y <- d$normalized[d$group == "Young"]
    o <- d$normalized[d$group == "Old"]
    r <- data.frame(n_cells_Young = length(y), n_cells_Old = length(o),
      n_donors_Young = n_distinct(d$donor[d$group == "Young"]),
      n_donors_Old = n_distinct(d$donor[d$group == "Old"]),
      mean_stored_Young = if (length(y)) mean(young) else NA_real_,
      mean_stored_Old = if (length(o)) mean(old) else NA_real_,
      mean_normalized_Young = if (length(y)) mean(y) else NA_real_,
      mean_normalized_Old = if (length(o)) mean(o) else NA_real_,
      W = NA_real_, p_value = NA_real_, log2FC_Old_vs_Young = NA_real_,
      fc_pseudocount = fc.pseudocount, status = "not_tested", reason = "", stringsAsFactors = FALSE)
    # FC is independent of whether Wilcoxon is available.
    if (length(y) && length(o)) r$log2FC_Old_vs_Young <-
      log2(mean(o) + fc.pseudocount) - log2(mean(y) + fc.pseudocount)
    if (!length(y) || !length(o)) r$reason <- "missing_group" else
      if (min(length(y), length(o)) < min.cells.per.group) r$reason <- "insufficient_cells" else
        if (length(unique(c(young, old))) < 2L) r$reason <- "all_values_identical"
    if (nzchar(r$reason)) return(r)
    fit <- tryCatch(stats::wilcox.test(young, old, alternative = "two.sided",
      paired = FALSE, exact = FALSE, correct = TRUE), error = function(e) e)
    if (inherits(fit, "error")) r$reason <- paste0("test_error: ", conditionMessage(fit)) else
      if (!is.finite(fit$p.value)) r$reason <- "nonfinite_p" else {
        r$W <- unname(fit$statistic)
        r$p_value <- fit$p.value
        r$status <- "tested"
      }
    r
  }
  effects <- donor.tables <- summaries <- coverages <- list()
  for (scope in names(selected)) {
    b <- selected[[scope]]
    if (!nrow(b)) stop("No profiles for scope: ", scope)
    s <- bind_rows(lapply(celltypes, function(ct) {
      r <- test.cells(b[b$celltype == ct, ])
      r$celltype <- ct
      r$scope <- scope
      r$fc_scale_note <- scale.note
      r$fc_log_base <- if (expression.scale == "linear") NA_real_ else log.info$base
      r$fc_scale_assumed <- scale.assumed
      r$fc_scale_evidence <- log.info$evidence
      r
    }))
    s$p_adj_BH <- NA_real_
    ok <- is.finite(s$p_value)
    s$p_adj_BH[ok] <- p.adjust(s$p_value[ok], method = "BH")
    s$test_unit <- "cell_or_nucleus_profile"
    effects[[scope]] <- s
    observed <- b %>% group_by(donor, celltype) %>% summarise(
      n_cells = n(), n_MFF_positive = sum(expression > 0),
      n_scRNA = sum(tech == "scRNA"), n_snRNA = sum(tech == "snRNA"),
      mean_stored_data = mean(expression), mean_normalized = mean(normalized),
      MFF_positive_percent = 100 * mean(expression > 0), .groups = "drop")
    grid <- merge(unique(b[, c("donor", "age", "Country", "group")]),
                  data.frame(celltype = celltypes), by = NULL, sort = FALSE)
    donors <- left_join(grid, observed, by = c("donor", "celltype"))
    for (nm in c("n_cells", "n_MFF_positive", "n_scRNA", "n_snRNA")) donors[[nm]][is.na(donors[[nm]])] <- 0L
    donors$scope <- scope
    donor.tables[[scope]] <- donors
    summaries[[scope]] <- b %>% group_by(group, celltype) %>% summarise(
      n_cells = n(), n_donors = n_distinct(donor), n_scRNA = sum(tech == "scRNA"),
      n_snRNA = sum(tech == "snRNA"), mean_stored_data = mean(expression),
      mean_normalized = mean(normalized), MFF_positive_percent = 100 * mean(expression > 0),
      .groups = "drop") %>% mutate(scope = scope)
    coverages[[scope]] <- b %>% group_by(tech, group, Country) %>% summarise(
      n_cells = n(), n_donors = n_distinct(donor), .groups = "drop") %>% mutate(scope = scope)
    stopifnot(sum(donors$n_cells) == nrow(b), !anyDuplicated(donors[, c("donor", "celltype")]),
      all(is.na(donors$mean_stored_data[donors$n_cells == 0])),
      sum(s$n_cells_Young + s$n_cells_Old) == nrow(b))
    if (any(s$n_donors_Young < 2 | s$n_donors_Old < 2)) warning(scope,
      ": some cell types have fewer than two donors in a group; see statistics CSV.", call. = FALSE)
  }
  statistics <- bind_rows(effects)
  print(bind_rows(coverages))

  # 4. Fixed output folder: reruns overwrite matching files.
  out <- file.path(output.root, analysis.folder)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)
  cat("Output directory (overwrite): ", out, "\n", sep = "")
  tables <- list("01.MFF.cell_statistics" = statistics,
    "02.MFF.donor_celltype" = bind_rows(donor.tables),
    "03.MFF.Young_Old_summary" = bind_rows(summaries), "04.MFF.coverage" = bind_rows(coverages))
  for (nm in names(tables)) {
    z <- as.data.frame(lapply(tables[[nm]], function(v) if (is.factor(v)) as.character(v) else v))
    rownames(z) <- NULL
    file <- file.path(out, paste0(nm, ".csv"))
    write.csv(z, file, row.names = FALSE, fileEncoding = "UTF-8")
    saved <- read.csv(file, check.names = FALSE, fileEncoding = "UTF-8")
    stopifnot(isTRUE(all.equal(saved, z, check.attributes = FALSE, tolerance = 1e-8)))
  }
  save.png <- function(p, file, width, height, res = plot.res) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(file, width = width, height = height, res = res, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return(paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 2)))
    paste0("p = ", formatC(p, format = if (p < 0.001) "e" else "f",
                          digits = if (p < 0.001) 2 else 3))
  }

  # 5. Per-type violin + white box + bracket + raw Wilcoxon p.
  pngs <- character()
  for (scope in names(selected)) {
    b <- selected[[scope]]
    folder <- file.path(out, scope, "cell_level_violin")
    plots <- vector("list", length(celltypes))
    for (i in seq_along(celltypes)) {
      ct <- celltypes[i]
      d <- b[b$celltype == ct, ]
      top <- max(c(1, d$expression))
      regular <- d %>% group_by(group) %>% filter(n() >= 2L, n_distinct(expression) > 1L) %>% ungroup()
      special <- d[!d$group %in% regular$group, ] %>% group_by(group) %>%
        summarise(value = first(expression), n = n(), .groups = "drop")
      p <- ggplot(d, aes(group, expression)) + geom_blank() +
        geom_violin(data = regular, aes(fill = group), trim = TRUE, scale = "width", width = 0.8) +
        geom_boxplot(data = regular, width = 0.18, fill = "white", alpha = 0.85,
                     linewidth = 0.45, outlier.shape = NA)
      if (any(special$n > 1)) p <- p + geom_errorbar(data = special[special$n > 1, ],
        aes(group, ymin = value, ymax = value), inherit.aes = FALSE, width = 0.45, colour = "grey20")
      if (any(special$n == 1)) p <- p + geom_point(data = special[special$n == 1, ],
        aes(group, value), inherit.aes = FALSE, size = 2)
      title <- paste(strwrap(paste(scope.labels[scope], ct, sep = " | "), width = 40), collapse = "\n")
      p <- p + scale_fill_manual(values = group.colors, guide = "none") +
        annotate("segment", x = 1, xend = 2, y = top * 1.10, yend = top * 1.10, linewidth = 0.4) +
        annotate("segment", x = c(1, 2), xend = c(1, 2), y = top * 1.065, yend = top * 1.10, linewidth = 0.4) +
        annotate("text", x = 1.5, y = top * 1.26, label = format.p(effects[[scope]]$p_value[i]), size = 10) +
        scale_x_discrete(limits = group.order, drop = FALSE) + coord_cartesian(ylim = c(0, top * 1.42)) +
        labs(title = title, x = NULL, y = "MFF expression (stored RNA/data)") +
        theme_bw(base_size = 13) + theme(legend.position = "none",
          plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 13),
          axis.text = element_text(size = 30), panel.grid.minor = element_blank())
      plots[[i]] <- p
      if (save.individual.violins) {
        file <- file.path(folder, "violin", sprintf("%02d.%s.Young_Old.violin.png", i,
          gsub("[^A-Za-z0-9_-]+", "_", ct)))
        save.png(p, file, panel.width, panel.height)
        pngs <- c(pngs, file)
      }
    }
    p <- wrap_plots(plots, ncol = n.columns, nrow = n.rows) + plot_annotation(
      caption = "Cell-level Wilcoxon; unadjusted p; donor dependence and technology are not adjusted.",
      theme = theme(plot.caption = element_text(size = 10, hjust = 0)))
    file <- file.path(folder, "overview", "MFF.Young_Old.violin.overview.png")
    save.png(p, file, panel.width * n.columns, panel.height * n.rows)
    pngs <- c(pngs, file)
  }

  # 6. Cell-mean log2FC (Old / Young); shared symmetric X limits.
  fc <- statistics$log2FC_Old_vs_Young
  x.limit <- max(0.25, max(c(0, abs(fc[is.finite(fc)]))) * 1.15)
  breaks <- pretty(c(-x.limit, x.limit), n = 5)
  breaks <- sort(unique(c(0, breaks[abs(breaks) <= x.limit])))
  significance.labels <- c("Below cutoff" = paste0("Unadjusted p < ", p.cutoff),
    "At or above cutoff" = paste0("Unadjusted p >= ", p.cutoff), "Not tested" = "p = NA")
  for (scope in names(effects)) {
    d <- effects[[scope]]
    d <- d[order(!is.finite(d$log2FC_Old_vs_Young), d$log2FC_Old_vs_Young, d$celltype, na.last = TRUE), ]
    d$celltype_axis <- factor(d$celltype, levels = rev(d$celltype))
    d$direction <- ifelse(d$log2FC_Old_vs_Young < 0, "Old lower",
                         ifelse(d$log2FC_Old_vs_Young > 0, "Old higher", "No change"))
    d$significance <- ifelse(!is.finite(d$p_value), "Not tested",
                            ifelse(d$p_value < p.cutoff, "Below cutoff", "At or above cutoff"))
    available <- d[is.finite(d$log2FC_Old_vs_Young), ]
    absent <- d[!is.finite(d$log2FC_Old_vs_Young), ]
    p <- ggplot(d, aes(y = celltype_axis)) +
      geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.6, linetype = "dashed") +
      geom_point(data = available, aes(x = log2FC_Old_vs_Young, colour = direction, shape = significance),
                 size = 4.5, stroke = 1.2) +
      scale_colour_manual(values = direction.colors, guide = "none") +
      scale_shape_manual(values = c("Below cutoff" = 16, "At or above cutoff" = 1, "Not tested" = 4),
        limits = names(significance.labels), labels = unname(significance.labels), drop = FALSE, name = NULL) +
      guides(shape = guide_legend(nrow = 1, override.aes = list(colour = "#444444", size = 4))) +
      geom_text(data = absent, x = x.limit * 0.96, label = "NA", hjust = 1, colour = "#888888", size = 4.5) +
      scale_x_continuous(breaks = breaks, limits = c(-x.limit, x.limit), expand = expansion(mult = 0.02)) +
      scale_y_discrete(drop = FALSE, expand = expansion(add = 0.7)) +
      labs(title = paste(scope.labels[scope], gene, sep = " | "),
        x = "Cell-mean log2 fold change (Old / Young)", y = NULL,
        caption = paste("Cell-level Wilcoxon; unadjusted p; not donor-adjusted.", scale.note, "| mean + 1.")) +
      theme_classic(base_size = 14) + theme(
        plot.title = element_text(size = 22, face = "bold", margin = margin(b = 16)),
        axis.title.x = element_text(size = 18, margin = margin(t = 12)),
        axis.text.x = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black", margin = margin(r = 10)),
        axis.ticks.y = element_blank(), panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.35),
        legend.position = "bottom", legend.text = element_text(size = 13),
        legend.key.width = grid::unit(1, "cm"), plot.caption = element_text(size = 9, hjust = 0),
        plot.margin = margin(15, 25, 12, 12))
    file <- file.path(out, scope, "cell_level_violin", "log2FC_dotplot", "MFF.Young_vs_Old.cell_mean.log2FC.png")
    save.png(p, file, width = 3150L, height = round(max(5.5, 2 + nrow(d) * 0.48) * 300), res = 300L)
    pngs <- c(pngs, file)
  }
  expected.pngs <- 6L + if (save.individual.violins) 45L else 0L
  stopifnot(length(pngs) == expected.pngs, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  writeLines(c("Human Young/Old: revised mouse 08.1 cell-expression plotting/statistics only.",
    "Scopes: combined (scRNA+snRNA), scRNA, snRNA. Same original cells and age_pop; no new embedding.",
    paste0("Overview: ", n.columns, " columns x ", n.rows, " rows; ", panel.width * n.columns,
      " x ", panel.height * n.rows, " pixels at ", plot.res, " dpi. Individual panels: ", panel.width, " x ", panel.height, "."),
    paste0("Fixed output folder: ", out, ". Reruns overwrite matching files."),
    "Violin: unchanged stored RNA/data, zeros retained; flat/single groups shown as lines/points.",
    "Test: two-sided unpaired cell-level Wilcoxon rank-sum, exact=FALSE, continuity correction=TRUE.",
    "At least two profiles per group; identical values or missing groups yield p=NA with a reason.",
    "Plot p is unadjusted. BH across 15 cell types separately per scope is also saved in CSV.",
    "Cell p does not adjust within-donor dependence, technology/country/batch, or paired donors.",
    "Filled dots describe cell distributions, not proof of donor-level ageing effects.",
    "FC = log2(mean inverse-transformed Old + 1) - log2(mean inverse-transformed Young + 1).",
    "FC is pooled-cell weighted, not equal-donor weighted, not counts pseudobulk, not a mean-log ratio.",
    paste0("Expression scale: ", scale.note, ". Evidence: ", log.info$evidence),
    paste0("Source h5ad for log metadata: ", source.h5ad),
    "Log metadata documents a transform; it does not validate counts matching or remove technical differences.",
    "If ln1p is assumed, FC is provisional until the transform of X is verified.",
    "Dot: Old lower=orange, Old higher=black; raw p<0.05=filled, p>=0.05=open, p=NA with FC=cross.",
    "Missing FC=text NA. All three dotplots share X limits; cell types sorted by FC within scope.",
    "Combined donor table pools profiles within donor/celltype and counts each donor once.",
    "n_scRNA/n_snRNA show the mixture. Donor/type without cells has expression NA, not zero.",
    "Three scopes overlap and are not independent validations. No counts/pseudobulk.",
    "File checks do not replace visual inspection of server-rendered figures."), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nMFF_YOUNG_OLD_SUMMARY_COMPLETE\nCells:", nrow(cells), "; PNGs:", length(pngs),
      "; main PNGs: 6\n", out, "\n")
})


