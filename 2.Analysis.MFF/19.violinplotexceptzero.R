# Paste ALL code into SERVER R, or source a server copy.
# Full Young/Old violin: positive-cell density is red; the remaining zero-cell
# contribution is blue. Both use the SAME bandwidth, grid and full-input scale.
# No mclust cutoff is used. Exact zeros remain in the shape, boxplot, n and p.
local({
  # 1. Settings. Outputs go directly into the requested folder.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/3.Mffexpression")
  profiles.file <- file.path(out, "mclust_cutoffs/05.MFF.profile_classification.combined.csv.gz")
  rds.file <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  scopes <- c("combined", "scRNA", "snRNA")  # Use "scRNA" here for the attached view only.
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  group.order <- c("Young", "Old")
  region.colors <- c(Zero = "#2878BD", Positive = "#C54045")
  n.columns <- 5L
  panel.width <- 900L
  panel.height <- 850L
  plot.res <- 160L
  violin.width <- 0.8
  min.cells.per.group <- 2L
  stopifnot(length(scopes) > 0L, !anyDuplicated(scopes), all(scopes %in% names(scope.labels)))
  needed <- c("ggplot2", "patchwork", "Cairo")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

  # 2. Input checks and counts. MFF = 0 and MFF > 0 partition every profile.
  validate.profiles <- function(b) {
    required <- c("cell_id", "donor", "celltype", "tech", "group", "expression")
    if (!is.data.frame(b) || !all(required %in% names(b))) stop("Missing profile columns.")
    b <- b[, required, drop = FALSE]
    for (nm in setdiff(required, "expression")) {
      b[[nm]] <- as.character(b[[nm]])
      if (anyNA(b[[nm]]) || any(!nzchar(trimws(b[[nm]])))) stop("Missing profile metadata: ", nm)
    }
    if (anyDuplicated(b$cell_id)) stop("Duplicated profile IDs.")
    if (!is.numeric(b$expression) || any(!is.finite(b$expression)) || any(b$expression < 0))
      stop("MFF expression must be finite and nonnegative.")
    if (!all(b$tech %in% c("scRNA", "snRNA")) || !all(b$group %in% group.order))
      stop("Unrecognized technology or Young/Old label.")
    if (anyDuplicated(unique(b[, c("donor", "group")])$donor))
      stop("A donor has inconsistent Young/Old labels.")
    b
  }
  prepare.panel <- function(b, scope, ct) {
    keep <- b$celltype == ct
    if (scope != "combined") keep <- keep & b$tech == scope
    d <- b[keep, , drop = FALSE]
    young <- d$expression[d$group == "Young"]
    old <- d$expression[d$group == "Old"]
    p.value <- NA_real_; reason <- ""
    if (!length(young) || !length(old)) reason <- "missing_age_group" else
      if (min(length(young), length(old)) < min.cells.per.group) reason <- "insufficient_profiles" else
        if (length(unique(c(young, old))) < 2L) reason <- "identical_values"
    if (!nzchar(reason)) {
      test <- stats::wilcox.test(young, old, exact = FALSE, correct = TRUE,
        alternative = "two.sided", paired = FALSE)
      if (is.finite(test$p.value)) p.value <- test$p.value else reason <- "nonfinite_p"
    }
    counts <- do.call(rbind, lapply(group.order, function(g) {
      z <- d[d$group == g, , drop = FALSE]; positive <- z$expression > 0
      data.frame(scope = scope, celltype = ct, group = g,
        n_profiles = nrow(z), n_zero = sum(z$expression == 0), n_positive = sum(positive),
        n_donors = length(unique(z$donor)), n_positive_donors = length(unique(z$donor[positive])),
        zero_percent = if (nrow(z)) 100 * mean(z$expression == 0) else NA_real_,
        positive_percent = if (nrow(z)) 100 * mean(positive) else NA_real_,
        mean_stored_expression = if (nrow(z)) mean(z$expression) else NA_real_,
        mean_positive = if (any(positive)) mean(z$expression[positive]) else NA_real_,
        p_Young_vs_Old = p.value, test_reason = reason, stringsAsFactors = FALSE)
    }))
    stopifnot(all(counts$n_zero + counts$n_positive == counts$n_profiles))
    list(data = d, positive = d[d$expression > 0, , drop = FALSE], counts = counts,
      p_value = p.value, reason = reason, celltype = ct)
  }

  # The positive subdensity has total mass n_positive / n_all, not mass 1.
  # Drawing it over the complete blue shape leaves exactly the zero contribution.
  density.parts <- function(v, centre, width = 0.8) {
    empty.polygons <- data.frame(x = numeric(), y = numeric(), polygon = character(), part = character())
    empty.special <- data.frame(x = numeric(), value = numeric(), n = integer(), region = character())
    empty.curves <- data.frame(y = numeric(), full = numeric(), positive = numeric(), zero = numeric(),
      full_width = numeric(), positive_width = numeric())
    answer <- list(curves = empty.curves, polygons = empty.polygons, special = empty.special)
    if (!is.numeric(v) || any(!is.finite(v)) || any(v < 0)) stop("Invalid density input.")
    if (!length(v)) return(answer)
    if (length(unique(v)) == 1L) {
      answer$special <- data.frame(x = centre, value = v[1], n = length(v),
        region = if (v[1] == 0) "Zero" else "Positive")
      return(answer)
    }
    bw <- stats::bw.nrd0(v)
    if (!is.finite(bw) || bw <= 0) stop("Invalid full-input density bandwidth.")
    full <- stats::density(v, bw = bw, kernel = "gaussian", n = 512L, from = min(v), to = max(v))
    full.y <- full$y
    n.positive <- sum(v > 0)
    if (!n.positive) positive.y <- rep(0, length(full.y)) else
      if (n.positive == length(v)) positive.y <- full.y else {
        weights <- as.numeric(v > 0) / n.positive
        pos <- stats::density(v, weights = weights, bw = bw, kernel = "gaussian", n = 512L,
          from = min(v), to = max(v))
        # FFT rounding can produce tiny negative residuals: constrain them to the full shape.
        positive.y <- pmin(full.y, pmax(0, pos$y * n.positive / length(v)))
      }
    half.full <- width / 2 * full.y / max(full.y)
    half.positive <- width / 2 * positive.y / max(full.y)
    answer$curves <- data.frame(y = full$x, full = full.y, positive = positive.y,
      zero = full.y - positive.y, full_width = half.full, positive_width = half.positive)
    polygon <- function(half, part) data.frame(x = c(centre - half, rev(centre + half)),
      y = c(full$x, rev(full$x)), polygon = paste(centre, part, sep = "_"), part = part)
    answer$polygons <- polygon(half.full, "Full")
    if (n.positive) answer$polygons <- rbind(answer$polygons, polygon(half.positive, "Positive"))
    answer
  }
  write.checked.csv <- function(z, path) {
    z <- as.data.frame(lapply(z, function(v) if (is.factor(v)) as.character(v) else v))
    rownames(z) <- NULL
    token <- ".__MFF_MISSING__."
    strings <- unlist(z[vapply(z, is.character, logical(1))], use.names = FALSE)
    while (any(strings == token, na.rm = TRUE)) token <- paste0(token, "_")
    write.csv(z, path, row.names = FALSE, na = token)
    saved <- read.csv(path, check.names = FALSE, na.strings = token,
      colClasses = vapply(z, function(v) class(v)[1], character(1)))
    stopifnot(identical(names(saved), names(z)), nrow(saved) == nrow(z))
    for (nm in names(z)) stopifnot(isTRUE(all.equal(saved[[nm]], z[[nm]],
      check.attributes = FALSE, tolerance = 1e-8)))
  }
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p < 2.2e-16) return("p < 2.2e-16")
    paste0("p = ", format(signif(p, 3), scientific = p < 0.001, trim = TRUE))
  }
  plot.panel <- function(panel, show.legend = TRUE) {
    d <- panel$data; d$position <- match(d$group, group.order)
    parts <- lapply(seq_along(group.order), function(i)
      density.parts(d$expression[d$group == group.order[i]], i, violin.width))
    polygons <- do.call(rbind, lapply(parts, `[[`, "polygons"))
    special <- do.call(rbind, lapply(parts, `[[`, "special"))
    outline <- polygons[polygons$part == "Full", , drop = FALSE]
    positive <- polygons[polygons$part == "Positive", , drop = FALSE]
    top <- max(c(1, d$expression))
    labels <- paste0(group.order, "\n(n=", format(panel$counts$n_profiles, big.mark = ",", trim = TRUE), ")")
    zero.labels <- ifelse(is.finite(panel$counts$zero_percent),
      paste0(format(round(panel$counts$zero_percent, 1), trim = TRUE), "%"), "NA")
    p <- ggplot(d, aes(position, expression)) + geom_blank() +
      geom_polygon(data = outline, aes(x, y, group = polygon, fill = "Zero"),
        inherit.aes = FALSE, colour = NA, show.legend = show.legend) +
      geom_polygon(data = positive, aes(x, y, group = polygon, fill = "Positive"),
        inherit.aes = FALSE, colour = NA, show.legend = show.legend) +
      geom_polygon(data = outline, aes(x, y, group = polygon), inherit.aes = FALSE,
        fill = NA, colour = "grey20", linewidth = 0.45)
    if (nrow(d)) p <- p + geom_boxplot(aes(group = group), width = 0.18,
      fill = "white", alpha = 0.85, colour = "grey20", linewidth = 0.45, outlier.shape = NA)
    if (nrow(special)) p <- p + geom_errorbar(data = special,
      aes(x, ymin = value, ymax = value, colour = region), inherit.aes = FALSE, width = 0.45, linewidth = 0.8)
    if (any(special$n == 1L)) p <- p + geom_point(data = special[special$n == 1L, , drop = FALSE],
      aes(x, value, colour = region), inherit.aes = FALSE, size = 2)
    p <- p + scale_fill_manual(name = NULL, values = region.colors, limits = names(region.colors),
      labels = c("MFF = 0 contribution", "MFF > 0 contribution"), drop = FALSE) +
      scale_colour_manual(values = region.colors, guide = "none") +
      scale_x_continuous(breaks = 1:2, labels = labels, limits = c(0.5, 2.5)) +
      annotate("segment", x = 1, xend = 2, y = top * 1.1, yend = top * 1.1, linewidth = 0.4) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = top * 1.065, yend = top * 1.1, linewidth = 0.4) +
      annotate("text", x = 1.5, y = top * 1.26, label = format.p(panel$p_value), size = 7) +
      coord_cartesian(ylim = c(0, top * 1.42)) +
      labs(title = panel$celltype,
        subtitle = paste0("MFF = 0: Young ", zero.labels[1], " | Old ", zero.labels[2]),
        x = NULL, y = "MFF expression (stored RNA/data)") + theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 20), panel.grid.minor = element_blank())
    if (!nrow(d)) p <- p + annotate("text", x = 1.5, y = top * 0.5, label = "No input profiles", size = 4)
    if (!show.legend) p <- p + guides(fill = "none")
    p
  }

  # 3. Reuse the exact profiles behind the existing figure; original RDS is a fallback.
  if (file.exists(profiles.file)) {
    con <- gzfile(profiles.file, "rt")
    b <- tryCatch(read.csv(con, colClasses = "character", na.strings = character(), check.names = FALSE),
      finally = close(con))
    if (!"expression" %in% names(b)) stop("Saved profiles lack expression.")
    b$expression <- suppressWarnings(as.numeric(b$expression))
    input.source <- profiles.file
  } else {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Saved profiles absent; install Seurat for RDS fallback.")
    x <- if (exists("hlma.umap", envir = .GlobalEnv, inherits = FALSE)) get("hlma.umap", envir = .GlobalEnv) else {
      if (!file.exists(rds.file)) stop("Neither saved profiles nor original RDS found.")
      readRDS(rds.file)
    }
    if (!inherits(x, "Seurat")) stop("Input is not a Seurat object.")
    meta <- x[[]]
    required <- c("sample", "Annotation", "tech", "age_pop")
    if (!all(required %in% names(meta)) || !identical(rownames(meta), colnames(x))) stop("Invalid RDS metadata/IDs.")
    if (anyNA(meta$age_pop) || !all(as.character(meta$age_pop) %in% c("young_pop", "old_pop")))
      stop("Invalid RDS age_pop labels.")
    rna <- Seurat::GetAssayData(x, assay = "RNA", layer = "data")
    if (sum(rownames(rna) == "MFF") != 1L || !identical(colnames(rna), colnames(x))) stop("MFF/RNA IDs do not match.")
    b <- data.frame(cell_id = colnames(x), donor = as.character(meta$sample),
      celltype = as.character(meta$Annotation), tech = as.character(meta$tech),
      group = ifelse(meta$age_pop == "young_pop", "Young", "Old"), expression = as.numeric(rna["MFF", ]))
    rm(rna)
    input.source <- "Existing hlma.umap or original Seurat RDS"
  }
  cells <- validate.profiles(b)
  if (!nrow(cells)) stop("No profiles to plot.")
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE))
  n.rows <- ceiling(length(celltypes) / n.columns)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output folder: ", out)

  # 4. Full-input Young/Old statistics and one 5-column overview per scope.
  count.rows <- list(); pngs <- character()
  for (scope in scopes) {
    panels <- lapply(celltypes, function(ct) prepare.panel(cells, scope, ct))
    counts <- do.call(rbind, lapply(panels, `[[`, "counts"))
    pv <- vapply(panels, `[[`, numeric(1), "p_value")
    counts$p_adj_BH <- rep(stats::p.adjust(pv, method = "BH"), each = length(group.order))
    count.rows[[scope]] <- counts
    plots <- lapply(seq_along(panels), function(i) plot.panel(panels[[i]], show.legend = i == 1L))
    overview <- (wrap_plots(plots, ncol = n.columns, nrow = n.rows, guides = "collect") +
      plot_annotation(title = paste(scope.labels[scope], "| MFF = 0 / MFF > 0 | Young / Old"),
        caption = paste0("Full violin, boxplot, n and p include all profiles. Red: positive-cell density; blue: remaining zero-cell density. Same bandwidth and full-input width scale.\n",
          "KDE spreads the zero contribution around 0; colors indicate density contributions, not a cutoff on the y-axis. Cell Wilcoxon: raw p, not donor-adjusted."),
        theme = theme(plot.title = element_text(face = "bold", size = 26),
          plot.caption = element_text(size = 10, hjust = 0)))) & theme(legend.position = "bottom")
    filename <- file.path(out, paste0("MFF.zero_positive.young_old.violin.overview.", scope, ".png"))
    Cairo::CairoPNG(filename = filename, width = panel.width * n.columns,
      height = panel.height * n.rows, res = plot.res, bg = "white")
    tryCatch(print(overview), finally = grDevices::dev.off())
    stopifnot(file.exists(filename), file.info(filename)$size > 1000)
    pngs <- c(pngs, filename)
  }
  statistics <- do.call(rbind, count.rows)
  write.checked.csv(statistics, file.path(out, "MFF.zero_positive.Young_Old.counts.csv"))
  writeLines(c("Blue = zero-cell contribution; red = positive-cell contribution to the complete age-specific KDE.",
    "Positive density is weighted by n_positive/n_all and uses the full-input bandwidth/grid/width scale.",
    "No new expression threshold, mclust fit, normalization or zero removal.",
    "Boxplots and raw cell-level Young/Old Wilcoxon use all profiles, including exact zeros.",
    "Constant/singleton inputs use lines/points; missing age/type combinations remain n=0 and p=NA.",
    "BH is applied across cell types separately within each scope; original and adjusted p are in CSV.",
    "These p-values do not account for donor dependence, country, technology or batch.",
    paste("Input:", input.source), paste("Scopes:", paste(scopes, collapse = ", "))),
    file.path(out, "MFF.zero_positive.README.txt"))
  cat("\nMFF_ZERO_POSITIVE_VIOLINS_COMPLETE\nProfiles:", nrow(cells), "| PNGs:", length(pngs),
    "\nOutput:", out, "\n", paste(pngs, collapse = "\n"), "\n")
})
