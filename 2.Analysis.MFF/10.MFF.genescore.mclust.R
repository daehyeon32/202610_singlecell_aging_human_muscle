# SERVER R: run after 33 (mclust annotation) and 34 (two scores), in the SAME session.
# Reuse original all/positive High/Low labels and full-reference scores. No RDS save.
local({
  packages <- c("Seurat", "ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/mclust_subsets_SenMayo_IFNalpha")
  score.object <- "hlma.scores"
  label.object <- "hlma.mff"  # Fallback if hlma.scores lacks the four classification columns.
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  class.columns <- c("MFF_mclust_all", "MFF_mclust_all_cutoff", "MFF_mclust_positive", "MFF_mclust_positive_cutoff")
  conditions <- c("all_high", "all_low", "positive_high", "positive_low")
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  group.order <- c("Young", "Old"); group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L; n.columns <- 5L

  align.labels <- function(meta, mff, labels, label.mff, ids = rownames(meta)) {
    shared <- c("Annotation", "sample", "tech", "age_pop", "Country")
    classification <- c("MFF_mclust_all", "MFF_mclust_all_cutoff", "MFF_mclust_positive", "MFF_mclust_positive_cutoff")
    lid <- rownames(labels)
    if (!identical(rownames(meta), ids) || anyNA(ids) || anyDuplicated(ids) || any(ids == "") ||
        is.null(lid) || anyNA(lid) || anyDuplicated(lid) || !setequal(ids, lid)) stop("Score/label cell IDs are not identical one-to-one sets.")
    if (!all(shared %in% names(meta)) || !all(c(shared, classification) %in% names(labels)))
      stop("Missing mclust labels/cutoffs or identity metadata; run script 33 first.")
    if (!identical(names(mff), ids) || !identical(names(label.mff), lid)) stop("MFF vector order mismatch.")
    labels <- labels[ids, , drop = FALSE]; label.mff <- label.mff[ids]
    for (nm in shared) if (!identical(as.character(meta[[nm]]), as.character(labels[[nm]])))
      stop("Score/label metadata differs: ", nm)
    if (any(!is.finite(mff)) || any(!is.finite(label.mff)) ||
        any(abs(mff - label.mff) > 1e-10 * pmax(1, abs(mff)))) stop("Score/label objects have different MFF expression.")
    meta[, classification] <- labels[, classification, drop = FALSE]
    meta
  }
  prepare.subsets <- function(meta, mff, columns, ids = rownames(meta)) {
    required <- c("Annotation", "sample", "tech", "age_pop", "Country", unname(columns),
      "MFF_mclust_all", "MFF_mclust_all_cutoff", "MFF_mclust_positive", "MFF_mclust_positive_cutoff")
    if (!all(required %in% names(meta))) stop("Missing metadata: ", paste(setdiff(required, names(meta)), collapse = ", "))
    if (!nrow(meta) || !identical(rownames(meta), ids) || anyNA(ids) || anyDuplicated(ids) || any(ids == "") ||
        !identical(names(mff), ids) || !is.numeric(mff) || any(!is.finite(mff)) || any(mff < 0)) stop("Invalid cell IDs/order or MFF values.")
    for (nm in c("Annotation", "sample", "tech", "age_pop", "Country"))
      if (anyNA(meta[[nm]]) || any(!nzchar(trimws(as.character(meta[[nm]]))))) stop("Missing metadata: ", nm)
    if (!all(meta$tech %in% c("scRNA", "snRNA")) || !all(meta$age_pop %in% c("young_pop", "old_pop"))) stop("Unknown tech/age_pop.")
    d <- data.frame(cell_id = ids, donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
      tech = as.character(meta$tech), Country = as.character(meta$Country),
      group = ifelse(meta$age_pop == "young_pop", "Young", "Old"), mff = unname(mff), stringsAsFactors = FALSE)
    for (nm in unname(columns)) {
      if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]]))) stop("Nonfinite/nonnumeric score: ", nm)
      d[[nm]] <- meta[[nm]]
    }
    if (anyDuplicated(unique(d[, c("donor", "group", "Country")])$donor)) stop("Inconsistent donor group/Country.")
    for (mode in c("all", "positive")) {
      column <- paste0("MFF_mclust_", mode); cut.column <- paste0(column, "_cutoff")
      label <- as.character(meta[[column]]); cutoff <- meta[[cut.column]]
      if (!is.numeric(cutoff) || any(!is.finite(cutoff) & !is.na(cutoff)) || any(cutoff < 0, na.rm = TRUE)) stop("Invalid cutoff: ", mode)
      for (ct in unique(d$celltype)) if (length(unique(cutoff[d$celltype == ct])) != 1L)
        stop("Cutoff varies within a cell type; expected the original combined fit: ", ct, " / ", mode)
      eligible <- if (mode == "all") rep(TRUE, nrow(d)) else d$mff > 0
      classified <- eligible & is.finite(cutoff)
      if (any(is.na(label[classified])) || any(!label[classified] %in% c("High", "Low")) ||
          any(!is.na(label[eligible & !is.finite(cutoff)]))) stop("Classification/cutoff mismatch: ", mode)
      if (mode == "positive" && any(is.na(label[!eligible]) | label[!eligible] != "Excluded_zero"))
        stop("Positive-mode MFF=0 must remain Excluded_zero.")
      delta <- d$mff[classified] - cutoff[classified]; tol <- 1e-10 * pmax(1, abs(cutoff[classified]))
      if (any((label[classified] == "High" & delta < -tol) | (label[classified] == "Low" & delta > tol)) ||
          any(label[classified] == "High" & d$mff[classified] == 0)) stop("Original High/Low labels disagree with MFF/cutoff: ", mode)
      d[[column]] <- label; d[[cut.column]] <- cutoff
    }
    masks <- list(all_high = !is.na(d$MFF_mclust_all) & d$MFF_mclust_all == "High",
      all_low = !is.na(d$MFF_mclust_all) & d$MFF_mclust_all == "Low",
      positive_high = !is.na(d$MFF_mclust_positive) & d$MFF_mclust_positive == "High",
      positive_low = !is.na(d$MFF_mclust_positive) & d$MFF_mclust_positive == "Low")
    list(cells = d, masks = masks)
  }
  subset.summaries <- function(cells, masks, columns, types, scope.names) {
    if (any(vapply(masks, function(v) length(v) != nrow(cells) || anyNA(v), logical(1)))) stop("Invalid subset masks.")
    statistics <- retention <- list()
    donor <- data.frame(condition = character(), scope = character(), score = character(), celltype = character(),
      donor = character(), group = character(), Country = character(), n_profiles = integer(), mean_score = numeric())
    for (condition in names(masks)) for (scope in scope.names) for (ct in types) {
      index <- cells$celltype == ct & (scope == "combined" | cells$tech == scope)
      full <- cells[index, , drop = FALSE]; d <- cells[index & masks[[condition]], , drop = FALSE]
      positive <- startsWith(condition, "positive")
      column <- if (positive) "MFF_mclust_positive" else "MFF_mclust_all"
      for (g in c("Young", "Old")) {
        z <- full[full$group == g, , drop = FALSE]; chosen <- d[d$group == g, , drop = FALSE]
        eligible <- if (positive) z$mff > 0 else rep(TRUE, nrow(z))
        retention[[length(retention) + 1L]] <- data.frame(condition = condition, scope = scope, celltype = ct, group = g,
          n_total = nrow(z), n_input_eligible = sum(eligible), n_selected = nrow(chosen),
          n_unclassified = sum(eligible & is.na(z[[column]])), n_excluded_zero = if (positive) sum(z$mff == 0) else 0L,
          n_donors_total = length(unique(z$donor)), n_donors_selected = length(unique(chosen$donor)),
          fraction_total = if (nrow(z)) nrow(chosen) / nrow(z) else NA_real_,
          fraction_eligible = if (sum(eligible)) nrow(chosen) / sum(eligible) else NA_real_)
      }
      for (nm in names(columns)) {
        d$score <- d[[columns[[nm]]]]
        ds <- donor[FALSE, , drop = FALSE]
        if (nrow(d)) for (idx in split(seq_len(nrow(d)), d$donor)) {
          z <- d[idx, , drop = FALSE]
          ds <- rbind(ds, data.frame(condition = condition, scope = scope, score = nm, celltype = ct,
            donor = z$donor[1], group = z$group[1], Country = z$Country[1], n_profiles = nrow(z), mean_score = mean(z$score)))
        }
        donor <- rbind(donor, ds)
        for (unit in c("cell", "donor")) {
          y <- if (unit == "cell") d$score[d$group == "Young"] else ds$mean_score[ds$group == "Young"]
          o <- if (unit == "cell") d$score[d$group == "Old"] else ds$mean_score[ds$group == "Old"]
          reason <- if (!length(y) || !length(o)) "missing_group" else if (min(length(y), length(o)) < 2L)
            "too_few_observations" else if (length(unique(c(y, o))) < 2L) "constant_score" else ""
          p <- W <- NA_real_
          if (reason == "") { w <- wilcox.test(y, o, alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE); p <- w$p.value; W <- unname(w$statistic) }
          my <- if (length(y)) mean(y) else NA_real_; mo <- if (length(o)) mean(o) else NA_real_
          statistics[[length(statistics) + 1L]] <- data.frame(condition = condition, scope = scope, celltype = ct,
            score = nm, unit = unit, n_Young = length(y), n_Old = length(o), n_profiles_Young = sum(d$group == "Young"),
            n_profiles_Old = sum(d$group == "Old"), n_donors_Young = sum(ds$group == "Young"), n_donors_Old = sum(ds$group == "Old"),
            mean_Young = my, mean_Old = mo, difference_Old_minus_Young = mo - my, W = W, p_value = p, p_adj_BH = NA_real_, reason = reason)
        }
      }
    }
    s <- do.call(rbind, statistics); k <- do.call(rbind, retention)
    for (condition in names(masks)) for (scope in scope.names) for (unit in c("cell", "donor")) {
      idx <- which(s$condition == condition & s$scope == scope & s$unit == unit)
      s$p_adj_BH[idx] <- p.adjust(s$p_value[idx], method = "BH")
    }
    rownames(s) <- rownames(k) <- rownames(donor) <- NULL
    list(statistics = s, retention = k, donor = donor)
  }
  write.checked.csv <- function(d, file) {
    token <- ".__HLMA_SUBSET_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    con <- if (grepl("[.]gz$", file)) gzfile(file, "wt", encoding = "UTF-8") else file(file, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, con, row.names = FALSE, na = token), finally = close(con))
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else if (is.logical(v)) "logical" else "character", character(1))
    con <- if (grepl("[.]gz$", file)) gzfile(file, "rt", encoding = "UTF-8") else file(file, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(con, colClasses = classes, check.names = FALSE, na.strings = token, stringsAsFactors = FALSE), finally = close(con))
    stopifnot(identical(dim(d), dim(restored)), identical(names(d), names(restored)))
    for (nm in names(d)) if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(d[[nm]], restored[[nm]], tolerance = 1e-10))) else stopifnot(identical(d[[nm]], restored[[nm]]))
    invisible(file)
  }

  # Reuse the existing full-reference scores. If needed, transfer labels by exact cell ID.
  if (!exists(score.object, envir = .GlobalEnv, inherits = FALSE)) stop("Run script 34 in this R session first: missing ", score.object)
  x <- get(score.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  ids <- colnames(x); meta <- x[[]]
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == "MFF") == 1L)
  mff <- setNames(as.numeric(rna["MFF", ]), ids); rm(rna)
  labels.from <- score.object
  if (!all(class.columns %in% names(meta))) {
    if (!exists(label.object, envir = .GlobalEnv, inherits = FALSE)) stop("Run script 33 first: missing classification columns and ", label.object)
    labelled <- get(label.object, envir = .GlobalEnv)
    stopifnot(inherits(labelled, "Seurat"), "RNA" %in% Assays(labelled), "data" %in% Layers(labelled[["RNA"]]))
    lrna <- GetAssayData(labelled, assay = "RNA", layer = "data"); lid <- colnames(labelled)
    stopifnot(identical(colnames(lrna), lid), sum(rownames(lrna) == "MFF") == 1L)
    lmff <- setNames(as.numeric(lrna["MFF", ]), lid); rm(lrna)
    meta <- align.labels(meta, mff, labelled[[]], lmff, ids)
    labels.from <- label.object
  }
  prepared <- prepare.subsets(meta, mff, score.columns, ids)
  cells <- prepared$cells; masks <- prepared$masks
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE)); n.rows <- ceiling(length(celltypes) / n.columns)
  cat("Reusing scores from ", score.object, "; mclust labels from ", labels.from, "; profiles: ", nrow(cells), "\n", sep = "")
  result <- subset.summaries(cells, masks, score.columns, celltypes, scopes)
  statistics <- result$statistics; donor <- result$donor; retention <- result$retention
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(statistics[statistics$unit == "cell", ], file.path(out, "01.mclust_scores.cell_statistics.csv"))
  write.checked.csv(statistics[statistics$unit == "donor", ], file.path(out, "02.mclust_scores.donor_statistics.csv"))
  write.checked.csv(donor, file.path(out, "03.mclust_scores.donor_means.csv"))
  write.checked.csv(retention, file.path(out, "04.mclust_scores.subset_retention.csv"))
  write.checked.csv(cells, file.path(out, "05.mclust_scores.profile_scores_labels.csv.gz"))

  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  score.violin <- function(d, full, stat, title, label) {
    d$group <- factor(d$group, levels = group.order)
    regular <- d[d$group %in% group.order[vapply(group.order, function(g) {
      v <- d$score[d$group == g]; length(v) >= 2L && length(unique(v)) > 1L
    }, logical(1))], , drop = FALSE]
    special <- d[!d$group %in% regular$group, , drop = FALSE]
    bounds <- if (nrow(full)) range(full$score) else c(-.1, .1); span <- max(diff(bounds), .05)
    bracket <- bounds[2] + .12 * span
    p <- ggplot(d, aes(group, score)) + geom_blank() +
      geom_violin(data = regular, aes(fill = group), trim = TRUE, scale = "width", width = .8)
    for (g in group.order) {
      z <- special[special$group == g, , drop = FALSE]
      if (nrow(z) == 1L) p <- p + geom_point(data = z, size = 2)
      if (nrow(z) > 1L) p <- p + geom_errorbar(data = z[1, , drop = FALSE], aes(ymin = score, ymax = score), width = .45)
    }
    axis.labels <- setNames(paste0(group.order, "\n(", c(stat$n_profiles_Young, stat$n_profiles_Old), ")"), group.order)
    p + scale_fill_manual(values = group.colors, guide = "none") +
      scale_x_discrete(limits = group.order, labels = axis.labels, drop = FALSE) +
      annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .4) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket - .04 * span, yend = bracket, linewidth = .4) +
      annotate("text", x = 1.5, y = bounds[2] + .28 * span, label = format.p(stat$p_value), size = 10) +
      coord_cartesian(ylim = c(bounds[1] - .08 * span, bounds[2] + .45 * span)) +
      labs(title = paste(strwrap(title, width = 40), collapse = "\n"), x = NULL, y = paste0(label, " score")) +
      theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 13), axis.text.y = element_text(size = 30), axis.text.x = element_text(size = 18), panel.grid.minor = element_blank())
  }
  score.difference <- function(s, limit, title, caption) {
    s <- s[order(is.na(s$difference_Old_minus_Young), s$difference_Old_minus_Young, s$celltype), ]
    s$celltype_axis <- factor(s$celltype, levels = rev(s$celltype))
    s$direction <- ifelse(s$difference_Old_minus_Young > 0, "Old higher", ifelse(s$difference_Old_minus_Young < 0, "Old lower", "No change"))
    s$significance <- ifelse(!is.finite(s$p_value), "Not tested", ifelse(s$p_value < .05, "p < 0.05", "p >= 0.05"))
    available <- s[is.finite(s$difference_Old_minus_Young), ]; absent <- s[!is.finite(s$difference_Old_minus_Young), ]
    ggplot(s, aes(y = celltype_axis)) + geom_vline(xintercept = 0, colour = "#999999", linetype = "dashed") +
      geom_point(data = available, aes(x = difference_Old_minus_Young, colour = direction, shape = significance), size = 4.5, stroke = 1.2) +
      geom_text(data = absent, x = limit * .95, label = "NA", hjust = 1, colour = "#888888") +
      scale_colour_manual(values = c("Old higher" = "#D73027", "Old lower" = "#000000", "No change" = "#777777"),
        limits = c("Old higher", "Old lower", "No change"), drop = FALSE, name = NULL) +
      scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1, "Not tested" = 4),
        limits = c("p < 0.05", "p >= 0.05", "Not tested"), drop = FALSE, name = "Unadjusted cell p") +
      scale_x_continuous(limits = c(-limit, limit)) + scale_y_discrete(drop = FALSE) +
      labs(title = title, x = "Cell-mean score difference (Old - Young)", y = NULL, caption = caption) +
      theme_classic(base_size = 14) + theme(plot.title = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 16), legend.position = "bottom",
        legend.text = element_text(size = 11), plot.caption = element_text(size = 9, hjust = 0), plot.margin = margin(15, 25, 12, 12))
  }
  save.png <- function(p, file, width, height, res = plot.res) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(file, width = width, height = height, res = res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  caption <- "Conditional cell Wilcoxon after MFF selection; unadjusted p; donor dependence, Country and technology are not adjusted."
  pngs <- character()
  for (condition in conditions) for (scope in scopes) for (nm in names(score.columns)) {
    scope.index <- if (scope == "combined") rep(TRUE, nrow(cells)) else cells$tech == scope
    b <- cells[scope.index & masks[[condition]], , drop = FALSE]
    reference <- cells[scope.index, , drop = FALSE]
    s <- statistics[statistics$condition == condition & statistics$scope == scope & statistics$score == nm & statistics$unit == "cell", ]
    plots <- list()
    for (ct in celltypes) {
      d <- b[b$celltype == ct, , drop = FALSE]; d$score <- d[[score.columns[[nm]]]]
      full <- reference[reference$celltype == ct, , drop = FALSE]; full$score <- full[[score.columns[[nm]]]]
      stat <- s[s$celltype == ct, ]; stopifnot(nrow(stat) == 1L)
      plots[[ct]] <- score.violin(d, full, stat, paste(scope.labels[scope], ct, sep = " | "), score.labels[nm])
    }
    folder <- file.path(out, condition, scope, nm)
    overview <- wrap_plots(plots, ncol = n.columns, nrow = n.rows) +
      plot_annotation(title = paste(score.labels[nm], condition, sep = " | "), caption = caption,
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 10, hjust = 0)))
    file <- file.path(folder, paste0(nm, ".", condition, ".Young_Old.violin.overview.png"))
    save.png(overview, file, panel.width * n.columns, panel.height * n.rows); pngs <- c(pngs, file)
    delta <- statistics$difference_Old_minus_Young[statistics$score == nm & statistics$unit == "cell"]
    limit <- max(.05, max(c(0, abs(delta[is.finite(delta)]))) * 1.15)
    p <- score.difference(s, limit, paste(scope.labels[scope], score.labels[nm], condition, sep = " | "), caption)
    file <- file.path(folder, paste0(nm, ".", condition, ".Young_Old.score_difference.png"))
    save.png(p, file, 3150L, round(max(5.5, 2 + length(celltypes) * .48) * 300), res = 300L); pngs <- c(pngs, file)
  }
  expected.pngs <- length(conditions) * length(scopes) * length(score.columns) * 2L
  stopifnot(length(pngs) == expected.pngs,
            all(file.exists(pngs)),
            all(file.info(pngs)$size > 8L))
  assign("hlma.mclust.scores", list(profiles = cells, masks = masks, statistics = statistics,
    donor_means = donor, retention = retention, label_source = labels.from), envir = .GlobalEnv)
  writeLines(c("Two PREVIOUSLY CALCULATED scores within four existing MFF mclust subsets, comparing Young/Old separately.",
    "all_high / all_low use MFF_mclust_all; positive_high / positive_low use MFF_mclust_positive.",
    "Original combined, cell-type-specific mclust cutoffs/labels are reused for all three technology scopes.",
    "Cutoffs, High/Low and scores are NOT recalculated within subset, age group or technology.",
    "Positive-mode zeros remain Excluded_zero; unavailable cutoff remains NA, never Low.",
    "Scores remain signed, including zero/negative scores; effect is arithmetic mean Old minus Young.",
    "5 columns x 3 rows for 15 types, 4500 x 2550 px at 160dpi. No boxplots.",
    "48 main PNGs = 4 subsets x 3 scopes x 2 scores x (violin overview + score-difference dotplot).",
    "Each cell-type/score/scope shares violin Y limits across the four subsets; each score's dotplots share X limits.",
    "Difference dot: positive red, negative black, zero grey. Raw cell p < .05 filled, >= .05 open, NA cross/label.",
    "Two-sided unpaired cell Wilcoxon exact=FALSE/correct=TRUE; at least 2 observations per age group.",
    "Missing/rare/constant groups have NA p with reason. These are conditional comparisons after selecting on MFF.",
    "Cell-level p ignores within-donor dependence, Country, technology and QC; not a test of the full population's age effect.",
    "Donor means and donor Wilcoxon are CSV only; observed donors have equal weight, absent donor/type is not a score zero.",
    "BH over both scores x all cell types separately within each condition x scope x unit.",
    "Retention includes original profiles, mode-eligible profiles, selected/unclassified/zero-excluded counts and donor coverage.",
    "The two all subsets partition classified all-mode profiles; the two positive subsets partition classified MFF-positive profiles.",
    "All and positive analyses overlap and are not independent validations. Scores use the same full-data reference from script 34.",
    "CSV/gzip outputs are readback-verified. No score calculation, mclust fit, normalization, correlation or RDS read/write.",
    paste0("Score object: ", score.object), paste0("Labels from: ", labels.from), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_MCLUST_TWO_SCORES_COMPLETE\nProfiles retained in reference: ", nrow(cells),
    "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
