# SERVER R: run ALL code in the SAME session after the GO-enabled script 34.
# OLD only. Compare stored MFF in bottom/top score quartiles, pooled AND per cell type.
# Quartiles are recalculated for each score and combined/scRNA/snRNA analysis separately.
# Keep score-boundary ties together and MFF=0. No score/mclust recalculation or RDS I/O.
local({
  packages <- c("Seurat", "ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. Fixed output files are overwritten on rerun.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/Old_score_quartiles_MFF")
  input.object <- "hlma.scores"
  groupings <- c("pooled", "celltype") # Both requested. Choose one here if needed.
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  group.order <- c("Bottom25", "Top25")
  group.colors <- c(Bottom25 = "#2878BD", Top25 = "#C54045")
  annotation.size <- 8
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L; n.columns <- 5L

  # 2. Pair stored scores and RNA/data MFF by exact original profile ID.
  make.pairs <- function(meta, mff, columns, ids = rownames(meta)) {
    required <- c("Annotation", "sample", "tech", "age_pop", "Country", unname(columns))
    if (!all(required %in% names(meta))) stop("Missing metadata/score columns: ", paste(setdiff(required, names(meta)), collapse = ", "))
    if (!nrow(meta) || !identical(rownames(meta), ids) || anyNA(ids) || anyDuplicated(ids) || any(ids == ""))
      stop("Invalid profile IDs/order.")
    if (!is.numeric(mff) || !identical(names(mff), ids) || length(mff) != length(ids) ||
        any(!is.finite(mff)) || any(mff < 0)) stop("Invalid MFF/profile pairing.")
    for (nm in c("Annotation", "sample", "tech", "age_pop", "Country"))
      if (anyNA(meta[[nm]]) || any(!nzchar(trimws(as.character(meta[[nm]]))))) stop("Missing metadata: ", nm)
    if (!all(meta$tech %in% c("scRNA", "snRNA")) || !all(meta$age_pop %in% c("young_pop", "old_pop")))
      stop("Unknown technology/age labels.")
    d <- data.frame(cell_id = ids, donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
      tech = as.character(meta$tech), Country = as.character(meta$Country),
      group = ifelse(meta$age_pop == "old_pop", "Old", "Young"), mff = unname(mff), stringsAsFactors = FALSE)
    for (nm in unname(columns)) {
      if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]]))) stop("Invalid score: ", nm)
      d[[nm]] <- meta[[nm]]
    }
    if (anyDuplicated(unique(d[, c("donor", "group", "Country")])$donor)) stop("Inconsistent donor age/Country.")
    d
  }
  quartile.groups <- function(scores, minimum = 4L) {
    if (!is.numeric(scores) || any(!is.finite(scores))) stop("Quartile scores must be finite/numeric.")
    n <- length(scores); group <- rep(NA_character_, n)
    q <- if (n >= minimum) unname(quantile(scores, c(.25, .75), type = 7)) else c(NA_real_, NA_real_)
    reason <- if (!n) "no_Old_profiles" else if (n < minimum) "fewer_than_four_Old_profiles" else
      if (q[1] == q[2]) "q25_equals_q75" else ""
    if (reason == "") {
      group[scores <= q[1]] <- "Bottom25"
      group[scores >= q[2]] <- "Top25"
      stopifnot(!any(scores <= q[1] & scores >= q[2]))
    }
    list(quartile = group, q25 = q[1], q75 = q[2], selection_reason = reason,
      n_ties_q25 = if (is.finite(q[1])) sum(scores == q[1]) else 0L,
      n_ties_q75 = if (is.finite(q[2])) sum(scores == q[2]) else 0L)
  }
  mff.test <- function(values, groups) {
    if (!is.numeric(values) || length(values) != length(groups) || any(!is.finite(values)) || any(values < 0) ||
        any(!is.na(groups) & !groups %in% c("Bottom25", "Top25"))) stop("Invalid MFF comparison.")
    b <- values[which(groups == "Bottom25")]; t <- values[which(groups == "Top25")]
    reason <- if (!length(b) || !length(t)) "missing_quartile_group" else
      if (min(length(b), length(t)) < 2L) "too_few_observations" else
        if (length(unique(c(b, t))) < 2L) "constant_MFF" else ""
    W <- p <- NA_real_
    if (reason == "") {
      w <- wilcox.test(b, t, alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)
      W <- unname(w$statistic); p <- w$p.value
    }
    data.frame(n_Bottom = length(b), n_Top = length(t), n_zero_Bottom = sum(b == 0), n_zero_Top = sum(t == 0),
      mean_Bottom = if (length(b)) mean(b) else NA_real_, mean_Top = if (length(t)) mean(t) else NA_real_,
      median_Bottom = if (length(b)) median(b) else NA_real_, median_Top = if (length(t)) median(t) else NA_real_,
      difference_Top_minus_Bottom = if (length(b) && length(t)) mean(t) - mean(b) else NA_real_,
      W = W, p_value = p, p_adj_BH = NA_real_, test_reason = reason)
  }
  quartile.tables <- function(cells, columns, scope.names, grouping.names, types) {
    if (!length(scope.names) || any(!scope.names %in% c("combined", "scRNA", "snRNA")) || anyDuplicated(scope.names) ||
        !length(grouping.names) || any(!grouping.names %in% c("pooled", "celltype")) || anyDuplicated(grouping.names))
      stop("Invalid scopes/groupings.")
    stats <- selected <- list()
    empty <- data.frame(scope = character(), grouping = character(), selection_celltype = character(), score = character(),
      quartile = character(), cell_id = character(), donor = character(), celltype = character(), tech = character(),
      Country = character(), group = character(), mff = numeric(), score_value = numeric(), stringsAsFactors = FALSE)
    for (scope in scope.names) for (grouping in grouping.names) {
      old <- cells[cells$group == "Old" & (scope == "combined" | cells$tech == scope), , drop = FALSE]
      targets <- if (grouping == "pooled") "ALL" else types
      for (ct in targets) for (nm in names(columns)) {
        d <- if (grouping == "pooled") old else old[old$celltype == ct, , drop = FALSE]
        v <- d[[columns[[nm]]]]; q <- quartile.groups(v); r <- mff.test(d$mff, q$quartile)
        n <- nrow(d)
        stats[[length(stats) + 1L]] <- cbind(data.frame(scope = scope, grouping = grouping,
          selection_celltype = ct, score = nm, n_Old = n, q25 = q$q25, q75 = q$q75,
          n_ties_q25 = q$n_ties_q25, n_ties_q75 = q$n_ties_q75,
          selection_reason = q$selection_reason), r,
          data.frame(n_Middle = sum(is.na(q$quartile)),
            fraction_Bottom = if (n) r$n_Bottom/n else NA_real_, fraction_Top = if (n) r$n_Top/n else NA_real_,
            n_donors_Bottom = length(unique(d$donor[which(q$quartile == "Bottom25")])),
            n_donors_Top = length(unique(d$donor[which(q$quartile == "Top25")]))))
        use <- which(!is.na(q$quartile))
        if (length(use)) selected[[length(selected) + 1L]] <- data.frame(scope = scope, grouping = grouping,
          selection_celltype = ct, score = nm, quartile = q$quartile[use],
          d[use, c("cell_id", "donor", "celltype", "tech", "Country", "group", "mff"), drop = FALSE],
          score_value = v[use], row.names = NULL)
      }
    }
    s <- do.call(rbind, stats)
    for (scope in scope.names) for (grouping in grouping.names) {
      idx <- which(s$scope == scope & s$grouping == grouping)
      s$p_adj_BH[idx] <- p.adjust(s$p_value[idx], method = "BH")
    }
    chosen <- if (length(selected)) do.call(rbind, selected) else empty
    rownames(s) <- rownames(chosen) <- NULL
    stopifnot(all(s$n_Bottom + s$n_Top + s$n_Middle == s$n_Old),
      !anyDuplicated(chosen[, c("scope", "grouping", "selection_celltype", "score", "cell_id")]))
    list(statistics = s, selected = chosen)
  }
  write.checked.csv <- function(d, path) {
    token <- ".__OLD_QUARTILE_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    con <- if (grepl("[.]gz$", path)) gzfile(path, "wt", encoding = "UTF-8") else file(path, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, con, row.names = FALSE, na = token), finally = close(con))
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else "character", character(1))
    con <- if (grepl("[.]gz$", path)) gzfile(path, "rt", encoding = "UTF-8") else file(path, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(con, colClasses = classes, na.strings = token, check.names = FALSE), finally = close(con))
    stopifnot(identical(dim(restored), dim(d)), identical(names(restored), names(d)))
    for (nm in names(d)) if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(restored[[nm]], d[[nm]], tolerance = 1e-10))) else
      stopifnot(identical(restored[[nm]], d[[nm]]))
    invisible(path)
  }

  # 3. Reuse script 34's scored object. Select Old only inside each analysis scope.
  if (!exists(input.object, envir = .GlobalEnv, inherits = FALSE)) stop("Run GO-enabled script 34 first: missing ", input.object)
  stopifnot(identical(names(score.columns), names(score.labels)))
  x <- get(input.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  ids <- colnames(x); rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == "MFF") == 1L)
  mff <- setNames(as.numeric(rna["MFF", ]), ids); rm(rna)
  cells <- make.pairs(x[[]], mff, score.columns, ids)
  types <- sort(unique(cells$celltype))
  info <- x@misc$HLMA_two_scores
  exclusion.confirmed <- !is.null(info) && "MFF" %in% info$pool_excludes &&
    all(names(score.columns) %in% names(info$matched_features)) &&
    !any(vapply(info$matched_features[names(score.columns)], function(v) "MFF" %in% v, logical(1)))
  if (!exclusion.confirmed) warning("MFF exclusion from score signatures/control pool is unconfirmed; verify script 34 provenance.")
  result <- quartile.tables(cells, score.columns, scopes, groupings, types)
  statistics <- result$statistics; selected <- result$selected
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(statistics, file.path(out, "01.Old_score_quartiles.MFF_statistics.csv"))
  write.checked.csv(selected, file.path(out, "02.Old_score_quartiles.selected_profiles.csv.gz"))
  write.checked.csv(cells[cells$group == "Old", ], file.path(out, "03.Old_all_profiles.reference.csv.gz"))

  # 4. Two violins per panel with narrow white boxes: bottom blue, top red.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  mff.violin <- function(d, stat, title) {
    d$quartile <- factor(d$quartile, levels = group.order)
    regular <- d[d$quartile %in% group.order[vapply(group.order, function(g) {
      v <- d$mff[d$quartile == g]; length(v) >= 2L && length(unique(v)) > 1L
    }, logical(1))], , drop = FALSE]
    special <- d[!d$quartile %in% regular$quartile, , drop = FALSE]
    bounds <- if (nrow(d)) range(d$mff) else c(0, .1); span <- max(diff(bounds), .05)
    bracket <- bounds[2] + .12 * span
    p <- ggplot(d, aes(quartile, mff)) + geom_blank() +
      geom_violin(data = regular, aes(fill = quartile), trim = TRUE, scale = "width", width = .8) +
      geom_boxplot(data = regular, width = .15, fill = "white", colour = "black",
        outlier.shape = NA, linewidth = .4)
    for (g in group.order) {
      z <- special[special$quartile == g, , drop = FALSE]
      if (nrow(z) == 1L) p <- p + geom_point(data = z, size = 2)
      if (nrow(z) > 1L) p <- p + geom_errorbar(data = z[1, , drop = FALSE], aes(ymin = mff, ymax = mff), width = .45)
    }
    labels <- setNames(paste0(c("Bottom 25%", "Top 25%"), "\n(n = ", c(stat$n_Bottom, stat$n_Top), ")"), group.order)
    p + scale_fill_manual(values = group.colors, guide = "none") +
      scale_x_discrete(limits = group.order, labels = labels, drop = FALSE) +
      annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .4) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket - .04 * span, yend = bracket, linewidth = .4) +
      annotate("text", x = 1.5, y = bounds[2] + .45 * span, label = format.p(stat$p_value), size = annotation.size) +
      coord_cartesian(ylim = c(max(0, bounds[1] - .08 * span), bounds[2] + .65 * span)) +
      labs(title = title, subtitle = if (stat$selection_reason != "") "Quartile groups unavailable" else NULL,
        x = NULL, y = "MFF expression (stored RNA/data)") +
      theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 13), axis.text = element_text(size = 18), panel.grid.minor = element_blank())
  }
  save.png <- function(p, path, width, height) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(path, width = width, height = height, res = plot.res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  panel <- function(scope, grouping, ct, score, title) {
    s <- statistics[statistics$scope == scope & statistics$grouping == grouping &
      statistics$selection_celltype == ct & statistics$score == score, , drop = FALSE]
    d <- selected[selected$scope == scope & selected$grouping == grouping &
      selected$selection_celltype == ct & selected$score == score, , drop = FALSE]
    stopifnot(nrow(s) == 1L)
    mff.violin(d, s, title)
  }
  caption <- paste("Old only. Bottom: score <= Q25; Top: score >= Q75. Boundary ties stay together, so n may exceed 25%.",
    "MFF zeros retained. Raw unpaired cell Wilcoxon p; donor/Country/composition not adjusted.")
  pngs <- character()
  for (scope in scopes) {
    if ("pooled" %in% groupings) {
      plots <- lapply(names(score.columns), function(nm) panel(scope, "pooled", "ALL", nm, score.labels[nm]))
      p <- wrap_plots(plots, ncol = length(score.columns), nrow = 1L) +
        plot_annotation(title = paste(scope.labels[scope], "Old | All cell types pooled | MFF by score quartiles", sep = " | "),
          caption = caption, theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "pooled", scope, "Old.MFF.by_three_score_quartiles.violin.png")
      save.png(p, path, panel.width * length(score.columns), panel.height); pngs <- c(pngs, path)
    }
    if ("celltype" %in% groupings) for (nm in names(score.columns)) {
      plots <- lapply(types, function(ct) panel(scope, "celltype", ct, nm, ct))
      cols <- min(n.columns, length(types))
      p <- wrap_plots(plots, ncol = cols) +
        plot_annotation(title = paste(scope.labels[scope], "Old", score.labels[nm], "MFF by within-cell-type score quartiles", sep = " | "),
          caption = caption, theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "celltype", scope, paste0("Old.MFF.by_", nm, ".quartiles.violin.overview.png"))
      save.png(p, path, panel.width * cols, panel.height * ceiling(length(types)/cols)); pngs <- c(pngs, path)
    }
  }
  expected <- length(scopes) * (as.integer("pooled" %in% groupings) + as.integer("celltype" %in% groupings) * length(score.columns))
  stopifnot(length(pngs) == expected, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  assign("hlma.old.score.quartiles", list(old_profiles = cells[cells$group == "Old", ],
    selected_profiles = selected, statistics = statistics, settings = list(groupings = groupings,
      scopes = scopes, score_columns = score.columns, quantile_type = 7L, boundary_ties = "included")), envir = .GlobalEnv)
  writeLines(c("Old-only MFF comparison by each previously calculated module score, both pooled and within cell type.",
    "Input: hlma.scores from GO-enabled script 34. No score/mclust recalculation, positive-MFF filtering or RDS I/O.",
    "Q25/Q75 use stats::quantile(type=7) on Old scores separately per score, scope and grouping/cell type.",
    "Bottom includes score <= Q25; Top includes score >= Q75. Ties remain together, so actual fractions may exceed 25%.",
    "If Q25==Q75, tails would overlap: both groups are unavailable and the selection reason is recorded.",
    "Fewer than 4 Old profiles: quartiles unavailable. Fewer than 2 profiles per tail/constant MFF: p=NA with reason.",
    "MFF is original RNA/data. Zero MFF and negative/zero module scores are preserved.",
    "Figures: Bottom25 blue versus Top25 red, with narrow white boxplots. Each panel has n and raw two-sided unpaired Wilcoxon p.",
    "Boxplots hide outlier markers only; all values remain in the violin/statistics. Zero IQR collapses the box to a line.",
    "Wilcoxon exact=FALSE/correct=TRUE. BH across all score/cell-type comparisons within each scope/grouping; CSV includes adjusted p.",
    "Pooled comparisons mix cell types. Cell p does not adjust donor dependence, Country, composition or technology.",
    "Combined and technical scopes, and pooled/cell-type analyses, overlap and are not independent replications.",
    "Default output: pooled 3 overview PNGs + celltype 9 overview PNGs = 12 PNGs.",
    "Three CSV/gzip tables are readback-verified: cutoffs/statistics, selected memberships and all Old reference profiles.",
    paste0("MFF scoring exclusion confirmed: ", exclusion.confirmed), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_OLD_SCORE_QUARTILES_MFF_COMPLETE\nOld profiles: ", sum(cells$group == "Old"),
    "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
