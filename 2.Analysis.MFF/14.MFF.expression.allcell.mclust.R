# SERVER R: run after 33 (mclust labels) and 34 (scores), in the SAME session.
# Four original mclust subsets; ALL cell types pooled within each subset.
# MFF + three-score Young/Old violins and three-score MFF correlations. No refitting or RDS I/O.
local({
  packages <- c("Seurat", "ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. Fixed outputs are overwritten on rerun.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/pooled_mclust_MFF_SenMayo_IFNalpha")
  input.object <- "hlma.scores"
  label.object <- "hlma.mff"  # Used only if hlma.scores lacks the four original mclust columns.
  class.columns <- c("MFF_mclust_all", "MFF_mclust_all_cutoff", "MFF_mclust_positive", "MFF_mclust_positive_cutoff")
  conditions <- c("all_high", "all_low", "positive_high", "positive_low")
  condition.labels <- c(all_high = "MFF all: High", all_low = "MFF all: Low",
    positive_high = "MFF positive: High", positive_low = "MFF positive: Low")
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  violin.columns <- c(MFF = "mff", score.columns)
  violin.labels <- c(MFF = "MFF expression", score.labels)
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  views <- c("All", "Young", "Old")
  plot.units <- "cell"  # c("cell", "donor") also plots donor means; both units' CSVs always saved.
  annotation.size <- 8  # Violin p and scatter rho / p / n.
  show.regression.line <- TRUE
  group.order <- c("Young", "Old"); group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L

  # 2. Pair stored MFF and all three existing scores by original profile ID.
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
  spearman.pair <- function(x, y, minimum = 3L) {
    if (length(x) != length(y) || any(!is.finite(x)) || any(!is.finite(y))) stop("Invalid correlation pairs.")
    n <- length(x); nx <- length(unique(x)); ny <- length(unique(y))
    reason <- if (!n) "no_observations" else if (n < minimum) "too_few_observations" else
      if (nx < 2L) "constant_MFF" else if (ny < 2L) "constant_score" else ""
    r <- data.frame(n_observations = n, n_unique_MFF = nx, n_unique_score = ny,
      rho = NA_real_, S = NA_real_, p_value = NA_real_, exact_used = NA, reason = reason)
    if (reason != "") return(r)
    exact <- n < 10L && !anyDuplicated(x) && !anyDuplicated(y)
    test <- cor.test(x, y, method = "spearman", alternative = "two.sided", exact = exact, continuity = FALSE)
    r$rho <- unname(test$estimate); r$S <- unname(test$statistic)
    r$p_value <- test$p.value; r$exact_used <- exact
    if (!is.finite(r$rho) || !is.finite(r$p_value)) r$reason <- "nonfinite_test"
    r
  }
  pooled.tables <- function(cells, columns, scope.names, view.names) {
    violin.columns <- c(MFF = "mff", columns)
    # Every profile contributes directly. Annotation never partitions or weights this analysis.
    donor <- data.frame(scope = character(), donor = character(), group = character(), Country = character(),
      n_profiles = integer(), n_scRNA = integer(), n_snRNA = integer(), mff = numeric(), stringsAsFactors = FALSE)
    for (nm in unname(columns)) donor[[nm]] <- numeric()
    for (scope in scope.names) {
      d <- cells[scope == "combined" | cells$tech == scope, , drop = FALSE]
      if (nrow(d)) for (idx in split(seq_len(nrow(d)), d$donor)) {
        z <- d[idx, , drop = FALSE]
        row <- data.frame(scope = scope, donor = z$donor[1], group = z$group[1], Country = z$Country[1],
          n_profiles = nrow(z), n_scRNA = sum(z$tech == "scRNA"), n_snRNA = sum(z$tech == "snRNA"), mff = mean(z$mff))
        for (nm in unname(columns)) row[[nm]] <- mean(z[[nm]])
        donor <- rbind(donor, row)
      }
    }
    violin <- correlation <- list()
    for (scope in scope.names) for (nm in names(violin.columns)) {
      d <- cells[scope == "combined" | cells$tech == scope, , drop = FALSE]
      ds <- donor[donor$scope == scope, , drop = FALSE]
      for (unit in c("cell", "donor")) {
        z <- if (unit == "cell") d else ds
        y <- z[[violin.columns[[nm]]]][z$group == "Young"]; o <- z[[violin.columns[[nm]]]][z$group == "Old"]
        reason <- if (!length(y) || !length(o)) "missing_group" else if (min(length(y), length(o)) < 2L)
          "too_few_observations" else if (length(unique(c(y, o))) < 2L) "constant_score" else ""
        p <- W <- NA_real_
        if (reason == "") {
          w <- wilcox.test(y, o, alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)
          p <- w$p.value; W <- unname(w$statistic)
        }
        my <- if (length(y)) mean(y) else NA_real_; mo <- if (length(o)) mean(o) else NA_real_
        violin[[length(violin) + 1L]] <- data.frame(scope = scope, score = nm, unit = unit,
          n_Young = length(y), n_Old = length(o), n_profiles_Young = sum(d$group == "Young"),
          n_profiles_Old = sum(d$group == "Old"), n_donors_Young = sum(ds$group == "Young"), n_donors_Old = sum(ds$group == "Old"),
          mean_Young = my, mean_Old = mo, difference_Old_minus_Young = mo - my,
          W = W, p_value = p, p_adj_BH = NA_real_, reason = reason)
        if (nm %in% names(columns)) for (view in view.names) {
          b <- if (view == "All") z else z[z$group == view, , drop = FALSE]
          dc <- if (view == "All") d else d[d$group == view, , drop = FALSE]
          dd <- if (view == "All") ds else ds[ds$group == view, , drop = FALSE]
          r <- spearman.pair(b$mff, b[[columns[[nm]]]])
          correlation[[length(correlation) + 1L]] <- cbind(data.frame(scope = scope, score = nm, view = view, unit = unit,
            n_profiles = nrow(dc), n_donors = nrow(dd), n_profiles_MFF_positive = sum(dc$mff > 0)), r)
        }
      }
    }
    v <- do.call(rbind, violin); c <- do.call(rbind, correlation); c$p_adj_BH <- NA_real_
    for (scope in scope.names) for (unit in c("cell", "donor")) {
      idx <- which(v$scope == scope & v$unit == unit); v$p_adj_BH[idx] <- p.adjust(v$p_value[idx], method = "BH")
      for (view in view.names) {
        idx <- which(c$scope == scope & c$unit == unit & c$view == view); c$p_adj_BH[idx] <- p.adjust(c$p_value[idx], method = "BH")
      }
    }
    rownames(v) <- rownames(c) <- rownames(donor) <- NULL
    list(violin = v, correlation = c, donor = donor)
  }
  mclust.pooled.tables <- function(cells, masks, columns, scope.names, view.names) {
    required <- c("all_high", "all_low", "positive_high", "positive_low")
    if (!identical(names(masks), required) || any(vapply(masks, function(v)
      !is.logical(v) || length(v) != nrow(cells) || anyNA(v), logical(1)))) stop("Invalid original subset masks.")
    violin <- correlation <- donor <- retention <- list()
    for (condition in required) {
      chosen <- cells[masks[[condition]], , drop = FALSE]
      r <- pooled.tables(chosen, columns, scope.names, view.names)
      tag <- function(d) cbind(data.frame(condition = rep(condition, nrow(d)), stringsAsFactors = FALSE), d)
      violin[[condition]] <- tag(r$violin); correlation[[condition]] <- tag(r$correlation)
      donor[[condition]] <- tag(r$donor)
      positive <- startsWith(condition, "positive")
      column <- if (positive) "MFF_mclust_positive" else "MFF_mclust_all"
      for (scope in scope.names) for (g in c("Young", "Old")) {
        index <- (scope == "combined" | cells$tech == scope) & cells$group == g
        z <- cells[index, , drop = FALSE]; selected <- cells[index & masks[[condition]], , drop = FALSE]
        eligible <- if (positive) z$mff > 0 else rep(TRUE, nrow(z))
        retention[[length(retention) + 1L]] <- data.frame(condition = condition, scope = scope, group = g,
          n_total = nrow(z), n_input_eligible = sum(eligible), n_selected = nrow(selected),
          n_unclassified = sum(eligible & is.na(z[[column]])), n_excluded_zero = if (positive) sum(z$mff == 0) else 0L,
          n_donors_total = length(unique(z$donor)), n_donors_selected = length(unique(selected$donor)),
          fraction_total = if (nrow(z)) nrow(selected) / nrow(z) else NA_real_,
          fraction_eligible = if (sum(eligible)) nrow(selected) / sum(eligible) else NA_real_)
      }
    }
    combine <- function(x) { d <- do.call(rbind, x); rownames(d) <- NULL; d }
    list(violin = combine(violin), correlation = combine(correlation), donor = combine(donor), retention = combine(retention))
  }
  write.checked.csv <- function(d, file) {
    token <- ".__HLMA_POOLED_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    con <- if (grepl("[.]gz$", file)) gzfile(file, "wt", encoding = "UTF-8") else file(file, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, con, row.names = FALSE, na = token), finally = close(con))
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else
      if (is.logical(v)) "logical" else "character", character(1))
    con <- if (grepl("[.]gz$", file)) gzfile(file, "rt", encoding = "UTF-8") else file(file, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(con, colClasses = classes, check.names = FALSE, na.strings = token,
      stringsAsFactors = FALSE), finally = close(con))
    stopifnot(identical(dim(restored), dim(d)), identical(names(restored), names(d)))
    for (nm in names(d)) if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(restored[[nm]], d[[nm]], tolerance = 1e-10))) else
      stopifnot(identical(restored[[nm]], d[[nm]]))
    invisible(file)
  }

  # 3. Reuse full-reference scores and original mclust labels, aligned by cell ID.
  if (!exists(input.object, envir = .GlobalEnv, inherits = FALSE))
    stop("Run the full script 34 in this R session first: missing ", input.object)
  if (!all(plot.units %in% c("cell", "donor")) || !length(plot.units) || anyDuplicated(plot.units)) stop("Invalid plot.units.")
  if (!is.numeric(annotation.size) || length(annotation.size) != 1L || !is.finite(annotation.size) || annotation.size <= 0)
    stop("annotation.size must be a finite positive number.")
  x <- get(input.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  ids <- colnames(x); rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == "MFF") == 1L)
  mff <- setNames(as.numeric(rna["MFF", ]), ids); rm(rna)
  meta <- x[[]]; labels.from <- input.object
  if (!all(class.columns %in% names(meta))) {
    if (!exists(label.object, envir = .GlobalEnv, inherits = FALSE))
      stop("Run script 33 in this session first: missing original labels and ", label.object)
    labelled <- get(label.object, envir = .GlobalEnv)
    stopifnot(inherits(labelled, "Seurat"), "RNA" %in% Assays(labelled), "data" %in% Layers(labelled[["RNA"]]))
    lrna <- GetAssayData(labelled, assay = "RNA", layer = "data"); lid <- colnames(labelled)
    stopifnot(identical(colnames(lrna), lid), sum(rownames(lrna) == "MFF") == 1L)
    lmff <- setNames(as.numeric(lrna["MFF", ]), lid); rm(lrna)
    meta <- align.labels(meta, mff, labelled[[]], lmff, ids); labels.from <- label.object
  }
  prepared <- prepare.subsets(meta, mff, score.columns, ids)
  cells <- prepared$cells; masks <- prepared$masks
  info <- x@misc$HLMA_two_scores
  exclusion.confirmed <- !is.null(info) && "MFF" %in% info$pool_excludes &&
    all(names(score.columns) %in% names(info$matched_features)) &&
    !any(vapply(info$matched_features[names(score.columns)], function(v) "MFF" %in% v, logical(1)))
  if (!exclusion.confirmed) warning("MFF exclusion from score signatures/control pool is unconfirmed; verify scoring provenance.")
  result <- mclust.pooled.tables(cells, masks, score.columns, scopes, views)
  violin.statistics <- result$violin; correlation.statistics <- result$correlation; donor <- result$donor
  retention <- result$retention
  for (condition in conditions) stopifnot(sum(donor$n_profiles[donor$condition == condition]) == 2L * sum(masks[[condition]]))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(violin.statistics[violin.statistics$unit == "cell", ], file.path(out, "01.pooled_scores.cell_young_old.csv"))
  write.checked.csv(violin.statistics[violin.statistics$unit == "donor", ], file.path(out, "02.pooled_scores.donor_young_old.csv"))
  write.checked.csv(correlation.statistics[correlation.statistics$unit == "cell", ], file.path(out, "03.pooled_MFF_score.cell_correlations.csv"))
  write.checked.csv(correlation.statistics[correlation.statistics$unit == "donor", ], file.path(out, "04.pooled_MFF_score.donor_correlations.csv"))
  write.checked.csv(donor, file.path(out, "05.pooled_scores.donor_means.csv"))
  write.checked.csv(cells, file.path(out, "06.pooled_scores.profile_pairs.csv.gz"))

  write.checked.csv(retention, file.path(out, "07.pooled_mclust.subset_retention.csv"))

  # 4. Violin: MFF first, then three scores. Correlation: three score panels.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  score.violin <- function(d, stat, unit, title, label, is.expression = FALSE) {
    d$group <- factor(d$group, levels = group.order)
    regular <- d[d$group %in% group.order[vapply(group.order, function(g) {
      v <- d$score[d$group == g]; length(v) >= 2L && length(unique(v)) > 1L
    }, logical(1))], , drop = FALSE]
    special <- d[!d$group %in% regular$group, , drop = FALSE]
    bounds <- if (nrow(d)) range(d$score) else c(-.1, .1); span <- max(diff(bounds), .05)
    bracket <- bounds[2] + .12 * span
    p <- ggplot(d, aes(group, score)) + geom_blank() +
      geom_violin(data = regular, aes(fill = group), trim = TRUE, scale = "width", width = .8)
    for (g in group.order) {
      z <- special[special$group == g, , drop = FALSE]
      if (nrow(z) == 1L) p <- p + geom_point(data = z, size = 2)
      if (nrow(z) > 1L) p <- p + geom_errorbar(data = z[1, , drop = FALSE], aes(ymin = score, ymax = score), width = .45)
    }
    axis.labels <- setNames(paste0(group.order, "\n(n = ", c(stat$n_Young, stat$n_Old), ")"), group.order)
    p + scale_fill_manual(values = group.colors, guide = "none") +
      scale_x_discrete(limits = group.order, labels = axis.labels, drop = FALSE) +
      annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .4) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket - .04 * span, yend = bracket, linewidth = .4) +
      annotate("text", x = 1.5, y = bounds[2] + .45 * span, label = format.p(stat$p_value), size = annotation.size) +
      coord_cartesian(ylim = c(bounds[1] - .08 * span, bounds[2] + .65 * span)) +
      labs(title = title, x = NULL, y = if (is.expression) {
        if (unit == "cell") "MFF expression (stored RNA/data)" else "Donor-mean MFF (stored RNA/data)"
      } else if (unit == "cell") paste0(label, " score") else paste0("Donor-mean ", label, " score")) +
      theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 13), axis.text = element_text(size = 18), panel.grid.minor = element_blank())
  }
  scatter.plot <- function(full, stat, unit, view, title, label) {
    d <- if (view == "All") full else full[full$group == view, , drop = FALSE]
    d$group <- factor(d$group, levels = group.order)
    bounds <- if (nrow(full)) range(full$score) else c(-.1, .1); span <- max(diff(bounds), .05)
    xb <- if (nrow(full)) range(full$mff) else c(0, 1)
    if (diff(xb) == 0) xb <- xb + c(-.05, .05) * max(1, abs(xb[1]))
    p <- ggplot(d, aes(mff, score)) + geom_point(aes(colour = group),
      size = if (unit == "cell") .55 else 3.5, alpha = if (unit == "cell") .70 else .95)
    if (show.regression.line && nrow(d) >= 3L && length(unique(d$mff)) > 1L && length(unique(d$score)) > 1L)
      p <- p + geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
        colour = "#222222", linewidth = 1, show.legend = FALSE)
    text <- paste0("rho = ", if (is.finite(stat$rho)) sprintf("%.2f", stat$rho) else "NA", "\n",
      format.p(stat$p_value), "\nn = ", stat$n_observations)
    p + scale_colour_manual(values = group.colors, limits = group.order,
      breaks = if (view == "All") group.order else view, drop = FALSE, name = NULL) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
      annotate("text", x = Inf, y = bounds[2] + .95 * span, label = text, hjust = 1.05, vjust = 1, size = annotation.size) +
      coord_cartesian(xlim = xb, ylim = c(bounds[1] - .08 * span, bounds[2] + span)) +
      labs(title = title, x = if (unit == "cell") "MFF expression (stored RNA/data)" else "Donor-mean MFF (stored RNA/data)",
        y = if (unit == "cell") paste0(label, " score") else paste0("Donor-mean ", label, " score")) +
      theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 13), axis.text = element_text(size = 18),
        panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 13))
  }
  save.png <- function(p, file, width, height, res = plot.res) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(file, width = width, height = height, res = res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  pngs <- character()
  for (condition in conditions) for (unit in plot.units) for (scope in scopes) {
    b <- if (unit == "cell") cells[masks[[condition]] & (scope == "combined" | cells$tech == scope), , drop = FALSE] else
      donor[donor$condition == condition & donor$scope == scope, , drop = FALSE]
    folder <- file.path(out, condition, scope, unit)
    caption <- if (unit == "cell")
      "All cell types pooled within the original MFF subset; conditional cell p after MFF selection; composition, donor, Country and tech not adjusted." else
      "All cell types pooled within each subset/donor; conditional donor p after MFF selection; composition, Country and tech not adjusted."
    plots <- list()
    for (nm in names(violin.columns)) {
      d <- b; d$score <- d[[violin.columns[[nm]]]]
      stat <- violin.statistics[violin.statistics$condition == condition & violin.statistics$scope == scope & violin.statistics$unit == unit & violin.statistics$score == nm, ]
      stopifnot(nrow(stat) == 1L)
      plots[[nm]] <- score.violin(d, stat, unit, violin.labels[nm], violin.labels[nm], is.expression = nm == "MFF")
    }
    overview <- wrap_plots(plots, ncol = length(violin.columns), nrow = 1L) +
      plot_annotation(title = paste(scope.labels[scope], condition.labels[condition], "All cell types | Young vs Old", sep = " | "), caption = caption,
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
    file <- file.path(folder, "Pooled.SenMayo_IFNalpha.Young_Old.violin.png")
    save.png(overview, file, panel.width * length(violin.columns), panel.height); pngs <- c(pngs, file)
    for (view in views) {
      plots <- list()
      for (nm in names(score.columns)) {
        full <- b; full$score <- full[[score.columns[[nm]]]]
        stat <- correlation.statistics[correlation.statistics$condition == condition & correlation.statistics$scope == scope & correlation.statistics$unit == unit &
          correlation.statistics$score == nm & correlation.statistics$view == view, ]
        stopifnot(nrow(stat) == 1L)
        plots[[nm]] <- scatter.plot(full, stat, unit, view, score.labels[nm], score.labels[nm])
      }
      overview <- (wrap_plots(plots, ncol = length(score.columns), nrow = 1L, guides = "collect") +
        plot_annotation(title = paste(scope.labels[scope], condition.labels[condition], "All cell types | MFF correlation", view, sep = " | "),
          caption = paste(caption, "Spearman rho/p; line: descriptive linear fit."),
          theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))) &
        theme(legend.position = "bottom")
      file <- file.path(folder, paste0("Pooled.SenMayo_IFNalpha.MFF_correlation.", view, ".png"))
      save.png(overview, file, panel.width * length(score.columns), panel.height); pngs <- c(pngs, file)
    }
  }
  stopifnot(length(pngs) == length(conditions) * length(plot.units) * length(scopes) * (1L + length(views)),
    all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  assign("hlma.pooled.mclust.scores", list(profiles = cells, masks = masks, donor_means = donor,
    violin_statistics = violin.statistics, correlation_statistics = correlation.statistics, retention = retention, label_source = labels.from), envir = .GlobalEnv)
  writeLines(c("All cell types pooled WITHIN each original MFF mclust subset.",
    "Run after scripts 33 and 34 in the same R session. Scripts 36/37/38 are not required.",
    "Input hlma.scores; missing classification columns are transferred from hlma.mff by exact cell IDs after metadata/MFF checks.",
    "Original combined, cell-type-specific labels/cutoffs are reused for combined, scRNA and snRNA. No global mclust refit.",
    "Subsets: all_high, all_low, positive_high, positive_low. They are overlapping selections, not four disjoint clusters.",
    "Classified MFF=0 cells remain in all_low. Positive subsets exclude zero. Unclassified profiles remain NA in the reference CSV.",
    "Violin panels, left to right: MFF expression (stored RNA/data), SaulSenMayo, Hallmark IFN-alpha, GO IFN-beta; 3600 x 850 px at 160dpi.",
    "Young/Old colours and raw two-sided unpaired Wilcoxon p/n match script 38. No boxplots; annotation.size=8.",
    "Three-score MFF Spearman scatter panels are also saved for All, Young and Old; 2700 x 850 px; optional descriptive LM line.",
    "Default cell output: 4 subsets x 3 scopes x (1 violin + 3 correlation) = 48 PNGs; cell+donor gives 96.",
    "Each subset/scope pools original profiles directly across all cell types; no equal weighting of cell-type means.",
    "Donor means average MFF and scores over identical selected profiles, with equal profile weights within donor and equal donor weights for tests.",
    "Wilcoxon needs at least 2 observations per age group; missing/rare/constant outcomes retain NA p and reason.",
    "Spearman needs at least 3 observations and two distinct values per axis; n<10 without ties uses exact permutation; other inputs are asymptotic.",
    "Violin BH: MFF + three scores within each condition/scope/unit; correlation BH: three scores within each condition/scope/unit/view.",
    "Scatter All/Young/Old share axis ranges within subset/scope/score/unit. Figures show raw p.",
    "MFF-based selection restricts MFF range; conditional age differences/correlations can reflect cell-type composition and donor/Country/tech effects.",
    "Seven CSV/gzip tables are readback-verified; retention includes selected, unclassified, eligible and excluded-zero profiles and donor coverage.",
    "No AddModuleScore, normalization, mclust fitting, RDS read/write or source Seurat modification.",
    paste0("MFF exclusion in scoring provenance confirmed: ", exclusion.confirmed),
    paste0("Input global object: ", input.object), paste0("Label source: ", labels.from), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_POOLED_MCLUST_MFF_THREE_SCORES_COMPLETE\nProfiles: ", nrow(cells), "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
