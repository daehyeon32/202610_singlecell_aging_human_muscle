# SERVER R: run after GO-enabled 34 in the SAME R session, with original script-32 CSVs available.
# Old only: filter cells FIRST, then select lower/upper MFF quartiles and compare three gene scores.
# Conditions are separate: MFF >= original zero-including mclust cutoff, or MFF > 0.
# No UMAP, score recalculation, mclust refitting, RDS I/O or dependency on scripts 43/44/45.
local({
  packages <- c("Seurat", "ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. Fixed output files are overwritten on rerun.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/Old_MFF_Top25_gene_scores")
  mclust.dir <- file.path(project, "analysis/3.Mffexpression/mclust_cutoffs")
  input.object <- "hlma.scores"
  conditions <- c("mclust_all_ge", "positive")
  condition.labels <- c(mclust_all_ge = "MFF >= original all-mode mclust cutoff", positive = "MFF > 0")
  groupings <- c("pooled", "celltype"); scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  group.order <- c("Bottom25", "Top25")
  group.colors <- c(Bottom25 = "#2878BD", Top25 = "#C54045")
  annotation.size <- 8; panel.width <- 900L; panel.height <- 850L; plot.res <- 160L; n.columns <- 5L

  # 2. Pure helpers. Scores stay signed; original cutoff sources are checked by exact profile ID.
  missing.value <- function(v) is.na(v) | v == "NA" | v == "" | grepl("^\\.__MFF_MISSING__\\._*$", v)
  make.current <- function(meta, mff, columns, ids = rownames(meta)) {
    required <- c("Annotation", "sample", "tech", "age_pop", "Country", unname(columns))
    if (!all(required %in% names(meta)) || !nrow(meta) || !identical(rownames(meta), ids) ||
        anyNA(ids) || anyDuplicated(ids) || any(ids == "")) stop("Missing metadata/score columns or invalid IDs.")
    if (!is.numeric(mff) || !identical(names(mff), ids) || length(mff) != length(ids) ||
        any(!is.finite(mff)) || any(mff < 0)) stop("Invalid RNA/data MFF or profile order.")
    for (nm in c("Annotation", "sample", "tech", "age_pop", "Country"))
      if (anyNA(meta[[nm]]) || any(!nzchar(trimws(as.character(meta[[nm]]))))) stop("Missing metadata: ", nm)
    if (!all(meta$tech %in% c("scRNA", "snRNA")) || !all(meta$age_pop %in% c("young_pop", "old_pop")))
      stop("Unknown technology or age labels.")
    d <- data.frame(cell_id = ids, donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
      tech = as.character(meta$tech), Country = as.character(meta$Country),
      group = ifelse(meta$age_pop == "old_pop", "Old", "Young"), mff = unname(mff), stringsAsFactors = FALSE)
    for (nm in unname(columns)) {
      if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]]))) stop("Invalid score: ", nm)
      d[[nm]] <- meta[[nm]]
    }
    d
  }
  number <- function(v, field) {
    missing <- missing.value(v)
    z <- suppressWarnings(as.numeric(as.character(v))); z[missing] <- NA_real_
    if (any(!missing & !is.finite(z))) stop("Invalid numeric field: ", field)
    z
  }
  read.saved <- function(path) {
    if (!file.exists(path)) stop("Missing original mclust result: ", path,
      "\nRun script 32 with save.classifications=TRUE first.")
    con <- if (grepl("[.]gz$", path)) gzfile(path, "rt", encoding = "UTF-8") else file(path, "rt", encoding = "UTF-8")
    on.exit(close(con))
    read.csv(con, colClasses = "character", na.strings = character(), check.names = FALSE, stringsAsFactors = FALSE)
  }
  join.all.classification <- function(current, profiles, fits, scope.names) {
    p.fields <- c("cell_id", "donor", "celltype", "tech", "group", "expression", "cutoff_all", "region_all")
    f.fields <- c("scope", "celltype", "input_mode", "cutoff", "status", "selected_G", "selected_return_code", "n_total", "n_input")
    if (!all(f.fields %in% names(fits)) || !all(scope.names %in% names(profiles))) stop("Missing original mclust columns/scopes.")
    if (anyNA(current$cell_id) || anyDuplicated(current$cell_id) || any(current$cell_id == "") ||
        !is.numeric(current$mff) || any(!is.finite(current$mff)) || any(current$mff < 0)) stop("Invalid current profiles/MFF.")
    types <- sort(unique(current$celltype)); rows <- list()
    for (scope in scope.names) {
      d <- current[scope == "combined" | current$tech == scope, , drop = FALSE]
      p <- profiles[[scope]]
      if (!all(p.fields %in% names(p)) || anyNA(p$cell_id) || anyDuplicated(p$cell_id) ||
          nrow(p) != nrow(d) || !setequal(p$cell_id, d$cell_id)) stop("Classification IDs do not match complete scope: ", scope)
      p <- p[match(d$cell_id, p$cell_id), , drop = FALSE]
      for (nm in c("donor", "celltype", "tech", "group"))
        if (!identical(as.character(p[[nm]]), as.character(d[[nm]]))) stop("Classification metadata mismatch: ", scope, " / ", nm)
      p$expression <- number(p$expression, "expression"); p$cutoff_all <- number(p$cutoff_all, "cutoff_all")
      p$region_all <- as.character(p$region_all); p$region_all[missing.value(p$region_all)] <- NA_character_
      if (any(!is.finite(p$expression)) || any(p$expression < 0) ||
          any(abs(p$expression - d$mff) > 1e-10 * pmax(1, abs(d$mff))) ||
          any(!is.na(p$region_all) & !p$region_all %in% c("High", "Low"))) stop("Original MFF/classification differs: ", scope)
      f <- fits[which(fits$scope == scope & fits$input_mode == "all"), , drop = FALSE]
      if (anyNA(f$celltype) || anyDuplicated(f$celltype) || nrow(f) != length(types) ||
          !setequal(f$celltype, types) || any(missing.value(f$status))) stop("Incomplete all-mode fits: ", scope)
      for (nm in c("cutoff", "selected_G", "selected_return_code", "n_total", "n_input")) f[[nm]] <- number(f[[nm]], nm)
      for (ct in types) {
        use <- which(p$celltype == ct); st <- f[f$celltype == ct, , drop = FALSE]; cutoff <- st$cutoff
        if (!is.finite(st$n_total) || !is.finite(st$n_input) || st$n_total != length(use) || st$n_input != length(use))
          stop("ALL fit counts mismatch: ", scope, " / ", ct)
        if (is.finite(cutoff)) {
          if (cutoff < 0 || st$status != "ok" || is.na(st$selected_G) || st$selected_G != 2 ||
              is.na(st$selected_return_code) || st$selected_return_code != 0) stop("Invalid converged ALL cutoff.")
          z <- p[use, , drop = FALSE]; tol <- 1e-10 * pmax(1, abs(cutoff))
          if (anyNA(z$region_all) || anyNA(z$cutoff_all) || any(abs(z$cutoff_all - cutoff) > tol) ||
              any(z$region_all == "High" & z$expression - cutoff < -tol) ||
              any(z$region_all == "Low" & z$expression - cutoff > tol) ||
              any(z$expression == 0 & z$region_all != "Low")) stop("Saved ALL labels/cutoff disagree: ", scope, " / ", ct)
        } else if (st$status == "ok" || any(!is.na(p$region_all[use])) || any(!is.na(p$cutoff_all[use]))) {
          stop("Unavailable cutoff must stay unclassified: ", scope, " / ", ct)
        }
      }
      status <- as.character(f$status[match(p$celltype, f$celltype)])
      rows[[length(rows) + 1L]] <- data.frame(scope = rep(scope, nrow(d)), d[, c("cell_id", "donor", "celltype", "tech", "group", "mff")],
        cutoff_all = p$cutoff_all, region_all = p$region_all, fit_status = status, row.names = NULL)
    }
    result <- do.call(rbind, rows)
    stopifnot(!anyDuplicated(result[, c("scope", "cell_id")]))
    result
  }
  write.checked.csv <- function(d, path) {
    token <- ".__MFF_QUARTILE_SCORE_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    con <- if (grepl("[.]gz$", path)) gzfile(path, "wt", encoding = "UTF-8") else file(path, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, con, row.names = FALSE, na = token), finally = close(con))
    classes <- vapply(d, function(v) if (is.logical(v)) "logical" else if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else "character", character(1))
    con <- if (grepl("[.]gz$", path)) gzfile(path, "rt", encoding = "UTF-8") else file(path, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(con, colClasses = classes, na.strings = token, check.names = FALSE), finally = close(con))
    stopifnot(identical(dim(restored), dim(d)), identical(names(restored), names(d)))
    for (nm in names(d)) if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(restored[[nm]], d[[nm]], tolerance = 1e-10))) else
      stopifnot(identical(restored[[nm]], d[[nm]]))
    invisible(path)
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

  score.test <- function(values, groups) {
    if (!is.numeric(values) || length(values) != length(groups) || any(!is.finite(values)) ||
        any(!is.na(groups) & !groups %in% c("Bottom25", "Top25"))) stop("Invalid gene-score comparison.")
    b <- values[which(groups == "Bottom25")]; t <- values[which(groups == "Top25")]
    reason <- if (!length(b) || !length(t)) "missing_MFF_quartile_group" else
      if (min(length(b), length(t)) < 2L) "too_few_observations" else
        if (length(unique(c(b, t))) < 2L) "constant_gene_score" else ""
    W <- p <- NA_real_
    if (reason == "") {
      w <- wilcox.test(b, t, alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE)
      W <- unname(w$statistic); p <- w$p.value
    }
    data.frame(n_Bottom = length(b), n_Top = length(t),
      mean_Bottom = if (length(b)) mean(b) else NA_real_, mean_Top = if (length(t)) mean(t) else NA_real_,
      median_Bottom = if (length(b)) median(b) else NA_real_, median_Top = if (length(t)) median(t) else NA_real_,
      difference_Top_minus_Bottom = if (length(b) && length(t)) mean(t) - mean(b) else NA_real_,
      W = W, p_value = p, p_adj_BH = NA_real_, test_reason = reason)
  }
  reverse.tables <- function(current, classified, columns, scope.names, grouping.names, types,
      condition.names = c("mclust_all_ge", "positive")) {
    if (!length(scope.names) || any(!scope.names %in% c("combined", "scRNA", "snRNA")) || anyDuplicated(scope.names) ||
        !length(grouping.names) || any(!grouping.names %in% c("pooled", "celltype")) || anyDuplicated(grouping.names) ||
        !length(condition.names) || any(!condition.names %in% c("mclust_all_ge", "positive")) || anyDuplicated(condition.names))
      stop("Invalid scopes/groupings/conditions.")
    fields <- c("cell_id", "donor", "celltype", "tech", "Country", "group", "mff", unname(columns))
    if (!all(fields %in% names(current)) || anyNA(current$cell_id) || anyDuplicated(current$cell_id) ||
        any(current$cell_id == "") || !is.numeric(current$mff) || any(!is.finite(current$mff)) || any(current$mff < 0) ||
        !all(c("scope", "cell_id", "donor", "celltype", "tech", "group", "mff", "cutoff_all", "region_all", "fit_status") %in% names(classified)) ||
        anyNA(classified[, c("scope", "cell_id")]) || anyDuplicated(classified[, c("scope", "cell_id")]))
      stop("Invalid current/source profile keys.")
    for (nm in unname(columns)) if (!is.numeric(current[[nm]]) || any(!is.finite(current[[nm]]))) stop("Invalid stored score: ", nm)
    if (!setequal(types, unique(current$celltype)) || anyDuplicated(types)) stop("Cell-type grid differs from current object.")
    complete <- list()
    for (scope in scope.names) {
      d <- current[scope == "combined" | current$tech == scope, , drop = FALSE]
      c <- classified[classified$scope == scope, , drop = FALSE]
      if (nrow(c) != nrow(d) || !setequal(c$cell_id, d$cell_id)) stop("Source classification does not cover complete scope: ", scope)
      c <- c[match(d$cell_id, c$cell_id), , drop = FALSE]
      for (nm in c("donor", "celltype", "tech", "group"))
        if (!identical(as.character(d[[nm]]), as.character(c[[nm]]))) stop("Source metadata mismatch: ", scope, " / ", nm)
      if (!is.numeric(c$mff) || !isTRUE(all.equal(c$mff, d$mff, tolerance = 1e-10, check.attributes = FALSE)) ||
          !is.numeric(c$cutoff_all) || any(!is.na(c$cutoff_all) & (!is.finite(c$cutoff_all) | c$cutoff_all < 0)))
        stop("Source MFF/cutoff differs: ", scope)
      d$cutoff_all <- c$cutoff_all; d$region_all <- c$region_all; d$fit_status <- c$fit_status
      complete[[scope]] <- d
    }
    empty <- data.frame(condition = character(), scope = character(), grouping = character(),
      selection_celltype = character(), quartile = character(), current[FALSE, fields, drop = FALSE],
      cutoff_all = numeric(), region_all = character(), fit_status = character(), row.names = NULL)
    filters <- stats <- selected <- list()
    for (condition in condition.names) for (scope in scope.names) for (grouping in grouping.names) {
      old <- complete[[scope]][complete[[scope]]$group == "Old", , drop = FALSE]
      targets <- if (grouping == "pooled") "ALL" else types
      for (ct in targets) {
        d <- if (grouping == "pooled") old else old[old$celltype == ct, , drop = FALSE]
        if (condition == "mclust_all_ge") {
          missing <- is.na(d$cutoff_all)
          keep <- !missing & d$mff >= d$cutoff_all
          removed.cutoff <- sum(missing); removed.below <- sum(!missing & !keep); removed.zero <- 0L
        } else {
          keep <- d$mff > 0
          removed.cutoff <- removed.below <- 0L; removed.zero <- sum(!keep)
        }
        eligible <- d[which(keep), , drop = FALSE]; q <- quartile.groups(eligible$mff)
        reason <- q$selection_reason
        if (reason == "no_Old_profiles") reason <- "no_eligible_Old_profiles"
        if (reason == "fewer_than_four_Old_profiles") reason <- "fewer_than_four_eligible_Old_profiles"
        use <- which(!is.na(q$quartile)); z <- eligible[use, , drop = FALSE]; g <- q$quartile[use]
        n <- nrow(eligible); nb <- sum(g == "Bottom25"); nt <- sum(g == "Top25")
        keys <- data.frame(condition = condition, scope = scope, grouping = grouping, selection_celltype = ct)
        f <- cbind(keys, data.frame(n_Old = nrow(d), n_eligible = n,
          n_removed_below_cutoff = removed.below, n_removed_cutoff_unavailable = removed.cutoff,
          n_removed_zero = removed.zero, n_zero_eligible = sum(eligible$mff == 0),
          q25 = q$q25, q75 = q$q75, n_ties_q25 = q$n_ties_q25, n_ties_q75 = q$n_ties_q75,
          n_Bottom = nb, n_Top = nt, n_Middle = sum(is.na(q$quartile)),
          fraction_Bottom_of_eligible = if (n) nb/n else NA_real_, fraction_Top_of_eligible = if (n) nt/n else NA_real_,
          n_zero_Bottom = sum(z$mff[g == "Bottom25"] == 0), n_zero_Top = sum(z$mff[g == "Top25"] == 0),
          n_donors_Bottom = length(unique(z$donor[g == "Bottom25"])), n_donors_Top = length(unique(z$donor[g == "Top25"])),
          selection_reason = reason))
        filters[[length(filters) + 1L]] <- f
        if (length(use)) selected[[length(selected) + 1L]] <- data.frame(keys, quartile = g, z, row.names = NULL)
        for (nm in names(columns)) {
          r <- score.test(z[[columns[[nm]]]], g)
          stopifnot(r$n_Bottom == nb, r$n_Top == nt)
          stats[[length(stats) + 1L]] <- cbind(f[, c(names(keys), "n_Old", "n_eligible", "q25", "q75",
            "selection_reason", "n_donors_Bottom", "n_donors_Top"), drop = FALSE],
            data.frame(score = nm, score_column = columns[[nm]]), r)
        }
      }
    }
    filters <- do.call(rbind, filters); stats <- do.call(rbind, stats)
    chosen <- if (length(selected)) do.call(rbind, selected) else empty
    for (condition in condition.names) for (scope in scope.names) for (grouping in grouping.names) {
      use <- which(stats$condition == condition & stats$scope == scope & stats$grouping == grouping)
      stats$p_adj_BH[use] <- p.adjust(stats$p_value[use], method = "BH", n = length(use))
    }
    rownames(filters) <- rownames(stats) <- rownames(chosen) <- NULL
    stopifnot(all(filters$n_eligible + filters$n_removed_below_cutoff + filters$n_removed_cutoff_unavailable +
      filters$n_removed_zero == filters$n_Old), all(filters$n_Bottom + filters$n_Top + filters$n_Middle == filters$n_eligible),
      !anyDuplicated(chosen[, c("condition", "scope", "grouping", "selection_celltype", "cell_id")]))
    list(filters = filters, statistics = stats, selected = chosen)
  }
  # 3. Reuse GO-enabled script 34's stored scores and RNA/data MFF.
  if (!exists(input.object, envir = .GlobalEnv, inherits = FALSE)) stop("Run GO-enabled script 34 first: missing ", input.object)
  stopifnot(identical(names(score.columns), names(score.labels)), n.columns >= 1L)
  x <- get(input.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  ids <- colnames(x); rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == "MFF") == 1L)
  mff <- setNames(as.numeric(rna["MFF", ]), ids); rm(rna)
  current <- make.current(x[[]], mff, score.columns, ids); types <- sort(unique(current$celltype))
  fits <- read.saved(file.path(mclust.dir, "02.MFF.mclust_fit_details.csv"))
  profiles <- setNames(lapply(scopes, function(scope) read.saved(file.path(mclust.dir,
    paste0("05.MFF.profile_classification.", scope, ".csv.gz")))), scopes)
  classified <- join.all.classification(current, profiles, fits, scopes)
  result <- reverse.tables(current, classified, score.columns, scopes, groupings, types, conditions)
  filters <- result$filters; statistics <- result$statistics; selected <- result$selected
  info <- x@misc$HLMA_two_scores
  exclusion.confirmed <- !is.null(info) && "MFF" %in% info$pool_excludes &&
    all(names(score.columns) %in% names(info$matched_features)) &&
    !any(vapply(info$matched_features[names(score.columns)], function(v) "MFF" %in% v, logical(1)))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(filters, file.path(out, "01.MFF_filter_quartile_counts.csv"))
  write.checked.csv(statistics, file.path(out, "02.gene_score_statistics.csv"))
  write.checked.csv(selected, file.path(out, "03.MFF_quartile_selected_profiles.csv.gz"))
  write.checked.csv(current[current$group == "Old", ], file.path(out, "04.Old_all_profiles.reference.csv.gz"))
  write.checked.csv(classified, file.path(out, "05.original_all_classification.audit.csv.gz"))

  # 4. Two violins per panel: SAME MFF lower/upper cells for all three score outcomes.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  score.violin <- function(d, stat, title) {
    d$score_value <- d[[stat$score_column]]
    d$quartile <- factor(d$quartile, levels = group.order)
    regular <- d[d$quartile %in% group.order[vapply(group.order, function(g) {
      v <- d$score_value[d$quartile == g]; length(v) >= 2L && length(unique(v)) > 1L
    }, logical(1))], , drop = FALSE]
    special <- d[!d$quartile %in% regular$quartile, , drop = FALSE]
    bounds <- if (nrow(d)) range(d$score_value) else c(-.05, .05); span <- max(diff(bounds), .05)
    bracket <- bounds[2] + .12 * span
    p <- ggplot(d, aes(quartile, score_value)) + geom_blank() +
      geom_violin(data = regular, aes(fill = quartile), trim = TRUE, scale = "width", width = .8) +
      geom_boxplot(data = regular, width = .15, fill = "white", colour = "black",
        outlier.shape = NA, linewidth = .4)
    for (g in group.order) {
      z <- special[special$quartile == g, , drop = FALSE]
      if (nrow(z) == 1L) p <- p + geom_point(data = z, size = 2)
      if (nrow(z) > 1L) p <- p + geom_errorbar(data = z[1, , drop = FALSE], aes(ymin = score_value, ymax = score_value), width = .45)
    }
    labels <- setNames(paste0(c("MFF lower 25%", "MFF upper 25%"), "\n(n = ", c(stat$n_Bottom, stat$n_Top), ")"), group.order)
    p + scale_fill_manual(values = group.colors, guide = "none") +
      scale_x_discrete(limits = group.order, labels = labels, drop = FALSE) +
      annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .4) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket - .04 * span, yend = bracket, linewidth = .4) +
      annotate("text", x = 1.5, y = bounds[2] + .45 * span, label = format.p(stat$p_value), size = annotation.size) +
      coord_cartesian(ylim = c(bounds[1] - .08 * span, bounds[2] + .65 * span)) +
      labs(title = title, subtitle = paste0("Old = ", stat$n_Old, "; after filter = ", stat$n_eligible,
        if (stat$selection_reason != "") " | MFF quartiles unavailable" else ""),
        x = NULL, y = paste0(score.labels[stat$score], " score")) +
      theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 13), axis.text = element_text(size = 18), panel.grid.minor = element_blank())
  }
  panel <- function(condition, scope, grouping, ct, score, title) {
    s <- statistics[statistics$condition == condition & statistics$scope == scope & statistics$grouping == grouping &
      statistics$selection_celltype == ct & statistics$score == score, , drop = FALSE]
    d <- selected[selected$condition == condition & selected$scope == scope & selected$grouping == grouping &
      selected$selection_celltype == ct, , drop = FALSE]
    stopifnot(nrow(s) == 1L); score.violin(d, s, title)
  }
  save.png <- function(p, path, width, height) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(path, width = width, height = height, res = plot.res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  caption <- paste("Old only. Apply the condition first, then MFF lower <= Q25 / upper >= Q75 within eligible cells.",
    "Boundary ties stay together; actual n may exceed 25%. Same MFF-selected cells across three signed gene scores.",
    "Raw unpaired cell Wilcoxon p; donor/Country/composition not adjusted. BH results in CSV.")
  manifest <- list(); pngs <- character()
  for (condition in conditions) for (scope in scopes) {
    if ("pooled" %in% groupings) {
      plots <- lapply(names(score.columns), function(nm) panel(condition, scope, "pooled", "ALL", nm, score.labels[nm]))
      heading <- paste(scope.labels[scope], "Old | All cell types pooled", condition.labels[condition], sep = " | ")
      p <- wrap_plots(plots, ncol = length(score.columns)) + plot_annotation(
        title = paste(strwrap(heading, width = 95), collapse = "\n"),
        caption = paste(strwrap(caption, width = 175), collapse = "\n"),
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, condition, "pooled", scope, "Old.three_scores.by_MFF_quartiles.violin.png")
      w <- panel.width * length(score.columns); h <- panel.height
      save.png(p, path, w, h); pngs <- c(pngs, path)
      manifest[[length(manifest) + 1L]] <- data.frame(condition = condition, scope = scope, grouping = "pooled",
        score = "ALL", width = w, height = h, path = path)
    }
    if ("celltype" %in% groupings) for (nm in names(score.columns)) {
      plots <- lapply(types, function(ct) panel(condition, scope, "celltype", ct, nm, ct))
      cols <- min(n.columns, length(types)); rows <- ceiling(length(types)/cols)
      heading <- paste(scope.labels[scope], "Old", score.labels[nm], "Within-cell-type MFF quartiles", condition.labels[condition], sep = " | ")
      p <- wrap_plots(plots, ncol = cols) + plot_annotation(
        title = paste(strwrap(heading, width = 32L * cols), collapse = "\n"),
        caption = paste(strwrap(caption, width = 45L * cols), collapse = "\n"),
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, condition, "celltype", scope, paste0("Old.", nm, ".by_MFF_quartiles.violin.overview.png"))
      w <- panel.width * cols; h <- panel.height * rows
      save.png(p, path, w, h); pngs <- c(pngs, path)
      manifest[[length(manifest) + 1L]] <- data.frame(condition = condition, scope = scope, grouping = "celltype_overview",
        score = nm, width = w, height = h, path = path)
    }
  }
  expected <- length(conditions) * length(scopes) * (as.integer("pooled" %in% groupings) +
    as.integer("celltype" %in% groupings) * length(score.columns))
  stopifnot(length(pngs) == expected, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  manifest <- do.call(rbind, manifest); write.checked.csv(manifest, file.path(out, "06.PNG_manifest.csv"))
  writeLines(c("Old-only inverse analysis: gene-score outcomes by MFF expression lower/upper quartiles.",
    "Input: hlma.scores from GO-enabled 34 in this R session, plus original 32 all-mode fit/classification CSVs.",
    "Scripts 43/44/45 do not need to be rerun or cached. No UMAP, score/normalization/mclust recalculation or RDS I/O.",
    "Apply two SEPARATE conditions before calculating any MFF quartile:",
    "1) mclust_all_ge: MFF >= numeric original all-mode cutoff, separately for combined/scRNA/snRNA and original Annotation.",
    "The original all-mode cutoffs include zero in the Young+Old fits; no Old-only fitting or positive-mode cutoff replacement.",
    "Inclusive >= was explicitly requested: equality is retained, including original Low labels at the boundary. This differs from original High (> cutoff).",
    "If cutoff=0, condition 1 can retain MFF=0; condition 2 independently excludes exact zeros. Unavailable original cutoffs exclude condition 1.",
    "2) positive: MFF > 0. No mclust label/cutoff requirement for selecting these cells; no positive-mode mclust threshold used.",
    "Complete original CSV IDs, metadata, RNA/data MFF, all-fit counts/convergence and labels are validated before selection.",
    "Each condition/scope/pooled-or-celltype unit uses its own eligible Old MFF Q25/Q75 from quantile(type=7).",
    "Lower includes MFF <= Q25; upper includes MFF >= Q75. Boundary ties stay together, so actual fractions can exceed 25%.",
    "Fewer than 4 eligible Old cells or Q25=Q75 gives unavailable, non-overlapping tails with reason recorded; never split ties randomly.",
    "The same lower/upper cell membership is reused for all three score outcomes. Middle MFF cells are omitted from violin comparison only.",
    "Module scores are stored AddModuleScores; negative and zero values stay in the data and visible y-axis range.",
    "Violin: blue MFF lower25 vs red MFF upper25 with narrow white boxplots. Singleton/constant score groups use a point/line; empty groups keep labeled axes.",
    "Boxplots hide outlier markers only; all values remain in the violin/statistics. Zero IQR collapses the box to a line.",
    "Each panel labels original Old n, eligible n, actual lower/upper n and raw two-sided unpaired Wilcoxon p.",
    "Wilcoxon exact=FALSE/correct=TRUE; missing/small/constant score comparison yields p=NA with reason.",
    "BH across all score/cell-type comparisons within condition x scope x grouping, explicit n includes unavailable tests in the family.",
    "Cell p does not adjust donor dependence, Country, technology or composition. Pooled types/combined technologies mix composition; analyses overlap.",
    "Counts CSV records filter removal reasons, eligible zeros, quartile limits/ties/fractions, middle cells, donors and unavailable selection reasons.",
    "Selected memberships are one row per cell/condition/scope/grouping, with all three original gene scores; never duplicated just to select per score.",
    "Default output: 2 conditions x (3 pooled three-score PNGs + 9 celltype overview PNGs) = 24 violin PNGs.",
    paste0("MFF exclusion from score signatures/control pool confirmed: ", exclusion.confirmed),
    paste0("PNGs: ", length(pngs)), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  assign("hlma.old.mff.quartile.scores", list(old_profiles = current[current$group == "Old", ],
    selected_profiles = selected, filters = filters, statistics = statistics, manifest = manifest,
    settings = list(conditions = conditions, groupings = groupings, scopes = scopes, score_columns = score.columns,
      age_group = "Old", mclust_mode = "all", cutoff_comparison = ">=", quantile_type = 7L, boundary_ties = "included")), envir = .GlobalEnv)
  cat("\nHLMA_OLD_MFF_QUARTILES_GENE_SCORES_VIOLIN_COMPLETE\nOld profiles: ", sum(current$group == "Old"),
    "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
