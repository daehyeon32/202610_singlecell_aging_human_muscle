# SERVER R: run after 43 and 44 have completed in the SAME R session.
# Default: intersect BOTH saved Old score quartiles with original zero-including MFF High.
# Reuse original score quartiles, all-mode mclust labels and saved UMAP; no refitting.
local({
  packages <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "scales", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. The comparison mode gets its own output folder.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  comparison.mode <- "both_high" # Filter BOTH score quartiles to original MFF all-High.
  # bottom_all: keep Bottom25 whole. top_only: draw only Top25 intersect all-High.
  out <- file.path(project, "analysis/4.genescore/Old_score_Top25_MFF_allHigh", comparison.mode)
  mclust.dir <- file.path(project, "analysis/3.Mffexpression/mclust_cutoffs")
  input.object <- "hlma.scores"
  cache.object <- "hlma.old.score.quartiles"; umap.cache.object <- "hlma.old.score.quartile.umap"
  assay <- "RNA"; gene <- "MFF"; reduction <- "umap"
  groupings <- c("pooled", "celltype"); scopes <- c("combined", "scRNA", "snRNA")
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  group.order <- if (comparison.mode == "top_only") "Top25" else c("Bottom25", "Top25")
  group.labels <- c(Bottom25 = if (comparison.mode == "both_high") "Bottom 25% AND MFF High" else "Bottom 25% (all)",
    Top25 = "Top 25% AND MFF High")
  group.colors <- c(Bottom25 = "#2878BD", Top25 = "#C54045")
  selection.label <- if (comparison.mode == "both_high") "Both score quartiles intersect MFF High" else "Top25 intersect MFF High"
  expression.colors <- c("#D9D9D9", "#FEE8C8", "#FDBB84", "#E34A33", "#990000")
  point.size <- .50; annotation.size <- 8
  font.base <- 16; font.title <- 20; font.axis <- 16; font.legend <- 13
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L
  n.columns <- 5L; n.pair.columns <- 5L
  save.individual.celltype.pairs <- FALSE

  # 2. Exact-ID validation of script 43 and original script-32 all-mode results.
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
  verify.cache <- function(cache, current, embedding, columns, scope.names, grouping.names, types) {
    fields <- c("cell_id", "donor", "celltype", "tech", "Country", "group", "mff", unname(columns))
    if (!all(fields %in% names(current)) || anyNA(current$cell_id) || anyDuplicated(current$cell_id) ||
        any(current$cell_id == "")) stop("Invalid current profile table.")
    if (!is.matrix(embedding) || !is.numeric(embedding) || ncol(embedding) < 2L ||
        is.null(rownames(embedding)) || anyNA(rownames(embedding)) || anyDuplicated(rownames(embedding)) ||
        !setequal(rownames(embedding), current$cell_id) || any(!is.finite(embedding[, 1:2, drop = FALSE])))
      stop("Saved UMAP must cover all original profiles with finite, unique coordinates.")
    if (is.null(cache$settings) || !identical(cache$settings$score_columns, columns) ||
        !identical(cache$settings$quantile_type, 7L) || !identical(cache$settings$boundary_ties, "included") ||
        !all(scope.names %in% cache$settings$scopes) || !all(grouping.names %in% cache$settings$groupings))
      stop("Script-43 settings mismatch; rerun script 43.")
    old <- cache$old_profiles; s <- cache$statistics; chosen <- cache$selected_profiles
    now <- current[current$group == "Old", , drop = FALSE]
    if (is.null(old) || !all(fields %in% names(old)) || anyNA(old$cell_id) || anyDuplicated(old$cell_id) ||
        nrow(old) != nrow(now) || !setequal(old$cell_id, now$cell_id)) stop("Stale/missing Old reference; rerun 43.")
    old <- old[match(now$cell_id, old$cell_id), , drop = FALSE]
    for (nm in fields) {
      ok <- if (nm %in% c("mff", unname(columns))) is.numeric(old[[nm]]) && all(is.finite(old[[nm]])) &&
        isTRUE(all.equal(old[[nm]], now[[nm]], tolerance = 1e-10, check.attributes = FALSE)) else
          identical(as.character(old[[nm]]), as.character(now[[nm]]))
      if (!ok) stop("Current object differs from script-43 Old reference: ", nm)
    }
    old$UMAP_1 <- embedding[old$cell_id, 1]; old$UMAP_2 <- embedding[old$cell_id, 2]
    stat.fields <- c("scope", "grouping", "selection_celltype", "score", "n_Old", "q25", "q75",
      "n_Bottom", "n_Top", "n_Middle", "selection_reason")
    chosen.fields <- c("scope", "grouping", "selection_celltype", "score", "quartile", "cell_id",
      "donor", "celltype", "tech", "Country", "group", "mff", "score_value")
    if (is.null(s) || is.null(chosen) || !all(stat.fields %in% names(s)) || !all(chosen.fields %in% names(chosen)))
      stop("Missing script-43 statistics/membership columns.")
    grid <- list()
    for (scope in scope.names) for (grouping in grouping.names) for (ct in if (grouping == "pooled") "ALL" else types)
      for (nm in names(columns)) grid[[length(grid) + 1L]] <- data.frame(scope = scope, grouping = grouping,
        selection_celltype = ct, score = nm)
    grid <- do.call(rbind, grid)
    keys <- function(d) paste(d$scope, d$grouping, d$selection_celltype, d$score, sep = "\034")
    s <- s[s$scope %in% scope.names & s$grouping %in% grouping.names, , drop = FALSE]
    chosen <- chosen[chosen$scope %in% scope.names & chosen$grouping %in% grouping.names, , drop = FALSE]
    if (anyNA(s[, stat.fields[1:4]]) || anyDuplicated(keys(s)) || nrow(s) != nrow(grid) ||
        !setequal(keys(s), keys(grid))) stop("Incomplete/duplicate statistic keys; rerun 43.")
    s <- s[match(keys(grid), keys(s)), , drop = FALSE]
    for (nm in c("n_Old", "n_Bottom", "n_Top", "n_Middle"))
      if (!is.numeric(s[[nm]]) || any(!is.finite(s[[nm]])) || any(s[[nm]] < 0 | s[[nm]] != floor(s[[nm]])))
        stop("Invalid membership counts: ", nm)
    if (anyNA(chosen[, chosen.fields[1:11]]) || any(!chosen$quartile %in% c("Bottom25", "Top25")) ||
        anyDuplicated(chosen[, c("scope", "grouping", "selection_celltype", "score", "cell_id")]) ||
        any(!keys(chosen) %in% keys(s))) stop("Invalid/duplicate selected memberships.")
    index <- match(chosen$cell_id, old$cell_id)
    if (anyNA(index)) stop("Selected cells are missing from the Old reference.")
    for (nm in c("donor", "celltype", "tech", "Country", "group"))
      if (!identical(as.character(chosen[[nm]]), as.character(old[[nm]][index]))) stop("Selected metadata differs: ", nm)
    if (!is.numeric(chosen$mff) || any(!is.finite(chosen$mff)) ||
        !isTRUE(all.equal(chosen$mff, old$mff[index], tolerance = 1e-10, check.attributes = FALSE)))
      stop("Selected MFF differs from original RNA/data.")
    for (nm in names(columns)) {
      use <- which(chosen$score == nm)
      if (!is.numeric(chosen$score_value) || any(!is.finite(chosen$score_value[use])) ||
          !isTRUE(all.equal(chosen$score_value[use], old[[columns[[nm]]]][index[use]],
            tolerance = 1e-10, check.attributes = FALSE))) stop("Selected score differs: ", nm)
    }
    chosen.index <- match(keys(chosen), keys(s))
    for (i in seq_len(nrow(s))) {
      st <- s[i, ]; d <- old[st$scope == "combined" | old$tech == st$scope, , drop = FALSE]
      if (st$grouping == "celltype") d <- d[d$celltype == st$selection_celltype, , drop = FALSE]
      if (st$n_Old != nrow(d) || st$n_Bottom + st$n_Top + st$n_Middle != st$n_Old || is.na(st$selection_reason))
        stop("Old/group counts mismatch at statistic row ", i)
      reason.ok <- if (st$selection_reason == "") nrow(d) >= 4L else
        if (st$selection_reason == "no_Old_profiles") nrow(d) == 0L else
          if (st$selection_reason == "fewer_than_four_Old_profiles") nrow(d) > 0L && nrow(d) < 4L else
            if (st$selection_reason == "q25_equals_q75") nrow(d) >= 4L &&
              is.finite(st$q25) && is.finite(st$q75) && st$q25 == st$q75 else FALSE
      if (!isTRUE(reason.ok)) stop("Saved unavailable-selection reason disagrees with Old profiles: row ", i)
      wanted <- list(Bottom25 = character(), Top25 = character())
      if (st$selection_reason == "") {
        if (!is.finite(st$q25) || !is.finite(st$q75) || st$q25 >= st$q75) stop("Invalid saved quartile limits.")
        v <- d[[columns[[st$score]]]]
        wanted$Bottom25 <- d$cell_id[v <= st$q25]; wanted$Top25 <- d$cell_id[v >= st$q75]
      } else if (!st$selection_reason %in% c("no_Old_profiles", "fewer_than_four_Old_profiles", "q25_equals_q75")) {
        stop("Unknown quartile selection reason.")
      }
      for (g in c("Bottom25", "Top25")) {
        actual <- chosen$cell_id[which(chosen.index == i & chosen$quartile == g)]
        count <- if (g == "Bottom25") st$n_Bottom else st$n_Top
        if (count != length(wanted[[g]]) || !setequal(actual, wanted[[g]]))
          stop("Membership differs from script-43 saved cutoff: row ", i, " / ", g)
      }
    }
    chosen$UMAP_1 <- embedding[chosen$cell_id, 1]; chosen$UMAP_2 <- embedding[chosen$cell_id, 2]
    rownames(old) <- rownames(s) <- rownames(chosen) <- NULL
    list(old = old, selected = chosen, statistics = s)
  }
  shared.limits <- function(old, modes, probability = .95) {
    if (!length(modes) || any(!modes %in% c("actual_max", "positive_quantile")) || anyDuplicated(modes) ||
        length(probability) != 1L || !is.finite(probability) || probability <= 0 || probability > 1 ||
        anyDuplicated(old$cell_id) || !is.numeric(old$mff) || any(!is.finite(old$mff)) || any(old$mff < 0))
      stop("Invalid unique Old color reference or color modes.")
    positive <- old$mff[old$mff > 0]
    caps <- c(actual_max = if (nrow(old)) max(old$mff) else 0,
      positive_quantile = if (length(positive)) unname(quantile(positive, probability, type = 7)) else 0)
    data.frame(mode = modes, reference_upper = unname(caps[modes]),
      display_upper = ifelse(unname(caps[modes]) > 0, unname(caps[modes]), 1),
      n_reference_Old = nrow(old), n_positive = length(positive),
      quantile_probability = ifelse(modes == "positive_quantile", probability, NA_real_))
  }
  ordered.panel <- function(selected, scope, grouping, ct, score, quartile) {
    z <- selected[selected$scope == scope & selected$grouping == grouping & selected$selection_celltype == ct &
      selected$score == score & selected$quartile == quartile, , drop = FALSE]
    z[order(z$mff, z$cell_id), , drop = FALSE]
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
  write.checked.csv <- function(d, path) {
    token <- ".__OLD_ALLHIGH_MISSING__."
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

  missing.value <- function(v) is.na(v) | v == "NA" | v == "" | grepl("^\\.__MFF_MISSING__\\._*$", v)
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
  filter.top.high <- function(selected, classified) {
    key <- function(d) paste(d$scope, d$cell_id, sep = "\034")
    if (anyDuplicated(key(classified)) || anyNA(selected$quartile) ||
        any(!selected$quartile %in% c("Bottom25", "Top25")) || any(selected$group != "Old"))
      stop("Invalid Old quartile/classification keys.")
    top <- selected[selected$quartile == "Top25", , drop = FALSE]
    idx <- match(key(top), key(classified))
    if (anyNA(idx)) stop("Top25 profiles missing from original scope classification.")
    if (!isTRUE(all.equal(top$mff, classified$mff[idx], tolerance = 1e-10, check.attributes = FALSE)))
      stop("Top25 MFF differs from original cutoff input.")
    top$cutoff_all <- classified$cutoff_all[idx]; top$region_all <- classified$region_all[idx]
    top$fit_status <- classified$fit_status[idx]
    top[which(top$region_all == "High"), , drop = FALSE]
  }
  analysis.tables <- function(selected, parent.statistics, classified, mode = "bottom_all") {
    if (length(mode) != 1L || is.na(mode) || !mode %in% c("bottom_all", "both_high", "top_only"))
      stop("comparison.mode must be bottom_all, both_high or top_only.")
    keys <- function(d) paste(d$scope, d$grouping, d$selection_celltype, d$score, sep = "\034")
    pk <- function(d) paste(d$scope, d$cell_id, sep = "\034")
    s <- parent.statistics
    if (!nrow(s) || anyDuplicated(keys(s)) || anyNA(selected$quartile) ||
        any(!selected$quartile %in% c("Bottom25", "Top25")) || any(selected$group != "Old") ||
        anyDuplicated(selected[, c("scope", "grouping", "selection_celltype", "score", "cell_id")]) ||
        any(!keys(selected) %in% keys(s)) || anyDuplicated(pk(classified))) stop("Invalid parent selection keys.")
    idx <- match(pk(selected), pk(classified))
    if (anyNA(idx)) stop("Parent quartiles missing from original scope classification.")
    for (nm in c("donor", "celltype", "tech", "group"))
      if (!identical(as.character(selected[[nm]]), as.character(classified[[nm]][idx]))) stop("Parent metadata differs: ", nm)
    if (!isTRUE(all.equal(selected$mff, classified$mff[idx], tolerance = 1e-10, check.attributes = FALSE)))
      stop("Parent MFF differs from original cutoff input.")
    d <- selected
    d$cutoff_all <- classified$cutoff_all[idx]; d$region_all <- classified$region_all[idx]
    d$fit_status <- classified$fit_status[idx]
    if (any(!is.na(d$region_all) & !d$region_all %in% c("High", "Low"))) stop("Invalid original region_all.")
    keep <- d$quartile == "Top25" & !is.na(d$region_all) & d$region_all == "High"
    if (mode == "bottom_all") keep <- keep | d$quartile == "Bottom25"
    if (mode == "both_high") keep <- keep | (d$quartile == "Bottom25" & !is.na(d$region_all) & d$region_all == "High")
    chosen <- d[which(keep), , drop = FALSE]
    stats <- retention <- list()
    parent.fields <- c("scope", "grouping", "selection_celltype", "score", "n_Old", "q25", "q75",
      "n_ties_q25", "n_ties_q75", "selection_reason")
    for (i in seq_len(nrow(s))) {
      st <- s[i, , drop = FALSE]; original <- d[keys(d) == keys(st), , drop = FALSE]
      z <- chosen[keys(chosen) == keys(st), , drop = FALSE]
      if (sum(original$quartile == "Bottom25") != st$n_Bottom || sum(original$quartile == "Top25") != st$n_Top)
        stop("Parent quartile counts differ: row ", i)
      r <- mff.test(z$mff, z$quartile)
      if (mode == "top_only") { r$W <- r$p_value <- NA_real_; r$test_reason <- "top_only_no_comparison" }
      stats[[i]] <- cbind(st[, parent.fields, drop = FALSE],
        data.frame(comparison_mode = mode, n_parent_Bottom = st$n_Bottom, n_parent_Top = st$n_Top), r,
        data.frame(fraction_Bottom_of_Old = if (st$n_Old) r$n_Bottom/st$n_Old else NA_real_,
          fraction_Top_of_Old = if (st$n_Old) r$n_Top/st$n_Old else NA_real_,
          n_donors_Bottom = length(unique(z$donor[z$quartile == "Bottom25"])),
          n_donors_Top = length(unique(z$donor[z$quartile == "Top25"]))))
      for (g in c("Bottom25", "Top25")) {
        p <- original[original$quartile == g, , drop = FALSE]
        excluded <- mode == "top_only" && g == "Bottom25"
        filtered <- g == "Top25" || mode == "both_high"
        n <- sum(z$quartile == g)
        retention[[length(retention) + 1L]] <- data.frame(st[, parent.fields[1:4], drop = FALSE],
          comparison_mode = mode, quartile = g, filter_applied = filtered,
          n_parent = nrow(p), n_parent_High = sum(p$region_all == "High", na.rm = TRUE),
          n_parent_Low = sum(p$region_all == "Low", na.rm = TRUE), n_parent_Unclassified = sum(is.na(p$region_all)),
          n_selected = n, n_removed_Low = if (filtered) sum(p$region_all == "Low", na.rm = TRUE) else 0L,
          n_removed_Unclassified = if (filtered) sum(is.na(p$region_all)) else 0L,
          n_removed_for_mode = if (excluded) nrow(p) else 0L,
          fraction_retained = if (nrow(p)) n/nrow(p) else NA_real_)
      }
    }
    stats <- do.call(rbind, stats); retention <- do.call(rbind, retention)
    for (scope in unique(stats$scope)) for (grouping in unique(stats$grouping)) {
      use <- which(stats$scope == scope & stats$grouping == grouping)
      stats$p_adj_BH[use] <- p.adjust(stats$p_value[use], method = "BH", n = length(use))
    }
    rownames(stats) <- rownames(retention) <- rownames(chosen) <- NULL
    stopifnot(all(retention$n_selected + retention$n_removed_Low + retention$n_removed_Unclassified +
      retention$n_removed_for_mode == retention$n_parent))
    list(statistics = stats, selected = chosen, retention = retention)
  }
  same.table <- function(saved, expected, key.fields, fields, what) {
    key <- function(d) do.call(paste, c(d[key.fields], sep = "\034"))
    if (!is.data.frame(saved) || !all(c(key.fields, fields) %in% names(saved)) ||
        anyNA(saved[key.fields]) || anyDuplicated(key(saved)) || nrow(saved) != nrow(expected) ||
        !setequal(key(saved), key(expected))) stop("Stale/incomplete script-44 ", what, "; rerun 44.")
    saved <- saved[match(key(expected), key(saved)), , drop = FALSE]
    for (nm in fields) {
      ok <- if (is.numeric(expected[[nm]])) is.numeric(saved[[nm]]) &&
        isTRUE(all.equal(saved[[nm]], expected[[nm]], tolerance = 1e-10, check.attributes = FALSE)) else
        identical(as.character(saved[[nm]]), as.character(expected[[nm]]))
      if (!ok) stop("Script-44 ", what, " differs: ", nm, "; rerun 44.")
    }
    invisible(TRUE)
  }
  verify.umap <- function(cache, verified, current, embedding, types) {
    if (!is.list(cache)) stop("Run script 44 first.")
    same.table(cache$old_profiles, verified$old, "cell_id", names(verified$old), "Old reference/coordinates")
    same.table(cache$selected_profiles, verified$selected,
      c("scope", "grouping", "selection_celltype", "score", "cell_id"), names(verified$selected), "memberships/coordinates")
    limits <- cache$color_limits; bounds <- cache$bounds
    if (!is.data.frame(limits) || !all(c("mode", "quantile_probability", "reference_upper", "display_upper") %in% names(limits)) ||
        !nrow(limits) || anyNA(limits$mode) || anyDuplicated(limits$mode) ||
        any(!limits$mode %in% c("actual_max", "positive_quantile"))) stop("Invalid script-44 color limits.")
    probability <- if ("positive_quantile" %in% limits$mode) limits$quantile_probability[limits$mode == "positive_quantile"] else .95
    expected <- shared.limits(verified$old, limits$mode, probability)
    same.table(limits, expected, "mode", names(expected), "shared color limits")
    wanted <- data.frame(grouping = c("pooled", rep("celltype", length(types))), selection_celltype = c("ALL", types))
    fields <- c("xmin", "xmax", "ymin", "ymax")
    if (!is.data.frame(bounds) || !all(c(names(wanted), fields) %in% names(bounds))) stop("Invalid script-44 crop windows.")
    same.table(bounds, wanted, names(wanted), names(wanted), "crop keys")
    if (any(!vapply(bounds[fields], is.numeric, logical(1))) || any(!is.finite(as.matrix(bounds[fields]))) ||
        any(bounds$xmin >= bounds$xmax | bounds$ymin >= bounds$ymax)) stop("Invalid saved crop coordinates.")
    for (i in seq_len(nrow(bounds))) {
      w <- bounds[i, ]; z <- verified$old
      if (w$grouping == "celltype") z <- z[z$celltype == w$selection_celltype, , drop = FALSE]
      if (!nrow(z)) {
        z <- current
        if (w$grouping == "celltype") z <- z[z$celltype == w$selection_celltype, , drop = FALSE]
        z$UMAP_1 <- embedding[z$cell_id, 1]; z$UMAP_2 <- embedding[z$cell_id, 2]
      }
      if (any(z$UMAP_1 < w$xmin | z$UMAP_1 > w$xmax | z$UMAP_2 < w$ymin | z$UMAP_2 > w$ymax))
        stop("Script-44 crop no longer covers original reference coordinates.")
    }
    list(limits = limits, bounds = bounds)
  }
  panel.summary <- function(selected, statistics, limits, groups) {
    rows <- list()
    for (i in seq_len(nrow(statistics))) for (mode in limits$mode) for (g in groups) {
      st <- statistics[i, ]; limit <- limits[limits$mode == mode, ]
      z <- ordered.panel(selected, st$scope, st$grouping, st$selection_celltype, st$score, g)
      n.parent <- if (g == "Top25") st$n_parent_Top else st$n_parent_Bottom
      rows[[length(rows) + 1L]] <- data.frame(st[, c("scope", "grouping", "selection_celltype", "score")],
        comparison_mode = st$comparison_mode, mode = mode, quartile = g, n_parent = n.parent,
        n_profiles = nrow(z), n_positive = sum(z$mff > 0), n_zero = sum(z$mff == 0),
        panel_max = if (nrow(z)) max(z$mff) else NA_real_, shared_color_upper = limit$reference_upper,
        display_upper = limit$display_upper, n_above_color_upper = sum(z$mff > limit$reference_upper),
        selection_reason = st$selection_reason)
    }
    do.call(rbind, rows)
  }

  # 3. Require both previous steps. Verify their cells and inherit 44's coordinates/scales/windows.
  required.objects <- c(input.object, cache.object, umap.cache.object)
  if (!all(vapply(required.objects, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)))
    stop("Complete GO-enabled 34, then 43 and 44, in this R session first.")
  if (!comparison.mode %in% c("bottom_all", "both_high", "top_only")) stop("Invalid comparison.mode.")
  stopifnot(identical(names(score.columns), names(score.labels)), n.columns >= 1L, n.pair.columns >= 1L, point.size > 0)
  x <- get(input.object, envir = .GlobalEnv)
  cache <- get(cache.object, envir = .GlobalEnv); umap.cache <- get(umap.cache.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), assay %in% Assays(x), "data" %in% Layers(x[[assay]]), reduction %in% Reductions(x))
  if (gene %in% names(x[[]])) stop("Metadata column MFF shadows the RNA feature; resolve it before FeaturePlot.")
  ids <- colnames(x); rna <- GetAssayData(x, assay = assay, layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == gene) == 1L)
  mff <- setNames(as.numeric(rna[gene, ]), ids); rm(rna)
  current <- make.current(x[[]], mff, score.columns, ids); types <- sort(unique(current$celltype))
  embedding <- Embeddings(x[[reduction]])
  verified <- verify.cache(cache, current, embedding, score.columns, scopes, groupings, types)
  inherited <- verify.umap(umap.cache, verified, current, embedding, types)
  old <- verified$old; limits <- inherited$limits; bounds <- inherited$bounds; color.modes <- limits$mode
  fits <- read.saved(file.path(mclust.dir, "02.MFF.mclust_fit_details.csv"))
  profiles <- setNames(lapply(scopes, function(scope) read.saved(file.path(mclust.dir,
    paste0("05.MFF.profile_classification.", scope, ".csv.gz")))), scopes)
  classified <- join.all.classification(current, profiles, fits, scopes)
  result <- analysis.tables(verified$selected, verified$statistics, classified, comparison.mode)
  statistics <- result$statistics; selected <- result$selected; retention <- result$retention
  summary <- panel.summary(selected, statistics, limits, group.order)
  file.map <- data.frame(celltype = types, file_stem = sprintf("%03d", seq_along(types)))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(statistics, file.path(out, "01.filtered_MFF_statistics.csv"))
  write.checked.csv(retention, file.path(out, "02.quartile_filter_retention.csv"))
  write.checked.csv(selected, file.path(out, "03.selected_profiles.csv.gz"))
  write.checked.csv(classified, file.path(out, "04.original_all_classification.audit.csv.gz"))
  write.checked.csv(limits, file.path(out, "05.inherited_color_limits.csv"))
  write.checked.csv(bounds, file.path(out, "06.inherited_UMAP_bounds.csv"))
  write.checked.csv(summary, file.path(out, "07.UMAP_panel_summary.csv"))
  write.checked.csv(file.map, file.path(out, "08.celltype_file_map.csv"))

  # 4. Conditional violin comparison. The score quartiles are never reselected.
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
    p <- ggplot(d, aes(quartile, mff)) + geom_blank() +
      geom_violin(data = regular, aes(fill = quartile), trim = TRUE, scale = "width", width = .8) +
      geom_boxplot(data = regular, width = .15, fill = "white", colour = "black",
        outlier.shape = NA, linewidth = .4)
    for (g in group.order) {
      z <- special[special$quartile == g, , drop = FALSE]
      if (nrow(z) == 1L) p <- p + geom_point(data = z, size = 2)
      if (nrow(z) > 1L) p <- p + geom_errorbar(data = z[1, , drop = FALSE], aes(ymin = mff, ymax = mff), width = .45)
    }
    count <- c(Bottom25 = stat$n_Bottom, Top25 = stat$n_Top)
    parent <- c(Bottom25 = stat$n_parent_Bottom, Top25 = stat$n_parent_Top)
    labels <- setNames(paste0(ifelse(group.order == "Top25", "Top 25%\nMFF High",
      if (comparison.mode == "both_high") "Bottom 25%\nMFF High" else "Bottom 25%\nAll"),
      "\n(n = ", count[group.order], "; parent = ", parent[group.order], ")"), group.order)
    if (length(group.order) == 2L) {
      bracket <- bounds[2] + .12 * span
      p <- p + annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .4) +
        annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket - .04 * span, yend = bracket, linewidth = .4) +
        annotate("text", x = 1.5, y = bounds[2] + .45 * span, label = format.p(stat$p_value), size = annotation.size)
    }
    p + scale_fill_manual(values = group.colors, guide = "none") +
      scale_x_discrete(limits = group.order, labels = labels, drop = FALSE) +
      coord_cartesian(ylim = c(max(0, bounds[1] - .08 * span), bounds[2] + .65 * span)) +
      labs(title = title, subtitle = if (stat$selection_reason != "") "Original quartile groups unavailable" else NULL,
        x = NULL, y = "MFF expression (stored RNA/data)") +
      theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 13), axis.text = element_text(size = 15), panel.grid.minor = element_blank())
  }
  save.png <- function(p, path, width, height) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(path, width = width, height = height, res = plot.res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  violin.panel <- function(scope, grouping, ct, score, title) {
    s <- statistics[statistics$scope == scope & statistics$grouping == grouping &
      statistics$selection_celltype == ct & statistics$score == score, , drop = FALSE]
    d <- selected[selected$scope == scope & selected$grouping == grouping &
      selected$selection_celltype == ct & selected$score == score, , drop = FALSE]
    stopifnot(nrow(s) == 1L); mff.violin(d, s, title)
  }
  selection.caption <- paste0("Old only; original script-43 score quartiles. ",
    if (comparison.mode == "top_only") "Top25 intersect MFF all-High only." else
      paste0("Blue: ", group.labels["Bottom25"], "; red: ", group.labels["Top25"], "."),
    " High uses original Young+Old, cell-type-specific ALL-mode fits separately for combined/scRNA/snRNA.")
  violin.caption <- paste(selection.caption, if (comparison.mode == "top_only") "No between-group test." else
    "Raw unpaired cell Wilcoxon p; groups selected on MFF, so the comparison is conditional on MFF selection. Donor/Country/composition unadjusted.")
  manifest <- list(); pngs <- character()
  record <- function(path, kind, grouping, scope, score, mode, ct, width, height) {
    pngs <<- c(pngs, path)
    manifest[[length(manifest) + 1L]] <<- data.frame(kind = kind, grouping = grouping, scope = scope,
      score = score, color_mode = mode, selection_celltype = ct, comparison_mode = comparison.mode,
      width = width, height = height, path = path)
  }
  for (scope in scopes) {
    if ("pooled" %in% groupings) {
      plots <- lapply(names(score.columns), function(nm) violin.panel(scope, "pooled", "ALL", nm, score.labels[nm]))
      p <- wrap_plots(plots, ncol = length(score.columns)) + plot_annotation(
        title = paste(scope.labels[scope], "Old | All cell types pooled", selection.label, sep = " | "),
        caption = paste(strwrap(violin.caption, width = 175), collapse = "\n"),
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "violin/pooled", scope, "Old.three_scores.Top25_allHigh.violin.png")
      w <- panel.width * length(score.columns); h <- panel.height
      save.png(p, path, w, h); record(path, "violin", "pooled", scope, "ALL", "", "ALL", w, h)
    }
    if ("celltype" %in% groupings) for (nm in names(score.columns)) {
      plots <- lapply(types, function(ct) violin.panel(scope, "celltype", ct, nm, ct))
      cols <- min(n.columns, length(types))
      p <- wrap_plots(plots, ncol = cols) + plot_annotation(
        title = paste(scope.labels[scope], "Old", score.labels[nm], "Within each cell type", selection.label, sep = " | "),
        caption = paste(strwrap(violin.caption, width = 220), collapse = "\n"),
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "violin/celltype", scope, paste0("Old.", nm, ".Top25_allHigh.violin.overview.png"))
      w <- panel.width * cols; h <- panel.height * ceiling(length(types)/cols)
      save.png(p, path, w, h); record(path, "violin", "celltype_overview", scope, nm, "", "ALL", w, h)
    }
  }

  # 5. FeaturePlot style from 44; retain its complete Old color reference and axes.
  DefaultAssay(x) <- assay
  template <- Seurat::FeaturePlot(object = x, features = gene, cells = ids[1], reduction = reduction,
    slot = "data", order = TRUE, min.cutoff = 0, max.cutoff = max(limits$display_upper),
    keep.scale = NULL, cols = c(expression.colors[1], tail(expression.colors, 1)), pt.size = point.size,
    label = FALSE, raster = FALSE, coord.fixed = FALSE, combine = FALSE)[[1]]
  umap.panel <- function(z, stat, limit, window, title, quartile) {
    cap <- limit$display_upper
    if (nrow(z)) {
      p <- Seurat::FeaturePlot(object = x, features = gene, cells = z$cell_id, reduction = reduction,
        slot = "data", order = TRUE, min.cutoff = 0, max.cutoff = cap, keep.scale = NULL,
        cols = c(expression.colors[1], tail(expression.colors, 1)), pt.size = point.size,
        label = FALSE, raster = FALSE, coord.fixed = FALSE, combine = FALSE)[[1]]
    } else p <- ggplot(z, aes(UMAP_1, UMAP_2, color = mff)) + geom_blank() + template$theme
    mode.label <- if (limit$mode == "actual_max") "Shared actual maximum" else paste0("Shared positive MFF Q", 100 * limit$quantile_probability)
    n.parent <- if (quartile == "Top25") stat$n_parent_Top else stat$n_parent_Bottom
    # Quiet the expected scale-replacement message; warnings/errors remain visible.
    p <- suppressMessages(p + scale_color_gradientn(colours = expression.colors, limits = c(0, cap), oob = scales::squish,
      breaks = c(0, cap/2, cap), labels = function(v) format(signif(v, 3), trim = TRUE), name = "MFF\nRNA/data",
      guide = guide_colorbar(barheight = grid::unit(4, "cm")))) +
      coord_cartesian(xlim = unname(c(window$xmin, window$xmax)), ylim = unname(c(window$ymin, window$ymax)), expand = FALSE) +
      labs(title = paste(paste(strwrap(title, width = 30), collapse = "\n"),
        paste(strwrap(group.labels[quartile], width = 25), collapse = "\n"), sep = "\n"),
        subtitle = paste0(mode.label, " = ", signif(limit$reference_upper, 4)),
        caption = paste0("n = ", nrow(z), " / parent quartile = ", n.parent, "; MFF > 0: ", sum(z$mff > 0),
          if (limit$mode == "positive_quantile") "\nAbove-limit values use the maximum color." else ""), x = "umap_1", y = "umap_2") +
      theme(text = element_text(size = font.base), plot.title = element_text(size = font.title, face = "bold", hjust = .5),
        axis.title = element_text(size = font.axis), axis.text = element_text(size = font.axis),
        legend.title = element_text(size = font.legend), legend.text = element_text(size = font.legend), plot.margin = margin(12, 12, 12, 12))
    if (!nrow(z)) p <- p + annotate("text", x = mean(c(window$xmin, window$xmax)), y = mean(c(window$ymin, window$ymax)),
      label = "No selected cells", color = "#666666", size = 5)
    p
  }
  paired.umap <- function(scope, grouping, ct, score, limit) {
    stat <- statistics[statistics$scope == scope & statistics$grouping == grouping &
      statistics$selection_celltype == ct & statistics$score == score, , drop = FALSE]
    window <- bounds[bounds$grouping == grouping & bounds$selection_celltype == ct, , drop = FALSE]
    stopifnot(nrow(stat) == 1L, nrow(window) == 1L)
    title <- if (grouping == "pooled") "Old | All cell types" else paste("Old", ct, sep = " | ")
    plots <- lapply(group.order, function(g) {
      z <- ordered.panel(selected, scope, grouping, ct, score, g)
      n <- if (g == "Bottom25") stat$n_Bottom else stat$n_Top
      stopifnot(nrow(z) == n); umap.panel(z, stat, limit, window, title, g)
    })
    wrap_plots(plots, ncol = length(group.order), guides = "collect")
  }
  umap.caption <- paste(selection.caption, "Saved UMAP; only each selected group's cells drawn. Color limits and crop windows inherited from 44.")
  for (mode in color.modes) for (scope in scopes) for (nm in names(score.columns)) {
    limit <- limits[limits$mode == mode, , drop = FALSE]
    heading <- paste(scope.labels[scope], score.labels[nm], "Old", selection.label, sep = " | ")
    if ("pooled" %in% groupings) {
      p <- paired.umap(scope, "pooled", "ALL", nm, limit) + plot_annotation(
        title = paste(strwrap(heading, width = 35L * length(group.order)), collapse = "\n"),
        caption = paste(strwrap(umap.caption, width = 120), collapse = "\n"),
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "UMAP/pooled", scope, mode, paste0("Old.", nm, ".Top25_allHigh.MFF_UMAP.png"))
      w <- length(group.order) * panel.width; h <- panel.height
      save.png(p, path, w, h); record(path, "UMAP", "pooled", scope, nm, mode, "ALL", w, h)
    }
    if ("celltype" %in% groupings) {
      pairs <- lapply(types, function(ct) paired.umap(scope, "celltype", ct, nm, limit))
      cols <- min(n.pair.columns, length(types)); rows <- ceiling(length(types)/cols)
      p <- wrap_plots(pairs, ncol = cols, guides = "collect") + plot_annotation(
        title = paste(heading, "Within each cell type", sep = " | "),
        caption = paste(strwrap(umap.caption, width = 180), collapse = "\n"),
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "UMAP/celltype", scope, mode, paste0("Old.", nm, ".Top25_allHigh.MFF_UMAP.overview.png"))
      w <- length(group.order) * panel.width * cols; h <- panel.height * rows
      save.png(p, path, w, h); record(path, "UMAP", "celltype_overview", scope, nm, mode, "ALL", w, h)
      if (save.individual.celltype.pairs) for (i in seq_along(types)) {
        path <- file.path(out, "UMAP/celltype", scope, mode, file.map$file_stem[i], paste0("Old.", nm, ".Top25_allHigh.MFF_UMAP.png"))
        p <- pairs[[i]] + plot_annotation(title = paste(strwrap(heading, width = 35L * length(group.order)), collapse = "\n"),
          theme = theme(plot.title = element_text(size = 22, face = "bold")))
        w <- length(group.order) * panel.width; h <- panel.height
        save.png(p, path, w, h); record(path, "UMAP", "celltype_pair", scope, nm, mode, types[i], w, h)
      }
    }
  }
  expected.violin <- length(scopes) * (as.integer("pooled" %in% groupings) +
    as.integer("celltype" %in% groupings) * length(score.columns))
  expected.umap <- length(color.modes) * length(scopes) * length(score.columns) *
    (as.integer("pooled" %in% groupings) + as.integer("celltype" %in% groupings) *
      (1L + as.integer(save.individual.celltype.pairs) * length(types)))
  stopifnot(length(pngs) == expected.violin + expected.umap, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  manifest <- do.call(rbind, manifest); write.checked.csv(manifest, file.path(out, "09.PNG_manifest.csv"))
  writeLines(c("Run after successful 43 and 44 in the same server R session; hlma.scores plus both caches required.",
    "Reuse original score quartile memberships without reranking after the MFF filter.",
    "Top25 is intersected with original MFF all-mode High, separately for combined/scRNA/snRNA and original Annotation.",
    "Original script-32 cutoffs used Young+Old together per cell type and included MFF=0. No new mclust fitting.",
    "Read 02.MFF.mclust_fit_details.csv and all three 05.MFF.profile_classification.<scope>.csv.gz files.",
    "Validate complete original profiles, IDs, metadata, expression, all-fit counts, convergence and saved labels.",
    "High means original full-precision MFF > cutoff. Low includes equality/zero. Saved labels retained at CSV-rounded boundaries.",
    "Unavailable original cutoff stays unclassified and cannot enter a filtered High group. Positive-mode results are ignored.",
    paste0("comparison.mode = ", comparison.mode),
    "bottom_all: unchanged Bottom25 whole versus Top25 intersect High; both_high: both tails intersect High; top_only: only Top25 intersect High.",
    "Retention CSV records parent quartile sizes, original High/Low/unclassified, selected/removed counts and retention fractions.",
    "Inclusive score-boundary ties and unavailable score quartiles are inherited from 43; no filtering of negative/zero scores.",
    "Violin: blue Bottom versus red Top with narrow white boxplots. Singleton/constant data use a point/line; empty groups keep labeled axes.",
    "Boxplots hide outlier markers only; all values remain in the violin/statistics. Zero IQR collapses the box to a line.",
    "Two-sided unpaired cell Wilcoxon exact=FALSE/correct=TRUE; insufficient/constant groups p=NA. Top-only has no test.",
    "BH across score/cell-type comparisons within each scope/grouping, including unavailable tests in the family.",
    "MFF is used to select High and is also the plotted outcome. These are conditional comparisons, not independent evidence of a score-MFF association.",
    "Cell tests do not adjust donor dependence, Country, technology or composition. Combined/technical and pooled/celltype results overlap.",
    "UMAP: exactly script-44 coordinates and common color limits/crop windows, with its grey-to-orange-to-dark-red continuous MFF palette.",
    "Only selected cells are shown; high MFF draws last; no downsampling, PCA, RunUMAP, score/normalization/mclust recalculation or RDS I/O.",
    "Each panel labels selected n / parent quartile n. Color Q95 saturates display only; original MFF values remain in CSV.",
    paste0("PNG outputs: violin=", expected.violin, "; UMAP=", expected.umap, "; total=", length(pngs)),
    "Default both groupings, three scopes/scores and both color modes: 12 violin + 36 UMAP = 48 PNGs.",
    "Figures are written into a separate comparison-mode folder; individual celltype UMAP pairs are optional.",
    paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  assign("hlma.old.score.Top25.allHigh", list(old_profiles = old, selected_profiles = selected,
    statistics = statistics, retention = retention, classification = classified, color_limits = limits,
    bounds = bounds, panel_summary = summary, manifest = manifest,
    settings = list(comparison_mode = comparison.mode, groupings = groupings, scopes = scopes,
      score_columns = score.columns, mclust_mode = "all", original_cutoff_directory = mclust.dir)), envir = .GlobalEnv)
  cat("\nHLMA_OLD_SCORE_TOP25_ALLHIGH_VIOLIN_UMAP_COMPLETE\nPNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
