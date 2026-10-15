# Run ALL code in the SAME SERVER R session after script 34 completes.
# Reuses hlma.scores: MFF versus SaulSenMayo / Hallmark IFN-alpha.
# No score recalculation or RDS read/write. All MFF zeros and signed scores retained.
local({
  packages <- c("Seurat", "ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. Reruns overwrite these correlation outputs.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/MFF_correlation_SenMayo_IFNalpha")
  input.object <- "hlma.scores"
  gene <- "MFF"
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  views <- c("All", "Young", "Old")
  plot.units <- "cell"  # c("cell", "donor") also plots donor-mean scatter. Both units' CSVs always saved.
  save.individual.scatter <- FALSE
  show.regression.line <- TRUE
  group.order <- c("Young", "Old"); group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L; n.columns <- 5L
  annotation.size <- 8

  # 2. Align each score and MFF value to the SAME original cell/nucleus.
  make.pairs <- function(meta, mff, columns, ids = rownames(meta)) {
    required <- c("Annotation", "sample", "tech", "age_pop", "Country", unname(columns))
    if (!all(required %in% names(meta))) stop("Missing metadata/score columns: ", paste(setdiff(required, names(meta)), collapse = ", "))
    if (!nrow(meta) || !identical(rownames(meta), ids) || anyNA(ids) || anyDuplicated(ids) || any(ids == ""))
      stop("Invalid cell IDs or metadata cell order.")
    if (!is.numeric(mff) || !identical(names(mff), ids) || length(mff) != length(ids) ||
        any(!is.finite(mff)) || any(mff < 0)) stop("Invalid MFF values or cell order.")
    for (nm in c("Annotation", "sample", "tech", "age_pop", "Country"))
      if (anyNA(meta[[nm]]) || any(!nzchar(trimws(as.character(meta[[nm]]))))) stop("Missing metadata: ", nm)
    if (!all(meta$tech %in% c("scRNA", "snRNA")) || !all(meta$age_pop %in% c("young_pop", "old_pop")))
      stop("Unexpected technology or age_pop labels.")
    d <- data.frame(cell_id = ids, donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
      tech = as.character(meta$tech), Country = as.character(meta$Country),
      group = ifelse(meta$age_pop == "young_pop", "Young", "Old"), mff = unname(mff), stringsAsFactors = FALSE)
    for (nm in unname(columns)) {
      if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]]))) stop("Nonfinite/nonnumeric score: ", nm)
      d[[nm]] <- meta[[nm]]
    }
    info <- unique(d[, c("donor", "group", "Country")])
    if (anyDuplicated(info$donor)) stop("Inconsistent group or Country within donor.")
    d
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
  correlation.tables <- function(cells, columns, types, scope.names, view.names) {
    donor <- data.frame(scope = character(), celltype = character(), donor = character(),
      group = character(), Country = character(), n_profiles = integer(), n_scRNA = integer(),
      n_snRNA = integer(), mff = numeric(), stringsAsFactors = FALSE)
    for (nm in unname(columns)) donor[[nm]] <- numeric()
    for (scope in scope.names) {
      b <- if (scope == "combined") cells else cells[cells$tech == scope, , drop = FALSE]
      for (ct in types) {
        d <- b[b$celltype == ct, , drop = FALSE]
        if (nrow(d)) for (idx in split(seq_len(nrow(d)), d$donor)) {
          z <- d[idx, , drop = FALSE]
          row <- data.frame(scope = scope, celltype = ct, donor = z$donor[1], group = z$group[1],
            Country = z$Country[1], n_profiles = nrow(z), n_scRNA = sum(z$tech == "scRNA"),
            n_snRNA = sum(z$tech == "snRNA"), mff = mean(z$mff))
          for (nm in unname(columns)) row[[nm]] <- mean(z[[nm]])
          donor <- rbind(donor, row)
        }
      }
    }
    statistics <- list()
    for (scope in scope.names) for (ct in types) for (nm in names(columns)) for (view in view.names) {
      d <- cells[cells$celltype == ct & (scope == "combined" | cells$tech == scope), , drop = FALSE]
      ds <- donor[donor$scope == scope & donor$celltype == ct, , drop = FALSE]
      if (view != "All") { d <- d[d$group == view, , drop = FALSE]; ds <- ds[ds$group == view, , drop = FALSE] }
      for (unit in c("cell", "donor")) {
        z <- if (unit == "cell") d else ds
        r <- spearman.pair(z$mff, z[[columns[[nm]]]])
        statistics[[length(statistics) + 1L]] <- cbind(data.frame(scope = scope, celltype = ct,
          score = nm, view = view, unit = unit, n_profiles = nrow(d), n_donors = nrow(ds),
          n_profiles_MFF_positive = sum(d$mff > 0)), r)
      }
    }
    s <- do.call(rbind, statistics); s$p_adj_BH <- NA_real_
    for (scope in scope.names) for (unit in c("cell", "donor")) for (view in view.names) {
      idx <- which(s$scope == scope & s$unit == unit & s$view == view)
      s$p_adj_BH[idx] <- p.adjust(s$p_value[idx], method = "BH")
    }
    rownames(donor) <- rownames(s) <- NULL
    list(donor = donor, statistics = s)
  }
  write.checked.csv <- function(d, file) {
    token <- ".__HLMA_COR_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    con <- if (grepl("[.]gz$", file)) gzfile(file, "wt", encoding = "UTF-8") else file(file, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, con, row.names = FALSE, na = token), finally = close(con))
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else
      if (is.logical(v)) "logical" else "character", character(1))
    con <- if (grepl("[.]gz$", file)) gzfile(file, "rt", encoding = "UTF-8") else file(file, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(con, check.names = FALSE, colClasses = classes, na.strings = token,
      stringsAsFactors = FALSE), finally = close(con))
    stopifnot(identical(dim(restored), dim(d)), identical(names(restored), names(d)))
    for (nm in names(d)) {
      if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(restored[[nm]], d[[nm]], tolerance = 1e-10))) else
        stopifnot(identical(restored[[nm]], d[[nm]]))
    }
    invisible(file)
  }

  # 3. Reuse already calculated scores, with no source()/AddModuleScore()/saveRDS().
  if (!exists(input.object, envir = .GlobalEnv, inherits = FALSE))
    stop("Run the full script 34 in this R session first: missing ", input.object)
  x <- get(input.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  if (!all(plot.units %in% c("cell", "donor")) || !length(plot.units) || anyDuplicated(plot.units)) stop("Invalid plot.units.")
  ids <- colnames(x); meta <- x[[]]
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == gene) == 1L)
  mff <- setNames(as.numeric(rna[gene, ]), ids); rm(rna)
  cells <- make.pairs(meta, mff, score.columns, ids)
  info <- x@misc$HLMA_two_scores
  exclusion.confirmed <- !is.null(info) && "MFF" %in% info$pool_excludes &&
    all(names(score.columns) %in% names(info$matched_features)) &&
    !any(vapply(info$matched_features[names(score.columns)], function(v) "MFF" %in% v, logical(1)))
  if (!exclusion.confirmed) warning("MFF exclusion from score signatures/control pool is unconfirmed; verify score provenance.")
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE)); n.rows <- ceiling(length(celltypes) / n.columns)
  cat("Reusing ", input.object, ": ", nrow(cells), " profiles; computing MFF-score Spearman correlations.\n", sep = "")
  result <- correlation.tables(cells, score.columns, celltypes, scopes, views)
  statistics <- result$statistics; donor <- result$donor
  stopifnot(sum(donor$n_profiles) == 2L * nrow(cells))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(statistics[statistics$unit == "cell", ], file.path(out, "01.MFF_score.cell_correlations.csv"))
  write.checked.csv(statistics[statistics$unit == "donor", ], file.path(out, "02.MFF_score.donor_correlations.csv"))
  write.checked.csv(donor, file.path(out, "03.MFF_score.donor_means.csv"))
  write.checked.csv(cells, file.path(out, "04.MFF_score.profile_pairs.csv.gz"))

  # 4. Scatter panels. Linear lines are descriptive; rho/p come from Spearman ranks.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  scatter.plot <- function(full, stat, unit, view, title, label) {
    d <- if (view == "All") full else full[full$group == view, , drop = FALSE]
    d$group <- factor(d$group, levels = group.order)
    bounds <- if (nrow(full)) range(full$score) else c(-.1, .1)
    span <- max(diff(bounds), .05)
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
      annotate("text", x = Inf, y = bounds[2] + .60 * span, label = text,
        hjust = 1.05, vjust = 1, size = annotation.size) +
      coord_cartesian(xlim = xb, ylim = c(bounds[1] - .08 * span, bounds[2] + .65 * span)) +
      labs(title = paste(strwrap(title, width = 40), collapse = "\n"),
        x = if (unit == "cell") "MFF expression (stored RNA/data)" else "Donor-mean MFF (stored RNA/data)",
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
  for (unit in plot.units) for (scope in scopes) for (nm in names(score.columns)) for (view in views) {
    b <- if (unit == "cell") {
      if (scope == "combined") cells else cells[cells$tech == scope, , drop = FALSE]
    } else donor[donor$scope == scope, , drop = FALSE]
    folder <- file.path(out, scope, unit, nm)
    plots <- list()
    for (ct in celltypes) {
      full <- b[b$celltype == ct, , drop = FALSE]; full$score <- full[[score.columns[[nm]]]]
      stat <- statistics[statistics$scope == scope & statistics$unit == unit & statistics$score == nm &
        statistics$celltype == ct & statistics$view == view, ]
      stopifnot(nrow(stat) == 1L)
      p <- scatter.plot(full, stat, unit, view, paste(scope.labels[scope], ct, sep = " | "), score.labels[nm])
      plots[[ct]] <- p
      if (save.individual.scatter) {
        file <- file.path(folder, "individual", paste0(gsub("[^A-Za-z0-9_-]+", "_", ct), ".", view, ".png"))
        save.png(p, file, panel.width, panel.height); pngs <- c(pngs, file)
      }
    }
    caption <- if (unit == "cell")
      "Spearman; unadjusted cell p; same-donor dependence, Country and technology are not adjusted. Line: descriptive linear fit." else
      "Spearman of same-profile donor means; unadjusted p; Country/age/technology not adjusted. Line: descriptive linear fit."
    overview <- wrap_plots(plots, ncol = n.columns, nrow = n.rows, guides = "collect") +
      plot_annotation(title = paste(score.labels[nm], "vs MFF", view, sep = " | "), caption = caption,
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 10, hjust = 0))) &
      theme(legend.position = "bottom")
    file <- file.path(folder, paste0(nm, ".MFF_correlation.", view, ".overview.png"))
    save.png(overview, file, panel.width * n.columns, panel.height * n.rows); pngs <- c(pngs, file)
  }
  expected <- length(plot.units) * length(scopes) * length(score.columns) * length(views) *
    (1L + if (save.individual.scatter) length(celltypes) else 0L)
  stopifnot(length(pngs) == expected, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  assign("hlma.correlations", list(profiles = cells, donor_means = donor, statistics = statistics), envir = .GlobalEnv)
  writeLines(c("MFF versus exactly two PREVIOUSLY CALCULATED module scores. No AddModuleScore or RDS read/write.",
    "X: stored RNA/data MFF. Y: SaulSenMayo or Hallmark IFN-alpha AddModuleScore from hlma.scores.",
    "MFF=0 and negative/zero scores are retained. No positive-only or High/Low filtering, normalization or pseudobulk.",
    "Views: All, Young, Old using original age_pop; scopes: combined, scRNA, snRNA.",
    "Correlation: two-sided Spearman rho/p; at least three observations and two distinct values on EACH axis.",
    "Exact permutation p only for fewer than 10 observations with no ties; otherwise asymptotic Spearman p.",
    "No/too few observations or constant axes yield NA rho/p with a reason; these panels remain in the figures.",
    "BH over both scores x all cell types separately for each scope x unit x view; figures display raw p.",
    "Default scatter unit is cell, 18 overview PNGs; optional plot.units=c('cell','donor') makes 36.",
    "Scatter keeps every profile; individual PNGs off by default. n is the observations actually used in that view.",
    "All/Young/Old share axis limits within the same scope/type/score/unit; scores are not recalculated by group.",
    "A straight regression line is descriptive, without CI; displayed correlation is rank-based Spearman.",
    "Cell p ignores within-donor dependence; use cell results as exploratory distributions, not donor-level evidence.",
    "Donor CSV uses ONE equally weighted observed donor per type/scope, averaging MFF and scores over identical profiles.",
    "Donor MFF is a mean of stored expression, not raw-count pseudobulk; no absent donor/type zero rows are invented.",
    "All-view correlations can reflect age, Country, technology or depth; separate age views do not remove these confounders.",
    "Correlations are unadjusted associations, not evidence that MFF causes either transcriptional programme.",
    "Combined and technical scopes overlap; All and age subsets overlap; these are not independent replications.",
    paste0("MFF exclusion in scoring provenance confirmed: ", exclusion.confirmed),
    paste0("Input global object: ", input.object), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_MFF_SCORE_CORRELATION_COMPLETE\nProfiles: ", nrow(cells), "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
