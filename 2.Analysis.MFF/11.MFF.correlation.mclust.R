# SERVER R: run ALL code in the SAME session after script 36 completes.
# Uses hlma.mclust.scores: four existing mclust subsets, MFF versus two existing scores.
# No Seurat/RDS reload, score calculation, mclust refit or RDS save.
local({
  packages <- c("ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

  # 1. Settings. Reruns overwrite these correlation outputs.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/mclust_MFF_correlation_SenMayo_IFNalpha")
  input.object <- "hlma.mclust.scores"
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  conditions <- c("all_high", "all_low", "positive_high", "positive_low")
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  views <- c("All", "Young", "Old")
  plot.units <- "cell"  # c("cell", "donor") also plots donor means; both units' CSVs always saved.
  save.individual.scatter <- FALSE
  show.regression.line <- TRUE
  group.order <- c("Young", "Old"); group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L; n.columns <- 5L
  annotation.size <- 8  # Size of rho / p / n text inside each panel.

  # 2. Verify the paired profiles and the original classification masks from script 36.
  prepare.correlations <- function(cache, columns) {
    modes <- c("all_high", "all_low", "positive_high", "positive_low")
    d <- cache$profiles
    required <- c("cell_id", "donor", "celltype", "tech", "Country", "group", "mff", unname(columns),
      "MFF_mclust_all", "MFF_mclust_all_cutoff", "MFF_mclust_positive", "MFF_mclust_positive_cutoff")
    if (!is.data.frame(d) || !nrow(d) || !all(required %in% names(d)))
      stop("Invalid cached profiles; run the full script 36 in this session first.")
    for (nm in c("cell_id", "donor", "celltype", "tech", "Country", "group"))
      if (!is.character(d[[nm]]) || anyNA(d[[nm]]) || any(!nzchar(trimws(d[[nm]])))) stop("Invalid profile metadata: ", nm)
    if (anyDuplicated(d$cell_id) || !all(d$tech %in% c("scRNA", "snRNA")) || !all(d$group %in% c("Young", "Old")))
      stop("Duplicate cell IDs or unknown technology/age group.")
    for (nm in c("mff", unname(columns)))
      if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]]))) stop("Nonfinite/nonnumeric paired value: ", nm)
    if (any(d$mff < 0)) stop("MFF must be nonnegative stored RNA/data expression.")
    if (anyDuplicated(unique(d[, c("donor", "group", "Country")])$donor)) stop("Inconsistent donor group/Country.")
    for (mode in c("all", "positive")) {
      column <- paste0("MFF_mclust_", mode); cut.column <- paste0(column, "_cutoff")
      label <- d[[column]]; cutoff <- d[[cut.column]]
      if (!is.character(label) || !is.numeric(cutoff) || any(!is.finite(cutoff) & !is.na(cutoff)) ||
          any(cutoff < 0, na.rm = TRUE)) stop("Invalid original labels/cutoff: ", mode)
      for (ct in unique(d$celltype)) if (length(unique(cutoff[d$celltype == ct])) != 1L)
        stop("Expected the original combined cell-type cutoff: ", ct, " / ", mode)
      eligible <- if (mode == "all") rep(TRUE, nrow(d)) else d$mff > 0
      classified <- eligible & is.finite(cutoff)
      if (any(is.na(label[classified])) || any(!label[classified] %in% c("High", "Low")) ||
          any(!is.na(label[eligible & !is.finite(cutoff)]))) stop("Original classification/cutoff mismatch: ", mode)
      if (mode == "positive" && any(is.na(label[!eligible]) | label[!eligible] != "Excluded_zero"))
        stop("Positive-mode MFF=0 must remain Excluded_zero.")
      delta <- d$mff[classified] - cutoff[classified]; tol <- 1e-10 * pmax(1, abs(cutoff[classified]))
      if (any((label[classified] == "High" & delta < -tol) | (label[classified] == "Low" & delta > tol)) ||
          any(label[classified] == "High" & d$mff[classified] == 0)) stop("Original classification disagrees with MFF: ", mode)
    }
    masks <- list(all_high = !is.na(d$MFF_mclust_all) & d$MFF_mclust_all == "High",
      all_low = !is.na(d$MFF_mclust_all) & d$MFF_mclust_all == "Low",
      positive_high = !is.na(d$MFF_mclust_positive) & d$MFF_mclust_positive == "High",
      positive_low = !is.na(d$MFF_mclust_positive) & d$MFF_mclust_positive == "Low")
    if (!is.list(cache$masks) || !identical(names(cache$masks), modes)) stop("Missing original subset masks; rerun script 36.")
    for (nm in modes) {
      v <- cache$masks[[nm]]
      if (!is.logical(v) || anyNA(v) || !identical(unname(v), masks[[nm]]) ||
          (!is.null(names(v)) && !identical(names(v), d$cell_id))) stop("Cached mask/profile order disagrees: ", nm)
    }
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
  correlation.tables <- function(cells, masks, columns, types, scope.names, view.names) {
    if (any(vapply(masks, function(v) !is.logical(v) || length(v) != nrow(cells) || anyNA(v), logical(1)))) stop("Invalid masks.")
    donor <- data.frame(condition = character(), scope = character(), celltype = character(), donor = character(),
      group = character(), Country = character(), n_profiles = integer(), n_scRNA = integer(),
      n_snRNA = integer(), mff = numeric(), stringsAsFactors = FALSE)
    for (nm in unname(columns)) donor[[nm]] <- numeric()
    for (condition in names(masks)) for (scope in scope.names) {
      b <- cells[masks[[condition]] & (scope == "combined" | cells$tech == scope), , drop = FALSE]
      for (ct in types) {
        d <- b[b$celltype == ct, , drop = FALSE]
        if (nrow(d)) for (idx in split(seq_len(nrow(d)), d$donor)) {
          z <- d[idx, , drop = FALSE]
          row <- data.frame(condition = condition, scope = scope, celltype = ct, donor = z$donor[1], group = z$group[1],
            Country = z$Country[1], n_profiles = nrow(z), n_scRNA = sum(z$tech == "scRNA"),
            n_snRNA = sum(z$tech == "snRNA"), mff = mean(z$mff))
          for (nm in unname(columns)) row[[nm]] <- mean(z[[nm]])
          donor <- rbind(donor, row)
        }
      }
    }
    statistics <- list()
    for (condition in names(masks)) for (scope in scope.names) for (ct in types) for (nm in names(columns)) for (view in view.names) {
      d <- cells[masks[[condition]] & cells$celltype == ct & (scope == "combined" | cells$tech == scope), , drop = FALSE]
      ds <- donor[donor$condition == condition & donor$scope == scope & donor$celltype == ct, , drop = FALSE]
      if (view != "All") { d <- d[d$group == view, , drop = FALSE]; ds <- ds[ds$group == view, , drop = FALSE] }
      for (unit in c("cell", "donor")) {
        z <- if (unit == "cell") d else ds
        r <- spearman.pair(z$mff, z[[columns[[nm]]]])
        statistics[[length(statistics) + 1L]] <- cbind(data.frame(condition = condition, scope = scope, celltype = ct,
          score = nm, view = view, unit = unit, n_profiles = nrow(d), n_donors = nrow(ds),
          n_profiles_MFF_positive = sum(d$mff > 0)), r)
      }
    }
    s <- do.call(rbind, statistics); s$p_adj_BH <- NA_real_
    for (condition in names(masks)) for (scope in scope.names) for (unit in c("cell", "donor")) for (view in view.names) {
      idx <- which(s$condition == condition & s$scope == scope & s$unit == unit & s$view == view)
      s$p_adj_BH[idx] <- p.adjust(s$p_value[idx], method = "BH")
    }
    rownames(donor) <- rownames(s) <- NULL
    list(donor = donor, statistics = s)
  }
  write.checked.csv <- function(d, file) {
    token <- ".__HLMA_MCLUST_COR_MISSING__."
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

  # 3. Reuse script 36's paired values and classifications without running that script again.
  if (!exists(input.object, envir = .GlobalEnv, inherits = FALSE))
    stop("Run the full script 36 in this R session first: missing ", input.object)
  cache <- get(input.object, envir = .GlobalEnv)
  prepared <- prepare.correlations(cache, score.columns)
  cells <- prepared$cells; masks <- prepared$masks
  if (!all(plot.units %in% c("cell", "donor")) || !length(plot.units) || anyDuplicated(plot.units)) stop("Invalid plot.units.")
  if (!is.numeric(annotation.size) || length(annotation.size) != 1L || !is.finite(annotation.size) || annotation.size <= 0)
    stop("annotation.size must be a finite positive number.")
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE)); n.rows <- ceiling(length(celltypes) / n.columns)
  cat("Reusing ", input.object, ": ", nrow(cells), " reference profiles; computing four-subset MFF-score Spearman correlations.\n", sep = "")
  result <- correlation.tables(cells, masks, score.columns, celltypes, scopes, views)
  statistics <- result$statistics; donor <- result$donor
  for (condition in conditions) stopifnot(sum(donor$n_profiles[donor$condition == condition]) == 2L * sum(masks[[condition]]))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(statistics[statistics$unit == "cell", ], file.path(out, "01.mclust_MFF_score.cell_correlations.csv"))
  write.checked.csv(statistics[statistics$unit == "donor", ], file.path(out, "02.mclust_MFF_score.donor_correlations.csv"))
  write.checked.csv(donor, file.path(out, "03.mclust_MFF_score.donor_means.csv"))
  write.checked.csv(cells, file.path(out, "04.mclust_MFF_score.profile_pairs_labels.csv.gz"))

  # 4. Scatter panels: signed scores, common All/Young/Old axes, rho / p / n at size 8.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
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
      annotate("text", x = Inf, y = bounds[2] + .95 * span, label = text,
        hjust = 1.05, vjust = 1, size = annotation.size) +
      coord_cartesian(xlim = xb, ylim = c(bounds[1] - .08 * span, bounds[2] + span)) +
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
  for (condition in conditions) for (unit in plot.units) for (scope in scopes) for (nm in names(score.columns)) for (view in views) {
    b <- if (unit == "cell") cells[masks[[condition]] & (scope == "combined" | cells$tech == scope), , drop = FALSE] else
      donor[donor$condition == condition & donor$scope == scope, , drop = FALSE]
    folder <- file.path(out, condition, scope, unit, nm)
    plots <- list()
    for (ct in celltypes) {
      full <- b[b$celltype == ct, , drop = FALSE]; full$score <- full[[score.columns[[nm]]]]
      stat <- statistics[statistics$condition == condition & statistics$scope == scope & statistics$unit == unit &
        statistics$score == nm & statistics$celltype == ct & statistics$view == view, ]
      stopifnot(nrow(stat) == 1L)
      p <- scatter.plot(full, stat, unit, view, paste(scope.labels[scope], ct, sep = " | "), score.labels[nm])
      plots[[ct]] <- p
      if (save.individual.scatter) {
        file <- file.path(folder, "individual", paste0(gsub("[^A-Za-z0-9_-]+", "_", ct), ".", view, ".png"))
        save.png(p, file, panel.width, panel.height); pngs <- c(pngs, file)
      }
    }
    caption <- if (unit == "cell")
      "Conditional Spearman after MFF selection; unadjusted cell p; donor, Country and technology not adjusted. Line: descriptive linear fit." else
      "Conditional Spearman of same-subset donor means; unadjusted p; Country/age/technology not adjusted. Line: descriptive linear fit."
    overview <- (wrap_plots(plots, ncol = n.columns, nrow = n.rows, guides = "collect") +
      plot_annotation(title = paste(score.labels[nm], "vs MFF", condition, view, sep = " | "), caption = caption,
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 10, hjust = 0)))) &
      theme(legend.position = "bottom")
    file <- file.path(folder, paste0(nm, ".", condition, ".MFF_correlation.", view, ".overview.png"))
    save.png(overview, file, panel.width * n.columns, panel.height * n.rows); pngs <- c(pngs, file)
  }
  expected <- length(conditions) * length(plot.units) * length(scopes) * length(score.columns) * length(views) *
    (1L + if (save.individual.scatter) length(celltypes) else 0L)
  stopifnot(length(pngs) == expected, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  assign("hlma.mclust.correlations", list(profiles = cells, masks = masks, donor_means = donor, statistics = statistics), envir = .GlobalEnv)
  writeLines(c("MFF versus exactly two previously calculated scores within four ORIGINAL mclust subsets from script 36.",
    "all_high/all_low use original all-mode labels; positive_high/positive_low use original positive-mode labels.",
    "Original combined cell-type labels/cutoffs are reused for all scopes. Neither mclust nor scores are recalculated.",
    "No RDS read/write, normalization, or pseudobulk. Source cache is left unchanged.",
    "Classified all-mode MFF=0 profiles remain in all_low; positive subsets exclude zeros. Unavailable labels remain NA.",
    "Negative/zero scores are retained. X is stored RNA/data MFF; Y is the corresponding full-reference AddModuleScore.",
    "Views: All, Young, Old using original age_pop; scopes: combined, scRNA, snRNA.",
    "Two-sided Spearman rho/p; at least three observations and two distinct values on EACH axis.",
    "Exact permutation p only for fewer than 10 observations without ties; otherwise asymptotic Spearman p.",
    "Empty/rare/constant inputs remain as NA rho/p panels with a reason in CSV. n is the observations in that view.",
    "BH over both scores x all cell types separately for each condition x scope x unit x view; figures show raw p.",
    "Default: 72 cell overview PNGs = 4 subsets x 3 scopes x 2 scores x 3 views.",
    "5 columns x 3 rows for 15 types, 4500 x 2550 px at 160dpi. rho/p/n annotation.size=8.",
    "Optional plot.units=c('cell','donor') gives 144 overview PNGs; individual scatter PNGs default off.",
    "All/Young/Old share axes within each condition/scope/type/score/unit. Different subsets may have different axes.",
    "The straight line is a descriptive LM fit without CI; rho and p come from Spearman ranks.",
    "Donor means use only selected profiles, averaging MFF and each score over the IDENTICAL profiles.",
    "Observed donors have equal weight; no missing donor/type score or expression zero is invented.",
    "Cell-level p ignores same-donor dependence, Country, technology and depth; these are exploratory conditional associations.",
    "Selecting on MFF restricts its range and can change correlation; subset rho is not the full-population association.",
    "Donor correlations and separate age views also do not adjust Country/technology or other confounders.",
    "All/positive, combined/technical and All/age analyses overlap and are not independent replications.",
    "The profile-pair CSV contains the FULL reference with original labels, including excluded/unclassified profiles.",
    "CSV/gzip outputs are readback-verified. Cell/donor statistical tables and donor means are always saved.",
    paste0("Input global cache: ", input.object), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_MCLUST_MFF_SCORE_CORRELATION_COMPLETE\nReference profiles: ", nrow(cells), "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
