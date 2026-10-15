# Paste ALL code in SERVER R, or source a server copy.
# Two mclust inputs per cell type: exact zeros included vs MFF > 0 only.
# Age-blind cutoffs; four conditional violin/FC sets plus paired Lower/Upper violins.
# Existing plots and original RDS stay unchanged; six additional overviews are saved.
# If needed, run install.packages("mclust") once before running this script.
local({
  needed <- c("Seurat", "mclust", "ggplot2", "patchwork", "Cairo")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(Seurat); library(ggplot2); library(patchwork)
  })

  # 1. Settings. Edit G.candidates <- 2L to FORCE two components.
  gene <- "MFF"
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  rds.file <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  source.h5ad <- "/BiO/Store/UNIST-Aging-MFF-2026-09/HLMA_scRNA/OMIX004308-03.h5ad"
  out <- file.path(project, "analysis/3.Mffexpression/mclust_cutoffs")
  G.candidates <- 1:2             # Compare 1 vs 2 by BIC; G=1 gives cutoff NA.
  model.names <- c("E", "V")      # Equal vs unequal univariate Gaussian variances.
  min.fit.n <- 20L                # Computational guard, not a biological QC threshold.
  seed <- 20261006L
  expression.scale <- "auto"      # auto, ln1p, log2p1, linear; clustering uses stored values.
  allow.assumed.ln1p <- TRUE      # Keep earlier FC assumption explicitly marked.
  fc.pseudocount <- 1
  min.cells.per.group <- 2L
  save.classifications <- TRUE
  group.order <- c("Young", "Old")
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  region.order <- c("Low", "High")
  region.labels <- c(Low = "Lower", High = "Upper")
  region.colors <- c(Low = "#2878BD", High = "#C54045")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  conditions <- c(all_high = "0 included fit | MFF > T_all",
    all_low = "0 included fit | MFF <= T_all",
    positive_high = "0 excluded fit | MFF > T_positive",
    positive_low = "0 excluded fit | 0 < MFF <= T_positive")
  n.columns <- 5L
  panel.width <- 900L
  panel.height <- 850L
  plot.res <- 160L
  stopifnot(length(G.candidates) > 0L, all(G.candidates %in% 1:2),
    length(model.names) > 0L, all(model.names %in% c("E", "V")),
    min.fit.n >= 3L, min.cells.per.group >= 2L, n.columns >= 1L,
    is.finite(fc.pseudocount), fc.pseudocount > 0)

  # 2. Pure helpers: no jitter, imputation, rounding or random subsampling.
  fit.input <- function(v, mode) {
    if (!is.numeric(v) || any(!is.finite(v)) || any(v < 0)) stop("MFF must be finite and nonnegative.")
    if (!mode %in% c("all", "positive")) stop("Unknown input mode.")
    if (mode == "all") v else v[v > 0]
  }
  posterior.high <- function(v, mu, variance, pro) {
    ord <- order(mu)
    mu <- mu[ord]; variance <- rep(variance, length.out = 2L)[ord]; pro <- pro[ord]
    high <- log(pro[2]) + dnorm(v, mu[2], sqrt(variance[2]), log = TRUE)
    low <- log(pro[1]) + dnorm(v, mu[1], sqrt(variance[1]), log = TRUE)
    plogis(high - low)
  }
  single.boundary <- function(mu, variance, pro, support) {
    result <- list(cutoff = NA_real_, status = "invalid_parameters", n_roots = 0L)
    if (length(mu) != 2L || !length(variance) %in% c(1L, 2L) || length(pro) != 2L ||
        length(support) != 2L || any(!is.finite(c(mu, variance, pro, support))) ||
        any(pro <= 0) || support[1] >= support[2]) return(result)
    ord <- order(mu); mu <- mu[ord]
    variance <- rep(variance, length.out = 2L)[ord]; pro <- pro[ord] / sum(pro)
    floor <- .Machine$double.eps * max(1, diff(support)^2)
    if (any(variance <= floor)) { result$status <- "degenerate_variance"; return(result) }
    if (diff(mu) <= sqrt(.Machine$double.eps) * max(1, abs(mu))) {
      result$status <- "indistinguishable_means"; return(result)
    }
    # Log weighted density(high) - log weighted density(low) = a*x^2+b*x+c.
    a <- 1 / (2 * variance[1]) - 1 / (2 * variance[2])
    b <- mu[2] / variance[2] - mu[1] / variance[1]
    cc <- log(pro[2] / pro[1]) + 0.5 * log(variance[1] / variance[2]) +
      mu[1]^2 / (2 * variance[1]) - mu[2]^2 / (2 * variance[2])
    if (a == 0) roots <- if (b == 0) numeric() else -cc / b else {
      z <- polyroot(c(cc, b, a))
      roots <- Re(z[abs(Im(z)) <= 1e-8 * pmax(1, abs(Re(z)))])
    }
    roots <- sort(unique(roots[is.finite(roots) & roots > support[1] & roots < support[2]]))
    result$n_roots <- length(roots)
    if (!length(roots)) { result$status <- "no_crossing_in_support"; return(result) }
    if (length(roots) != 1L) { result$status <- "multiple_crossings"; return(result) }
    endpoint.p <- posterior.high(support, mu, variance, pro)
    if (any(!is.finite(endpoint.p)) || endpoint.p[1] >= 0.5 || endpoint.p[2] <= 0.5) {
      result$status <- "not_low_to_high"; return(result)
    }
    if (abs(posterior.high(roots, mu, variance, pro) - 0.5) > 1e-5) {
      result$status <- "boundary_numerical_failure"; return(result)
    }
    result$cutoff <- roots[1]; result$status <- "ok"
    result
  }
  evaluate.fit <- function(fit, input, mode) {
    r <- data.frame(selected_G = NA_integer_, model_name = "", selected_return_code = NA_integer_, BIC_best = NA_real_,
      BIC_G1 = NA_real_, BIC_G2 = NA_real_, BIC_delta_G2_minus_G1 = NA_real_,
      cutoff = NA_real_, status = "fit_failed", n_roots_in_support = 0L,
      mean_low = NA_real_, mean_high = NA_real_, sd_low = NA_real_, sd_high = NA_real_,
      weight_low = NA_real_, weight_high = NA_real_,
      fraction_max_p_below_0_8 = NA_real_, stringsAsFactors = FALSE)
    if (is.null(fit) || is.null(fit$G)) return(r)
    r$selected_G <- as.integer(fit$G); r$model_name <- as.character(fit$modelName)
    r$BIC_best <- as.numeric(fit$bic)
    if (is.matrix(fit$BIC)) for (g in c("1", "2")) {
      if (g %in% rownames(fit$BIC)) {
        v <- as.numeric(fit$BIC[g, ]); v <- v[is.finite(v)]
        if (length(v)) r[[paste0("BIC_G", g)]] <- max(v)
      }
    }
    r$BIC_delta_G2_minus_G1 <- r$BIC_G2 - r$BIC_G1
    codes <- attr(fit$BIC, "returnCodes")
    if (is.matrix(codes) && as.character(r$selected_G) %in% rownames(codes) &&
        r$model_name %in% colnames(codes))
      r$selected_return_code <- as.integer(codes[as.character(r$selected_G), r$model_name])
    if (r$selected_G != 2L) { r$status <- "one_component_selected"; return(r) }
    if (!is.finite(r$selected_return_code)) { r$status <- "selected_convergence_unknown"; return(r) }
    if (r$selected_return_code != 0L) { r$status <- "selected_fit_not_converged"; return(r) }
    mu <- as.numeric(fit$parameters$mean)
    variance <- as.numeric(fit$parameters$variance$sigmasq)
    pro <- as.numeric(fit$parameters$pro)
    if (!length(pro)) pro <- c(0.5, 0.5)
    boundary <- single.boundary(mu, variance, pro, range(input))
    r$cutoff <- boundary$cutoff; r$status <- boundary$status
    r$n_roots_in_support <- boundary$n_roots
    if (length(mu) == 2L && length(variance) %in% c(1L, 2L) && length(pro) == 2L) {
      ord <- order(mu); variance <- rep(variance, length.out = 2L)[ord]
      r$mean_low <- mu[ord][1]; r$mean_high <- mu[ord][2]
      r$sd_low <- sqrt(pmax(0, variance[1])); r$sd_high <- sqrt(pmax(0, variance[2]))
      r$weight_low <- pro[ord][1]; r$weight_high <- pro[ord][2]
      if (all(is.finite(c(mu, variance, pro))) && all(variance > 0) && all(pro > 0)) {
        p <- posterior.high(input, mu, as.numeric(fit$parameters$variance$sigmasq), pro)
        r$fraction_max_p_below_0_8 <- mean(pmax(p, 1 - p) < 0.8)
      }
    }
    r
  }
  classify.cutoff <- function(v, cutoff, mode) {
    fit.input(v, mode) # Validate without changing v/order.
    region <- rep(NA_character_, length(v))
    eligible <- if (mode == "all") rep(TRUE, length(v)) else v > 0
    if (mode == "positive") region[!eligible] <- "Excluded_zero"
    if (length(cutoff) == 1L && is.finite(cutoff))
      region[eligible] <- ifelse(v[eligible] > cutoff, "High", "Low")
    region
  }
  test.subset <- function(d, cutoff.available, min.n, pseudocount) {
    y <- d$expression[d$group == "Young"]; o <- d$expression[d$group == "Old"]
    yn <- d$normalized[d$group == "Young"]; on <- d$normalized[d$group == "Old"]
    r <- data.frame(n_Young = length(y), n_Old = length(o),
      n_donors_Young = length(unique(d$donor[d$group == "Young"])),
      n_donors_Old = length(unique(d$donor[d$group == "Old"])),
      mean_stored_Young = if (length(y)) mean(y) else NA_real_,
      mean_stored_Old = if (length(o)) mean(o) else NA_real_,
      mean_normalized_Young = if (length(yn)) mean(yn) else NA_real_,
      mean_normalized_Old = if (length(on)) mean(on) else NA_real_,
      log2FC_Old_vs_Young = NA_real_, W = NA_real_, p_value = NA_real_, reason = "",
      stringsAsFactors = FALSE)
    if (!cutoff.available) { r$reason <- "cutoff_unavailable"; return(r) }
    if (length(yn) && length(on)) r$log2FC_Old_vs_Young <-
      log2(mean(on) + pseudocount) - log2(mean(yn) + pseudocount)
    if (!length(y) || !length(o)) r$reason <- "missing_age_group" else
      if (min(length(y), length(o)) < min.n) r$reason <- "insufficient_profiles" else
        if (length(unique(c(y, o))) < 2L) r$reason <- "identical_values"
    if (nzchar(r$reason)) return(r)
    fit <- tryCatch(stats::wilcox.test(y, o, exact = FALSE, correct = TRUE,
      alternative = "two.sided", paired = FALSE), error = function(e) e)
    if (inherits(fit, "error")) r$reason <- conditionMessage(fit) else
      if (!is.finite(fit$p.value)) r$reason <- "nonfinite_p" else {
        r$p_value <- fit$p.value; r$W <- unname(fit$statistic)
      }
    r
  }
  write.checked.csv <- function(z, path) {
    z <- as.data.frame(lapply(z, function(v) if (is.factor(v)) as.character(v) else v))
    rownames(z) <- NULL
    token <- ".__MFF_MISSING__."
    strings <- unlist(z[vapply(z, is.character, logical(1))], use.names = FALSE)
    while (any(strings == token, na.rm = TRUE)) token <- paste0(token, "_")
    opener <- if (grepl("\\.gz$", path)) gzfile else file
    con <- opener(path, "wt", encoding = "UTF-8")
    tryCatch(write.csv(z, con, row.names = FALSE, na = token), finally = close(con))
    con <- opener(path, "rt", encoding = "UTF-8")
    saved <- tryCatch(read.csv(con, check.names = FALSE, na.strings = token,
      colClasses = vapply(z, function(v) class(v)[1], character(1))), finally = close(con))
    stopifnot(identical(names(saved), names(z)), nrow(saved) == nrow(z))
    for (field in names(z)) {
      if (is.numeric(z[[field]])) stopifnot(isTRUE(all.equal(saved[[field]], z[[field]],
        check.attributes = FALSE, tolerance = 1e-8))) else stopifnot(identical(saved[[field]], z[[field]]))
    }
    invisible(path)
  }

  # Additive Lower/Upper views reuse the saved cutoff classifications, pooling Young+Old.
  prepare.lower.upper <- function(b, fits, scope, mode, ct) {
    if (!mode %in% c("all", "positive")) stop("Unknown Lower/Upper input mode.")
    column <- paste0("region_", mode)
    required <- c("cell_id", "donor", "celltype", "group", "expression", column)
    if (!is.data.frame(b) || !all(required %in% names(b)) ||
        !is.numeric(b$expression) || any(!is.finite(b$expression)) || any(b$expression < 0) ||
        anyNA(b$cell_id) || anyDuplicated(b$cell_id) || any(b$cell_id == "") ||
        anyNA(b$celltype) || anyNA(b$group) || !all(b$group %in% c("Young", "Old"))) {
      stop("Invalid profiles for Lower/Upper violin.")
    }
    if (!is.data.frame(fits) || !all(c("scope", "celltype", "input_mode", "cutoff", "status") %in% names(fits)))
      stop("Missing Lower/Upper cutoff table fields.")
    s <- fits[which(fits$scope == scope & fits$celltype == ct & fits$input_mode == mode), , drop = FALSE]
    if (nrow(s) != 1L) stop("Expected one fit row for ", scope, "/", ct, "/", mode)
    part <- b[which(b$celltype == ct), , drop = FALSE]
    eligible <- if (mode == "all") rep(TRUE, nrow(part)) else part$expression > 0
    n.input <- sum(eligible)
    if ("n_input" %in% names(s) && (!is.finite(s$n_input) || s$n_input != n.input))
      stop("Profile count does not match the original fit input: ", scope, "/", ct, "/", mode)
    available <- is.numeric(s$cutoff) && is.finite(s$cutoff) && identical(as.character(s$status), "ok")
    raw.region <- as.character(part[[column]])
    if (available && (anyNA(raw.region[eligible]) ||
        any(!raw.region[eligible] %in% c("Low", "High")) ||
        (mode == "positive" && any(is.na(raw.region[!eligible]) | raw.region[!eligible] != "Excluded_zero")))) {
      stop("Classification/zero inclusion does not match the original fit input.")
    }
    d <- if (available) part[which(eligible & !is.na(raw.region) & raw.region %in% c("Low", "High")), , drop = FALSE] else
      part[FALSE, , drop = FALSE]
    d$region <- factor(as.character(d[[column]]), levels = c("Low", "High"))
    regular <- d[FALSE, , drop = FALSE]
    special <- data.frame(region = character(), value = numeric(), n = integer(), stringsAsFactors = FALSE)
    counts <- list()
    for (g in c("Low", "High")) {
      dd <- d[which(d$region == g), , drop = FALSE]
      if (nrow(dd) >= 2L && length(unique(dd$expression)) > 1L) regular <- rbind(regular, dd) else
        if (nrow(dd)) special <- rbind(special, data.frame(region = g, value = dd$expression[1L], n = nrow(dd)))
      counts[[g]] <- data.frame(scope = scope, celltype = ct, input_mode = mode,
        region = g, display_group = if (g == "Low") "Lower" else "Upper",
        cutoff = s$cutoff, cutoff_status = as.character(s$status),
        n_input = as.integer(n.input), n_zero_input = as.integer(sum(part$expression[eligible] == 0)),
        n_excluded_zero = if (mode == "positive") as.integer(sum(part$expression == 0)) else 0L,
        n_profiles = if (available) nrow(dd) else NA_integer_,
        n_zero = if (available) as.integer(sum(dd$expression == 0)) else NA_integer_,
        n_donors = if (available) length(unique(dd$donor)) else NA_integer_,
        n_Young = if (available) as.integer(sum(dd$group == "Young")) else NA_integer_,
        n_Old = if (available) as.integer(sum(dd$group == "Old")) else NA_integer_,
        mean_stored_expression = if (nrow(dd)) mean(dd$expression) else NA_real_,
        median_stored_expression = if (nrow(dd)) stats::median(dd$expression) else NA_real_,
        stringsAsFactors = FALSE)
    }
    list(data = d, regular = regular, special = special, counts = do.call(rbind, counts),
      available = available, cutoff = s$cutoff, status = as.character(s$status),
      n_input = n.input, celltype = ct, mode = mode)
  }

  lower.upper.violin <- function(panel) {
    d <- panel$data; regular <- panel$regular; special <- panel$special
    top <- max(c(0.2, d$expression, if (panel$available) panel$cutoff else 0))
    labels <- setNames(paste0(region.labels[region.order], "\n(n=",
      ifelse(is.na(panel$counts$n_profiles), "NA", format(panel$counts$n_profiles, big.mark = ",", trim = TRUE)), ")"), region.order)
    p <- ggplot2::ggplot(d, ggplot2::aes(region, expression)) + ggplot2::geom_blank() +
      ggplot2::geom_violin(data = regular, ggplot2::aes(fill = region),
        trim = TRUE, scale = "width", width = 0.8) +
      ggplot2::geom_boxplot(data = regular, width = 0.18, fill = "white", alpha = 0.85, outlier.shape = NA)
    if (nrow(special)) p <- p + ggplot2::geom_errorbar(data = special,
      ggplot2::aes(region, ymin = value, ymax = value, colour = region), inherit.aes = FALSE, width = 0.45)
    if (any(special$n == 1L)) p <- p + ggplot2::geom_point(data = special[special$n == 1L, , drop = FALSE],
      ggplot2::aes(region, value, colour = region), inherit.aes = FALSE, size = 2.5)
    if (panel$available) p <- p + ggplot2::geom_hline(yintercept = panel$cutoff,
      linetype = "dashed", colour = "grey35", linewidth = 0.6)
    p <- p + ggplot2::scale_fill_manual(values = region.colors, guide = "none") +
      ggplot2::scale_colour_manual(values = region.colors, guide = "none") +
      ggplot2::scale_x_discrete(limits = region.order, drop = FALSE, labels = labels) +
      ggplot2::coord_cartesian(ylim = c(0, top * 1.15)) +
      ggplot2::labs(title = panel$celltype,
        subtitle = if (panel$available) paste0("T = ", signif(panel$cutoff, 4), " | Young + Old pooled") else panel$status,
        x = NULL, y = "MFF expression (stored RNA/data)") + ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face = "bold", size = 20),
        axis.text.x = ggplot2::element_text(size = 17), panel.grid.minor = ggplot2::element_blank())
    if (!nrow(d)) p <- p + ggplot2::annotate("text", x = 1.5, y = top * 0.5,
      label = if (panel$available) "No classified profiles" else "Cutoff unavailable", size = 4)
    p
  }

  # 3. Input and FC scale; only a single expression vector is extracted.
  x <- if (exists("hlma.umap", envir = .GlobalEnv, inherits = FALSE))
    get("hlma.umap", envir = .GlobalEnv) else {
      if (!file.exists(rds.file)) stop("RDS not found: ", rds.file)
      readRDS(rds.file)
    }
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  meta <- x[[]]
  required <- c("Annotation", "sample", "tech", "age_pop", "age", "Country")
  stopifnot(all(required %in% names(meta)), identical(rownames(meta), colnames(x)),
    !anyDuplicated(colnames(x)))
  for (nm in required) if (anyNA(meta[[nm]]) || any(trimws(as.character(meta[[nm]])) == ""))
    stop("Missing metadata: ", nm)
  stopifnot(all(as.character(meta$tech) %in% c("scRNA", "snRNA")),
    all(as.character(meta$age_pop) %in% c("young_pop", "old_pop")))
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(sum(rownames(rna) == gene) == 1L, identical(colnames(rna), colnames(x)))
  cells <- data.frame(cell_id = colnames(x), donor = as.character(meta$sample),
    celltype = as.character(meta$Annotation), tech = as.character(meta$tech),
    group = ifelse(meta$age_pop == "young_pop", "Young", "Old"),
    age = suppressWarnings(as.numeric(as.character(meta$age))), Country = as.character(meta$Country),
    expression = as.numeric(rna[gene, ]), stringsAsFactors = FALSE)
  rm(rna)
  fit.input(cells$expression, "all")
  stopifnot(all(is.finite(cells$age)))
  donor.info <- unique(cells[, c("donor", "group", "age", "Country")])
  if (anyDuplicated(donor.info$donor)) stop("Inconsistent donor age/group/Country.")
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE))
  n.rows <- ceiling(length(celltypes) / n.columns)
  import <- x@misc$HLMA_import
  if (!is.null(import$source_h5ad)) source.h5ad <- import$source_h5ad
  read.log.base <- function(path) {
    absent <- list(base = NA_real_, evidence = "no_usable_h5ad_log1p_metadata")
    if (!file.exists(path) || !requireNamespace("hdf5r", quietly = TRUE)) return(absent)
    h <- hdf5r::H5File$new(path, mode = "r"); on.exit(h$close_all(), add = TRUE)
    if (!h$exists("uns/log1p/base")) return(absent)
    b <- h[["uns/log1p/base"]]
    encoding <- tryCatch(hdf5r::h5attr(b, "encoding-type"), error = function(e) "")
    if (identical(as.character(encoding), "null")) return(list(base = exp(1), evidence = "h5ad_log1p_base_null"))
    value <- suppressWarnings(as.numeric(b$read()))
    if (length(value) == 1L && is.finite(value) && value > 1)
      list(base = value, evidence = "h5ad_log1p_numeric_base") else absent
  }
  stopifnot(expression.scale %in% c("auto", "ln1p", "log2p1", "linear"))
  log.info <- if (expression.scale == "auto") tryCatch(read.log.base(source.h5ad),
    error = function(e) list(base = NA_real_, evidence = conditionMessage(e))) else
      list(base = switch(expression.scale, ln1p = exp(1), log2p1 = 2, linear = NA_real_),
        evidence = paste0("manual_setting_", expression.scale))
  scale.assumed <- expression.scale != "auto"
  if (expression.scale == "auto" && !is.finite(log.info$base)) {
    if (!allow.assumed.ln1p) stop("Unconfirmed log scale; set documented expression.scale.")
    log.info$base <- exp(1); scale.assumed <- TRUE
    log.info$evidence <- paste0(log.info$evidence, "; assumed_ln1p")
  }
  cells$normalized <- if (expression.scale == "linear") cells$expression else
    expm1(cells$expression * log(log.info$base))
  stopifnot(all(is.finite(cells$normalized)), all(cells$normalized >= 0))
  scale.note <- if (scale.assumed) paste("FC scale:",
    if (expression.scale == "auto") "ln1p assumed" else paste(expression.scale, "manual")) else
      paste("FC log base:", signif(log.info$base, 5))
  cat(scale.note, "\n", log.info$evidence, "\n")
  warning("Threshold-selected distributions are descriptive. Cell-level p is not donor-adjusted; cutoffs do not correct detection/QC differences.", call. = FALSE)
  if (scale.assumed) warning("LOG2FC_SCALE_ASSUMED: ", scale.note, call. = FALSE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)

  # 4. Fit each scope/type twice. G/model choice uses no age label.
  fit.mff <- function(v, mode, fit.seed) {
    input <- fit.input(v, mode)
    info <- data.frame(input_mode = mode, n_total = length(v), n_zero_total = sum(v == 0),
      n_positive_total = sum(v > 0), n_input = length(input),
      n_distinct_input = length(unique(input)), zero_fraction_total = if (length(v)) mean(v == 0) else NA_real_,
      input_min = if (length(input)) min(input) else NA_real_,
      input_max = if (length(input)) max(input) else NA_real_,
      G_candidates = paste(G.candidates, collapse = ","), model_candidates = paste(model.names, collapse = ","),
      seed = fit.seed, fit_warnings = "", stringsAsFactors = FALSE)
    model <- NULL; early <- ""
    if (!length(input)) early <- "no_input_profiles" else
      if (length(input) < min.fit.n) early <- "insufficient_input_profiles" else
        if (length(unique(input)) < 3L) early <- "fewer_than_three_distinct_values"
    messages <- character()
    if (!nzchar(early)) {
      set.seed(fit.seed)
      model <- tryCatch(withCallingHandlers(mclust::Mclust(input,
        G = G.candidates, modelNames = model.names,
        initialization = list(subset = seq_along(input)),
        control = mclust::emControl(itmax = c(1000L, 1000L)), warn = TRUE, verbose = FALSE),
        warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error = function(e) { messages <<- c(messages, paste0("ERROR: ", conditionMessage(e))); NULL })
    }
    result <- evaluate.fit(model, input, mode)
    if (nzchar(early)) result$status <- early
    info$fit_warnings <- paste(unique(messages), collapse = " | ")
    list(summary = cbind(info, result), model = model,
      region = classify.cutoff(v, result$cutoff, mode))
  }
  cutoff.rows <- stats.rows <- coverage.rows <- list()
  models <- assignments <- list()
  counter <- 0L
  for (scope in names(scope.labels)) {
    b <- if (scope == "combined") cells else cells[cells$tech == scope, ]
    assignments[[scope]] <- b[, c("cell_id", "donor", "celltype", "tech", "group", "expression")]
    for (mode in c("all", "positive")) {
      assignments[[scope]][[paste0("cutoff_", mode)]] <- rep(NA_real_, nrow(b))
      assignments[[scope]][[paste0("region_", mode)]] <- rep(NA_character_, nrow(b))
    }
    models[[scope]] <- list()
    for (ct in celltypes) {
      idx <- which(b$celltype == ct); d <- b[idx, ]
      models[[scope]][[ct]] <- list()
      for (mode in c("all", "positive")) {
        counter <- counter + 1L
        cat(sprintf("FIT %d/%d | %s | %s | %s\n", counter, 6L * length(celltypes), scope, ct, mode))
        fit <- fit.mff(d$expression, mode, seed + counter)
        row <- fit$summary; row$scope <- scope; row$celltype <- ct
        cutoff.rows[[counter]] <- row
        models[[scope]][[ct]][[mode]] <- fit$model
        assignments[[scope]][[paste0("cutoff_", mode)]][idx] <- row$cutoff
        assignments[[scope]][[paste0("region_", mode)]][idx] <- fit$region
        available <- is.finite(row$cutoff)
        if (available) stopifnot(sum(fit$region %in% c("Low", "High")) == row$n_input)
        for (region in c("High", "Low")) {
          condition <- paste(mode, tolower(region), sep = "_")
          use <- which(!is.na(fit$region) & fit$region == region)
          part <- d[use, ]
          r <- test.subset(part, available, min.cells.per.group, fc.pseudocount)
          r$scope <- scope; r$celltype <- ct; r$condition <- condition
          r$input_mode <- mode; r$cutoff <- row$cutoff; r$cutoff_status <- row$status
          r$fc_scale_note <- scale.note; r$fc_log_base <- log.info$base
          r$fc_scale_assumed <- scale.assumed; r$fc_scale_evidence <- log.info$evidence
          r$fc_pseudocount <- fc.pseudocount
          stats.rows[[length(stats.rows) + 1L]] <- r
          # One coverage row per observed donor/type; unavailable thresholds stay NA.
          for (donor in unique(d$donor)) {
            dd <- d[d$donor == donor, ]; ss <- part[part$donor == donor, ]
            coverage.rows[[length(coverage.rows) + 1L]] <- data.frame(scope = scope,
              celltype = ct, condition = condition, donor = donor, group = dd$group[1],
              Country = dd$Country[1], n_total = nrow(dd), n_positive = sum(dd$expression > 0),
              n_selected = if (available) nrow(ss) else NA_integer_,
              selected_fraction_of_all = if (available) nrow(ss) / nrow(dd) else NA_real_,
              mean_stored_selected = if (nrow(ss)) mean(ss$expression) else NA_real_,
              cutoff_status = row$status, stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
  cutoffs <- do.call(rbind, cutoff.rows); rownames(cutoffs) <- NULL
  statistics <- do.call(rbind, stats.rows); rownames(statistics) <- NULL
  statistics$p_adj_BH <- NA_real_
  for (scope in names(scope.labels)) for (condition in names(conditions)) {
    idx <- which(statistics$scope == scope & statistics$condition == condition & is.finite(statistics$p_value))
    statistics$p_adj_BH[idx] <- p.adjust(statistics$p_value[idx], "BH")
  }
  wide <- merge(cutoffs[cutoffs$input_mode == "all", c("scope", "celltype", "cutoff", "status")],
    cutoffs[cutoffs$input_mode == "positive", c("scope", "celltype", "cutoff", "status")],
    by = c("scope", "celltype"), suffixes = c("_all", "_positive"), sort = FALSE)
  stopifnot(nrow(cutoffs) == 6L * length(celltypes), nrow(wide) == 3L * length(celltypes),
    nrow(statistics) == 12L * length(celltypes))
  write.checked.csv(wide, file.path(out, "01.MFF.two_cutoffs_by_celltype.csv"))
  write.checked.csv(cutoffs, file.path(out, "02.MFF.mclust_fit_details.csv"))
  write.checked.csv(statistics, file.path(out, "03.MFF.four_subset_statistics.csv"))
  write.checked.csv(do.call(rbind, coverage.rows), file.path(out, "04.MFF.donor_subset_coverage.csv"))
  if (save.classifications) for (scope in names(assignments))
    write.checked.csv(assignments[[scope]], file.path(out, paste0("05.MFF.profile_classification.", scope, ".csv.gz")))
  saveRDS(models, file.path(out, "MFF.mclust_models.RDS"))
  print(wide)

  # 5. Mixture diagnostics: actual histograms, weighted component densities, cutoff.
  pngs <- character()
  save.png <- function(p, file, width, height, res = plot.res) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(file, width = width, height = height, res = res, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE); print(p)
    pngs <<- c(pngs, file)
  }
  for (scope in names(scope.labels)) for (mode in c("all", "positive")) {
    b <- assignments[[scope]]; plots <- list()
    for (ct in celltypes) {
      v <- fit.input(b$expression[b$celltype == ct], mode)
      s <- cutoffs[cutoffs$scope == scope & cutoffs$celltype == ct & cutoffs$input_mode == mode, ]
      p <- ggplot(data.frame(expression = v), aes(expression)) + theme_bw(base_size = 13)
      if (length(v)) p <- p + geom_histogram(aes(y = after_stat(density)), bins = 60,
        fill = "#AFC3D3", colour = "white", linewidth = 0.1)
      fit <- models[[scope]][[ct]][[mode]]
      if (!is.null(fit) && is.finite(s$selected_G) && s$selected_G %in% 1:2) {
        mu <- as.numeric(fit$parameters$mean)
        vv <- rep(as.numeric(fit$parameters$variance$sigmasq), length.out = length(mu))
        pro <- as.numeric(fit$parameters$pro)
        if (!length(pro)) pro <- rep(1 / length(mu), length(mu))
        if (length(v) && all(is.finite(c(mu, vv, pro))) && all(vv > 0)) {
          xx <- seq(min(v), max(v), length.out = 500L)
          curves <- do.call(rbind, lapply(order(mu), function(k) data.frame(x = xx,
            density = pro[k] * dnorm(xx, mu[k], sqrt(vv[k])), component = paste0("Mean ", signif(mu[k], 3)))))
          p <- p + geom_line(data = curves, aes(x, density, colour = component), inherit.aes = FALSE, linewidth = 0.8)
        }
      }
      if (is.finite(s$cutoff)) p <- p + geom_vline(xintercept = s$cutoff, linetype = "dashed", linewidth = 0.8)
      subtitle <- paste0("n=", format(length(v), big.mark = ","), " | G=", s$selected_G,
        " | ", if (is.finite(s$cutoff)) paste0("T=", signif(s$cutoff, 4)) else s$status,
        if (nzchar(s$fit_warnings)) " | fit warning: see CSV" else "")
      plots[[ct]] <- p + labs(title = ct, subtitle = subtitle,
        x = "Stored MFF expression", y = "Density") + theme(legend.position = "bottom",
        legend.title = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9), panel.grid.minor = element_blank())
    }
    p <- wrap_plots(plots, ncol = n.columns) + plot_annotation(
      title = paste(scope.labels[scope], gene, mode, "mclust", sep = " | "),
      caption = "Young + Old pooled; cell-weighted Gaussian fits. Exact zeros are a point mass, not a Gaussian population; inspect warnings/BIC. Dashed line = posterior-0.5 single boundary.",
      theme = theme(plot.title = element_text(size = 26, face = "bold"),
        plot.caption = element_text(size = 10, hjust = 0)))
    save.png(p, file.path(out, scope, paste0("MFF.mclust.", mode, ".overview.png")),
      panel.width * n.columns, panel.height * n.rows)
  }

  # 6. Four conditional violin/FC sets per scope. Missing cutoffs/types remain visible.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.23e-308")
    paste0("p = ", formatC(p, format = if (p < 0.001) "e" else "f", digits = if (p < 0.001) 2 else 3))
  }
  fc <- statistics$log2FC_Old_vs_Young
  x.limit <- max(0.25, max(c(0, abs(fc[is.finite(fc)]))) * 1.15)
  for (scope in names(scope.labels)) for (condition in names(conditions)) {
    b <- assignments[[scope]]
    mode <- if (startsWith(condition, "all_")) "all" else "positive"
    region <- if (endsWith(condition, "_high")) "High" else "Low"
    region.column <- paste0("region_", mode)
    st <- statistics[statistics$scope == scope & statistics$condition == condition, ]
    plots <- list()
    for (ct in celltypes) {
      s <- st[st$celltype == ct, ]
      idx <- which(b$celltype == ct & !is.na(b[[region.column]]) & b[[region.column]] == region)
      d <- b[idx, ]; d$group <- factor(d$group, levels = group.order)
      regular <- d[FALSE, ]; special <- data.frame(group = character(), value = numeric(), n = integer())
      for (g in group.order) {
        dd <- d[d$group == g, ]
        if (nrow(dd) >= 2L && length(unique(dd$expression)) > 1L) regular <- rbind(regular, dd) else
          if (nrow(dd)) special <- rbind(special, data.frame(group = g, value = dd$expression[1], n = nrow(dd)))
      }
      top <- max(c(0.2, d$expression))
      p <- ggplot(d, aes(group, expression)) + geom_blank() +
        geom_violin(data = regular, aes(fill = group), trim = TRUE, scale = "width", width = 0.8) +
        geom_boxplot(data = regular, width = 0.18, fill = "white", alpha = 0.85, outlier.shape = NA)
      if (nrow(special)) p <- p + geom_errorbar(data = special,
        aes(group, ymin = value, ymax = value), inherit.aes = FALSE, width = 0.45)
      if (any(special$n == 1L)) p <- p + geom_point(data = special[special$n == 1L, ],
        aes(group, value), inherit.aes = FALSE, size = 2)
      p <- p + annotate("text", x = 1.5, y = top * 1.25, label = format.p(s$p_value), size = 7) +
        annotate("segment", x = 1, xend = 2, y = top * 1.1, yend = top * 1.1) +
        scale_fill_manual(values = group.colors, guide = "none") +
        scale_x_discrete(limits = group.order, drop = FALSE,
          labels = setNames(paste0(group.order, "\n(n=", format(c(s$n_Young, s$n_Old), big.mark = ",", trim = TRUE), ")"), group.order)) +
        coord_cartesian(ylim = c(0, top * 1.42)) +
        labs(title = ct, subtitle = if (is.finite(s$cutoff)) paste0("T = ", signif(s$cutoff, 4)) else s$cutoff_status,
          x = NULL, y = "MFF expression (stored RNA/data)") + theme_bw(base_size = 13) +
        theme(legend.position = "none", plot.title = element_text(face = "bold", size = 20),
          axis.text.x = element_text(size = 17), panel.grid.minor = element_blank())
      if (!nrow(d)) p <- p + annotate("text", x = 1.5, y = top * 0.5,
        label = if (!is.finite(s$cutoff)) "Cutoff unavailable" else "No selected profiles", size = 4)
      plots[[ct]] <- p
    }
    caption <- "Conditional on selected MFF range; raw cell-level Wilcoxon p, not donor-adjusted. Thresholding does not correct QC/detection bias."
    p <- wrap_plots(plots, ncol = n.columns) + plot_annotation(
      title = paste(scope.labels[scope], conditions[condition], sep = " | "), caption = caption,
      theme = theme(plot.title = element_text(size = 26, face = "bold"),
        plot.caption = element_text(size = 10, hjust = 0)))
    folder <- file.path(out, scope, condition)
    save.png(p, file.path(folder, "MFF.Young_Old.violin.overview.png"),
      panel.width * n.columns, panel.height * n.rows)
    st <- st[order(!is.finite(st$log2FC_Old_vs_Young), st$log2FC_Old_vs_Young, st$celltype, na.last = TRUE), ]
    st$celltype_axis <- factor(st$celltype, levels = rev(st$celltype))
    st$direction <- ifelse(st$log2FC_Old_vs_Young < 0, "Old lower", "Old higher")
    st$significance <- ifelse(!is.finite(st$p_value), "p = NA", ifelse(st$p_value < 0.05, "p < 0.05", "p >= 0.05"))
    available <- st[is.finite(st$log2FC_Old_vs_Young), ]; absent <- st[!is.finite(st$log2FC_Old_vs_Young), ]
    p <- ggplot(st, aes(y = celltype_axis)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(data = available, aes(x = log2FC_Old_vs_Young, colour = direction, shape = significance),
        size = 4.5, stroke = 1.2, show.legend = TRUE) +
      geom_text(data = absent, x = x.limit * 0.95, label = "NA", hjust = 1, colour = "grey60") +
      scale_colour_manual(values = c("Old lower" = "#D95F02", "Old higher" = "black"), guide = "none") +
      scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1, "p = NA" = 4),
        limits = c("p < 0.05", "p >= 0.05", "p = NA"), drop = FALSE, name = "Raw cell p") +
      guides(shape = guide_legend(nrow = 1, override.aes = list(colour = "#444444"))) +
      scale_x_continuous(limits = c(-x.limit, x.limit), expand = expansion(mult = 0.02)) +
      scale_y_discrete(drop = FALSE) + labs(title = paste(strwrap(
        paste(scope.labels[scope], conditions[condition], sep = " | "), width = 75), collapse = "\n"),
        x = "Selected-profile mean log2 fold change (Old / Young)", y = NULL,
        caption = paste(strwrap(paste(caption, scale.note, "| inverse-transformed mean +", fc.pseudocount),
          width = 110), collapse = "\n")) +
      theme_classic(base_size = 15) + theme(legend.position = "bottom", legend.title = element_text(size = 11),
        plot.title = element_text(size = 19, face = "bold"), axis.text.y = element_text(size = 16),
        panel.grid.major.y = element_line(colour = "#EEEEEE"), plot.caption = element_text(size = 9, hjust = 0))
    save.png(p, file.path(folder, "MFF.Young_vs_Old.cell_mean.log2FC.png"), 3150L,
      round(max(5.5, 2 + length(celltypes) * 0.48) * 300), res = 300L)
  }
  # 7. Additional paired Lower/Upper violins: six overviews, with the original fits/labels.
  lower.upper.rows <- list()
  for (scope in names(scope.labels)) for (mode in c("all", "positive")) {
    b <- assignments[[scope]]; plots <- list()
    for (ct in celltypes) {
      panel <- prepare.lower.upper(b, cutoffs, scope, mode, ct)
      plots[[ct]] <- lower.upper.violin(panel)
      lower.upper.rows[[length(lower.upper.rows) + 1L]] <- panel$counts
    }
    input.label <- if (mode == "all") "0 included fit" else "0 excluded fit (MFF > 0)"
    p <- patchwork::wrap_plots(plots, ncol = n.columns) + patchwork::plot_annotation(
      title = paste(scope.labels[scope], input.label, "Lower / Upper", sep = " | "),
      caption = "Young + Old pooled within each cell type. Blue: Lower (MFF <= T); red: Upper (MFF > T). Dashed line: original mclust cutoff. Positive mode excludes exact zeros. Descriptive distributions; no Lower/Upper hypothesis test.",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 26, face = "bold"),
        plot.caption = ggplot2::element_text(size = 10, hjust = 0)))
    save.png(p, file.path(out, scope, paste0("MFF.mclust.", mode, ".lower_upper.violin.overview.png")),
      panel.width * n.columns, panel.height * n.rows)
  }
  lower.upper.counts <- do.call(rbind, lower.upper.rows)
  write.checked.csv(lower.upper.counts, file.path(out, "06.MFF.lower_upper_violin_counts.csv"))
  stopifnot(nrow(lower.upper.counts) == 12L * length(celltypes),
    length(pngs) == 36L, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  writeLines(c("MFF mclust: two input populations and four conditional expression comparisons.",
    "Input: existing nonnegative RNA/data; original Seurat/RDS untouched. No jitter or normalization.",
    "Full input is explicitly used for initialization as well as EM; no implicit random initialization subset.",
    "All fit includes exact zeros; positive fit uses expression > 0. Young/Old pooled within scope/type.",
    paste0("G candidates: ", paste(G.candidates, collapse = ","), "; variance models: ", paste(model.names, collapse = ",")),
    "BIC_G1/G2 is the best finite BIC for that G; positive delta favors G=2. NA means untested/failed.",
    "A Gaussian mixture is misspecified for an exact zero point mass; all-fit thresholds are exploratory.",
    "Warnings can concern candidate fits; inspect BIC and selected model. G=2 does not establish two biological states.",
    "Cutoff: weighted Gaussian densities equal (P(high)=0.5), one low-to-high crossing in observed support.",
    "Multiple/no crossings, G=1, unknown/failed convergence, invalid/degenerate parameters, sparse/constant inputs -> cutoff NA, explicit status.",
    "selected_return_code is from mclust BIC returnCodes for the selected G/model; 0=converged, 1=iteration cap. Other candidate failures remain warnings/BIC diagnostics.",
    "No median fallback. High uses > cutoff; low uses <= cutoff. Positive-mode subsets always exclude zero.",
    "min.fit.n=20 and at least 3 distinct values are computational guards, not QC filtering of the original object.",
    "01 CSV has two cutoff/status columns per scope/type; 02 CSV includes BIC/components/uncertainty/warnings.",
    "03 CSV contains four-subset conditional FC/p/coverage; BH across available types separately per scope/condition.",
    "04 CSV counts observed donors/types only; n_selected=NA for unavailable cutoff, 0 for an empty valid subset.",
    "05 compressed CSVs retain every profile with both region labels; Excluded_zero is specific to positive fit.",
    "CSV missing token is .__MFF_MISSING__. (additional underscores only if needed to avoid literal-string collision).",
    "Models saved in MFF.mclust_models.RDS include mclust fits. Failed fits may be absent from nested lists.",
    "36 PNGs: original 30 + 6 paired Lower/Upper violin overviews (3 scopes x 2 input modes).",
    "Lower/Upper overviews reuse original cutoff classifications; Young + Old are pooled, not split on their x-axis.",
    "Blue Lower versus red Upper, n per region, dashed original cutoff, violin plus boxplot; constant/singleton groups use lines/points.",
    "Cutoff-unavailable types remain visible with NA counts; positive mode never includes exact zeros.",
    "06 CSV audits Lower/Upper counts, zero counts, donors and age composition; no new fit, p-value or fold change.",
    "Subset selection conditions on the tested expression; results do not estimate the total-population ageing effect.",
    "FC = log2(mean inverse-transformed Old + pseudocount) - log2(mean inverse-transformed Young + pseudocount).",
    paste0("Pseudocount: ", fc.pseudocount, "; ", scale.note, "; evidence: ", log.info$evidence),
    "FC is pooled-profile weighted; cell p does not account for donor, paired donors, country, technology or batch.",
    "Combined/scRNA/snRNA overlap; combined mixes technologies. These are not independent validations.",
    "Exact zeros can reflect nondetection and biology. Thresholding is not detection-bias correction.",
    "Reruns overwrite matching files in this fixed folder. Completion/file checks do not replace visual review."),
    file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nMFF_MCLUST_CUTOFFS_COMPLETE\nProfiles:", nrow(cells), "; cutoff rows:", nrow(cutoffs),
    "; valid cutoffs:", sum(is.finite(cutoffs$cutoff)), "; PNGs:", length(pngs), "\n", out, "\n")
})
