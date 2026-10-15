# Donor-level mean MFF expression: boxplots + two-sided Wilcoxon.
# Paste ALL code in SERVER R, or source a server copy. No counts required.
local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(Cairo)
  })

  # 1. Settings. Reruns overwrite matching outputs.
  gene <- "MFF"
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  rds.file <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  out <- file.path(project, "analysis/3.Mffexpression/young_old_donors")
  scope.labels <- c(integrated = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  group.order <- c("Young", "Old")
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  min.profiles.per.donor <- 1L       # All observed donor/type pairs retained by default.
  low.profile.warning <- 20L       # Warning only; does NOT exclude profiles/donors.
  min.donors.per.group <- 2L
  n.columns <- 3L
  n.rows <- 5L                     # For 5 columns x 3 rows: change these two settings.
  panel.width <- 900L
  panel.height <- 850L
  plot.res <- 160L
  x.axis.text.size <- 18L
  save.individual.boxplots <- FALSE
  stopifnot(min.profiles.per.donor >= 1L, min.donors.per.group >= 2L)

  # 2. Stored RNA/data is used directly; no inverse transform or new normalization.
  x <- if (exists("hlma.umap", envir = .GlobalEnv, inherits = FALSE)) {
    get("hlma.umap", envir = .GlobalEnv)
  } else {
    if (!file.exists(rds.file)) stop("RDS not found: ", rds.file)
    readRDS(rds.file)
  }
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  meta <- x[[]]
  required <- c("Annotation", "sample", "tech", "age_pop", "Country")
  stopifnot(all(required %in% names(meta)), identical(rownames(meta), colnames(x)),
            !anyDuplicated(colnames(x)))
  for (nm in required) {
    z <- as.character(meta[[nm]])
    if (anyNA(z) || any(trimws(z) == "")) stop("Missing metadata: ", nm)
  }
  stopifnot(all(meta$tech %in% c("scRNA", "snRNA")), all(meta$age_pop %in% c("young_pop", "old_pop")))
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(sum(rownames(rna) == gene) == 1L, identical(colnames(rna), colnames(x)))
  cells <- data.frame(donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
    tech = as.character(meta$tech), Country = as.character(meta$Country),
    group = ifelse(meta$age_pop == "young_pop", "Young", "Old"), expression = as.numeric(rna[gene, ]))
  rm(rna)
  stopifnot(nrow(cells) > 0, all(is.finite(cells$expression)), all(cells$expression >= 0))
  donor.info <- unique(cells[, c("donor", "group", "Country")])
  if (anyDuplicated(donor.info$donor)) stop("Inconsistent age_pop/Country within donor.")
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE))
  stopifnot(length(celltypes) == n.columns * n.rows)
  selected <- list(integrated = cells, scRNA = cells[cells$tech == "scRNA", ],
                   snRNA = cells[cells$tech == "snRNA", ])

  # 3. One row per observed donor/type, including measured zero expression.
  summarize.donors <- function(d) {
    groups <- split(seq_len(nrow(d)), interaction(d$donor, d$celltype, drop = TRUE, sep = "\r"))
    z <- do.call(rbind, lapply(groups, function(ii) {
      b <- d[ii, , drop = FALSE]
      data.frame(donor = b$donor[1], celltype = b$celltype[1], group = b$group[1], Country = b$Country[1],
        n_profiles = nrow(b), n_scRNA = sum(b$tech == "scRNA"), n_snRNA = sum(b$tech == "snRNA"),
        mean_stored_data = mean(b$expression), n_MFF_positive = sum(b$expression > 0),
        MFF_positive_percent = 100 * mean(b$expression > 0))
    }))
    rownames(z) <- NULL
    z
  }
  test.donors <- function(d) {
    y <- d$mean_stored_data[d$group == "Young"]
    o <- d$mean_stored_data[d$group == "Old"]
    r <- data.frame(n_donors_Young = length(y), n_donors_Old = length(o),
      mean_donor_Young = if (length(y)) mean(y) else NA_real_,
      mean_donor_Old = if (length(o)) mean(o) else NA_real_,
      median_donor_Young = if (length(y)) median(y) else NA_real_,
      median_donor_Old = if (length(o)) median(o) else NA_real_,
      W = NA_real_, p_value = NA_real_, test_method = "not_tested", reason = "")
    if (!length(y) || !length(o)) r$reason <- "missing_group" else
      if (min(length(y), length(o)) < min.donors.per.group) r$reason <- "insufficient_donors" else
        if (length(unique(c(y, o))) < 2L) r$reason <- "all_values_identical"
    if (nzchar(r$reason)) return(r)
    exact <- !anyDuplicated(c(y, o)) && max(length(y), length(o)) < 50L
    fit <- tryCatch(stats::wilcox.test(y, o, alternative = "two.sided", paired = FALSE,
      exact = exact, correct = TRUE), error = function(e) e)
    if (inherits(fit, "error")) r$reason <- conditionMessage(fit) else
      if (!is.finite(fit$p.value)) r$reason <- "nonfinite_p" else {
        r$W <- unname(fit$statistic)
        r$p_value <- fit$p.value
        r$test_method <- if (exact) "Wilcoxon_exact" else "Wilcoxon_normal_approximation_tie_corrected"
      }
    r
  }
  donors <- effects <- list()
  for (scope in names(selected)) {
    b <- selected[[scope]]
    if (!nrow(b)) stop("No profiles: ", scope)
    z <- summarize.donors(b)
    z$scope <- scope
    z$included <- z$n_profiles >= min.profiles.per.donor
    stopifnot(sum(z$n_profiles) == nrow(b), all(z$n_profiles == z$n_scRNA + z$n_snRNA),
      !anyDuplicated(z[, c("donor", "celltype")]),
      isTRUE(all.equal(sum(z$n_profiles * z$mean_stored_data), sum(b$expression), tolerance = 1e-10)))
    donors[[scope]] <- z
    s <- do.call(rbind, lapply(celltypes, function(ct) {
      r <- test.donors(z[z$celltype == ct & z$included, ])
      r$celltype <- ct; r$scope <- scope
      r
    }))
    # Adjustment family: all cell types within each scope.
    s$p_adj_BH <- p.adjust(s$p_value, method = "BH", n = length(celltypes))
    effects[[scope]] <- s
    if (any(z$n_profiles[z$included] < low.profile.warning)) warning(scope,
      ": some donor/type means use fewer than ", low.profile.warning,
      " profiles; retained by default, inspect 01.MFF.donor_celltype.csv.", call. = FALSE)
    if (any(s$n_donors_Young < 3 | s$n_donors_Old < 3)) warning(scope,
      ": some cell types have fewer than three donors in a group; inspect test counts/reasons.", call. = FALSE)
  }
  warning("Donor-level tests do not adjust country/batch. Integrated donor means can depend on the scRNA/snRNA mixture; interpret alongside separate scopes.", call. = FALSE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)
  write.csv(do.call(rbind, donors), file.path(out, "01.MFF.donor_celltype.csv"), row.names = FALSE)
  write.csv(do.call(rbind, effects), file.path(out, "02.MFF.donor_Wilcoxon.csv"), row.names = FALSE)

  # 4. Each dot = one donor; X labels count INCLUDED donors, not cells.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return(paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 2)))
    paste0("p = ", formatC(p, format = if (p < 0.001) "e" else "f", digits = if (p < 0.001) 2 else 3))
  }
  save.png <- function(p, file, width, height) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(file, width = width, height = height, res = plot.res, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
  pngs <- character()
  for (scope in names(donors)) {
    plots <- vector("list", length(celltypes))
    for (i in seq_along(celltypes)) {
      ct <- celltypes[i]
      d <- donors[[scope]]
      d <- d[d$celltype == ct & d$included, ]
      d$group <- factor(d$group, levels = group.order)
      s <- effects[[scope]][i, ]
      labels <- setNames(paste0(group.order, "\n(", c(s$n_donors_Young, s$n_donors_Old), " donors)"), group.order)
      top <- max(c(0.05, d$mean_stored_data))
      title <- paste(strwrap(paste(scope.labels[scope], ct, sep = " | "), width = 40), collapse = "\n")
      p <- ggplot(d, aes(group, mean_stored_data, fill = group)) +
        geom_boxplot(width = 0.55, alpha = 0.35, linewidth = 0.5, outlier.shape = NA) +
        geom_point(shape = 21, colour = "black", stroke = 0.3, size = 3,
          position = position_jitter(width = 0.12, height = 0, seed = 42)) +
        scale_fill_manual(values = group.colors, guide = "none") +
        scale_x_discrete(limits = group.order, labels = labels, drop = FALSE) +
        annotate("segment", x = 1, xend = 2, y = top * 1.10, yend = top * 1.10, linewidth = 0.4) +
        annotate("segment", x = c(1, 2), xend = c(1, 2), y = top * 1.065, yend = top * 1.10, linewidth = 0.4) +
        annotate("text", x = 1.5, y = top * 1.26, label = format.p(s$p_value), size = 8) +
        coord_cartesian(ylim = c(0, top * 1.42)) +
        labs(title = title, x = NULL, y = "Donor mean MFF (stored RNA/data)") +
        theme_bw(base_size = 13) + theme(legend.position = "none",
          plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 13),
          axis.text.y = element_text(size = 22), axis.text.x = element_text(size = x.axis.text.size),
          panel.grid.minor = element_blank())
      plots[[i]] <- p
      if (save.individual.boxplots) {
        file <- file.path(out, scope, "individual", sprintf("%02d.%s.Young_Old.donor_boxplot.png", i,
          gsub("[^A-Za-z0-9_-]+", "_", ct)))
        save.png(p, file, panel.width, panel.height)
        pngs <- c(pngs, file)
      }
    }
    p <- wrap_plots(plots, ncol = n.columns, nrow = n.rows) + plot_annotation(
      caption = "Each point = one donor; mean includes zeros. Donor-level Wilcoxon: unadjusted p; country/batch not adjusted.",
      theme = theme(plot.caption = element_text(size = 10, hjust = 0)))
    file <- file.path(out, paste0(scope, ".MFF.Young_Old.donor_boxplot.png"))
    save.png(p, file, panel.width * n.columns, panel.height * n.rows)
    pngs <- c(pngs, file)
  }
  stopifnot(length(pngs) == 3L + if (save.individual.boxplots) 45L else 0L,
            all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  writeLines(c("Donor-level Young/Old MFF boxplots: integrated, scRNA, snRNA.",
    "Input is unchanged stored RNA/data; arithmetic means include zeros. No counts, inverse transform, log2FC or regression.",
    "Each observed donor/celltype has one mean per scope. Absent types are omitted, not assigned mean zero.",
    "Integrated pools profiles within each donor/type; the same donor is counted once, irrespective of technology.",
    "Donors contribute equally to boxplots/tests; integrated means within each donor are weighted by profiles.",
    paste0("Minimum profiles per donor/type: ", min.profiles.per.donor, "; low-profile warning threshold: ", low.profile.warning, "."),
    "CSV retains observed donor/type pairs and marks inclusion. X labels show included donor numbers.",
    "Two-sided unpaired Wilcoxon: exact for small samples without ties; tie-corrected normal approximation otherwise.",
    "At least two donors/group are required. All identical values or insufficient/missing donors yield p=NA with a CSV reason.",
    "Plots show raw p; CSV additionally saves BH across all cell types separately per scope.",
    "No adjustment for country, technology or batch; donor aggregation does not remove confounding.",
    "Scopes overlap and are not independent replications. Cell numbers affect precision of donor means.",
    "Fixed folder and filenames: reruns overwrite matching outputs."), file.path(out, "README.txt"))
  cat("\nMFF_DONOR_BOXPLOTS_COMPLETE\n", out, "\n", sep = "")
})
