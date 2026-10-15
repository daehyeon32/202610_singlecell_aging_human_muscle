# Metadata-only composition check before adapting 4.Analysis_Mff/01.
# Paste into the same SERVER R console, or source a server copy.
# Uses hlma.umap; does not read expression, exclude cells, or change the RDS.

local({
  packages <- c("ggplot2", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  library(ggplot2)
  if (!exists("hlma.umap", envir = .GlobalEnv, inherits = FALSE)) {
    stop("Run the UMAP input check in this R session first to create hlma.umap.")
  }
  x <- get("hlma.umap", envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"))
  meta <- x[[]]
  required <- c("Annotation", "sample", "tech", "age", "age_pop", "Country")
  stopifnot(all(required %in% names(meta)), nrow(meta) == ncol(x),
            identical(rownames(meta), colnames(x)), !anyDuplicated(rownames(meta)))
  meta <- meta[, required, drop = FALSE]
  for (nm in required) {
    meta[[nm]] <- as.character(meta[[nm]])
    if (anyNA(meta[[nm]]) || any(trimws(meta[[nm]]) == "")) stop("Missing metadata: ", nm)
  }
  stopifnot(all(meta$tech %in% c("scRNA", "snRNA")),
            all(meta$age_pop %in% c("young_pop", "old_pop")))
  donor.info <- unique(meta[, c("sample", "age", "age_pop", "Country")])
  if (anyDuplicated(donor.info$sample)) stop("A donor has inconsistent age, age_pop, or Country.")
  donor.info$age <- suppressWarnings(as.numeric(donor.info$age))
  stopifnot(all(is.finite(donor.info$age)))
  donor.info$group <- ifelse(donor.info$age_pop == "young_pop", "Young", "Old")
  donor.info <- donor.info[order(match(donor.info$group, c("Young", "Old")),
                                donor.info$age, donor.info$Country, donor.info$sample), ]
  rownames(donor.info) <- NULL
  celltypes <- names(sort(table(meta$Annotation), decreasing = TRUE))
  color.names <- sort(celltypes)
  celltype.colors <- setNames(grDevices::hcl.colors(length(color.names), "Dark 3"), color.names)
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")

  long.tables <- list()
  summaries <- list()
  for (technology in c("scRNA", "snRNA")) {
    m <- meta[meta$tech == technology, , drop = FALSE]
    if (!nrow(m)) next
    donors <- donor.info$sample[donor.info$sample %in% m$sample]
    tab <- table(celltype = factor(m$Annotation, levels = celltypes),
                 donor = factor(m$sample, levels = donors))
    totals <- colSums(tab)
    z <- as.data.frame(tab, responseName = "n_cells", stringsAsFactors = FALSE)
    z$celltype <- as.character(z$celltype)
    z$donor <- as.character(z$donor)
    z$n_cells <- as.integer(z$n_cells)
    z$tech <- technology
    z$n_total_cells <- as.integer(totals[z$donor])
    z$percent <- 100 * z$n_cells / z$n_total_cells
    info <- donor.info[match(z$donor, donor.info$sample), ]
    for (nm in c("age", "age_pop", "Country", "group")) z[[nm]] <- info[[nm]]
    stopifnot(sum(z$n_cells) == nrow(m), nrow(z) == length(celltypes) * length(donors),
              all(z$n_total_cells > 0), all(is.finite(z$percent)),
              all(abs(tapply(z$percent, z$donor, sum) - 100) < 1e-8))
    long.tables[[technology]] <- z
    for (group in c("Young", "Old")) {
      for (ct in celltypes) {
        a <- z[z$group == group & z$celltype == ct, , drop = FALSE]
        if (!nrow(a)) next
        summaries[[length(summaries) + 1L]] <- data.frame(
          tech = technology, group = group, celltype = ct,
          n_donors = nrow(a), n_donors_with_celltype = sum(a$n_cells > 0),
          n_cells = sum(a$n_cells), n_total_cells_in_group = sum(a$n_total_cells),
          pooled_cell_percent = 100 * sum(a$n_cells) / sum(a$n_total_cells),
          donor_mean_percent = mean(a$percent), donor_sd_percent = stats::sd(a$percent)
        )
      }
    }
  }
  if (!length(long.tables)) stop("No observed donor/technology combinations")
  long <- do.call(rbind, long.tables)
  summary <- do.call(rbind, summaries)
  rownames(long) <- rownames(summary) <- NULL
  stopifnot(sum(long$n_cells) == nrow(meta))
  coverage <- expand.grid(tech = c("scRNA", "snRNA"), group = c("Young", "Old"),
                          Country = sort(unique(meta$Country)), stringsAsFactors = FALSE)
  coverage$n_cells <- coverage$n_donors <- 0L
  for (i in seq_len(nrow(coverage))) {
    group.label <- if (coverage$group[i] == "Young") "young_pop" else "old_pop"
    a <- meta[meta$tech == coverage$tech[i] & meta$age_pop == group.label &
                meta$Country == coverage$Country[i], , drop = FALSE]
    coverage$n_cells[i] <- nrow(a)
    coverage$n_donors[i] <- length(unique(a$sample))
  }
  stopifnot(sum(coverage$n_cells) == nrow(meta))

  root <- paste0("/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle/",
                 "analysis/1.processing/celltype_counts")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root)) stop("Cannot create output root: ", root)
  out <- tempfile(pattern = paste0("composition_",
    format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Seoul"), "_"), tmpdir = root)
  if (!dir.create(out)) stop("Cannot create new output directory: ", out)
  cat("Output directory: ", out, "\n", sep = "")
  tables <- list("01.donor_metadata" = donor.info,
                 "02.donor_celltype_composition.long" = long,
                 "03.Young_Old_celltype_summary" = summary,
                 "04.tech_group_country_coverage" = coverage)
  for (nm in names(tables)) {
    filename <- file.path(out, paste0(nm, ".csv"))
    utils::write.csv(tables[[nm]], filename, row.names = FALSE, fileEncoding = "UTF-8")
    saved <- utils::read.csv(filename, check.names = FALSE, fileEncoding = "UTF-8")
    stopifnot(identical(names(saved), names(tables[[nm]])), nrow(saved) == nrow(tables[[nm]]))
    if ("n_cells" %in% names(saved)) stopifnot(sum(saved$n_cells) == sum(tables[[nm]]$n_cells))
    cat("Saved: ", filename, "\n", sep = "")
  }
  save.png <- function(p, filename, width = 2300L, height = 1700L) {
    Cairo::CairoPNG(filename = filename, width = width, height = height, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
  theme.shared <- theme_classic(base_size = 22) +
    theme(plot.title = element_text(size = 28, face = "bold"),
          axis.title = element_text(size = 26), axis.text = element_text(size = 20),
          legend.text = element_text(size = 20), plot.margin = margin(15, 20, 15, 15))
  pngs <- character()
  for (technology in names(long.tables)) {
    z <- long.tables[[technology]]
    donor.order <- donor.info$sample[donor.info$sample %in% z$donor]
    info <- donor.info[match(donor.order, donor.info$sample), ]
    donor.labels <- setNames(paste0(info$sample, " (", info$group, ", ",
                                   info$Country, ", ", info$age, "y)"), info$sample)
    z$donor_plot <- factor(z$donor, levels = rev(donor.order))
    z$celltype <- factor(z$celltype, levels = celltypes)
    z$group <- factor(z$group, levels = c("Young", "Old"))
    for (metric in c("n_cells", "percent")) {
      z$.value <- z[[metric]]
      p <- ggplot(z, aes(.value, donor_plot, fill = celltype)) +
        geom_col(width = 0.8, orientation = "y") +
        scale_fill_manual(values = celltype.colors, drop = FALSE) +
        scale_y_discrete(labels = donor.labels) +
        labs(title = paste(technology, "Donor cell-type composition", sep = " | "),
             x = if (metric == "percent") "Profiles within donor (%)" else "Number of cell/nucleus profiles",
             y = NULL, fill = "Cell type") + theme.shared
      filename <- file.path(out, paste0(technology, ".donor_composition.", metric, ".png"))
      save.png(p, filename, height = max(1500L, 300L + 55L * length(donor.order)))
      pngs <- c(pngs, filename)
    }
    p <- ggplot(z, aes(group, percent)) +
      geom_boxplot(width = 0.55, fill = "grey95", outlier.shape = NA) +
      geom_point(aes(colour = group), size = 2.8,
                 position = position_jitter(width = 0.1, height = 0, seed = 1234)) +
      scale_colour_manual(values = group.colors, guide = "none") +
      facet_wrap(~celltype, ncol = 5, scales = "free_y", drop = FALSE) +
      scale_y_continuous(limits = c(0, NA)) +
      labs(title = paste(technology, "Young / Old donor composition", sep = " | "),
           subtitle = "Each point = one donor; cell-type absence within a measured donor = 0%",
           x = NULL, y = "Profiles within donor (%)") + theme.shared +
      theme(strip.text = element_text(size = 19), axis.text = element_text(size = 17),
            plot.subtitle = element_text(size = 19))
    filename <- file.path(out, paste0(technology, ".Young_Old.donor_percent.png"))
    save.png(p, filename, width = 2800L, height = 1800L)
    pngs <- c(pngs, filename)
  }
  stopifnot(length(pngs) == 3L * length(long.tables), all(file.exists(pngs)),
            all(file.info(pngs)$size > 8L))
  writeLines(c("Metadata-only composition analysis; source hlma.umap is unchanged.",
    "All cells/nuclei retained; no expression-count matrix is used.",
    "scRNA and snRNA are analyzed separately. sample is donor; orig.ident is not donor.",
    "Only observed donor/technology combinations are included. An absent technology is not a measured zero.",
    "Every cell type is included within each observed donor/technology, including zero-cell combinations.",
    "Percent denominator is total cells/nuclei for the same donor and technology.",
    "Young/Old donor_mean_percent weights donors equally; pooled_cell_percent weights cell profiles.",
    "No minimum-cell threshold, statistical test, pseudobulk, or age effect inference is applied.",
    "Completion verifies output files exist; figure layout still needs visual inspection."),
    file.path(out, "README.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  print(coverage, row.names = FALSE)
  cat("\nCells:", nrow(meta), "; donors:", nrow(donor.info),
      "; cell types:", length(celltypes), "; PNGs:", length(pngs), "\n")
  cat("CELLTYPE_COMPOSITION_COMPLETE\n", out, "\n")
})
