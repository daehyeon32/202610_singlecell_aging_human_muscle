# MFF expression/detection overview: 3 technologies x 3 age subsets.
# Paste ALL code in SERVER R, or source a server copy.
# MFF-positive = stored expression > 0. >= 0 is also exported as requested.
local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(Cairo)
  })

  # 1. Settings: fixed output folder; matching files are overwritten.
  gene <- "MFF"
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  rds.file <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  out <- file.path(project, "analysis/3.Mffexpression/expression_detection")
  scope.labels <- c(integrated = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  age.labels <- c(Young_Old = "Young + Old", Young = "Young", Old = "Old")
  width <- 3000L
  height <- 2200L
  resolution <- 180L
  comma <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)

  # 2. Read the unchanged stored MFF vector and metadata.
  x <- if (exists("hlma.umap", envir = .GlobalEnv, inherits = FALSE)) {
    get("hlma.umap", envir = .GlobalEnv)
  } else {
    if (!file.exists(rds.file)) stop("RDS not found: ", rds.file)
    readRDS(rds.file)
  }
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  meta <- x[[]]
  stopifnot(all(c("Annotation", "sample", "tech", "age_pop") %in% names(meta)),
            identical(rownames(meta), colnames(x)), !anyDuplicated(colnames(x)))
  for (nm in c("Annotation", "sample", "tech", "age_pop")) {
    z <- as.character(meta[[nm]])
    if (anyNA(z) || any(trimws(z) == "")) stop("Missing metadata: ", nm)
  }
  stopifnot(all(as.character(meta$tech) %in% c("scRNA", "snRNA")))
  stopifnot(all(as.character(meta$age_pop) %in% c("young_pop", "old_pop")))
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(sum(rownames(rna) == gene) == 1L, identical(colnames(rna), colnames(x)))
  cells <- data.frame(celltype = as.character(meta$Annotation), donor = as.character(meta$sample),
                      tech = as.character(meta$tech), age_group = ifelse(meta$age_pop == "young_pop", "Young", "Old"),
                      expression = as.numeric(rna[gene, ]))
  rm(rna)
  stopifnot(nrow(cells) > 0, all(is.finite(cells$expression)), all(cells$expression >= 0))
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE))
  select.profiles <- function(d) {
    selected <- list()
    for (scope in names(scope.labels)) for (age in names(age.labels)) {
      keep <- (scope == "integrated" | d$tech == scope) & (age == "Young_Old" | d$age_group == age)
      selected[[paste(scope, age, sep = ".")]] <- d[keep, , drop = FALSE]
    }
    selected
  }
  selected <- select.profiles(cells)

  # 3. Mean includes zeros; absent types have NA expression/percentage.
  summary.of <- function(d) {
    n <- nrow(d)
    positive <- d$expression > 0
    data.frame(n_profiles = n, n_donors = length(unique(d$donor)),
      n_MFF_positive = sum(positive), n_MFF_zero = sum(d$expression == 0),
      n_MFF_ge0 = sum(d$expression >= 0),
      MFF_positive_percent = if (n) 100 * mean(positive) else NA_real_,
      mean_stored_data = if (n) mean(d$expression) else NA_real_,
      mean_positive_stored_data = if (any(positive)) mean(d$expression[positive]) else NA_real_)
  }
  overall <- by.type <- list()
  for (key in names(selected)) {
    parts <- strsplit(key, ".", fixed = TRUE)[[1]]
    scope <- parts[1]; age <- parts[2]
    b <- selected[[key]]
    if (!nrow(b)) stop("No profiles: ", key)
    s <- do.call(rbind, lapply(celltypes, function(ct) {
      z <- summary.of(b[b$celltype == ct, ])
      z$celltype <- ct
      z$scope <- scope
      z$age_group <- age
      z
    }))
    total <- summary.of(b)
    total$scope <- scope
    total$age_group <- age
    by.type[[key]] <- s
    overall[[key]] <- total
    stopifnot(sum(s$n_profiles) == nrow(b), sum(s$n_MFF_positive) == total$n_MFF_positive,
      all(s$n_MFF_positive + s$n_MFF_zero == s$n_profiles), all(s$n_MFF_ge0 == s$n_profiles),
      isTRUE(all.equal(sum(s$n_profiles * s$mean_stored_data, na.rm = TRUE),
                      sum(b$expression), tolerance = 1e-10)))
  }
  for (age in names(age.labels)) for (field in c("n_profiles", "n_MFF_positive")) {
    stopifnot(all(by.type[[paste0("integrated.", age)]][[field]] ==
      by.type[[paste0("scRNA.", age)]][[field]] + by.type[[paste0("snRNA.", age)]][[field]]))
  }
  for (scope in names(scope.labels)) for (field in c("n_profiles", "n_MFF_positive")) {
    stopifnot(all(by.type[[paste0(scope, ".Young_Old")]][[field]] ==
      by.type[[paste0(scope, ".Young")]][[field]] + by.type[[paste0(scope, ".Old")]][[field]]))
  }
  overall.table <- do.call(rbind, overall)
  type.table <- do.call(rbind, by.type)
  rownames(overall.table) <- rownames(type.table) <- NULL
  print(overall.table)
  message("MFF > 0 means detected. MFF >= 0 includes zeros and equals all profiles in this input.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)
  write.csv(overall.table, file.path(out, "01.MFF.overall.csv"), row.names = FALSE)
  write.csv(type.table, file.path(out, "02.MFF.by_celltype.csv"), row.names = FALSE)

  # 4. Same panels, type order and axes across all nine figures.
  common.theme <- theme_classic(base_size = 15) + theme(
    plot.title = element_text(size = 20, face = "bold"), axis.text = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)), axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    plot.margin = margin(10, 15, 10, 10), legend.position = "none")
  mean.limit <- max(0.01, type.table$mean_stored_data, na.rm = TRUE) * 1.30
  count.limit <- max(type.table$n_profiles) * 1.55
  pngs <- character()
  for (key in names(selected)) {
    total <- overall[[key]]
    scope <- total$scope; age <- total$age_group
    s <- by.type[[key]]
    s$celltype_axis <- factor(s$celltype, levels = rev(celltypes))
    available <- s[s$n_profiles > 0, ]
    absent <- s[s$n_profiles == 0, ]
    p.mean <- ggplot(s, aes(y = celltype_axis)) +
      geom_col(data = available, aes(x = mean_stored_data), fill = "#A50F15", width = 0.65, orientation = "y") +
      geom_text(data = available, aes(x = mean_stored_data, label = sprintf("%.3f", mean_stored_data)),
                hjust = -0.15, size = 4) +
      geom_text(data = absent, x = 0, label = "NA", hjust = -0.15, size = 4, colour = "grey50") +
      scale_x_continuous(limits = c(0, mean.limit), expand = expansion(mult = c(0, 0))) +
      scale_y_discrete(drop = FALSE) + labs(title = "Mean MFF expression", x = "Mean stored RNA/data (zeros included)") + common.theme
    p.count <- ggplot(s, aes(y = celltype_axis)) +
      geom_col(aes(x = n_profiles), fill = "grey88", width = 0.65, orientation = "y") +
      geom_col(aes(x = n_MFF_positive), fill = "#377EB8", width = 0.65, orientation = "y") +
      geom_text(aes(x = n_profiles, label = paste0(comma(n_MFF_positive), " / ", comma(n_profiles))),
                hjust = -0.10, size = 4) +
      scale_x_continuous(limits = c(0, count.limit), labels = comma, expand = expansion(mult = c(0, 0))) +
      scale_y_discrete(drop = FALSE) + labs(title = "MFF-positive profiles", x = "Profiles: blue = MFF > 0; grey = total") +
      common.theme + theme(axis.text.y = element_blank())
    p.percent <- ggplot(s, aes(y = celltype_axis)) +
      geom_col(data = available, aes(x = MFF_positive_percent), fill = "#1B9E77", width = 0.65, orientation = "y") +
      geom_text(data = available, aes(x = MFF_positive_percent, label = sprintf("%.1f%%", MFF_positive_percent)),
                hjust = -0.15, size = 4) +
      geom_text(data = absent, x = 0, label = "NA", hjust = -0.15, size = 4, colour = "grey50") +
      scale_x_continuous(limits = c(0, 115), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
      scale_y_discrete(drop = FALSE) + labs(title = "MFF-positive percentage", x = "MFF > 0 / all profiles of this type (%)") +
      common.theme + theme(axis.text.y = element_blank())
    p <- wrap_plots(p.mean, p.count, p.percent, ncol = 3, widths = c(1.2, 1.3, 1)) +
      plot_annotation(title = paste(scope.labels[scope], "MFF expression and detection", sep = " | "),
        subtitle = paste0(age.labels[age], " profiles | total = ", comma(total$n_profiles),
          " | MFF > 0 = ", comma(total$n_MFF_positive), " (", sprintf("%.2f", total$MFF_positive_percent), "%)"),
        caption = "Within each cell type: labels in the middle panel = MFF-positive / total. Profiles represent cells or nuclei; means are pooled-cell weighted.",
        theme = theme(plot.title = element_text(size = 26, face = "bold"), plot.subtitle = element_text(size = 15),
                      plot.caption = element_text(size = 11, hjust = 0)))
    # Keep the original Young+Old filenames; append age for the six added figures.
    prefix <- if (age == "Young_Old") scope else paste(scope, age, sep = ".")
    file <- file.path(out, paste0(prefix, ".MFF.expression_detection.png"))
    Cairo::CairoPNG(file, width = width, height = height, res = resolution, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    pngs <- c(pngs, file)
    cat("Saved: ", file, "\n", sep = "")
  }
  stopifnot(length(pngs) == 9L, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  writeLines(c("MFF expression/detection overview: integrated, scRNA, snRNA x Young+Old, Young, Old (9 figures).",
    "All original Young/Old profiles retained. Integrated pools scRNA/snRNA; no new embedding or expression correction.",
    "Age subsets use author age_pop: young_pop=Young, old_pop=Old; no new age cutoffs.",
    "MFF-positive is stored RNA/data > 0. Literal >= 0 is exported in n_MFF_ge0 and includes all zeros.",
    "Means include zeros and remain on the stored scale; no inverse transform, log2FC, counts, or pseudobulk.",
    "The positive percentage denominator is all measured profiles of the same cell type within the same technology and age subset.",
    "Absent types have n_profiles=0 and mean/percentage=NA. All-zero measured types have mean=0 and positive percentage=0.",
    "Means and percentages weight profiles, not donors equally; no significance testing or age effect claim.",
    "The same donor observed in both technologies is counted once in integrated n_donors.",
    "Fixed folder and filenames: reruns overwrite matching outputs. All nine figures share the same axis limits.",
    "Young+Old retains original scope.MFF.expression_detection.png names; age-specific files add .Young or .Old after scope."),
    file.path(out, "README.txt"))
  cat("\nMFF_DETECTION_FIGURES_COMPLETE\n", out, "\n")
})
