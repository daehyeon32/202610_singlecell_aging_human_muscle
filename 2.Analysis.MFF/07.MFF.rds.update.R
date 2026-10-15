local({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  suppressPackageStartupMessages(library(Seurat))
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  source.rds <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  result.folder <- file.path(project, "analysis/3.Mffexpression/mclust_cutoffs")
  profile.file <- file.path(result.folder, "05.MFF.profile_classification.combined.csv.gz")
  fit.file <- file.path(result.folder, "02.MFF.mclust_fit_details.csv")
  output.rds <- file.path(result.folder, "HLMA.scsn.seurat.MFF_mclust.RDS")
  verify.saved.rds <- TRUE  # Read back saved object; requires memory for an additional object.
  assign.hlma.mff <- TRUE  # Original hlma.umap remains unchanged.

  # Numeric missing values are decoded only in numeric/label fields, preserving literal cell IDs.
  is.saved.missing <- function(v) {
    is.na(v) | v == "NA" | v == "" | grepl("^\\.__MFF_MISSING__\\._*$", v)
  }
  saved.numeric <- function(v, field) {
    if (is.numeric(v)) {
      if (any(!is.finite(v) & !is.na(v))) stop("Nonfinite numeric field: ", field)
      return(as.numeric(v))
    }
    v <- as.character(v); missing <- is.saved.missing(v)
    r <- suppressWarnings(as.numeric(v)); r[missing] <- NA_real_
    if (any(!missing & !is.finite(r))) stop("Invalid numeric field: ", field)
    r
  }
  read.saved.table <- function(path) {
    if (!file.exists(path)) stop("Missing mclust result file: ", path,
      "\nRun the full mclust script first with save.classifications=TRUE.")
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt", encoding = "UTF-8") else
      file(path, "rt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    read.csv(con, check.names = FALSE, colClasses = "character", na.strings = character(),
      stringsAsFactors = FALSE)
  }
  make.mff.metadata <- function(ids, meta, expression, profiles, fits) {
    required.meta <- c("Annotation", "sample", "tech", "age_pop")
    required.profiles <- c("cell_id", "donor", "celltype", "tech", "group", "expression",
      "cutoff_all", "region_all", "cutoff_positive", "region_positive")
    required.fits <- c("scope", "celltype", "input_mode", "cutoff", "status", "selected_G",
      "model_name", "selected_return_code", "n_total", "n_input")
    if (!all(required.meta %in% names(meta)) || !all(required.profiles %in% names(profiles)) ||
        !all(required.fits %in% names(fits))) stop("Required metadata/result columns missing.")
    if (length(ids) != nrow(meta) || !identical(ids, rownames(meta)) || anyNA(ids) ||
        any(ids == "") || anyDuplicated(ids)) stop("Invalid original cell IDs/order.")
    if (length(expression) != length(ids) || any(!is.finite(expression)) || any(expression < 0))
      stop("Original MFF expression must be finite/nonnegative.")
    for (field in required.meta) if (anyNA(meta[[field]]) || any(as.character(meta[[field]]) == ""))
      stop("Missing original metadata: ", field)
    if (!all(as.character(meta$age_pop) %in% c("young_pop", "old_pop"))) stop("Unknown age_pop label.")
    if (!all(as.character(meta$tech) %in% c("scRNA", "snRNA"))) stop("Unknown technology label.")
    if (anyNA(profiles$cell_id) || any(profiles$cell_id == "") || anyDuplicated(profiles$cell_id) ||
        nrow(profiles) != length(ids) || !setequal(profiles$cell_id, ids))
      stop("Classification CSV cell IDs do not match ALL object profiles one-to-one.")
    profiles <- profiles[match(ids, profiles$cell_id), , drop = FALSE]
    comparisons <- list(celltype = as.character(meta$Annotation), donor = as.character(meta$sample),
      tech = as.character(meta$tech), group = ifelse(meta$age_pop == "young_pop", "Young", "Old"))
    for (field in names(comparisons)) if (!identical(as.character(profiles[[field]]), comparisons[[field]]))
      stop("Object/classification metadata mismatch: ", field)
    original.stored <- saved.numeric(profiles$expression, "expression")
    if (anyNA(original.stored) || any(abs(original.stored - expression) > 1e-10 * pmax(1, abs(expression))))
      stop("Current RNA/data MFF differs from the classified expression. Use the original input object.")
    # All cells are matched to COMBINED fits, estimated separately within their original cell type.
    f <- fits[!is.na(fits$scope) & fits$scope == "combined", , drop = FALSE]
    if (anyDuplicated(f[, c("celltype", "input_mode")]) || anyNA(f$celltype) || anyNA(f$input_mode))
      stop("Duplicate/missing combined celltype/input_mode fit rows.")
    types <- unique(as.character(meta$Annotation))
    if (nrow(f) != 2L * length(types) || !setequal(f$celltype, types) ||
        !setequal(f$input_mode, c("all", "positive"))) stop("Incomplete combined fit coverage.")
    added <- data.frame(row.names = ids)
    for (mode in c("all", "positive")) {
      fm <- f[f$input_mode == mode, , drop = FALSE]
      if (nrow(fm) != length(types) || !setequal(fm$celltype, types)) stop("Missing fit rows: ", mode)
      fm <- fm[match(as.character(meta$Annotation), fm$celltype), , drop = FALSE]
      cutoff <- saved.numeric(fm$cutoff, paste0(mode, " cutoff"))
      g <- saved.numeric(fm$selected_G, "selected_G")
      code <- saved.numeric(fm$selected_return_code, "selected_return_code")
      total <- saved.numeric(fm$n_total, "n_total"); count <- saved.numeric(fm$n_input, "n_input")
      expected.total <- as.numeric(table(as.character(meta$Annotation))[as.character(meta$Annotation)])
      expected.input <- if (mode == "all") expected.total else {
        counts <- tapply(expression > 0, as.character(meta$Annotation), sum)
        as.numeric(counts[as.character(meta$Annotation)])
      }
      if (anyNA(total) || anyNA(count) || any(total != expected.total) || any(count != expected.input))
        stop("Fit input counts do not match original profiles: ", mode)
      status <- as.character(fm$status)
      if (anyNA(status) || any(is.saved.missing(status))) stop("Missing fit status.")
      available <- is.finite(cutoff)
      if (any(available & (cutoff < 0 | status != "ok" | is.na(g) | g != 2 | is.na(code) | code != 0)))
        stop("Usable cutoff lacks a converged two-component fit: ", mode)
      if (any(!available & status == "ok")) stop("Fit status ok without cutoff.")
      if (any(!is.na(g) & !g %in% c(1, 2))) stop("Unexpected selected component count.")
      source.cutoff <- saved.numeric(profiles[[paste0("cutoff_", mode)]], paste0("profile cutoff_", mode))
      if (!identical(is.na(source.cutoff), is.na(cutoff)) ||
          any(abs(source.cutoff[available] - cutoff[available]) > 1e-10 * pmax(1, abs(cutoff[available]))))
        stop("Profile/fit cutoff mismatch: ", mode)
      label <- as.character(profiles[[paste0("region_", mode)]])
      label[is.saved.missing(label)] <- NA_character_
      eligible <- if (mode == "all") rep(TRUE, length(ids)) else expression > 0
      classified <- available & eligible
      if (any(is.na(label[classified])) || any(!label[classified] %in% c("High", "Low")))
        stop("Missing/unknown classification at a usable cutoff: ", mode)
      if (any(!is.na(label[!available & eligible]))) stop("Unclassified fit was assigned High/Low.")
      if (mode == "positive" && any(is.na(label[!eligible]) | label[!eligible] != "Excluded_zero"))
        stop("Zero-expression profiles must be Excluded_zero in positive mode.")
      # Retain the original full-precision classification; CSV cutoffs may be rounded near a boundary.
      delta <- expression[classified] - cutoff[classified]
      tol <- 1e-10 * pmax(1, abs(cutoff[classified]))
      if (any((label[classified] == "High" & delta < -tol) |
              (label[classified] == "Low" & delta > tol))) stop("High/Low label disagrees with expression/cutoff.")
      prefix <- paste0("MFF_mclust_", mode)
      added[[prefix]] <- label
      added[[paste0(prefix, "_cutoff")]] <- cutoff
    }
    added
  }
  label.summary <- function(meta, added) {
    rows <- list()
    for (mode in c("all", "positive")) {
      label <- added[[paste0("MFF_mclust_", mode)]]
      label[is.na(label)] <- "Unclassified"
      d <- data.frame(celltype = as.character(meta$Annotation),
        group = ifelse(meta$age_pop == "young_pop", "Young", "Old"), label = label,
        n_profiles = 1L, stringsAsFactors = FALSE)
      s <- aggregate(n_profiles ~ celltype + group + label, data = d, FUN = sum)
      s$input_mode <- mode; rows[[mode]] <- s
    }
    do.call(rbind, rows)
  }
  write.annotated.rds <- function(object, path, source.path, verify = NULL) {
    if (identical(normalizePath(path, mustWork = FALSE), normalizePath(source.path, mustWork = FALSE)))
      stop("Output RDS must differ from original source RDS.")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(pattern = "MFF_mclust_annotated_", tmpdir = dirname(path), fileext = ".RDS")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(object, temporary, compress = "gzip")
    if (!is.null(verify)) verify(readRDS(temporary))
    if (!file.rename(temporary, path)) stop("Cannot replace output RDS: ", path)
    stopifnot(file.exists(path), file.info(path)$size > 0)
    invisible(path)
  }

  # Start from the original FILE every time, including when hlma.umap/hlma.mff already exists.
  loaded.from <- source.rds
  x <- {
    if (!file.exists(source.rds)) stop("Source RDS missing: ", source.rds)
    readRDS(source.rds)
  }
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  ids <- colnames(x); original.meta <- x[[]]
  if (any(grepl("^MFF_mclust_", names(original.meta))))
    stop("Source RDS already contains MFF_mclust metadata. Set source.rds to the ORIGINAL unannotated RDS.")
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == "MFF") == 1L)
  expression <- as.numeric(rna["MFF", ]); rm(rna)
  profiles <- read.saved.table(profile.file); fits <- read.saved.table(fit.file)
  added <- make.mff.metadata(ids, original.meta, expression, profiles, fits)
  annotation.columns <- c("MFF_mclust_all", "MFF_mclust_all_cutoff",
    "MFF_mclust_positive", "MFF_mclust_positive_cutoff")
  stopifnot(identical(names(added), annotation.columns), ncol(added) == 4L)
  annotated <- SeuratObject::AddMetaData(x, metadata = added)
  stopifnot(identical(colnames(annotated), ids), identical(annotated[[]][, names(added), drop = FALSE], added),
    ncol(annotated[[]]) == ncol(original.meta) + 4L,
    setequal(setdiff(names(annotated[[]]), names(original.meta)), annotation.columns))
  retained.columns <- setdiff(names(original.meta), names(added))
  stopifnot(identical(annotated[[]][, retained.columns, drop = FALSE], original.meta[, retained.columns, drop = FALSE]))
  # Assays, reductions/UMAP, identities, graphs and other slots must stay unchanged in memory.
  for (slot.name in setdiff(slotNames(x), c("meta.data", "misc")))
    if (!identical(slot(x, slot.name), slot(annotated, slot.name))) stop("Unexpected object change: ", slot.name)
  if (!identical(x[[]], original.meta)) stop("Original object's metadata changed.")
  info <- list(annotation_scope = "combined", gene = "MFF", assay = "RNA", layer = "data",
    grouping = "Young + Old pooled, separately per original Annotation; scRNA + snRNA combined",
    label_definition = "Original cutoff labels: High > T; Low <= T; positive mode excludes zeros",
    unavailable = "NA label when no usable cutoff; fit status kept in cutoff_details, not cell metadata",
    provenance_files = normalizePath(c(profile.file, fit.file), mustWork = TRUE),
    provenance_md5 = tools::md5sum(c(profile.file, fit.file)),
    input_object = loaded.from, n_profiles = length(ids),
    cutoff_details = fits[fits$scope == "combined", , drop = FALSE],
    annotation_columns = names(added), timestamp = as.character(Sys.time()))
  annotated@misc$MFF_mclust <- info

  # Save to a NEW file. Existing matching annotated output is overwritten on rerun.
  cat("Saving annotated object (all profiles): ", length(ids), "\n", output.rds, "\n", sep = "")
  verify <- if (verify.saved.rds) function(restored) {
    stopifnot(inherits(restored, "Seurat"), identical(colnames(restored), ids),
      identical(restored[[]], annotated[[]]), identical(restored@misc$MFF_mclust, info))
    rr <- GetAssayData(restored, assay = "RNA", layer = "data")
    stopifnot(identical(as.numeric(rr["MFF", ]), expression))
    for (reduction in Reductions(x))
      stopifnot(identical(Embeddings(restored[[reduction]]), Embeddings(x[[reduction]])))
    invisible(TRUE)
  } else NULL
  write.annotated.rds(annotated, output.rds, source.rds, verify = verify)
  invisible(gc())
  summary <- label.summary(original.meta, added)
  stopifnot(sum(summary$n_profiles) == 2L * length(ids))
  write.csv(summary, file.path(result.folder, "06.MFF.metadata_annotation_summary.csv"), row.names = FALSE)
  writeLines(c("Annotated RDS retains every original profile, assay, reduction and identity.",
    "MFF_mclust_all is exactly the combined all_high/all_low cutoff classification from the saved analysis.",
    "Fits include zero expression and pool Young+Old separately within each original cell type.",
    "MFF_mclust_positive stores the separate positive-only fit; zero profiles remain Excluded_zero.",
    "NA group means unavailable cutoff. Do not treat NA as Low.",
    "Exactly four new metadata columns: MFF_mclust_all, MFF_mclust_all_cutoff, MFF_mclust_positive, MFF_mclust_positive_cutoff.",
    "G, model and status are retained only in the existing fit CSV and @misc$MFF_mclust$cutoff_details.",
    "These are relative expression annotations; original Annotation, RNA/data and UMAP are unchanged.",
    "Near-boundary CSV rounding is tolerated for label validation; original CSV classification is retained.",
    "Fit provenance and file checksums are in @misc$MFF_mclust.",
    paste0("Loaded input: ", loaded.from), paste0("Output RDS: ", output.rds),
    paste0("Saved object readback verification: ", verify.saved.rds),
    "The original source.rds FILE is read on every run, even when global objects exist.",
    "Reruns overwrite this annotated output using the original RDS and add these four columns.",
    "Global hlma.mff is assigned if enabled; global hlma.umap and original RDS are preserved."),
    file.path(result.folder, "README.MFF_annotated_RDS.txt"))
  if (assign.hlma.mff) assign("hlma.mff", annotated, envir = .GlobalEnv)
  print(head(added))
  print(table(added$MFF_mclust_all, useNA = "ifany"))
  print(table(added$MFF_mclust_positive, useNA = "ifany"))
  cat("\nMFF_MCLUST_ANNOTATED_RDS_COMPLETE\nProfiles: ", length(ids),
    "\nNew metadata columns: ", paste(names(added), collapse = ", "), "\n", output.rds, "\n", sep = "")
})
