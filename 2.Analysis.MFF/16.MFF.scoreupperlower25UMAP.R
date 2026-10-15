# SERVER R: run ALL code after 43 in the SAME R session.
# Reuse EXACT Old bottom/top memberships from hlma.old.score.quartiles.
# Saved hlma.scores UMAP, MFF color, pooled AND cell-type crops. No refitting or RDS I/O.
# Style: GitHub 7.FeaturePlot.tissue.Mff.expression.132430.R and preserved script 8 crops.
local({
  packages <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "scales", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. Fixed output files are overwritten on rerun.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  out <- file.path(project, "analysis/4.genescore/Old_score_quartiles_MFF_UMAP")
  input.object <- "hlma.scores"; cache.object <- "hlma.old.score.quartiles"
  assay <- "RNA"; gene <- "MFF"; reduction <- "umap"
  groupings <- c("pooled", "celltype")
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  group.order <- c("Bottom25", "Top25")
  group.labels <- c(Bottom25 = "Bottom 25%", Top25 = "Top 25%")
  color.modes <- c("actual_max", "positive_quantile") # Use only "positive_quantile" for 18 main PNGs.
  robust.quantile <- .95
  expression.colors <- c("#D9D9D9", "#FEE8C8", "#FDBB84", "#E34A33", "#990000")
  point.size <- .50
  font.base <- 16; font.title <- 20; font.axis <- 16; font.legend <- 13
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L
  n.pair.columns <- 5L # Cell-type overview: 5 bottom/top pairs per row.
  save.individual.celltype.pairs <- FALSE # TRUE also writes a separate pair per cell type.
  min.crop.fraction <- .02; crop.padding <- .04

  # 2. Check current values and the complete cached memberships by exact original ID.
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
  shared.bounds <- function(old, current, embedding, types, minimum.fraction = .02, padding = .04) {
    all <- current
    all$UMAP_1 <- embedding[all$cell_id, 1]; all$UMAP_2 <- embedding[all$cell_id, 2]
    reference <- if (nrow(old)) old else all
    global.span <- max(diff(range(reference$UMAP_1)), diff(range(reference$UMAP_2)), 1e-6)
    pool <- function(v) {
      limits <- range(v); span <- max(diff(limits), global.span * minimum.fraction, 1e-6)
      mean(limits) + c(-1, 1) * span * (.5 + padding)
    }
    xr <- pool(reference$UMAP_1); yr <- pool(reference$UMAP_2)
    rows <- list(data.frame(grouping = "pooled", selection_celltype = "ALL", xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2]))
    for (ct in types) {
      d <- old[old$celltype == ct, , drop = FALSE]
      if (!nrow(d)) d <- all[all$celltype == ct, , drop = FALSE] # Empty Old panel still has valid axes.
      x <- range(d$UMAP_1); y <- range(d$UMAP_2)
      half <- max(diff(x), diff(y), global.span * minimum.fraction, 1e-6) * (.5 + padding)
      rows[[length(rows) + 1L]] <- data.frame(grouping = "celltype", selection_celltype = ct,
        xmin = mean(x) - half, xmax = mean(x) + half, ymin = mean(y) - half, ymax = mean(y) + half)
    }
    do.call(rbind, rows)
  }
  ordered.panel <- function(selected, scope, grouping, ct, score, quartile) {
    z <- selected[selected$scope == scope & selected$grouping == grouping & selected$selection_celltype == ct &
      selected$score == score & selected$quartile == quartile, , drop = FALSE]
    z[order(z$mff, z$cell_id), , drop = FALSE]
  }
  panel.summary <- function(selected, statistics, limits) {
    rows <- list()
    for (i in seq_len(nrow(statistics))) for (mode in limits$mode) for (g in c("Bottom25", "Top25")) {
      st <- statistics[i, ]; limit <- limits[limits$mode == mode, ]
      z <- ordered.panel(selected, st$scope, st$grouping, st$selection_celltype, st$score, g)
      rows[[length(rows) + 1L]] <- data.frame(st[, c("scope", "grouping", "selection_celltype", "score")],
        mode = mode, quartile = g, n_profiles = nrow(z), n_positive = sum(z$mff > 0), n_zero = sum(z$mff == 0),
        panel_max = if (nrow(z)) max(z$mff) else NA_real_,
        shared_color_upper = limit$reference_upper, display_upper = limit$display_upper,
        n_above_color_upper = sum(z$mff > limit$reference_upper), selection_reason = st$selection_reason)
    }
    do.call(rbind, rows)
  }
  write.checked.csv <- function(d, path) {
    token <- ".__OLD_UMAP_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    con <- if (grepl("[.]gz$", path)) gzfile(path, "wt", encoding = "UTF-8") else file(path, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, con, row.names = FALSE, na = token), finally = close(con))
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else "character", character(1))
    con <- if (grepl("[.]gz$", path)) gzfile(path, "rt", encoding = "UTF-8") else file(path, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(con, colClasses = classes, na.strings = token, check.names = FALSE), finally = close(con))
    stopifnot(identical(dim(restored), dim(d)), identical(names(restored), names(d)))
    for (nm in names(d)) if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(restored[[nm]], d[[nm]], tolerance = 1e-10))) else
      stopifnot(identical(restored[[nm]], d[[nm]]))
    invisible(path)
  }

  # 3. Use the scored object and script-43 cache from this session.
  if (!exists(input.object, envir = .GlobalEnv, inherits = FALSE) || !exists(cache.object, envir = .GlobalEnv, inherits = FALSE))
    stop("Run GO-enabled script 34 and script 43 in this R session first.")
  stopifnot(identical(names(score.columns), names(score.labels)), n.pair.columns >= 1L,
    point.size > 0, min.crop.fraction > 0, crop.padding >= 0)
  x <- get(input.object, envir = .GlobalEnv); cache <- get(cache.object, envir = .GlobalEnv)
  stopifnot(inherits(x, "Seurat"), assay %in% Assays(x), "data" %in% Layers(x[[assay]]), reduction %in% Reductions(x))
  if (gene %in% names(x[[]])) stop("Metadata column MFF shadows the RNA feature; resolve that name before FeaturePlot.")
  ids <- colnames(x); rna <- GetAssayData(x, assay = assay, layer = "data")
  stopifnot(identical(colnames(rna), ids), sum(rownames(rna) == gene) == 1L)
  mff <- setNames(as.numeric(rna[gene, ]), ids); rm(rna)
  current <- make.current(x[[]], mff, score.columns, ids)
  types <- sort(unique(current$celltype))
  embedding <- Embeddings(x[[reduction]])
  verified <- verify.cache(cache, current, embedding, score.columns, scopes, groupings, types)
  old <- verified$old; selected <- verified$selected; statistics <- verified$statistics
  limits <- shared.limits(old, color.modes, robust.quantile)
  bounds <- shared.bounds(old, current, embedding, types, min.crop.fraction, crop.padding)
  summary <- panel.summary(selected, statistics, limits)
  file.map <- data.frame(celltype = types, file_stem = sprintf("%03d", seq_along(types)))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(old, file.path(out, "01.Old_UMAP.reference.csv.gz"))
  write.checked.csv(selected, file.path(out, "02.Old_quartile_UMAP.selected_profiles.csv.gz"))
  write.checked.csv(limits, file.path(out, "03.color_scale_limits.csv"))
  write.checked.csv(bounds, file.path(out, "04.shared_UMAP_bounds.csv"))
  write.checked.csv(summary, file.path(out, "05.panel_summary.csv"))
  write.checked.csv(file.map, file.path(out, "06.celltype_file_map.csv"))

  # 4. Match the GitHub FeaturePlot style; only each quartile's own cells are drawn.
  DefaultAssay(x) <- assay # Local copy only, for FeaturePlot's RNA feature lookup.
  template <- Seurat::FeaturePlot(object = x, features = gene, cells = ids[1], reduction = reduction,
    slot = "data", order = TRUE, min.cutoff = 0, max.cutoff = max(limits$display_upper),
    keep.scale = NULL, cols = c(expression.colors[1], tail(expression.colors, 1)),
    pt.size = point.size, label = FALSE, raster = FALSE, coord.fixed = FALSE, combine = FALSE)[[1]]
  umap.panel <- function(z, stat, limit, window, title, quartile) {
    cap <- limit$display_upper
    if (nrow(z)) {
      p <- Seurat::FeaturePlot(object = x, features = gene, cells = z$cell_id, reduction = reduction,
        slot = "data", order = TRUE, min.cutoff = 0, max.cutoff = cap, keep.scale = NULL,
        cols = c(expression.colors[1], tail(expression.colors, 1)), pt.size = point.size,
        label = FALSE, raster = FALSE, coord.fixed = FALSE, combine = FALSE)[[1]]
    } else {
      p <- ggplot(z, aes(UMAP_1, UMAP_2, color = mff)) + geom_blank() + template$theme
    }
    mode.label <- if (limit$mode == "actual_max") "Shared actual maximum" else paste0("Shared positive MFF Q", 100 * robust.quantile)
    p <- p + scale_color_gradientn(colours = expression.colors, limits = c(0, cap), oob = scales::squish,
      breaks = c(0, cap/2, cap), labels = function(v) format(signif(v, 3), trim = TRUE),
      name = "MFF\nRNA/data", guide = guide_colorbar(barheight = grid::unit(4, "cm"))) +
      coord_cartesian(xlim = unname(c(window$xmin, window$xmax)), ylim = unname(c(window$ymin, window$ymax)), expand = FALSE) +
      labs(title = paste(paste(strwrap(title, width = 30), collapse = "\n"), group.labels[quartile], sep = "\n"),
        subtitle = paste0(mode.label, " = ", signif(limit$reference_upper, 4)),
        caption = paste0("n = ", nrow(z), "; MFF > 0: ", sum(z$mff > 0),
          if (limit$mode == "positive_quantile") "\nAbove-limit values use the maximum color." else ""),
        x = "umap_1", y = "umap_2") +
      theme(text = element_text(size = font.base), plot.title = element_text(size = font.title, face = "bold", hjust = .5),
        axis.title = element_text(size = font.axis), axis.text = element_text(size = font.axis),
        legend.title = element_text(size = font.legend), legend.text = element_text(size = font.legend),
        plot.margin = margin(12, 12, 12, 12))
    if (!nrow(z)) p <- p + annotate("text", x = mean(c(window$xmin, window$xmax)),
      y = mean(c(window$ymin, window$ymax)), label = "No selected cells", color = "#666666", size = 5)
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
      stopifnot(nrow(z) == n)
      umap.panel(z, stat, limit, window, title, g)
    })
    wrap_plots(plots, ncol = 2L, guides = "collect")
  }
  save.png <- function(p, path, width, height) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(path, width = width, height = height, res = plot.res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  manifest <- list(); pngs <- character()
  for (mode in color.modes) for (scope in scopes) for (nm in names(score.columns)) {
    limit <- limits[limits$mode == mode, , drop = FALSE]
    heading <- paste(scope.labels[scope], score.labels[nm], "Old score quartiles", sep = " | ")
    caption <- "Exact script-43 memberships and saved UMAP. Left: bottom quartile. Right: top quartile. Color: stored MFF, shared across every panel."
    if ("pooled" %in% groupings) {
      p <- paired.umap(scope, "pooled", "ALL", nm, limit) +
        plot_annotation(title = heading, caption = paste(strwrap(caption, width = 120), collapse = "\n"),
          theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "pooled", scope, mode, paste0("Old.", nm, ".Bottom25_vs_Top25.MFF_UMAP.png"))
      save.png(p, path, 2L * panel.width, panel.height); pngs <- c(pngs, path)
      manifest[[length(manifest) + 1L]] <- data.frame(grouping = "pooled", scope = scope, score = nm, mode = mode,
        selection_celltype = "ALL", n_pairs = 1L, width = 2L * panel.width, height = panel.height, path = path)
    }
    if ("celltype" %in% groupings) {
      pairs <- lapply(types, function(ct) paired.umap(scope, "celltype", ct, nm, limit))
      cols <- min(n.pair.columns, length(types)); rows <- ceiling(length(types)/cols)
      p <- wrap_plots(pairs, ncol = cols, guides = "collect") +
        plot_annotation(title = paste(heading, "Within each cell type", sep = " | "),
          caption = paste(strwrap(caption, width = 180), collapse = "\n"),
          theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 9, hjust = 0)))
      path <- file.path(out, "celltype", scope, mode, paste0("Old.", nm, ".MFF_UMAP.celltype_overview.png"))
      save.png(p, path, 2L * panel.width * cols, panel.height * rows); pngs <- c(pngs, path)
      manifest[[length(manifest) + 1L]] <- data.frame(grouping = "celltype_overview", scope = scope, score = nm, mode = mode,
        selection_celltype = "ALL", n_pairs = length(types), width = 2L * panel.width * cols, height = panel.height * rows, path = path)
      if (save.individual.celltype.pairs) for (i in seq_along(types)) {
        path <- file.path(out, "celltype", scope, mode, file.map$file_stem[i], paste0("Old.", nm, ".Bottom25_vs_Top25.MFF_UMAP.png"))
        p <- pairs[[i]] + plot_annotation(title = heading,
          theme = theme(plot.title = element_text(size = 22, face = "bold")))
        save.png(p, path, 2L * panel.width, panel.height); pngs <- c(pngs, path)
        manifest[[length(manifest) + 1L]] <- data.frame(grouping = "celltype_pair", scope = scope, score = nm, mode = mode,
          selection_celltype = types[i], n_pairs = 1L, width = 2L * panel.width, height = panel.height, path = path)
      }
    }
  }
  expected <- length(color.modes) * length(scopes) * length(score.columns) *
    (as.integer("pooled" %in% groupings) + as.integer("celltype" %in% groupings) *
      (1L + as.integer(save.individual.celltype.pairs) * length(types)))
  stopifnot(length(pngs) == expected, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))
  manifest <- do.call(rbind, manifest)
  write.checked.csv(manifest, file.path(out, "07.PNG_manifest.csv"))
  assign("hlma.old.score.quartile.umap", list(old_profiles = old, selected_profiles = selected,
    color_limits = limits, bounds = bounds, panel_summary = summary, manifest = manifest), envir = .GlobalEnv)
  writeLines(c("Old-only score-quartile MFF UMAPs: exact script-43 bottom/top cell memberships reused.",
    "Input: hlma.scores + hlma.old.score.quartiles in the same R session. Cached IDs/metadata/MFF/all scores checked against current object.",
    "No quartile reselection, normalization, AddModuleScore, mclust, PCA, RunUMAP or RDS reading/writing.",
    "Only the selected group's own cells are drawn. MFF=0 stays grey; high expression uses orange to dark red, full opacity.",
    "Continuous palette: #D9D9D9, #FEE8C8, #FDBB84, #E34A33, #990000. High expression draws last; no cell downsampling.",
    "All genes/scopes/groupings share actual-max and positive-MFF-Q95 limits computed ONCE from unique Old reference profiles, including middle scores.",
    "Q95 is a visualization limit, not a selection cutoff. Above-limit colors saturate only; CSV retains original MFF values.",
    "All-zero/empty reference: actual limit=0, display limit=1. Empty quartiles remain labeled panels, not zero-expression fake cells.",
    "Pooled views share global Old UMAP axes. Cell-type crops share each type's complete Old window across scores/scopes/quartiles/color modes.",
    "If a type has no Old profiles, its full original type coordinates define an empty panel window. No Young points are drawn.",
    "FeaturePlot coord.fixed=FALSE and explicit coord_cartesian match the reference's flexible panel aspect.",
    "Default main outputs: 2 modes x 3 scopes x 3 scores x (pooled pair + celltype overview) = 36 PNGs.",
    "Celltype overview has 5 bottom/top pairs per row; each pair contains 2 UMAP panels. Optional individual pairs use numeric type folder mapping.",
    "UMAP is descriptive, with no statistical test or donor/composition adjustment added.",
    "Reference: https://github.com/daehyeon32/202609_singlecell_aging_MFF/blob/main/scripts/3.Preprocessing_seurat/7.FeaturePlot.tissue.Mff.expression.132430.R",
    "Verified reference blob SHA: 815913a2519327cc6e3431252ef79d7e6f1b44e7",
    "Cell-type crop follows the preserved 8.FeaturePlot.celltype.Mff.expression.132430.R snapshot.",
    paste0("Saved reduction: ", reduction), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_OLD_SCORE_QUARTILES_MFF_UMAP_COMPLETE\nOld reference profiles: ", nrow(old),
    "; PNGs: ", length(pngs), "\n", out, "\n", sep = "")
})
