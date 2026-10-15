# SERVER R: run ALL code after scripts 35 and 37 in the SAME R session.
# Ten heatmaps: basic + four mclust conditions, each combined and snRNA/scRNA.
# Reuse cell-type-specific Spearman results. No correlation/score/mclust recalculation or RDS I/O.
local({
  packages <- c("ggplot2", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(ggplot2))

  # 1. Settings. All ten figures share one signed color scale and one row order.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  basic.object <- "hlma.correlations"          # Script 35, cell-type-specific full population.
  mclust.object <- "hlma.mclust.correlations"  # Script 37, original four mclust conditions.
  analysis.unit <- "cell"  # Both source caches also contain "donor" statistics.
  analysis.view <- "All"   # All = Young+Old together. "Young" or "Old" selects that existing result.
  out <- file.path(project, "analysis/4.genescore/MFF_score_heatmaps", paste(analysis.unit, analysis.view, sep = "_"))
  conditions <- c("basic", "all_high", "all_low", "positive_high", "positive_low")
  condition.labels <- c(basic = "Basic", all_high = "All-mode High", all_low = "All-mode Low",
    positive_high = "Positive-mode High", positive_low = "Positive-mode Low")
  score.order <- c("Hallmark_IFN_alpha", "SaulSenMayo", "GO_IFN_beta")
  score.labels <- c(Hallmark_IFN_alpha = "Hallmark\nIFN-alpha", SaulSenMayo = "SaulSenMayo", GO_IFN_beta = "GO\nIFN-beta")
  celltype.order <- NULL  # NULL uses alphabetical order, identical across all ten heatmaps.
  expected.celltypes <- 15L
  rho.maximum <- .3
  significance.p.cutoff <- .05  # STRICT raw p < 0.05; not p_adj_BH and not <= 0.05.
  rho.colors <- c(negative = "#EF3B2C", zero = "#FFFFFF", positive = "#2166AC", missing = "#E5E5E5")
  border.colors <- c(negative = "#FF0000", positive = "#0000FF")
  border.linewidth <- 1.4; border.tile.size <- .92
  font.family <- "sans"; font.base.size <- 15; font.title.size <- 24
  font.axis.title.size <- 19; font.score.size <- 12; font.celltype.size <- 15
  font.legend.title.size <- 15; font.legend.text.size <- 13
  legend.bar.height.cm <- 8; legend.bar.width.cm <- .7
  tile.aspect <- .5
  plot.width.in <- c(combined = 8, technologies = 10)
  plot.height.minimum.in <- 6; plot.height.per.celltype.in <- .48; png.dpi <- 300L

  # 2. Check that baseline and subset results refer to the same full scored profiles.
  check.references <- function(basic, subset) {
    identity.fields <- c("cell_id", "celltype", "donor", "tech", "group", "Country")
    paired.fields <- c("mff", paste0("AMS_", score.order))
    fields <- c(identity.fields, paired.fields)
    for (d in list(basic, subset)) {
      if (!is.data.frame(d) || !nrow(d) || !all(fields %in% names(d))) stop("Missing full paired reference profiles; rerun scripts 35 and 37.")
      for (nm in identity.fields) if (anyNA(d[[nm]]) || any(!nzchar(trimws(as.character(d[[nm]]))))) stop("Invalid reference metadata: ", nm)
      if (anyDuplicated(d$cell_id)) stop("Duplicate reference cell IDs.")
      for (nm in paired.fields) if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]]))) stop("Invalid reference paired values: ", nm)
    }
    if (!setequal(basic$cell_id, subset$cell_id)) stop("Baseline/mclust reference cell IDs differ; rerun both on the same object.")
    subset <- subset[match(basic$cell_id, subset$cell_id), , drop = FALSE]
    for (nm in identity.fields) if (!identical(as.character(basic[[nm]]), as.character(subset[[nm]])))
      stop("Baseline/mclust reference metadata differs: ", nm)
    for (nm in paired.fields) if (any(abs(basic[[nm]] - subset[[nm]]) > 1e-10 * pmax(1, abs(basic[[nm]]))))
      stop("Baseline/mclust reference paired values differ: ", nm)
    invisible(TRUE)
  }
  prepare.heatmap.tables <- function(basic, subset, types, unit, view, maximum, p.cutoff) {
    cases <- c("basic", "all_high", "all_low", "positive_high", "positive_low")
    scores <- score.order; scopes <- c("combined", "snRNA", "scRNA")
    fields <- c("scope", "celltype", "score", "view", "unit", "rho", "p_value", "p_adj_BH", "n_observations", "n_profiles", "n_donors", "reason")
    if (length(unit) != 1L || !unit %in% c("cell", "donor") || length(view) != 1L || !view %in% c("All", "Young", "Old")) stop("Invalid analysis unit/view.")
    if (length(maximum) != 1L || !is.finite(maximum) || maximum <= 0 || maximum > 1 ||
        length(p.cutoff) != 1L || !is.finite(p.cutoff) || p.cutoff <= 0 || p.cutoff > 1) stop("Invalid rho maximum or p cutoff.")
    if (!length(types) || anyNA(types) || anyDuplicated(types) || any(types == "")) stop("Invalid cell-type order.")
    if (!is.data.frame(basic) || !is.data.frame(subset) || !all(fields %in% names(basic)) ||
        !all(c(fields, "condition") %in% names(subset))) stop("Missing correlation-statistic columns; use scripts 35 and 37.")
    if (anyNA(subset$condition) || any(!subset$condition %in% cases[-1])) stop("Unknown mclust condition in statistics.")
    b <- basic[basic$unit == unit & basic$view == view & basic$score %in% scores, fields, drop = FALSE]
    b$condition <- "basic"
    m <- subset[subset$unit == unit & subset$view == view & subset$score %in% scores, c(fields, "condition"), drop = FALSE]
    d <- rbind(b, m)
    for (nm in c("condition", "scope", "celltype", "score", "view", "unit", "reason")) {
      d[[nm]] <- as.character(d[[nm]])
      if (anyNA(d[[nm]])) stop("Missing statistic label: ", nm)
    }
    if (any(!d$scope %in% scopes) || any(!d$celltype %in% types)) stop("Unexpected scope/cell type in selected statistics.")
    for (nm in c("rho", "p_value", "p_adj_BH", "n_observations", "n_profiles", "n_donors"))
      if (!is.numeric(d[[nm]])) stop("Nonnumeric correlation statistic: ", nm)
    if (any(!is.finite(d$rho) & !is.na(d$rho)) || any(abs(d$rho) > 1, na.rm = TRUE) ||
        any(!is.finite(d$p_value) & !is.na(d$p_value)) || any(d$p_value < 0 | d$p_value > 1, na.rm = TRUE)) stop("Invalid raw rho/p values.")
    for (nm in c("n_observations", "n_profiles", "n_donors"))
      if (any(!is.finite(d[[nm]])) || any(d[[nm]] < 0 | d[[nm]] != floor(d[[nm]]))) stop("Invalid observation counts: ", nm)
    tested <- d$reason == ""
    if (any(tested & (!is.finite(d$rho) | !is.finite(d$p_value) | d$n_observations < 3L)))
      stop("A tested row lacks a valid rho/p or enough observations.")
    key <- c("condition", "scope", "celltype", "score")
    if (anyDuplicated(d[, key])) stop("Duplicated condition/scope/celltype/score correlation keys.")
    grid <- expand.grid(condition = cases, scope = scopes, celltype = types, score = scores, stringsAsFactors = FALSE)
    grid$.grid_order <- seq_len(nrow(grid)); d$.present <- TRUE
    if (nrow(d) != nrow(grid)) stop("Incomplete selected correlation grid; missing or extra results.")
    d <- merge(grid, d, by = key, all.x = TRUE, sort = FALSE)
    if (anyNA(d$.present)) stop("Missing condition/scope/celltype/score correlation results.")
    d <- d[order(d$.grid_order), ]; d$.present <- d$.grid_order <- NULL; rownames(d) <- NULL
    d$tested <- d$reason == "" & is.finite(d$rho) & is.finite(d$p_value)
    d$rho_display <- ifelse(d$tested, pmax(-maximum, pmin(maximum, d$rho)), NA_real_)
    d$color_clipped <- d$tested & is.finite(d$rho) & abs(d$rho) > maximum
    d$border <- ""
    significant <- d$tested & d$p_value < p.cutoff
    d$border[significant & d$rho < 0] <- "negative"
    d$border[significant & d$rho > 0] <- "positive"
    d
  }
  layout.data <- function(data, condition, layout, types) {
    scores <- score.order
    scopes <- if (layout == "combined") "combined" else if (layout == "technologies") c("snRNA", "scRNA") else stop("Invalid heatmap layout.")
    d <- data[data$condition == condition & data$scope %in% scopes, , drop = FALSE]
    if (nrow(d) != length(types) * length(scores) * length(scopes)) stop("Incorrect heatmap grid size: ", condition, " / ", layout)
    columns <- if (layout == "combined") scores else
      unlist(lapply(scopes, function(tech) paste(tech, scores, sep = ".")), use.names = FALSE)
    d$column_id <- if (layout == "combined") d$score else paste(d$scope, d$score, sep = ".")
    d$column_axis <- factor(d$column_id, levels = columns)
    d$celltype_axis <- factor(d$celltype, levels = rev(types))
    if (anyNA(d$column_axis) || anyNA(d$celltype_axis) || anyDuplicated(d[, c("celltype", "column_id")])) stop("Invalid heatmap positions.")
    d
  }
  write.checked.csv <- function(d, file) {
    token <- ".__HLMA_HEATMAP_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    write.csv(d, file, row.names = FALSE, na = token, fileEncoding = "UTF-8")
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else if (is.logical(v)) "logical" else "character", character(1))
    restored <- read.csv(file, colClasses = classes, check.names = FALSE, na.strings = token, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    stopifnot(identical(dim(d), dim(restored)), identical(names(d), names(restored)))
    for (nm in names(d)) if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(d[[nm]], restored[[nm]], tolerance = 1e-10))) else
      stopifnot(identical(d[[nm]], restored[[nm]]))
    invisible(file)
  }

  # 3. Select the existing cell-type results, not script 38's pooled statistics.
  for (nm in c(basic.object, mclust.object)) if (!exists(nm, envir = .GlobalEnv, inherits = FALSE))
    stop("Run scripts 35 and 37 in this same R session first: missing ", nm)
  basic <- get(basic.object, envir = .GlobalEnv); subsets <- get(mclust.object, envir = .GlobalEnv)
  if (!is.list(basic) || !is.list(subsets)) stop("Invalid correlation caches; rerun scripts 35 and 37.")
  check.references(basic$profiles, subsets$profiles)
  types <- if (is.null(celltype.order)) sort(unique(as.character(basic$profiles$celltype))) else celltype.order
  if (length(types) != expected.celltypes || anyDuplicated(types) || !setequal(types, basic$profiles$celltype))
    stop("Expected exactly ", expected.celltypes, " distinct cell types, covering the original full reference.")
  data <- prepare.heatmap.tables(basic$statistics, subsets$statistics, types, analysis.unit, analysis.view, rho.maximum, significance.p.cutoff)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.checked.csv(data, file.path(out, "01.heatmap_raw_and_display_values.csv"))

  # 4. Reference-style tiles: no numbers inside, white dividers, inset significant boxes.
  make.heatmap <- function(d, title, layout) {
    labels <- if (layout == "combined") score.labels[score.order] else
      unlist(lapply(c("snRNA", "scRNA"), function(tech) setNames(
        paste(tech, gsub("\n", " ", unname(score.labels[score.order]), fixed = TRUE), sep = "\n"),
        paste(tech, score.order, sep = "."))), use.names = TRUE)
    positive <- d[d$border == "positive", , drop = FALSE]; negative <- d[d$border == "negative", , drop = FALSE]
    caption <- paste0("Inset boxes: raw p < ", significance.p.cutoff, "; negative red, positive blue. Grey: unavailable.\n",
      if (analysis.unit == "cell") "Cell Spearman; donor dependence, Country and technology not adjusted." else
        "Donor-mean Spearman; Country, age and technology not adjusted.")
    ggplot(d, aes(column_axis, celltype_axis, fill = rho_display)) +
      geom_tile(colour = "white", linewidth = .8) +
      geom_tile(data = negative, aes(column_axis, celltype_axis), inherit.aes = FALSE, fill = NA,
        colour = border.colors["negative"], linewidth = border.linewidth, width = border.tile.size, height = border.tile.size, show.legend = FALSE) +
      geom_tile(data = positive, aes(column_axis, celltype_axis), inherit.aes = FALSE, fill = NA,
        colour = border.colors["positive"], linewidth = border.linewidth, width = border.tile.size, height = border.tile.size, show.legend = FALSE) +
      scale_fill_gradient2(low = rho.colors["negative"], mid = rho.colors["zero"], high = rho.colors["positive"],
        midpoint = 0, limits = c(-rho.maximum, rho.maximum), breaks = c(-rho.maximum, -rho.maximum / 2, 0, rho.maximum / 2, rho.maximum),
        na.value = rho.colors["missing"], name = "Spearman\nrho",
        guide = guide_colorbar(barheight = grid::unit(legend.bar.height.cm, "cm"), barwidth = grid::unit(legend.bar.width.cm, "cm"),
          frame.colour = "black", ticks.colour = "black")) +
      scale_x_discrete(limits = levels(d$column_axis), labels = labels, drop = FALSE) + scale_y_discrete(drop = FALSE) +
      coord_fixed(ratio = tile.aspect) + labs(title = title, x = "Score", y = "Cell type", caption = caption) +
      theme_bw(base_size = font.base.size, base_family = font.family) +
      theme(plot.title = element_text(size = font.title.size, face = "bold", hjust = .5, margin = margin(b = 15)),
        axis.title = element_text(size = font.axis.title.size), axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.text.x = element_text(size = font.score.size, colour = "black", angle = 35, hjust = 1, vjust = 1, lineheight = .9),
        axis.text.y = element_text(size = font.celltype.size, colour = "black"), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right", legend.title = element_text(size = font.legend.title.size), legend.text = element_text(size = font.legend.text.size),
        plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 12)), plot.margin = margin(15, 20, 15, 15))
  }
  save.png <- function(p, file, width.in, height.in) {
    Cairo::CairoPNG(file, width = round(width.in * png.dpi), height = round(height.in * png.dpi), res = png.dpi, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  files <- character(); manifest <- list()
  height.in <- max(plot.height.minimum.in, 2.7 + length(types) * plot.height.per.celltype.in)
  for (condition in conditions) for (layout in c("combined", "technologies")) {
    d <- layout.data(data, condition, layout, types)
    title <- paste("HLMA", condition.labels[condition], if (layout == "combined") "Combined" else "snRNA + scRNA", analysis.view, sep = " | ")
    p <- make.heatmap(d, title, layout)
    file <- file.path(out, paste0("HLMA.", condition, ".", layout, ".MFF_score_rho_heatmap.png"))
    save.png(p, file, plot.width.in[layout], height.in); files <- c(files, file)
    manifest[[length(manifest) + 1L]] <- data.frame(condition = condition, layout = layout, unit = analysis.unit, view = analysis.view,
      n_rows = length(types), n_columns = length(levels(d$column_axis)), n_tiles = nrow(d), n_tested = sum(d$tested),
      n_significant_negative = sum(d$border == "negative"), n_significant_positive = sum(d$border == "positive"),
      n_unavailable = sum(!d$tested), n_color_clipped = sum(d$color_clipped), file = file)
  }
  manifest <- do.call(rbind, manifest); rownames(manifest) <- NULL
  stopifnot(length(files) == 10L, all(file.exists(files)), all(file.info(files)$size > 8L),
    all(manifest$n_rows == 15L), all(manifest$n_columns == rep(c(length(score.order), 2L * length(score.order)), length(conditions))))
  write.checked.csv(manifest, file.path(out, "02.heatmap_manifest.csv"))
  assign("hlma.score.heatmaps", list(data = data, manifest = manifest, celltype_order = types), envir = .GlobalEnv)
  writeLines(c("Ten CELL-TYPE-SPECIFIC MFF-score Spearman heatmaps using existing script 35 and 37 statistics.",
    "Cases: basic (full population), all_high, all_low, positive_high, positive_low. No correlations, scores or mclust fits recalculated.",
    "Each case: combined 15 rows x 3 score columns; technical panel 15 rows x 6 columns.",
    "Combined columns: Hallmark IFN-alpha, SaulSenMayo, GO IFN-beta. Technical columns: snRNA Hallmark, snRNA SenMayo, snRNA GO IFN-beta, scRNA Hallmark, scRNA SenMayo, scRNA GO IFN-beta.",
    "Same alphabetical row order across ten figures (or the explicitly configured celltype.order). All 15 types remain, including unavailable results.",
    "Fill: negative red, zero white, positive blue; fixed -0.3 to +0.3. Beyond-range rho saturates at the endpoint color.",
    "The CSV keeps original rho and raw p; rho_display is only a clipped display value, not a recalculated correlation.",
    "Strict RAW p < 0.05: rho<0 red inset box; rho>0 blue inset box. p=0.05, rho=0 and unavailable tests have no colored box.",
    "No numbers inside tiles. Grey indicates an unavailable correlation; the original reason and n are retained in the CSV.",
    "Statistical families/BH values come from source analyses; box decisions use raw p_value, not p_adj_BH.",
    "Baseline and subset full-profile IDs, metadata, MFF and all three scores must agree; duplicate/incomplete statistic grids are rejected.",
    "Cell p ignores within-donor dependence, Country/technology and other confounders. Subset correlations are conditional on MFF selection.",
    "All and positive subsets overlap, combined and technical scopes overlap; these are not independent replications.",
    "No use of script 38 pooled statistics, RDS read/write, profile filtering, normalization or source-object modification.",
    "Ten PNGs, raw/display long CSV and manifest CSV. Both CSVs are readback-verified.",
    paste0("Analysis unit: ", analysis.unit), paste0("Analysis view: ", analysis.view), paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_TEN_SCORE_HEATMAPS_COMPLETE\nHeatmaps: ", length(files), "; cell types: ", length(types), "\n", out, "\n", sep = "")
})
