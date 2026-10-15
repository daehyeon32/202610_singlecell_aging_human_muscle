# Human HLMA: THREE AddModuleScores, cell-type-specific Young/Old comparisons.
# Paste ALL code into SERVER R, or source a server copy. Violin plots include boxplots.
# Reference: preserved GitHub 11.1.Mff.genescore.new.R (data-only score portion).
local({
  packages <- c("Seurat", "ggplot2", "patchwork", "Cairo")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

  # 1. Settings. Fixed folders/files are overwritten on rerun.
  project <- "/BiO/Live/dleogus32/202609Agineg_MFF_human_muscle"
  source.rds <- file.path(project, "analysis/1.processing/HLMA.scsn.seurat.RDS")
  out <- file.path(project, "analysis/4.genescore/young_old_SenMayo_IFNalpha")
  input.object <- "auto"  # Or an explicit global object name, e.g. "hlma" or "hlma.mff".
  save.individual.violins <- FALSE  # Main outputs: 9 violin overviews + 9 difference dotplots.
  min.matched.genes <- 1L          # Same computational minimum as mouse 11.1; coverage is reported.
  low.coverage.warning <- 0.50
  nbin <- 24L; ctrl <- 100L; score.seed <- 1L
  scopes <- c("combined", "scRNA", "snRNA")
  scope.labels <- c(combined = "scRNA + snRNA", scRNA = "scRNA", snRNA = "snRNA")
  score.labels <- c(SaulSenMayo = "SaulSenMayo", Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta")
  score.columns <- c(SaulSenMayo = "AMS_SaulSenMayo", Hallmark_IFN_alpha = "AMS_Hallmark_IFN_alpha", GO_IFN_beta = "AMS_GO_IFN_beta")
  group.order <- c("Young", "Old")
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  panel.width <- 900L; panel.height <- 850L; plot.res <- 160L; n.columns <- 5L

  # 2. HUMAN gene-symbol snapshots from official MSigDB: JSONs on 2026-10-06,
  # and the GO IFN-beta human gene-set member table on 2026-10-07.
  # SenMayo: 125 original source identifiers map to 124 current human symbols.
  # These are NOT the mouse .Mm.grp files and no synonym/ortholog search is performed.
  gene.sources <- c(
    SaulSenMayo = "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SAUL_SEN_MAYO",
    Hallmark_IFN_alpha = "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE",
    GO_IFN_beta = "https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_RESPONSE_TO_INTERFERON_BETA")
  source.dates <- c(SaulSenMayo = "2026-10-06", Hallmark_IFN_alpha = "2026-10-06", GO_IFN_beta = "2026-10-07")
  genesets <- list(
    SaulSenMayo = strsplit(paste(
      "ACVR1B ANG ANGPT1 ANGPTL4 AREG AXL BEX3 BMP2 BMP6 C3 CCL1 CCL13 CCL16 CCL2 CCL20 CCL24 CCL26 CCL3 CCL4 CCL5 CCL7 CCL8",
      "CD55 CD9 CSF1 CSF2 CSF2RB CST4 CTNNB1 CTSB CXCL1 CXCL10 CXCL12 CXCL16 CXCL2 CXCL3 CXCL8 CXCR2 DKK1 EDN1 EGF EGFR EREG ESM1 ETS2 FAS",
      "FGF1 FGF2 FGF7 GDF15 GEM GMFG HGF HMGB1 ICAM1 ICAM3 IGF1 IGFBP1 IGFBP2 IGFBP3 IGFBP4 IGFBP5 IGFBP6 IGFBP7 IL10 IL13 IL15 IL18 IL1A IL1B IL2 IL32 IL6 IL6ST IL7",
      "INHA IQGAP2 ITGA2 ITPKA JUN KITLG LCP1 MIF MMP1 MMP10 MMP12 MMP13 MMP14 MMP2 MMP3 MMP9 NAP1L4 NRG1 PAPPA PECAM1 PGF PIGF PLAT PLAU PLAUR PTBP1 PTGER2 PTGES",
      "RPS6KA5 SCAMP4 SELPLG SEMA3F SERPINB4 SERPINE1 SERPINE2 SPP1 SPX TIMP2 TNF TNFRSF10C TNFRSF11B TNFRSF1A TNFRSF1B TUBGCP2 VEGFA VEGFC VGF WNT16 WNT2"), "\\s+")[[1]],
    Hallmark_IFN_alpha = strsplit(paste(
      "ADAR B2M BATF2 BST2 C1S CASP1 CASP8 CCRL2 CD47 CD74 CMPK2 CMTR1 CNP CSF1 CXCL10 CXCL11 DDX60 DHX58 EIF2AK2 ELF1 EPSTI1 GBP2 GBP4 GMPR HELZ2 HERC6 HLA-C",
      "IFI27 IFI30 IFI35 IFI44 IFI44L IFIH1 IFIT2 IFIT3 IFITM1 IFITM2 IFITM3 IL15 IL4R IL7 IRF1 IRF2 IRF7 IRF9 ISG15 ISG20 LAMP3 LAP3 LGALS3BP LPAR6 LY6E MOV10 MVB12A MX1 NCOA7 NMI NUB1",
      "OAS1 OASL OGFR PARP12 PARP14 PARP9 PLSCR1 PNPT1 PROCR PSMA3 PSMB8 PSMB9 PSME1 PSME2 RIPK2 RNF31 RSAD2 RTP4 SAMD9 SAMD9L SELL SLC25A28 SP110 STAT2 TAP1 TDRD7 TENT5A TMEM140 TRAFD1",
      "TRIM14 TRIM21 TRIM25 TRIM26 TRIM5 TXNIP UBA7 UBE2L6 USP18 WARS1"), "\\s+")[[1]],
    GO_IFN_beta = strsplit(paste(
      "ACOD1 AIM2 BST2 CALM1 CAMK2A CAPN2 CDC34 DDX41 GPR146 HCN1",
      "HTRA2 IFI16 IFITM1 IFITM2 IFITM3 IFNAR1 IFNAR2 IFNB1 IKBKE IRF1",
      "IRGM MAVS MNDA NDUFA13 OAS1 PLSCR1 PNPT1 PYDC5 PYHIN1 SHFL",
      "STAT1 STING1 TRAF3IP3 TREX1 TRIM6 UBE2G2 UBE2K XAF1"), "\\s+")[[1]])
  stopifnot(identical(lengths(genesets), c(SaulSenMayo = 124L, Hallmark_IFN_alpha = 97L, GO_IFN_beta = 38L)),
    identical(names(genesets), names(score.columns)), identical(names(score.labels), names(score.columns)),
    identical(names(gene.sources), names(score.columns)), identical(names(source.dates), names(score.columns)))

  match.genesets <- function(sets, rna.genes, minimum = 1L) {
    if (!length(rna.genes) || anyNA(rna.genes) || anyDuplicated(rna.genes))
      stop("RNA/data feature names must be present and unique.")
    if (!length(sets) || is.null(names(sets)) || anyDuplicated(names(sets))) stop("Invalid gene sets.")
    features <- list(); coverage <- detail <- list()
    for (nm in names(sets)) {
      genes <- unique(trimws(sets[[nm]]))
      if (!length(genes) || anyNA(genes) || any(!nzchar(genes))) stop("Invalid genes in ", nm)
      eligible <- setdiff(genes, "MFF")
      matched <- intersect(eligible, rna.genes)
      if (length(matched) < minimum) stop("Insufficient matched human genes for ", nm, ": ", length(matched))
      features[[nm]] <- matched
      coverage[[nm]] <- data.frame(score = nm, n_unique_input = length(genes),
        n_excluded_MFF = sum(genes == "MFF"), n_matched = length(matched),
        n_missing = length(setdiff(eligible, rna.genes)),
        coverage = if (length(eligible)) length(matched) / length(eligible) else NA_real_)
      detail[[nm]] <- data.frame(score = nm, gene = genes, status = ifelse(genes == "MFF",
        "excluded_MFF", ifelse(genes %in% rna.genes, "matched", "not_in_RNA_data")))
    }
    list(features = features, coverage = do.call(rbind, coverage), detail = do.call(rbind, detail))
  }
  prepare.cells <- function(meta, columns, ids = rownames(meta)) {
    required <- c("Annotation", "sample", "tech", "age_pop", "Country", unname(columns))
    if (!all(required %in% names(meta))) stop("Missing metadata: ", paste(setdiff(required, names(meta)), collapse = ", "))
    if (!nrow(meta) || !identical(rownames(meta), ids) || anyNA(ids) || anyDuplicated(ids) || any(ids == ""))
      stop("Invalid cell IDs or metadata/cell-order mismatch.")
    for (nm in c("Annotation", "sample", "tech", "age_pop", "Country"))
      if (anyNA(meta[[nm]]) || any(!nzchar(trimws(as.character(meta[[nm]]))))) stop("Missing metadata: ", nm)
    if (!all(meta$tech %in% c("scRNA", "snRNA")) || !all(meta$age_pop %in% c("young_pop", "old_pop")))
      stop("Unexpected technology or age_pop labels.")
    d <- data.frame(cell_id = ids, donor = as.character(meta$sample), celltype = as.character(meta$Annotation),
      tech = as.character(meta$tech), Country = as.character(meta$Country),
      group = ifelse(meta$age_pop == "young_pop", "Young", "Old"), stringsAsFactors = FALSE)
    for (nm in unname(columns)) {
      if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]]))) stop("Invalid/nonfinite score: ", nm)
      d[[nm]] <- meta[[nm]]  # Negative and zero module scores are retained.
    }
    info <- unique(d[, c("donor", "group", "Country")])
    if (anyDuplicated(info$donor)) stop("Inconsistent Young/Old or Country within donor.")
    d
  }
  summarize.scores <- function(cells, columns, types, scope.names) {
    statistics <- list()
    donor <- data.frame(scope = character(), score = character(), celltype = character(), donor = character(),
      group = character(), Country = character(), n_profiles = integer(), n_scRNA = integer(),
      n_snRNA = integer(), mean_score = numeric(), stringsAsFactors = FALSE)
    for (scope in scope.names) {
      b <- if (scope == "combined") cells else cells[cells$tech == scope, , drop = FALSE]
      for (ct in types) for (nm in names(columns)) {
        d <- b[b$celltype == ct, , drop = FALSE]; d$score <- d[[columns[[nm]]]]
        ds <- donor[FALSE, , drop = FALSE]
        if (nrow(d)) for (idx in split(seq_len(nrow(d)), d$donor)) {
          z <- d[idx, , drop = FALSE]
          ds <- rbind(ds, data.frame(scope = scope, score = nm, celltype = ct, donor = z$donor[1],
            group = z$group[1], Country = z$Country[1], n_profiles = nrow(z),
            n_scRNA = sum(z$tech == "scRNA"), n_snRNA = sum(z$tech == "snRNA"), mean_score = mean(z$score)))
        }
        donor <- rbind(donor, ds)
        for (unit in c("cell", "donor")) {
          z <- if (unit == "cell") data.frame(group = d$group, value = d$score) else
            data.frame(group = ds$group, value = ds$mean_score)
          young <- z$value[z$group == "Young"]; old <- z$value[z$group == "Old"]
          my <- if (length(young)) mean(young) else NA_real_; mo <- if (length(old)) mean(old) else NA_real_
          reason <- if (!length(young) || !length(old)) "missing_group" else
            if (min(length(young), length(old)) < 2L) "too_few_observations" else
            if (length(unique(c(young, old))) < 2L) "constant_score" else ""
          p <- W <- NA_real_
          if (reason == "") {
            w <- wilcox.test(young, old, alternative = "two.sided", paired = FALSE,
              exact = FALSE, correct = TRUE)
            p <- w$p.value; W <- unname(w$statistic)
            if (!is.finite(p)) reason <- "nonfinite_p"
          }
          statistics[[length(statistics) + 1L]] <- data.frame(scope = scope, score = nm, celltype = ct, unit = unit,
            n_Young = length(young), n_Old = length(old),
            n_profiles_Young = sum(d$group == "Young"), n_profiles_Old = sum(d$group == "Old"),
            n_donors_Young = sum(ds$group == "Young"), n_donors_Old = sum(ds$group == "Old"),
            mean_Young = my, mean_Old = mo, difference_Old_minus_Young = mo - my,
            W = W, p_value = p, p_adj_BH = NA_real_, reason = reason)
        }
      }
    }
    s <- do.call(rbind, statistics); rownames(s) <- rownames(donor) <- NULL
    for (scope in scope.names) for (unit in c("cell", "donor")) {
      idx <- which(s$scope == scope & s$unit == unit)
      s$p_adj_BH[idx] <- p.adjust(s$p_value[idx], method = "BH")
    }
    list(statistics = s, donor = donor)
  }
  write.checked.csv <- function(d, file) {
    token <- ".__HLMA_SCORE_MISSING__."
    while (any(vapply(d, function(v) any(as.character(v) == token, na.rm = TRUE), logical(1)))) token <- paste0(token, "_")
    write.con <- if (grepl("[.]gz$", file)) gzfile(file, "wt", encoding = "UTF-8") else file(file, "wt", encoding = "UTF-8")
    tryCatch(write.csv(d, write.con, row.names = FALSE, na = token), finally = close(write.con))
    classes <- vapply(d, function(v) if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else
      if (is.logical(v)) "logical" else "character", character(1))
    read.con <- if (grepl("[.]gz$", file)) gzfile(file, "rt", encoding = "UTF-8") else file(file, "rt", encoding = "UTF-8")
    restored <- tryCatch(read.csv(read.con, check.names = FALSE, colClasses = classes, na.strings = token,
      stringsAsFactors = FALSE), finally = close(read.con))
    stopifnot(identical(dim(restored), dim(d)), identical(names(restored), names(d)))
    for (nm in names(d)) {
      if (is.numeric(d[[nm]])) stopifnot(isTRUE(all.equal(restored[[nm]], d[[nm]], tolerance = 1e-10))) else
        stopifnot(identical(restored[[nm]], d[[nm]]))
    }
    invisible(file)
  }

  # 3. Load the whole muscle object. If present, the MFF-annotated object is preserved.
  candidates <- if (input.object == "auto") c("hlma.mff", "hlma.umap", "hlma", "hlma.rds") else input.object
  available <- candidates[vapply(candidates, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)]
  selected <- if (length(available)) available[1] else NA_character_
  if (input.object != "auto" && is.na(selected)) stop("Global input object does not exist: ", input.object)
  x <- if (!is.na(selected)) get(selected, envir = .GlobalEnv) else {
    if (!file.exists(source.rds)) stop("RDS missing: ", source.rds)
    readRDS(source.rds)
  }
  loaded.from <- if (!is.na(selected)) paste0("global:", selected) else source.rds
  stopifnot(inherits(x, "Seurat"), "RNA" %in% Assays(x), "data" %in% Layers(x[["RNA"]]))
  ids <- colnames(x); original.meta <- x[[]]
  DefaultAssay(x) <- "RNA"
  rna <- GetAssayData(x, assay = "RNA", layer = "data")
  stopifnot(identical(colnames(rna), ids))
  matched <- match.genesets(genesets, rownames(rna), min.matched.genes)
  pool <- setdiff(rownames(rna), "MFF"); rm(rna)
  stopifnot(identical(names(matched$features), names(score.columns)), !"MFF" %in% pool,
    all(unlist(matched$features, use.names = FALSE) %in% pool))
  print(matched$coverage)
  if (any(matched$coverage$coverage < low.coverage.warning))
    warning("Low human gene-set coverage: see 01.gene_set_coverage.csv before interpreting scores.")
  # Validate all metadata before the potentially expensive scoring step.
  preview <- original.meta; preview[, unname(score.columns)] <- 0
  invisible(prepare.cells(preview, score.columns, ids)); rm(preview)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  matched$coverage$source <- unname(gene.sources[matched$coverage$score])
  matched$coverage$retrieved <- unname(source.dates[matched$coverage$score])
  write.checked.csv(matched$coverage, file.path(out, "01.gene_set_coverage.csv"))
  write.checked.csv(matched$detail, file.path(out, "02.gene_set_members.csv"))

  # 4. Compute ONCE using all original profiles, before Young/Old or tech splitting.
  # Same nbin/ctrl/seed as mouse 11.1; MFF excluded from signatures AND control pool.
  # Uses existing RNA/data as processed expression; no NormalizeData() or counts required.
  # Source X's precise transformation remains unverified; this is a data-only adaptation.
  prefix <- "HLMA_score34_TMP"
  while (any(startsWith(names(x[[]]), prefix))) prefix <- paste0(prefix, "_")
  cat("Scoring all ", length(ids), " profiles from ", loaded.from, "\n", sep = "")
  x <- AddModuleScore(x, features = unname(matched$features), pool = pool,
    nbin = nbin, ctrl = ctrl, assay = "RNA", name = prefix, seed = score.seed, search = FALSE, slot = "data")
  for (i in seq_along(score.columns)) {
    temporary <- paste0(prefix, i)
    values <- x[[]][ids, temporary]
    if (length(values) != length(ids) || any(!is.finite(values))) stop("Nonfinite/missing AddModuleScore.")
    x[[unname(score.columns[i])]] <- setNames(values, ids)
    x[[temporary]] <- NULL
  }
  stopifnot(identical(colnames(x), ids))
  retained <- setdiff(names(original.meta), unname(score.columns))
  stopifnot(identical(x[[]][, retained, drop = FALSE], original.meta[, retained, drop = FALSE]))
  cells <- prepare.cells(x[[]], score.columns, ids)
  celltypes <- names(sort(table(cells$celltype), decreasing = TRUE)); n.rows <- ceiling(length(celltypes) / n.columns)
  summary <- summarize.scores(cells, score.columns, celltypes, scopes)
  statistics <- summary$statistics; donor <- summary$donor
  write.checked.csv(statistics[statistics$unit == "cell", ], file.path(out, "03.score.Young_Old.cell_statistics.csv"))
  write.checked.csv(statistics[statistics$unit == "donor", ], file.path(out, "04.score.Young_Old.donor_statistics.csv"))
  write.checked.csv(donor, file.path(out, "05.score.donor_means.csv"))
  write.checked.csv(cells, file.path(out, "06.score.profile_scores.csv.gz"))

  # 5. Signed-score Young/Old violins, 5 columns x 3 rows for the current 15 types.
  format.p <- function(p) {
    if (!is.finite(p)) return("p = NA")
    if (p == 0) return("p < 2.2e-308")
    paste0("p = ", if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  save.png <- function(p, file, width, height, res = plot.res) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(file, width = width, height = height, res = res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
  }
  caption <- "Cell-level Wilcoxon; unadjusted p; donor dependence, technology and Country are not adjusted."
  pngs <- character()
  for (nm in names(score.columns)) for (scope in scopes) {
    b <- if (scope == "combined") cells else cells[cells$tech == scope, , drop = FALSE]
    folder <- file.path(out, scope, nm)
    s <- statistics[statistics$scope == scope & statistics$score == nm & statistics$unit == "cell", ]
    plots <- list()
    for (ct in celltypes) {
      d <- b[b$celltype == ct, , drop = FALSE]; d$score <- d[[score.columns[[nm]]]]
      d$group <- factor(d$group, levels = group.order)
      stat <- s[s$celltype == ct, ]
      regular <- d[d$group %in% group.order[vapply(group.order, function(g) {
        v <- d$score[d$group == g]; length(v) >= 2L && length(unique(v)) > 1L
      }, logical(1))], , drop = FALSE]
      special <- d[!d$group %in% regular$group, , drop = FALSE]
      bounds <- if (nrow(d)) range(d$score) else c(-.1, .1)
      span <- max(diff(bounds), .05); bracket <- bounds[2] + .12 * span
      axis.labels <- setNames(paste0(group.order, "\n(", c(stat$n_profiles_Young, stat$n_profiles_Old), ")"), group.order)
      p <- ggplot(d, aes(group, score)) + geom_blank() +
        geom_violin(data = regular, aes(fill = group), trim = TRUE, scale = "width", width = .8) +
        geom_boxplot(data = regular, width = .15, fill = "white", colour = "black",
          outlier.shape = NA, linewidth = .4)
      for (g in group.order) {
        z <- special[special$group == g, , drop = FALSE]
        if (nrow(z) == 1L) p <- p + geom_point(data = z, size = 2)
        if (nrow(z) > 1L) p <- p + geom_errorbar(data = z[1, , drop = FALSE],
          aes(ymin = score, ymax = score), width = .45, linewidth = .6)
      }
      p <- p + scale_fill_manual(values = group.colors, guide = "none") +
        scale_x_discrete(limits = group.order, labels = axis.labels, drop = FALSE) +
        annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .4) +
        annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket - .04 * span, yend = bracket, linewidth = .4) +
        annotate("text", x = 1.5, y = bounds[2] + .28 * span, label = format.p(stat$p_value), size = 10) +
        coord_cartesian(ylim = c(bounds[1] - .08 * span, bounds[2] + .45 * span)) +
        labs(title = paste(strwrap(paste(scope.labels[scope], ct, sep = " | "), width = 40), collapse = "\n"),
          x = NULL, y = paste0(score.labels[nm], " score")) +
        theme_bw(base_size = 13) + theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 13), axis.text.y = element_text(size = 30),
          axis.text.x = element_text(size = 18), panel.grid.minor = element_blank())
      plots[[ct]] <- p
      if (save.individual.violins) {
        file <- file.path(folder, "violin", paste0(gsub("[^A-Za-z0-9_-]+", "_", ct), ".Young_Old.violin.png"))
        save.png(p, file, panel.width, panel.height); pngs <- c(pngs, file)
      }
    }
    overview <- wrap_plots(plots, ncol = n.columns, nrow = n.rows) +
      plot_annotation(title = score.labels[nm], caption = caption,
        theme = theme(plot.title = element_text(size = 22, face = "bold"), plot.caption = element_text(size = 10, hjust = 0)))
    file <- file.path(folder, paste0(nm, ".Young_Old.violin.overview.png"))
    save.png(overview, file, panel.width * n.columns, panel.height * n.rows); pngs <- c(pngs, file)

    # 6. Arithmetic score difference. Never exponentiate or compute log2FC of scores.
    s <- s[order(is.na(s$difference_Old_minus_Young), s$difference_Old_minus_Young, s$celltype), ]
    s$celltype_axis <- factor(s$celltype, levels = rev(s$celltype))
    s$direction <- ifelse(s$difference_Old_minus_Young > 0, "Old higher",
      ifelse(s$difference_Old_minus_Young < 0, "Old lower", "No change"))
    s$significance <- ifelse(!is.finite(s$p_value), "Not tested", ifelse(s$p_value < .05, "p < 0.05", "p >= 0.05"))
    delta <- statistics$difference_Old_minus_Young[statistics$score == nm & statistics$unit == "cell"]
    limit <- max(.05, max(c(0, abs(delta[is.finite(delta)]))) * 1.15)
    colors <- c("Old higher" = "#D73027", "Old lower" = "#000000", "No change" = "#777777")
    available.delta <- s[is.finite(s$difference_Old_minus_Young), ]; absent <- s[!is.finite(s$difference_Old_minus_Young), ]
    p <- ggplot(s, aes(y = celltype_axis)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "#999999") +
      geom_point(data = available.delta, aes(x = difference_Old_minus_Young, colour = direction, shape = significance), size = 4.5, stroke = 1.2) +
      geom_text(data = absent, x = limit * .95, label = "NA", hjust = 1, colour = "#888888") +
      scale_colour_manual(values = colors, name = NULL) +
      scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1, "Not tested" = 4),
        limits = c("p < 0.05", "p >= 0.05", "Not tested"), drop = FALSE, name = "Unadjusted cell p") +
      scale_x_continuous(limits = c(-limit, limit)) + scale_y_discrete(drop = FALSE) +
      labs(title = paste(scope.labels[scope], score.labels[nm], sep = " | "),
        x = "Cell-mean score difference (Old - Young)", y = NULL, caption = caption) +
      theme_classic(base_size = 14) + theme(plot.title = element_text(size = 22, face = "bold"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 16),
        legend.position = "bottom", legend.text = element_text(size = 11), plot.caption = element_text(size = 9, hjust = 0),
        plot.margin = margin(15, 25, 12, 12))
    file <- file.path(folder, paste0(nm, ".Young_Old.score_difference.png"))
    save.png(p, file, 3150L, round(max(5.5, 2 + length(celltypes) * .48) * 300), res = 300L); pngs <- c(pngs, file)
  }
  expected.pngs <- length(score.columns) * length(scopes) *
    (2L + if (save.individual.violins) length(celltypes) else 0L)
  stopifnot(length(pngs) == expected.pngs, all(file.exists(pngs)), all(file.info(pngs)$size > 8L))

  # 7. Keep scores and provenance in memory as global hlma.scores.
  # Retain this provenance key for compatibility with existing correlation scripts.
  x@misc$HLMA_two_scores <- list(method = "Seurat::AddModuleScore", assay = "RNA", layer = "data",
    scope = "all original profiles pooled; same scores sliced for combined/scRNA/snRNA figures",
    sources = gene.sources, source_date = source.dates, genesets = genesets, matched_features = matched$features,
    coverage = matched$coverage, pool_excludes = "MFF", nbin = nbin, ctrl = ctrl, seed = score.seed,
    score_columns = score.columns, input = loaded.from, input_scale = "stored processed RNA/data; exact transformation unverified",
    timestamp = as.character(Sys.time()))
  assign("hlma.scores", x, envir = .GlobalEnv)
  writeLines(c("Human Young/Old comparison of exactly THREE scores; violin plots include boxplots.",
    "White boxes show the median and interquartile range inside each variable-group violin. Outlier markers are hidden; all scores remain in the violin and statistics.",
    "Constant/singleton groups retain line/point displays. No age regression or MFF-score correlation.",
    "Human gene sets: official MSigDB snapshots, 124 SaulSenMayo and 97 Hallmark IFN-alpha symbols (2026-10-06), and 38 GO IFN-beta symbols (2026-10-07).",
    "GO IFN-beta: GOBP_RESPONSE_TO_INTERFERON_BETA, GO:0035456; human members, not the mouse .Mm.grp file.",
    "Gene coverage/missing members are exported. MFF is excluded from signatures and the control pool as in mouse 11.1.",
    "AddModuleScore uses stored RNA/data: mean signature expression minus expression-matched control mean.",
    "Existing processed expression is assumed suitable; its exact source transformation has not been independently confirmed.",
    "Scores are computed once across ALL original muscle profiles with nbin=24, ctrl=100, seed=1, search=FALSE.",
    "Young/Old and technology subsets are taken AFTER scoring; they use a common reference, not independent scores.",
    "Negative/zero scores are valid and retained. A score of zero is not a senescent-cell cutoff.",
    "Effect = arithmetic mean score Old minus Young; no inverse transformation, ratio or log2FC.",
    "All original profiles and annotations are retained. No normalization, embedding, counts or pseudobulk is rerun.",
    "18 main PNGs: 3 scores x 3 scopes x (violin overview + score difference dotplot).",
    "Violin overview: 5 columns x 3 rows for 15 cell types, 4500 x 2550 pixels at 160 dpi.",
    "Violin and dotplot p: two-sided unpaired CELL Wilcoxon, exact=FALSE, correct=TRUE; exploratory and unadjusted.",
    "Cell-level p does not handle same-donor dependence or technical/Country differences.",
    "Donor means and separate donor-level Wilcoxon statistics are saved as CSV; no donor boxplots.",
    "Donor summaries use equally weighted observed donors. Absent donor/type combinations are not zero scores.",
    "BH is calculated over all three scores x all cell types separately within each scope and unit.",
    "At least 2 observations per group; missing/rare/constant groups get NA p with a reason.",
    "Combined/scRNA/snRNA scopes overlap and are not independent replications; scores are not batch-adjusted.",
    "Score magnitudes from different signatures need not be comparable.",
    paste0("Input: ", loaded.from),
    "Score metadata and provenance remain in memory as hlma.scores; score CSVs are readback-verified.",
    "The compatibility provenance key hlma.scores@misc$HLMA_two_scores contains all three signatures and matched feature lists.",
    "Global hlma.scores contains the scored object; input global objects retain their previous metadata.",
    paste0("Output: ", out)), file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("\nHLMA_THREE_SCORES_YOUNG_OLD_COMPLETE\nProfiles: ", nrow(cells), "; PNGs: ", length(pngs),
    "\nScore columns: ", paste(score.columns, collapse = ", "), "\n", out, "\n", sep = "")
})
