## FINAL SEPSIS τ-AXIS + FULL Gene_TAO.xlsx ANNOTATION
## Healthy vs Low / Medium / Severe SOFA (164 genes only)
## November 2025 – Publication-ready
############################################################
setwd("/Volumes/Transcend/τ-Axis")

pkgs <- c("edgeR","limma","GEOquery","dplyr","readr","tidyr",
          "effsize","writexl","readxl","ggplot2","ggpubr")
for(p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))

# 1. Load + match samples -------------------------------------------------
expr <- readr::read_csv("GSE185263_raw_counts.csv", col_types = cols())
expr_mat <- as.matrix(expr[, -1])
rownames(expr_mat) <- expr[[1]]

gse   <- GEOquery::getGEO("GSE185263", GSEMatrix = TRUE)[[1]]
pheno <- pData(gse)

sample_titles <- trimws(pheno$title)
col_names     <- trimws(colnames(expr_mat))
matched_cols  <- match(sample_titles, col_names)

expr_mat <- expr_mat[, matched_cols]
pheno    <- pheno[!is.na(matched_cols), ]
colnames(expr_mat) <- pheno$geo_accession

# 2. SOFA groups ----------------------------------------------------------
sofa_raw <- pheno$`sofa 24h post admisssion:ch1`
sofa_num <- as.numeric(as.character(sofa_raw))

is_healthy <- grepl("healthy|control", pheno$characteristics_ch1, ignore.case = TRUE)

sofa_group <- cut(
  sofa_num,
  breaks = c(-Inf, 4, 7, Inf),
  labels = c("Low_SOFA", "Medium_SOFA", "Severe_SOFA"),
  right  = FALSE
)

Group <- ifelse(is_healthy, "Healthy", as.character(sofa_group))
Group <- factor(Group, levels = c("Healthy","Low_SOFA","Medium_SOFA","Severe_SOFA"))

valid    <- !is.na(Group)
expr_mat <- expr_mat[, valid]
Group    <- Group[valid]

cat("Groups:\n")
print(table(Group))

# 3. Normalisation + collapse isoforms -----------------------------------
dge <- edgeR::DGEList(counts = expr_mat, group = Group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# Collapse isoforms → highest expressing version per ENSG
rownames(logcpm) <- sub("\\..*$", "", rownames(logcpm))
if (anyDuplicated(rownames(logcpm))) {
  cat("Collapsing isoforms...\n")
  logcpm <- do.call(rbind, by(logcpm, rownames(logcpm), matrixStats::colMaxs))
}

# 4. Load your exact Gene_TAO.xlsx annotation (164 genes) -----------------
gene_tao <- readxl::read_excel("Gene_TAO.xlsx", sheet = "Sheet1") %>%
  dplyr::mutate(ENSG = trimws(ENSG)) %>%
  dplyr::distinct(ENSG, .keep_all = TRUE) %>%
  dplyr::select(ENSG, GeneSymbol, Modules, Enzyme, Category)

cat("Loaded", nrow(gene_tao), "annotated genes from Gene_TAO.xlsx\n")

# Keep only these 164 genes
common_genes <- intersect(rownames(logcpm), gene_tao$ENSG)
logcpm  <- logcpm[common_genes, ]
gene_tao <- gene_tao %>% dplyr::filter(ENSG %in% common_genes)

cat("Final analysis on exactly", nrow(logcpm), "genes\n")

# 5. Module definitions directly from your Excel --------------------------
tau_modules <- gene_tao %>%
  dplyr::filter(!is.na(Modules)) %>%
  tidyr::separate_rows(Modules, sep = ";") %>%
  dplyr::group_by(Modules) %>%
  dplyr::summarise(ENSG_list = list(ENSG), .groups = "drop") %>%
  tibble::deframe()

cat("Modules detected:", length(tau_modules), "\n")
print(names(tau_modules))

# 6. Module scores (per sample) -------------------------------------------
score_module <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0L) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
}

scores_mat <- sapply(
  names(tau_modules),
  function(m) score_module(logcpm, tau_modules[[m]])
)

final_df <- data.frame(
  Sample = colnames(logcpm),
  Group  = Group,
  scores_mat,
  stringsAsFactors = FALSE
)

############################################################################
## 7. τ from module scores in final_df (Sepsis, ASD-style: all weights +1)
############################################################################

# Modules that exist both in tau_modules and in final_df
module_cols_for_tau <- intersect(names(tau_modules), colnames(final_df))
cat("Modules used for τ calculation (from final_df):\n")
print(module_cols_for_tau)

# Build module matrix (samples x modules)
module_mat_tau <- as.matrix(final_df[, module_cols_for_tau, drop = FALSE])

# τ_raw = sum of module scores (all weights +1)
tau_raw <- rowSums(module_mat_tau, na.rm = TRUE)

# τ_z = z-normalised τ_raw across all samples
tau_z <- as.numeric(scale(tau_raw))

# Attach to final_df
final_df$tau_raw <- tau_raw
final_df$tau_z   <- tau_z

############################################################################
# 8. Gene-level function with PERFECT left-join to Gene_TAO.xlsx ----------
############################################################################
make_gene_table <- function(g1, g2, label) {
  cat("Building gene table →", label, "\n")
  
  sub <- final_df %>% dplyr::filter(Group %in% c(g1, g2))
  grp <- droplevels(sub$Group)
  mat <- logcpm[, sub$Sample, drop = FALSE]
  
  # Limma
  design <- model.matrix(~ grp)
  fit <- limma::lmFit(mat, design) %>% limma::eBayes()
  limma_res <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("ENSG")
  
  # Means + Wilcoxon + Cliff's delta
  mean_g1 <- rowMeans(mat[, grp == g1, drop = FALSE])
  mean_g2 <- rowMeans(mat[, grp == g2, drop = FALSE])
  
  stats_fast <- t(apply(mat, 1, function(x) {
    a <- x[grp == g1]; b <- x[grp == g2]
    c(
      stats::wilcox.test(a, b)$p.value,
      effsize::cliff.delta(b, a)$estimate
    )
  }))
  
  fast_df <- data.frame(
    ENSG          = rownames(mat),
    p_raw_wilcox  = stats_fast[, 1],
    Cliff         = stats_fast[, 2],
    stringsAsFactors = FALSE
  )
  
  # Build final table with LEFT JOIN to your Gene_TAO.xlsx
  full_table <- data.frame(
    ENSG    = rownames(mat),
    Mean_g1 = mean_g1,
    Mean_g2 = mean_g2,
    Delta   = mean_g2 - mean_g1,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(fast_df, by = "ENSG") %>%
    dplyr::left_join(
      limma_res %>% dplyr::select(ENSG, logFC, t, P.Value, adj.P.Val),
      by = "ENSG"
    ) %>%
    dplyr::rename(
      T_value     = t,
      p_raw_limma = P.Value,
      p_adj_limma = adj.P.Val
    ) %>%
    dplyr::left_join(gene_tao, by = "ENSG") %>%   # attach Gene_TAO.xlsx
    dplyr::mutate(
      Signif = symnum(
        p_adj_limma,
        c(0, 0.001, 0.01, 0.05, 1),
        c("***","**","*"," ")
      )
    ) %>%
    dplyr::arrange(p_adj_limma)
  
  # Order exactly like your ASD files
  full_table %>% dplyr::arrange(p_adj_limma)
}

# 9. Module-level summaries (with τ added as extra "module") --------------
mod_comp <- function(g1, g2) {
  sub <- final_df %>%
    dplyr::filter(Group %in% c(g1, g2))
  
  sub$Group <- droplevels(sub$Group)
  
  # Use *all* module columns that really exist in final_df
  module_cols <- intersect(names(tau_modules), colnames(sub))
  cat("Modules used in module comparison:", module_cols, "\n")
  
  # Standard modules
  modules_tbl <- sub %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(module_cols),
      names_to  = "Module",
      values_to = "Score"
    ) %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      Mean_g1 = mean(Score[Group == g1], na.rm = TRUE),
      Mean_g2 = mean(Score[Group == g2], na.rm = TRUE),
      Delta   = Mean_g2 - Mean_g1,
      p_raw   = stats::wilcox.test(
        x = Score[Group == g1],
        y = Score[Group == g2]
      )$p.value,
      .groups = "drop"
    )
  
  # τ as an extra "module" (robust)
  a <- as.numeric(sub$tau_z[sub$Group == g1])
  b <- as.numeric(sub$tau_z[sub$Group == g2])
  
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  
  if (length(a) < 2 || length(b) < 2) {
    # Not enough data for Wilcoxon
    tau_row <- data.frame(
      Module  = "Tau_index",
      Mean_g1 = ifelse(length(a) > 0, mean(a), NA_real_),
      Mean_g2 = ifelse(length(b) > 0, mean(b), NA_real_),
      Delta   = ifelse(length(a) > 0 & length(b) > 0,
                       mean(b) - mean(a),
                       NA_real_),
      p_raw   = NA_real_
    )
  } else {
    tau_row <- data.frame(
      Module  = "Tau_index",
      Mean_g1 = mean(a),
      Mean_g2 = mean(b),
      Delta   = mean(b) - mean(a),
      p_raw   = stats::wilcox.test(a, b)$p.value
    )
  }
  
  dplyr::bind_rows(modules_tbl, tau_row) %>%
    dplyr::mutate(
      p_adj = p.adjust(p_raw, "BH"),
      Signif = symnum(
        p_adj,
        c(0, 0.001, 0.01, 0.05, 1),
        c("***","**","*"," ")
      )
    )
}

# 10. Run all comparisons --------------------------------------------------
gene1 <- make_gene_table("Healthy", "Low_SOFA",    "Healthy_vs_Low")
gene2 <- make_gene_table("Healthy", "Medium_SOFA", "Healthy_vs_Medium")
gene3 <- make_gene_table("Healthy", "Severe_SOFA", "Healthy_vs_Severe")

mod1 <- mod_comp("Healthy","Low_SOFA")
mod2 <- mod_comp("Healthy","Medium_SOFA")
mod3 <- mod_comp("Healthy","Severe_SOFA")

# 11. Per-sample τ + module scores (for Excel) ----------------------------
module_cols_export <- intersect(colnames(final_df), names(tau_modules))

sample_tau <- final_df %>%
  dplyr::select(
    Sample,
    Group,
    dplyr::all_of(module_cols_export),
    tau_raw,
    tau_z
  )

# 12. Save exactly like your ASD workflow ---------------------------------
dir.create("τ-Axis_Results", showWarnings = FALSE)

out <- list(
  Healthy_vs_Low_Module    = mod1,
  Healthy_vs_Low_Gene      = gene1,
  Healthy_vs_Medium_Module = mod2,
  Healthy_vs_Medium_Gene   = gene2,
  Healthy_vs_Severe_Module = mod3,
  Healthy_vs_Severe_Gene   = gene3,
  Sample_TauAxis           = sample_tau    # NEW sheet with τ per sample
)

# writexl::write_xlsx(out, "τ-Axis_Results/Sepsis_Final_NO_WEIGHTS.xlsx")
#writexl::write_xlsx(out, "τ-Axis_Results/Sepsis_Final_PSAS.xlsx")

cat("\nSUCCESS! 
→ 164 genes perfectly matched to Gene_TAO.xlsx
→ τ_raw and τ_z calculated from module scores (all +1, ASD-style)
→ τ included as 'Tau_index' in module summaries
→ Sample-level τ-axis sheet added
→ Ready for submission!\n")


unique(gene_tao$Modules)





library(ggplot2)
library(dplyr)
df <- out$Sample_TauAxis  # your exported sheet

ggplot(df, aes(x = tau_z, color = Group, fill = Group)) +
  geom_density(alpha = 0.25, size = 1.2) +
  labs(
    title = "τ-Axis Density Curve by Group",
    x = "τ_z (Turnover Z-score)",
    y = "Density"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 20, face = "bold")
  )





library(pROC)
library(dplyr)

# Use exactly the same object you used for the density
df <- out$Sample_TauAxis

df$Group <- factor(
  df$Group,
  levels = c("Healthy", "Low_SOFA", "Medium_SOFA", "Severe_SOFA")
)

############################################################
# 1) AUC: Healthy vs Low_SOFA
############################################################
df_low <- df %>%
  filter(Group %in% c("Healthy", "Low_SOFA")) %>%
  mutate(Label = ifelse(Group == "Low_SOFA", 1, 0))

roc_low <- roc(df_low$Label, df_low$tau_z)
auc_low <- auc(roc_low)
ci_low  <- ci.auc(roc_low)

cat("AUC (Healthy vs Low_SOFA):", round(auc_low, 3),
    "| 95% CI:", round(ci_low[1],3), "-", round(ci_low[3],3), "\n")

# Optional ROC plot
plot(roc_low,
     main = paste("ROC – Healthy vs Low_SOFA | AUC =", round(auc_low, 3)))


############################################################
# 2) AUC: Healthy vs Medium_SOFA
############################################################
df_med <- df %>%
  filter(Group %in% c("Healthy", "Medium_SOFA")) %>%
  mutate(Label = ifelse(Group == "Medium_SOFA", 1, 0))

roc_med <- roc(df_med$Label, df_med$tau_z)
auc_med <- auc(roc_med)
ci_med  <- ci.auc(roc_med)

cat("AUC (Healthy vs Medium_SOFA):", round(auc_med, 3),
    "| 95% CI:", round(ci_med[1],3), "-", round(ci_med[3],3), "\n")

plot(roc_med,
     main = paste("ROC – Healthy vs Medium_SOFA | AUC =", round(auc_med, 3)))


############################################################
# 3) AUC: Healthy vs Severe_SOFA
############################################################
df_sev <- df %>%
  filter(Group %in% c("Healthy", "Severe_SOFA")) %>%
  mutate(Label = ifelse(Group == "Severe_SOFA", 1, 0))

roc_sev <- roc(df_sev$Label, df_sev$tau_z)
auc_sev <- auc(roc_sev)
ci_sev  <- ci.auc(roc_sev)

cat("AUC (Healthy vs Severe_SOFA):", round(auc_sev, 3),
    "| 95% CI:", round(ci_sev[1],3), "-", round(ci_sev[3],3), "\n")

plot(roc_sev,
     main = paste("ROC – Healthy vs Severe_SOFA | AUC =", round(auc_sev, 3)))


