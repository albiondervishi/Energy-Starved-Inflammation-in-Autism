
##############################################
## τ-Axis – GSE18123 (Dec 2025, IL10 extremes)
## Comparisons:
##   1) Control vs ASD_Low_IL10
##   2) Control vs ASD_High_IL10
##############################################

## 0. Setup --------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/Volumes/Transcend/τ-Axis")

## 1. Libraries -----------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(matrixStats)
library(AnnotationDbi)
library(hgu133plus2.db)
library(edgeR)
library(limma)
library(effsize)
library(pROC)
library(writexl)

## 2. Load raw files ------------------------------------------------------
exprSet   <- read.csv("GSE18123_expression.csv",
                      row.names = 1, check.names = FALSE)
pheno_raw <- read.csv("pheno.csv",
                      row.names = 1, check.names = FALSE)

cat("Raw expression:", nrow(exprSet), "probes x", ncol(exprSet), "samples\n")

## 3. Build pheno: GSM + ASD vs Control -----------------------------------
pheno <- data.frame(
  GSM   = trimws(as.character(pheno_raw$geo_accession)),
  Group = ifelse(
    grepl("autism|ASD|PDD|Asperger",
          pheno_raw$`diagnosis:ch1`, ignore.case = TRUE),
    "ASD", "Control"
  ),
  stringsAsFactors = FALSE
) %>%
  filter(GSM != "", !is.na(GSM))

cat("Pheno parsed from geo_accession →", nrow(pheno), "samples\n")
head(pheno)

## 4. Align expression columns to GSM IDs ---------------------------------
colnames(exprSet) <- trimws(as.character(pheno_raw$geo_accession))

common <- intersect(colnames(exprSet), pheno$GSM)
cat("Overlap:", length(common), "samples\n")

expr_mat <- as.matrix(exprSet)[, common, drop = FALSE]
pheno    <- pheno %>%
  filter(GSM %in% common) %>%
  arrange(match(GSM, common))

group <- factor(pheno$Group, levels = c("Control", "ASD"))

cat("SUCCESS →", ncol(expr_mat), "samples aligned!\n")
cat("Controls:", sum(group == "Control"), "| ASD:", sum(group == "ASD"), "\n")

## 5. Probe → ENSG mapping and collapse -----------------------------------
cat("Mapping probes to ENSG via hgu133plus2.db...\n")

probe2ensg <- AnnotationDbi::select(
  hgu133plus2.db,
  keys    = rownames(expr_mat),
  columns = "ENSEMBL",
  keytype = "PROBEID"
) %>%
  filter(!is.na(ENSEMBL))

cat("Mapped", nrow(probe2ensg), "probes to ENSG IDs\n")

expr_mapped <- expr_mat[probe2ensg$PROBEID, , drop = FALSE]
rownames(expr_mapped) <- probe2ensg$ENSEMBL

cat("Collapsing probes → one row per ENSG (max mean expression)...\n")

ensg_ids  <- rownames(expr_mapped)
split_idx <- split(seq_len(nrow(expr_mapped)), ensg_ids)

collapsed_list <- lapply(split_idx, function(idx) {
  submat <- expr_mapped[idx, , drop = FALSE]
  if (nrow(submat) == 1L) return(submat)
  means <- rowMeans(submat, na.rm = TRUE)
  submat[which.max(means), , drop = FALSE]
})

expr_ensg <- do.call(rbind, collapsed_list)
rownames(expr_ensg) <- names(split_idx)

cat("→ Collapsed to", nrow(expr_ensg), "unique ENSGs\n")

## 6. TMM + logCPM --------------------------------------------------------
cat("Running edgeR TMM + logCPM...\n")
dge <- edgeR::DGEList(counts = expr_ensg, group = group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

cat("logCPM matrix:", nrow(logcpm), "genes x", ncol(logcpm), "samples\n")

## 7. Load Gene_TAO.xlsx and intersect genes ------------------------------
gene_tao <- readxl::read_excel("Gene_TAO.xlsx", sheet = "Sheet1") %>%
  mutate(ENSG = trimws(ENSG)) %>%
  distinct(ENSG, .keep_all = TRUE) %>%
  dplyr::select(ENSG, GeneSymbol, Modules, Enzyme, Category)

cat("Loaded", nrow(gene_tao), "genes from Gene_TAO.xlsx\n")

common_genes <- intersect(rownames(logcpm), gene_tao$ENSG)
logcpm       <- logcpm[common_genes, , drop = FALSE]
gene_tao     <- gene_tao %>% filter(ENSG %in% common_genes)

cat("Final τ-Axis ASD gene set:", nrow(logcpm), "genes\n")

## 8. Build τ-modules from Gene_TAO.xlsx ----------------------------------
tau_modules <- gene_tao %>%
  filter(!is.na(Modules)) %>%
  separate_rows(Modules, sep = ";") %>%
  group_by(Modules) %>%
  summarise(ENSG_list = list(unique(ENSG)), .groups = "drop") %>%
  tibble::deframe()

cat("Modules detected:", length(tau_modules), "\n")
print(names(tau_modules))

## 9. Module scoring (per sample) -----------------------------------------
score_module <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0L) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[g, , drop = FALSE])
}

# scores_mat: samples x modules (cols = modules, rows = GSMs)
scores_mat <- sapply(names(tau_modules),
                     function(m) score_module(logcpm, tau_modules[[m]]))

# Ensure samples in rows and modules in columns
scores_df <- as.data.frame(scores_mat) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(Sample = as.character(Sample)) %>%
  left_join(pheno %>% mutate(GSM = as.character(GSM)),
            by = c("Sample" = "GSM"))

cat("scores_df:", nrow(scores_df), "samples x", ncol(scores_df), "columns\n")
print(head(scores_df[, c("Sample", "Group")]))

# module names for τ
tau_module_names <- names(tau_modules)

## 10. τ-index (tau_raw, tau_z) over ALL samples --------------------------
missing <- setdiff(tau_module_names, colnames(scores_df))
if (length(missing) > 0) {
  cat("⚠️ WARNING: These modules are in tau_modules but NOT in scores_df:\n")
  print(missing)
}
tau_module_names <- intersect(tau_module_names, colnames(scores_df))
cat("Using", length(tau_module_names), "modules for τ calculation:\n")
print(tau_module_names)

module_mat <- as.matrix(scores_df[, tau_module_names, drop = FALSE])

tau_weights <- setNames(rep(1, length(tau_module_names)), tau_module_names)

tau_raw <- as.numeric(module_mat %*% tau_weights[tau_module_names])
tau_z   <- as.numeric(scale(tau_raw))

scores_df$tau_raw <- tau_raw
scores_df$tau_z   <- tau_z

## 11. Define IL-10 groups within ASD (LOW vs HIGH) -----------------------
if (!"IL10" %in% colnames(scores_df)) {
  stop("Module 'IL10' not found in scores_df. Check 'Modules' column in Gene_TAO.xlsx.")
}

asd_scores <- scores_df %>%
  filter(Group == "ASD" & !is.na(IL10))

q33 <- quantile(asd_scores$IL10, 1/3, na.rm = TRUE)
q66 <- quantile(asd_scores$IL10, 2/3, na.rm = TRUE)

cat("IL10 tertiles in ASD → 33%:", q33, " | 66%:", q66, "\n")

scores_df$ImmuneState <- NA_character_
scores_df$ImmuneState[scores_df$Group == "Control"] <- "Control"

scores_df$ImmuneState[
  scores_df$Group == "ASD" & scores_df$IL10 <= q33
] <- "ASD_Low_IL10"

scores_df$ImmuneState[
  scores_df$Group == "ASD" & scores_df$IL10 >= q66
] <- "ASD_High_IL10"

cat("ImmuneState distribution:\n")
print(table(scores_df$ImmuneState, useNA = "ifany"))

## 12. FINAL DF for analysis (3 states, but we use 2 comparisons) ----------
final_df <- scores_df %>%
  filter(ImmuneState %in% c("Control", "ASD_Low_IL10", "ASD_High_IL10")) %>%
  mutate(
    ImmuneState = factor(
      ImmuneState,
      levels = c("Control", "ASD_Low_IL10", "ASD_High_IL10")
    )
  )

cat("Final immune states:\n")
print(table(final_df$ImmuneState))

## 13. Gene-level comparison function -------------------------------------
make_gene_table <- function(g1, g2, label) {
  cat("Building gene table →", label, "\n")
  
  sub <- final_df %>% dplyr::filter(ImmuneState %in% c(g1, g2))
  grp <- droplevels(stats::relevel(sub$ImmuneState, ref = g1))
  mat <- logcpm[, sub$Sample, drop = FALSE]
  
  # limma model
  design <- model.matrix(~ grp)
  fit    <- lmFit(mat, design) %>% eBayes()
  
  limma_res <- topTable(fit, coef = 2, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("ENSG")
  
  # Robust p-value naming: P.Value or P
  if ("P.Value" %in% colnames(limma_res)) {
    limma_res <- limma_res %>%
      dplyr::rename(p_limma_raw = P.Value)
  } else if ("P" %in% colnames(limma_res)) {
    limma_res <- limma_res %>%
      dplyr::rename(p_limma_raw = P)
  } else {
    stop("No p-value column (P.Value or P) found in limma_res.")
  }
  
  # Adjusted p-value
  if ("adj.P.Val" %in% colnames(limma_res)) {
    limma_res <- limma_res %>%
      dplyr::rename(p_limma_adj = adj.P.Val)
  } else {
    limma_res$p_limma_adj <- p.adjust(limma_res$p_limma_raw, method = "BH")
  }
  
  # Means + stats (Wilcoxon + Cliff)
  means_g1 <- rowMeans(mat[, grp == g1, drop = FALSE])
  means_g2 <- rowMeans(mat[, grp == g2, drop = FALSE])
  
  fast_stats <- t(apply(mat, 1, function(x) {
    a <- x[grp == g1]; b <- x[grp == g2]
    c(
      suppressWarnings(wilcox.test(a, b)$p.value),
      suppressWarnings(effsize::cliff.delta(b, a)$estimate)
    )
  }))
  
  stats_df <- data.frame(
    ENSG         = rownames(mat),
    Mean_g1      = means_g1,
    Mean_g2      = means_g2,
    Delta        = means_g2 - means_g1,
    p_raw_wilcox = fast_stats[, 1],
    Cliff        = fast_stats[, 2],
    stringsAsFactors = FALSE
  ) %>% tibble::as_tibble()
  
  full <- stats_df %>%
    left_join(
      limma_res %>%
        dplyr::select(ENSG, logFC, t, p_limma_raw, p_limma_adj),
      by = "ENSG"
    ) %>%
    dplyr::rename(
      T_value     = t,
      p_raw_limma = p_limma_raw,
      p_adj_limma = p_limma_adj
    ) %>%
    left_join(gene_tao, by = "ENSG") %>%
    mutate(
      Signif = symnum(
        p_adj_limma,
        c(0, 0.001, 0.01, 0.05, 1),
        c("***", "**", "*", " ")
      )
    ) %>%
    dplyr::select(
      ENSG, GeneSymbol, Modules, Enzyme, Category,
      Mean_g1, Mean_g2, Delta, logFC, T_value,
      p_raw_limma, p_adj_limma, p_raw_wilcox, Cliff, Signif
    ) %>%
    arrange(p_adj_limma)
  
  full
}

## 14. Run ONLY the two comparisons you want -------------------------------
gene_ctrl_low  <- make_gene_table("Control", "ASD_Low_IL10",
                                  "Control_vs_ASD_Low_IL10")
gene_ctrl_high <- make_gene_table("Control", "ASD_High_IL10",
                                  "Control_vs_ASD_High_IL10")

## 15. Module-level + τ comparison ----------------------------------------
mod_comp <- function(g1, g2) {
  sub <- final_df %>% dplyr::filter(ImmuneState %in% c(g1, g2))
  
  modules_tbl <- sub %>%
    pivot_longer(
      cols      = all_of(tau_module_names),
      names_to  = "Module",
      values_to = "Score"
    ) %>%
    group_by(Module) %>%
    summarise(
      Mean_g1 = mean(Score[ImmuneState == g1], na.rm = TRUE),
      Mean_g2 = mean(Score[ImmuneState == g2], na.rm = TRUE),
      Delta   = Mean_g2 - Mean_g1,
      p_raw   = wilcox.test(Score ~ ImmuneState)$p.value,
      Cliff   = effsize::cliff.delta(
        Score[ImmuneState == g2],
        Score[ImmuneState == g1]
      )$estimate,
      .groups = "drop"
    )
  
  # τ_z as extra row
  a <- sub$tau_z[sub$ImmuneState == g1]
  b <- sub$tau_z[sub$ImmuneState == g2]
  
  tau_row <- data.frame(
    Module  = "Tau_index",
    Mean_g1 = mean(a, na.rm = TRUE),
    Mean_g2 = mean(b, na.rm = TRUE),
    Delta   = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE),
    p_raw   = wilcox.test(a, b)$p.value,
    Cliff   = effsize::cliff.delta(b, a)$estimate
  )
  
  bind_rows(modules_tbl, tau_row) %>%
    mutate(
      p_adj = p.adjust(p_raw, "BH"),
      Signif = symnum(
        p_adj,
        c(0, 0.001, 0.01, 0.05, 1),
        c("***", "**", "*", " ")
      )
    )
}

mod_ctrl_low  <- mod_comp("Control", "ASD_Low_IL10")
mod_ctrl_high <- mod_comp("Control", "ASD_High_IL10")

## 16. Sample-level τ-sheet (for plots, AUC etc.) --------------------------
sample_tau <- final_df %>%
  dplyr::select(
    Sample, Group, ImmuneState,
    dplyr::all_of(tau_module_names),
    tau_raw, tau_z
  )

## 17. Save everything to one Excel file -----------------------------------
dir.create("τ-Axis_Results", showWarnings = FALSE)

final_excel <- list(
  Ctrl_vs_LowIL10_Module  = mod_ctrl_low,
  Ctrl_vs_LowIL10_Gene    = gene_ctrl_low,
  Ctrl_vs_HighIL10_Module = mod_ctrl_high,
  Ctrl_vs_HighIL10_Gene   = gene_ctrl_high,
  Sample_TauAxis          = sample_tau
)

#writexl::write_xlsx(final_excel,
#                    "τ-Axis_Results/ASD_Control_IL10_Comparisons.xlsx")


#writexl::write_xlsx(final_excel, "τ-Axis_Results/ASD_Control_IL10_Comparisons_PSA.xlsx")

cat("✓ DONE: Results written to τ-Axis_Results/ASD_Control_IL10_Comparisons.xlsx\n")







library(ggplot2)
library(dplyr)
final_excel <- data.frame(final_excel)
df <- final_excel$Sample_TauAxis  # your exported sheet

ggplot(df, aes(x = tau_z, color = ImmuneState, fill = ImmuneState)) +
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

# Use exactly the same df as for your density plot
df <- final_excel$Sample_TauAxis  

# Ensure factors are correct
df$ImmuneState <- factor(
  df$ImmuneState,
  levels = c("Control", "ASD_Low_IL10", "ASD_High_IL10")
)

# ─────────────────────────────────────────────
# 1) AUC: Control vs ASD_Low_IL10
# ─────────────────────────────────────────────
df_low <- df %>%
  filter(ImmuneState %in% c("Control", "ASD_Low_IL10")) %>%
  mutate(Label = ifelse(ImmuneState == "ASD_Low_IL10", 1, 0))

roc_low <- roc(df_low$Label, df_low$tau_z)
auc_low <- auc(roc_low)

cat("AUC (Control vs ASD_Low_IL10):", round(auc_low, 3), "\n")

# Optional ROC plot
plot(roc_low,
     main = paste("ROC – Control vs ASD_Low_IL10 | AUC =", round(auc_low, 3)))


# ─────────────────────────────────────────────
# 2) AUC: Control vs ASD_High_IL10
# ─────────────────────────────────────────────
df_high <- df %>%
  filter(ImmuneState %in% c("Control", "ASD_High_IL10")) %>%
  mutate(Label = ifelse(ImmuneState == "ASD_High_IL10", 1, 0))

roc_high <- roc(df_high$Label, df_high$tau_z)
auc_high <- auc(roc_high)

cat("AUC (Control vs ASD_High_IL10):", round(auc_high, 3), "\n")

# Optional ROC plot
plot(roc_high,
     main = paste("ROC – Control vs ASD_High_IL10 | AUC =", round(auc_high, 3)))




cat("AUC (Control vs ASD_Low_IL10):", round(auc_low, 3), "\n")
cat("AUC (Control vs ASD_High_IL10):", round(auc_high, 3), "\n")


ci_low  <- ci.auc(roc_low)
ci_high <- ci.auc(roc_high)

cat("AUC Low IL10 95% CI:", round(ci_low[1],3), "-", round(ci_low[3],3), "\n")
cat("AUC High IL10 95% CI:", round(ci_high[1],3), "-", round(ci_high[3],3), "\n")


