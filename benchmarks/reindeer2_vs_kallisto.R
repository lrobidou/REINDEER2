# ------------------------------------
# Required Libraries
# ------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)

# ------------------------------------
# Command-line arguments
# ------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || !(args[1] %in% c("tpm", "raw", "seqc")) || !(args[2] %in% c("dense", "sparse"))) {
  stop("Usage: Rscript script.R [tpm|raw|seqc] [dense|sparse]")
}

abundance_type <- args[1]  # "tpm", "raw", or "seqc"
index_type <- args[2]      # "dense" or "sparse"

index_suffix <- if (index_type == "dense") "with dense index" else "without dense index"
index_filename <- index_type

# ------------------------------------
# Load GTF annotation
# ------------------------------------
gtf <- read.table(gzfile("Homo_sapiens.GRCh38.113.chr.gtf.gz"),
                  header = FALSE, sep = "\t", comment.char = "#",
                  stringsAsFactors = FALSE) %>%
  filter(V3 == "transcript") %>%
  mutate(
    transcript = str_extract(V9, "ENST[0-9]+"),
    gene = str_extract(V9, "ENSG[0-9]+"),
    name = str_extract(V9, "gene_name [^;]+") %>% str_remove("gene_name ")
  ) %>%
  select(transcript, gene, name) %>%
  distinct()

# ------------------------------------
# Load reference (Kallisto) file
# ------------------------------------
if (abundance_type == "tpm") {
  reference <- read.table("gene-tpm_10-CCLE_gene-symbols.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(name = .[[1]]) %>%
    left_join(gtf, by = "name")
  title1 <- "Kallisto TPM (log scale)"
  title2 <- paste0("reindeer_", index_filename, "_vs_kallisto_tpm_ccle.pdf")
  title3 <- paste0("reindeer_", index_filename, "_vs_kallisto_tpm_ccle_2.pdf")
  title4 <- "Estimated abundance in 10 CCLE samples"
} else if (abundance_type == "raw") {
  reference <- read.table("kallisto_raw_counts_10_ccle.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(gene = str_remove(gene, "\\.[0-9]+"))
  title1 <- "Kallisto raw counts (log scale)"
  title2 <- paste0("reindeer_", index_filename, "_vs_kallisto_raw_ccle.pdf")
  title3 <- paste0("reindeer_", index_filename, "_vs_kallisto_raw_ccle_2.pdf")
  title4 <- "Estimated abundance in 10 CCLE samples"
} else if (abundance_type == "seqc") {
  reference <- read.table("gene-tpm_SEQC_16-samples_gene-symbols.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(name = .[[1]]) %>%
    left_join(gtf, by = "name")
  title1 <- "Kallisto TPM (log scale)"
  title2 <- paste0("reindeer_", index_filename, "_vs_kallisto_tpm_seqc.pdf")
  title3 <- paste0("reindeer_", index_filename, "_vs_kallisto_tpm_seqc_2.pdf")
  title4 <- "Estimated abundance in 16 SEQC samples"
}

sample_cols <- grep("^SRR|^seqc", colnames(reference), value = TRUE)

# ------------------------------------
# Load REINDEER2 results
# ------------------------------------
file <- paste0("result_reindeer2_", index_type, "_raw_C50_genes.kallistoformat.tsv")
if (!file.exists(file)) stop("File not found: ", file)

rd2_results <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(NAME = str_remove(NAME, "^>"),
         transcript = str_extract(NAME, "ENST[0-9]+")) %>%
  left_join(gtf, by = "transcript")

if (any(is.na(rd2_results$gene))) warning("Some transcripts could not be mapped to genes")

# Aggregate by gene
rd2_results_agg <- rd2_results %>%
  group_by(gene) %>%
  summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Harmonize gene lists
common_genes <- intersect(rd2_results_agg$gene, if("gene.y" %in% colnames(reference)) reference$gene.y else reference$gene)

if("gene.y" %in% colnames(reference)) {
  results_common <- rd2_results_agg %>% filter(gene %in% common_genes) %>% arrange(gene)
  reference_common <- reference %>% filter(gene.y %in% common_genes) %>% arrange(gene.y) %>% rename(gene = gene.y)
} else {
  results_common <- rd2_results_agg %>% filter(gene %in% common_genes)
  reference_common <- reference %>% filter(gene %in% common_genes)
}

# Reshape to long format
results_long <- results_common %>% pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "result_count")
reference_long <- reference_common %>% pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "ref_count")
df_long <- left_join(results_long, reference_long, by = c("gene", "sample"))

# Global Pearson correlation
agg_corr <- cor(df_long$result_count, df_long$ref_count, use = "complete.obs")

# Plot per-sample faceted plot
p <- ggplot(df_long, aes(x = result_count, y = ref_count)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ sample, scales = "free") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = paste("Reindeer2 abundances (", index_suffix, ", log scale)", sep = ""),
    y = title1,
    title = title4
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 24, face = "bold")
  ) +
  annotate("text",
           x = max(df_long$result_count, na.rm = TRUE),
           y = max(df_long$ref_count, na.rm = TRUE),
           label = paste0("r = ", round(agg_corr, 4)),
           hjust = 0, vjust = 0.5,
           size = 6, color = "red", fontface = "bold")

# Save plots
pdf(title2, width = 12, height = 8)
print(p)
dev.off()

p_single <- ggplot(df_long, aes(x = result_count, y = ref_count)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = paste("Reindeer2 abundances (", index_suffix, ", log scale)", sep = ""),
    y = title1,
    title = title4
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  ) +
  annotate("text",
           x = max(df_long$result_count, na.rm = TRUE),
           y = max(df_long$ref_count, na.rm = TRUE),
           label = paste0("r = ", round(agg_corr, 4)),
           hjust = 5, vjust = 0.5,
           size = 6, color = "red", fontface = "bold")

pdf(title3, width = 12, height = 8)
print(p_single)
dev.off()
