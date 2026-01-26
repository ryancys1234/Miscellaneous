.libPaths("~/miniconda3/envs/epi293/lib/R/library") # Note: This line prevents errors with loading dplyr.
library(data.table); library(dplyr); library(fastman)
options(datatable.fread.datatable = F)

CWD <- "~/EPI293/"
phenoCol <- "predictScore" # Replace this with the desired variable
platforms <- c("AffymetrixData", "GlobalScreeningArrayData", "HumanCoreExData2", "IlluminaHumanHapData", "OmniExpressData", "OncoArrayData")

# Manhattan and QQ plots
all_df <- fread(paste0(CWD, "meta_analysis_chrall_", phenoCol, "_results1.txt"))
all_df <- all_df |>
  mutate(CHR = sub(":.*", "", MarkerName)) |>
  mutate(CHR = as.numeric(substr(CHR, 4, 5))) |>
  mutate(POS = sub("^[^:]*:([^:]*):.*", "\\1", MarkerName))
all_df$CHR <- factor(all_df$CHR, levels = 1:22)
all_df$POS <- as.numeric(all_df$POS)

png(paste0(CWD, "plots/METAL_manhattan_plot.png"), width = 14, height = 5, units = "in", res = 300)
fastman(all_df, chr = "CHR", bp = "POS", p = "P-value", col = "greys")
dev.off()

png(paste0(CWD, "plots/METAL_qq_plot.png"), width = 5, height = 5, units = "in", res = 300)
fastqq(all_df$`P-value`)
dev.off()

# MAF and imputation quality plots
for (p in platforms) {
  my_data <- fread(paste0(CWD, "platform_", p, "/regenie/step2/hpfs_step2_chrall_", phenoCol, "_withP.regenie"))

  png(file = paste0(CWD, "plots/", p, "_MAF.png"))
  hist(as.numeric(my_data$A1FREQ), main = paste("MAF distribution -", p), xlab = "MAF")
  dev.off()

  png(file = paste0(CWD, "plots/", p, "_quality.png"))
  hist(as.numeric(my_data$INFO), main = paste("Imputation quality distribution -", p), xlab = "Imputation quality")
  dev.off()
}

# Top 10 SNPs
df <- fread(paste0(CWD, "/meta_analysis_chrall_", phenoCol, "_results1.txt"))
sink(paste0(CWD, "top_hits_results.txt"))
head(df[order(df$`P-value`),], 10)
sink()
