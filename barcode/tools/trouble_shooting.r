#!/usr/bin/env Rscript

library(tidyverse)
library(ggseqlogo)
library(this.path)

# wt1 <- read_tsv("barcode/sxanalysis/D2-wt1-g1n-1.csv.wfilter", col_names = c("sgRNA", "del", "indel", "ins", "1bp", "<2bp", "<3bp", "<4bp", "<5bp", "<6bp", "<7bp", "<8bp", "<9bp", "<10bp")) |>
#   mutate(sample = "wt1")
# wt2 <- read_tsv("barcode/sxanalysis/D2-wt2-g1n-1.csv.wfilter", col_names = c("sgRNA", "del", "indel", "ins", "1bp", "<2bp", "<3bp", "<4bp", "<5bp", "<6bp", "<7bp", "<8bp", "<9bp", "<10bp")) |>
#   mutate(sample = "wt2")
# bind_rows(wt1, wt2) |>
#   mutate(minus4 = factor(substring(sgRNA, 17, 17), levels = c("A", "T", "C", "G"))) |>
#   filter(del >= 0 & del <= 1 & indel >= 0 & indel <= 1 & ins >= 0 & ins <= 1) |>
#   unite(sample4, sample, minus4, remove = FALSE) |>
#   ggplot(aes(minus4, indel + ins)) +
#   stat_summary(fun = "mean") +
#   scale_y_continuous(labels = scales::percent) -> ggfig
# ggsave("percentage.sx.png", path = "barcode/sxanalysis", width = 22, height = 12)

#######################################################
# Usage: trouble_shooting.r fqfile score_thres(default: -Inf)
#######################################################
args <- commandArgs(trailingOnly = TRUE)
fqfile <- args[1]
score_thres <- -Inf
if (length(args) > 1) {
  score_thres <- as.double(args[2])
}
csvfile <- case_match(
  str_to_upper(substring(strsplit(fqfile, "-")[[1]][2], 1, 2)),
  "A1" ~ "final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv",
  "A2" ~ "final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A2.csv",
  "A3" ~ "final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A3.csv",
  "G1" ~ "final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv",
  "G2" ~ "final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G2.csv",
  "G3" ~ "final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G3.csv",
)
csvfile <- file.path(dirname(this.path()), "../csvfiles", csvfile)


# read barcode and sgRNA
barcode_sgRNA <- read_csv(pipe(sprintf("rev %s | cut -c23-40 | tr 'ACGT' 'TGCA'", csvfile)), col_names = "barcode", col_types = "c") |>
  bind_cols(read_csv(pipe(sprintf("sed -r 's/[acgt]+.*$//' %s | rev | cut -c1-20 | rev", csvfile)), col_names = "sgRNA", col_types="c"))

# load table
idtable <- read_tsv(sprintf("%s.table", fqfile), col_types = "ciiiciiiiciiiiciic")

# all insertions
idtable |>
  filter(score >= score_thres) |> # filter score larger than score_thres
  left_join(barcode_sgRNA, by = "barcode") |>
  mutate(m6 = factor(substring(sgRNA, 15, 15), levels = c("A", "T", "C", "G"))) |>
  mutate(m5 = factor(substring(sgRNA, 16, 16), levels = c("A", "T", "C", "G"))) |>
  mutate(m4 = factor(substring(sgRNA, 17, 17), levels = c("A", "T", "C", "G"))) |>
  mutate(m3 = factor(substring(sgRNA, 18, 18), levels = c("A", "T", "C", "G"))) |>
  mutate(m2 = factor(substring(sgRNA, 19, 19), levels = c("A", "T", "C", "G"))) |>
  mutate(m1 = factor(substring(sgRNA, 20, 20), levels = c("A", "T", "C", "G"))) |>
  mutate(insertion = factor(ref_end1 > cut1 | !is.na(random_insertion) | ref_start2 < cut2, levels = c(FALSE, TRUE), labels = c("not_insertion", "insertion"))) |>
  mutate(templated = factor(ref_end1 > cut1 | ref_start2 < cut2, levels = c(FALSE, TRUE), labels = c("not_templated", "templated"))) |>
  pivot_longer(cols = c(m1, m2, m3, m4, m5, m6), names_to = "pos", values_to = "base") ->
  long_idtable

## insertion
# write tsv
long_idtable |>
  summarise(count = sum(count), .by = c(pos, base, insertion)) |>
  arrange(pos, base, insertion) |>
  write_tsv(sprintf("%s.insertion.tsv", fqfile))
# draw figure
long_idtable |>
  ggplot(aes(base, fill = insertion, weight = count)) +
  facet_wrap(~ pos) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) -> ggfig
ggsave(sprintf("%s.insertion.png", basename(fqfile)), path = dirname(fqfile), width = 22, height = 12)

## templated
# write tsv
long_idtable |>
  summarise(count = sum(count), .by = c(pos, base, templated)) |>
  arrange(pos, base, templated) |>
  write_tsv(sprintf("%s.templated.tsv", fqfile))
# draw figure
long_idtable |>
  ggplot(aes(base, fill = templated, weight = count)) +
  facet_wrap(~ pos) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent, , expand = c(0,0)) -> ggfig
ggsave(sprintf("%s.templated.png", basename(fqfile)), path = dirname(fqfile), width = 22, height = 12)

