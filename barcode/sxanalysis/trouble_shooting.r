library(tidyverse)
library(Biostrings)

wt1 <- read_tsv("barcode/sxanalysis/D2-wt1-g1n-1.csv.wfilter", col_names = c("sgRNA", "del", "indel", "ins", "1bp", "<2bp", "<3bp", "<4bp", "<5bp", "<6bp", "<7bp", "<8bp", "<9bp", "<10bp")) |>
  mutate(sample = "wt1")
wt2 <- read_tsv("barcode/sxanalysis/D2-wt2-g1n-1.csv.wfilter", col_names = c("sgRNA", "del", "indel", "ins", "1bp", "<2bp", "<3bp", "<4bp", "<5bp", "<6bp", "<7bp", "<8bp", "<9bp", "<10bp")) |>
  mutate(sample = "wt2")
bind_rows(wt1, wt2) |>
  mutate(minus4 = factor(substring(sgRNA, 17, 17), levels = c("A", "T", "C", "G"))) |>
  filter(del >= 0 & del <= 1 & indel >= 0 & indel <= 1 & ins >= 0 & ins <= 1) |>
  unite(sample4, sample, minus4, remove = FALSE) |>
  ggplot(aes(minus4, indel + ins)) +
  stat_summary(fun = "mean") +
  scale_y_continuous(labels = scales::percent) -> ggfig
ggsave("percentage.sx.png", path = "barcode/sxanalysis", width = 22, height = 12)

#######################################################

barcode_sgRNA <- read_csv(pipe("rev barcode/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv | cut -c23-40 | tr 'ACGT' 'TGCA'"), col_names = "barcode") |>
  bind_cols(read_csv(pipe("sed -r 's/[acgt]+.*$//' barcode/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv | rev | cut -c1-20 | rev"), col_names = "sgRNA"))

idtable <- read_tsv("barcode/test5/D2-g1n-1.fq.table")

idtable |>
  left_join(barcode_sgRNA, by = "barcode") |>
  mutate(minus4 = factor(substring(sgRNA, 17, 17), levels = c("A", "T", "C", "G"))) |>
  mutate(ins = factor(ref_end1 > cut1 | !is.na(random_insertion) | ref_start2 < cut2, levels = c(FALSE, TRUE), labels = c("not insertion", "insertion"))) |>
  ggplot(aes(minus4, fill = ins, weight = count)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) -> ggfig
ggsave("percentage.insertion.png", path = "barcode/sxanalysis", width = 22, height = 12)

idtable |>
  left_join(barcode_sgRNA, by = "barcode") |>
  mutate(minus4 = factor(substring(sgRNA, 17, 17), levels = c("A", "T", "C", "G"))) |>
  mutate(ins = factor(ref_end1 > cut1 | ref_start2 < cut2, levels = c(FALSE, TRUE), labels = c("not insertion", "insertion"))) |>
  ggplot(aes(minus4, fill = ins, weight = count)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) -> ggfig
ggsave("percentage.templated_insertion.png", path = "barcode/sxanalysis", width = 22, height = 12)
