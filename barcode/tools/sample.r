#!/usr/bin/env Rscript

library(tidyverse)

fqfile <- "barcode/test1/D2-g1n-1.fq"
idtable <- read_tsv(sprintf("%s.table", fqfile), col_types = "ciiiciiiiciiiicii")

idtable |>
  select(barcode, count, score, ref_end1, query_end1, ref_start2, query_start2, cut1, cut2) |>
  mutate(uptemp = ref_end1 > cut1, random = query_end1 < query_start2, downtemp = ref_start2 < cut2) |>
  arrange(desc(score), ref_end1) |> print(width = Inf) |>
  summarise(count = sum(count), .by = c(uptemp, downtemp)) |>
  ggplot(aes(uptemp, weight = count)) +
  geom_bar()

idtable |>
  select(barcode, count, score, ref_end1, query_end1, ref_start2, query_start2, cut1, cut2) |>
  arrange(desc(count))
