# browser()

library(tidyverse)

# args <- commandArgs(trailingOnly = TRUE)
filename <- "bench/runs/double.200.0.03/rearr/random.fq.alg.100.100.200"
fields <- strsplit(filename, ".", fixed = TRUE)[[1]]
ref1len <- strtoi(last(fields))
cut1 <- strtoi(fields[length(fields) - 2])
cut2 <- strtoi(fields[length(fields) - 1])
lineNum <- strtoi(system(sprintf("wc -l <%s", filename), intern = TRUE))
fd <- file(description = filename, open = "r")
counts <- reflines <- querylines <- rep(NA, lineNum / 3)
lposes <- rposes <- vector(mode = "list", length = lineNum / 3)
for (i in seq(lineNum / 3)){
    lines <- readLines(con = fd, n = 3)
    counts[i] <- strtoi(strsplit(lines[1], "\t")[[1]][2])
    reflines[i] <- lines[2] |> toupper()
    querylines[i] <- lines[3]
    rposes[[i]] <- cumsum(!(strsplit(reflines[i], "")[[1]] %in% c(" ", "-")))
    lposes[[i]] <- c(0, rposes[[i]][1:(length(rposes[[i]]) - 1)])
}
close(fd)
base_pos_df <- tibble(
    refbase = unlist(strsplit(reflines, NULL)),
    querybase = unlist(strsplit(querylines, NULL)),
    lpos = unlist(lposes),
    rpos = unlist(rposes),
    count = rep(counts, times = sapply(querylines, nchar, USE.NAMES = FALSE))
    )
ref1 <- gsub("[- ]", "", reflines[1])
ref2 <- ref1 |> substr(ref1len + 1, nchar(ref1))
ref1 <- ref1 |> substr(1, ref1len)

base_pos_df |>
  mutate(part = factor(ifelse(rpos == 0, "updangle", ifelse(lpos == nchar(ref1) + nchar(ref2), "downdangle", ifelse(lpos == nchar(ref1) & rpos == nchar(ref1), "randomins", ifelse(lpos < nchar(ref1), "ref1", "ref2")))), levels = c("updangle", "ref1", "randomins", "ref2", "downdangle"))) |>
  mutate(posstatus = factor(ifelse(refbase == "-", "ins", ifelse(querybase == "-", "del", ifelse(querybase == refbase, "match", "SNP"))), levels = c("ins", "del", "SNP", "match"))) ->
  base_pos_df

base_pos_df |>
  mutate(rel1pos = (lpos + rpos) / 2 - cut1) |>
  filter(part == "ref1") |>
  ggplot(aes(rel1pos, fill = posstatus, weight = count)) +
  stat_bin(data = ~ filter(.x, posstatus != "ins"), breaks = seq(-cut1, nchar(ref1) - cut1)) +
  stat_bin(mapping = aes(color = "black"), data = ~ filter(.x, posstatus == "ins"), geom = "step", direction = "mid", breaks = seq(-cut1 + 0.5, nchar(ref1) - cut1 - 0.5)) +
  scale_x_continuous(name = "position relative to cut1") +
  scale_fill_discrete(name = NULL, limits = c("del", "SNP", "match"), labels = c("deletion", "SNP", "match")) +
  scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion") ->
  ggfig
ggsave("indel.pos.ref1.pdf", plot = ggfig, path = "figures")

base_pos_df |>
  mutate(rel2pos = (lpos + rpos) / 2 - nchar(ref1) - cut2) |>
  filter(part == "ref2") |>
  ggplot(aes(rel2pos, fill = posstatus, weight = count)) +
  stat_bin(data = ~ filter(.x, posstatus != "ins"), breaks = seq(-cut1, nchar(ref1) - cut1)) +
  stat_bin(mapping = aes(color = "black"), data = ~ filter(.x, posstatus == "ins"), geom = "step", direction = "mid", breaks = seq(-cut1 + 0.5, nchar(ref1) - cut1 - 0.5)) +
  scale_x_continuous(name = "position relative to cut2") +
  scale_fill_discrete(name = NULL, limits = c("del", "SNP", "match"), labels = c("deletion", "SNP", "match")) +
  scale_color_identity(name = NULL, guide = guide_legend(), labels = "insertion") ->
  ggfig
ggsave("indel.pos.ref2.pdf", plot = ggfig, path = "figures")

base_pos_df |>
  ggplot(aes(part, weight = count)) +
  geom_bar() +
  scale_x_discrete(limits = c("updangle", "randomins", "downdangle")) +
  scale_y_log10() ->
  ggfig
ggsave("indel.pos.dangle.pdf", plot = ggfig, path = "figures")

base_pos_df |>
  mutate(trans = factor(sprintf("%s>%s", refbase, querybase),
    levels = c("A>A", "C>C", "G>G", "T>T",
    "A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G",
    "A>-", "C>-", "G>-", "T>-",
    "->A", "->C", "->G", "->T"
    ))) |>
  ggplot(aes(trans, weight = count)) +
  geom_bar() +
  scale_y_log10() ->
  ggfig
ggsave("indel.pos.trans.pdf", plot = ggfig, path = "figures")

ref1vec <- strsplit({{ref1}}, "")[[1]]
ref2vec <- strsplit({{ref2}}, "")[[1]]
microhomo <- matrix(as.integer(rep(ref1vec, time = length(ref2vec)) == rep(ref2vec, each = length(ref1vec))), nrow = length(ref1vec))
for (i in seq(nrow(microhomo))) {
  for (j in seq(ncol(microhomo))) {
    if (i > 1 && j > 1)
    {
      if (microhomo[i, j] > 0) {
        microhomo[i, j] <- microhomo[i - 1, j - 1] + microhomo[i, j]
      }
    }
  }
}
rc <- which(microhomo > 3, arr.ind = T)
pos12 <- rc[rep(seq(nrow(rc)), time = microhomo[rc]),] -
  rep(unlist(sapply(microhomo[rc], function(i) seq(i - 1, 0, by = -1))), time = 2)
tibble(pos1 = pos12[, 1], pos2 = pos12[, 2]) |>
  ggplot(aes(pos1, pos2)) +
  geom_point() ->
  ggfig
ggsave("ref1.ref2.micro.homo.pdf", plot = ggfig, path = "figures")


tablefile <- sub(".alg.", ".table.", filename)
indel_tsv <- read_tsv(tablefile, col_types = "iiiciiiiciiiiciid", na = "NA")
indel_tsv |>
  mutate(insertion = (ref_end1 > cut1 | ref_start2 < cut2 | nchar(random_insertion) > 0)) |>
  mutate(deletion = (ref_end1 < cut1 | ref_start2 > cut2)) |>
  mutate(indel_type = factor(ifelse(insertion & !deletion, "insertion",
                             ifelse(deletion & !insertion, "deletion",
                             ifelse(insertion & deletion, "indel", "WT"
                      ))), levels = c("WT", "deletion", "insertion", "indel"))) ->
  indel_tsv

indel_tsv |>
  summarise(count = n(), .by = indel_type) |>
  mutate(percent = count / sum(count), perlabel = scales::percent(percent, accuracy = 0.01)) |>
  mutate(type_count = sprintf("%s: %d", indel_type, count)) |>
  mutate(type_count = factor(type_count, levels = c(type_count[startsWith(type_count, "WT")], type_count[startsWith(type_count, "deletion")], type_count[startsWith(type_count, "insertion")], type_count[startsWith(type_count, "indel")]))) |>
  ggplot(aes(1, percent, fill = type_count, weight = count)) +
  geom_col() +
  geom_text(aes(label = perlabel), position = position_stack(vjust = 0.5)) +
  scale_x_discrete(name = NULL, breaks = NULL) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = scales::percent) +
  coord_polar(theta = "y") ->
  ggfig
ggsave("indel.type.pie.pdf", plot = ggfig, path = "figures")

indel_tsv |>
  mutate(uplen = nchar(updangle)) |>
  mutate(randlen = nchar(random_insertion)) |>
  mutate(downlen = nchar(downdangle)) |>
  pivot_longer(cols = c(uplen, randlen, downlen), names_to = "unmapped", values_to = "unmapped_length") |>
  ggplot(aes(unmapped_length, fill = unmapped)) +
  geom_bar() +
  scale_y_log10() +
  scale_fill_discrete(limits = c("uplen", "randlen", "downlen"), labels = c("updangle length", "random insertion length", "downdangle length")) ->
  ggfig
ggsave("unmapped.length.pdf", plot = ggfig, path = "figures")

