library(tidyverse)
library(patchwork)
library(extrafont)

loadfonts(device = "all")
mytheme <- theme_bw() + theme(text = element_text(family = "Arial", colour = "black", size = 15))
softwares <- c("rearr", "RESSO", "AMP", "CRVS", "ADIV", "CRGR", "ZhangFeng", "SelfTarget")
brdf <-
  read_delim("bench/benchresult", col_types = "cidicdddiciiii", na = "*")  

#############################################
# memory and speed
#############################################

brdf |>
  distinct(mode, reflen, probability, readnum, software, usertime, systime, realtime, memory) |>
  filter(software != "simulate") |>
  pivot_longer(cols = c(usertime, systime, realtime, memory), names_to = "performance", values_to = "cost") |>
  mutate(reflen = factor(reflen)) |>
  mutate(software = factor(software, levels = softwares)) |>
  mutate(mode = factor(mode, levels = c("single", "double"))) |>
  mutate(probability = factor(probability)) |>
  mutate(performance = factor(performance, levels = c("usertime", "systime", "realtime", "memory"))) |>
  ggplot(aes(reflen, cost, color = software, shape = software)) +
  geom_point(size = 3) +
  facet_grid(performance ~ mode + probability, scales = "free_y") +
  scale_color_manual(values = rep(c("black", "red", "blue", "purple"), times = 2)) +
  scale_shape_manual(values = rep(c(1, 2, 3, 4), each = 2)) +
  scale_y_log10() +
  mytheme -> ggfig
ggsave("performance.png", path = "bench", width = 22, height = 12)

############################################
# report rate
############################################

brdf |>
  filter(!is.na(refup) & !is.na(refdown) & !is.na(queryup) & !is.na(querydown) & software != "simulate") |>
  summarise(count = n(), .by = c(mode, reflen, probability, software)) |>
  mutate(reflen = factor(reflen)) |>
  mutate(software = factor(software, levels = softwares)) |>
  mutate(mode = factor(mode, levels = c("single", "double"))) |>
  mutate(probability = factor(probability)) |>
  ggplot(aes(reflen, software, size = count, color = software)) +
  geom_point() +
  facet_grid(mode ~ probability) +
  scale_size_area(limits = c(0, NA)) +
  mytheme -> ggfig
ggsave("report_rate.png", path = "bench", width = 22, height = 12)

###########################################
# accuracy
###########################################

brdf |>
  filter(software != "simulate") |>
  left_join(brdf |> filter(software == "simulate") |> select(!(readnum:memory)), by = c("mode", "reflen", "probability", "query"), suffix = c("", ".sim")) ->
  brdf_long

brdf_long |>
  pivot_longer(cols = c(refup.sim, refdown.sim, queryup.sim, querydown.sim), names_to = "end.sim", values_to = "pos.sim") |>
  pull(pos.sim) ->
  pos.sims

brdf_long |>
  select(!(refup.sim:querydown.sim)) |>
  pivot_longer(cols = c(refup, refdown, queryup, querydown), names_to = "end", values_to = "pos") |>
  mutate(pos.sim = pos.sims) |>
  mutate(diff = abs(pos - pos.sim)) ->
  brdf_tidy

brdf_tidy |>
  pivot_wider(id_cols = c(mode, reflen, probability, query, end), names_from = software, values_from = diff) |>
  mutate(diffmax = pmax(rearr, RESSO, AMP, CRVS, ADIV, CRGR, ZhangFeng, SelfTarget, na.rm = TRUE)) |>
  select(mode, reflen, probability, query, end, diffmax) |>
  left_join(brdf_tidy, y = _, by = c("mode", "reflen", "probability", "query", "end")) |>
  mutate(diff = replace(diff, which(is.na(diff)), diffmax[is.na(diff)])) ->
  brdf_tidy_diffnona

varcombs = list()
for (m in 0:4) {
    varcombs <- c(varcombs, combn(c("mode", "reflen", "probability", "end"), m, simplify = FALSE))
}
brdf_tidy_summary <- tibble(mode = character(), reflen = character(), probability = character(), software = character(), end = character(), diffmean = double())
for (varcomb in varcombs) {
    convarcomb <- c(setdiff(c("mode", "reflen", "probability", "end"), varcomb), "software")
    brdf_tidy_diffnona |> summarise(diffmean = mean(diff), .by = {convarcomb}) -> brdf_tidy_summary_
    for (var in convarcomb) {
        brdf_tidy_summary_[[var]] = as.character(brdf_tidy_summary_[[var]])
    }
    for (var in varcomb) {
        brdf_tidy_summary_[var] = "all"
    }
    print(brdf_tidy_summary_)
    brdf_tidy_summary <- bind_rows(brdf_tidy_summary, brdf_tidy_summary_)
}

brdf_tidy_summary |>
  mutate(reflen = factor(reflen, levels = c("100", "200", "300", "400", "500", "all"))) |>
  mutate(probability = factor(probability, levels = c("0.01", "0.02", "0.03", "0.04", "0.05", "all"))) |>
  mutate(software = factor(software, levels = softwares)) |>
  mutate(mode = factor(mode, levels = c("single", "double", "all"))) |>
  mutate(end = factor(end, levels = c("refup", "refdown", "queryup", "querydown", "all"))) |>
  ggplot(aes(reflen, diffmean, color = software, shape = software)) +
  geom_point(size = 3) +
  facet_grid(end ~ mode + probability) +
  scale_color_manual(values = rep(c("black", "red", "blue", "purple"), times = 2)) +
  scale_shape_manual(values = rep(c(1, 2, 3, 4), each = 2)) +
  scale_y_log10() +
  mytheme -> ggfig
ggsave("accuracy.png", path = "bench", width = 22, height = 12)
