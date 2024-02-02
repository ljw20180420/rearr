library(ggplot2)
library(patchwork)
library(ggtext)
library(Cairo)
library(reticulate)
library(tidyverse)
library(lintr)

##############
# theme
##############
theme(axis.title.x = ggtext::element_markdown())

##############
# read data
##############
# args <- commandArgs(trailingOnly = TRUE)
# filename <- args[1]
# filename <- 'barcode/test1/A2_TEST.fq.final_hgsgrna_libb_all_0811-NGG.csv.barcode.table'
filename <- 'barcode/test4/A7-G1-n-1.R2.fq.final_hgsgrna_libb_all_0811-NGG.csv.barcode.table'
csvfile <- 'barcode/final_hgsgrna_libb_all_0811-NGG.csv'
filename %>% sprintf('rearr_barcode_post_process.sh <%s', .) %>% pipe %>% read.table(header = TRUE, sep = "\t") -> idtable
csvfile %>% sprintf("tail -n+2 %s | cut -d',' -f12 | tr ACGT TGCA | rev | sort", .) %>% pipe %>% read.table(col.names="barcode") %>% .$barcode -> barcodes
idtable$barcode <- factor(idtable$barcode, levels = barcodes)
idtable$indel_type <- factor(idtable$indel_type, levels = c("WT", "del", "ins", "indel"))


###################
# stat_sum geom_raster version
###################

ggmain <- ggplot(idtable, aes(barcode, indel_type)) +
stat_sum(aes(fill = after_stat(n)), geom = "raster") +
scale_x_discrete(limits = barcodes, breaks = NULL) +
scale_y_discrete(limits  = c("WT", "del", "ins", "indel")) +
scale_fill_viridis_c(guide = guide_colorbar(barheight = unit(0.5, "inches")))

ggmargin_x <- ggplot(idtable, aes(barcode)) +
geom_bar() +
scale_x_discrete(name = NULL, limits = barcodes, breaks = NULL)

ggmargin_y <- ggplot(idtable, aes(y = indel_type)) +
geom_bar(orientation = "y") +
scale_y_discrete(name = NULL, limits  = c("WT", "del", "ins", "indel"), breaks = NULL)

ggtotal <- ggmain + ggmargin_x + ggmargin_y + guide_area() + plot_layout(guides = "collect", heights = c(1, 3), widths = c(10, 1), design = "BD\nAC")
ggtotal[[3]] <- ggtotal[[3]] + scale_x_continuous(guide = guide_axis(angle = 90))
ggsave(sprintf("%s.png", filename), width = 12, height = 4, ggtotal)

###################
# stat_count version
###################

ggmain <- ggplot(idtable %>% count(barcode, indel_type, .drop = FALSE), aes(barcode, n)) +
geom_point(shape = 1) +
facet_grid(rows = vars(indel_type)) +
# stat_smooth(method = "lm") +
scale_x_discrete(limits = barcodes, breaks = NULL) +
scale_y_continuous(limits = c(0, NA))

ggmargin_x <- ggplot(idtable, aes(barcode)) +
geom_bar(width = 1, color = "black") +
scale_x_discrete(name = NULL, limits = barcodes, breaks = NULL)

ggmargin_y <- ggplot(idtable, aes(y = indel_type)) +
geom_bar(width = 1, color = "black", orientation = "y") +
scale_y_discrete(name = NULL, limits  = c("WT", "del", "ins", "indel"), breaks = NULL) +
scale_x_continuous(guide = guide_axis(angle = -90))

ggtotal <- ggmain + ggmargin_x + ggmargin_y + plot_spacer() + plot_layout(guides = "collect", heights = c(1, 3), widths = c(10, 1), design = "BD\nAC")
ggsave(sprintf("%s.png", filename), width = 12, height = 4, ggtotal)

#######################
# test
#######################
testdf <- data.frame(x = rnorm(1000), y = rnorm(1000), z = c(rnorm(10000, 2), rnorm(10000, -2)))
testfig <- ggplot(testdf, aes(y = y)) +
# geom_point() +
# geom_smooth(method = "lm", formula = "y~1")
print(testfig)

library(EBImage)
otsu(array(testdf$z, dim = c(length(testdf$z), 1)), range = range(testdf$z))
