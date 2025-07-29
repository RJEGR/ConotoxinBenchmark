# From complexity.md, bind data with ConoServerDB

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)


# Complexity data here

# Correlate 
# data %>%
#   select_if(is.double) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   cor(method = "spearman") -> M
# 
# 
# testRes <- corrplot::cor.mtest(M, conf.level = 0.95)
# 
# corrplot::corrplot(M, p.mat = testRes$p ,method = "color", type="upper", order = "hclust", insig = "label_sig")
