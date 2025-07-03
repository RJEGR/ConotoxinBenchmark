

# Read transrate scores (to evals accuracy of assembly methods)
# Read detonate scores (optional)
#


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


dir <- "~/Documents/GitHub/ConotoxinBenchmark/transrate/"

subdirs <- list.files(dir, pattern = "_transrate")

read_transrate_scores <- function(path) {
  # Path where 
  
  subdir <- path[1]
  
  # basedir <- sapply(strsplit(subdir, "/"), `[`, 3)
  
  basedir <- list.files(file.path(dir, subdir), pattern = "trinity_ensamble.Trinity")
  
  f <- file.path(dir, subdir, basedir, "contigs.csv")
  
  read_csv(f) %>% mutate(Method = basedir)
  
}

transratedf <- lapply(subdirs, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>% count(Method)

transratedf <- transratedf %>% 
  separate(Method, into = c("xcov", "sf"), sep = "_", remove = F) 

data_text <- transratedf %>% 
  drop_na(hits) %>%
  count(sf, xcov) %>%
  mutate(xcov = gsub("x$","", xcov), xcov = as.numeric(xcov)) %>%
  filter(xcov == 500) %>% distinct(sf, xcov, n)

transratedf %>%
  drop_na(hits) %>%
  count(sf, xcov) %>%
  mutate(xcov = gsub("x$","", xcov), xcov = as.numeric(xcov)) %>%
  ggplot(aes(x = xcov, y = n, group = sf)) +
  geom_point() +
  geom_line() +
  geom_text(data = data_text, aes(label = sf), hjust = 2, vjust = -2)
  
