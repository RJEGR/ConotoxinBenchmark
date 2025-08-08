
#. ESTIMATES Effective sample size!!!
# The definition of Big Data is somewhat nebulous. Typically, this term implies a large number of data points (as opposed to variables) and it is worth noting that the effective sample size might be smaller than the actual data size. In other words, big data does not necessarily mean better data.

# For example (http://www.feat.engineering/important-concepts#model-bias-and-variance):
# --if there is a severe class imbalance or rare event rate, the number of events in the data might be fairly small OR 
# -- particular region of the predictor space is abundantly sampled. Suppose a data set had billions of records but most correspond to redundant within a certain factor range. The number of distinct samples might be low, resulting in a data set that is not diverse.
# -- there can be large amounts of unlabeled data where the outcome is unknown. For example, pharmaceutical companies have large databases of chemical compounds that have been designed but their important characteristics have not been measured (which can be expensive). Other examples include public governmental databases where there is an abundance of data that have not been connected to a specific outcome.


# splintting stratified resampling with randomeness to capture nature variance of data without systemic bias
# Set a random seed for reproducibility, or omit it for fresh random draws.
# Stratified random sampling helps minimize systematic bias by ensuring each group is represented proportionally/randomly, capturing natural variability.

# Resampling methods can generate different versions of our training set that can be used to simulate how well models would perform on new data (http://www.feat.engineering/resampling.html).

# some complex data sets can have multiple levels of data hierarchies which might need to be collapsed so that the modeling data set is consistent with the unit of prediction. There are good and bad ways of summarizing the data, but the tactics are usually subject-specific (http://www.feat.engineering/profile-data#profile-data).

# When the outcome is categorical, stratified splitting techniques can also be applied here to make sure that the analysis and assessment sets produce the same frequency distribution of the outcome. Again, this is a good idea when a continuous outcome is skewed or a categorical outcome is imbalanced, but is unlikely to be problematic otherwise 

## stratified 10-fold cross-validation may be used.

library(rsample)
library(tidyverse)

outdir <- "~/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = outdir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

df <- read_rds(f) %>% dplyr::rename("hits" = "entry_id")

# https://rsample.tidymodels.org



# Stratified bootstrap: resample by 'organismdiet'
boot <- bootstraps(df, strata = organismdiet, times = 100)

# Each bootstrap split captures the natural variance within each group
# Access the resampled data from the first split:
resampled_data <- analysis(boot$splits[[1]])


strata <- df %>%
  mutate(strata = interaction(organismdiet, genesuperfamily, drop = TRUE))

# sampled_df <- df %>%
#   group_by(strata) %>%
#   sample_n(size = 5, replace = TRUE) %>%
#   ungroup()

set.seed(123)

# folds <- rsample::vfold_cv(df, v = 10, strata = organismdiet)

folds <- rsample::vfold_cv(strata, v = 10, strata = strata)

# estimate and diagnose sample bias within and across your folds. 

ess_per_fold <- sapply(folds$splits, function(split) nrow(analysis(split)))
mean_ess <- mean(ess_per_fold)
print(mean_ess)

# continue w 3.4.6 Variance and Bias in Resampling

# f you were to conduct 10-fold cross-validation many times on the same data set, the variation in the resampling scheme could be measured by determining the spread of the resulting averages.This variation could be compared to a different scheme that was repeated the same number of times (and so on) to get relative comparisons of the amount of noise in each scheme.
