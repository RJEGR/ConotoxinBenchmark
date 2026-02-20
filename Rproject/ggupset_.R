

# Convert toy set to assembler set

library(tidyverse)
set.seed(32353)
schools_dat <- tibble(
  school_id = paste('School', LETTERS),
  below_in_math = sample(
    c(TRUE, FALSE), 
    length(school_id), 
    replace = TRUE
  ),
  below_in_reading = sample(
    c(TRUE, FALSE), 
    length(school_id), 
    replace = TRUE
  ),
  below_in_writing = sample(
    c(TRUE, FALSE), 
    length(school_id), 
    replace = TRUE
  ),
)

counts_combinations <- schools_dat |>  
  mutate(
    combination = pmap_chr(
      list(below_in_math, below_in_reading, below_in_writing),
      \(lgl1, lgl2, lgl3) {
        c('math', 'reading', 'writing')[c(lgl1, lgl2, lgl3)] |> 
          paste(collapse = ',')
      }
    )
  ) |> 
  count(combination)

counts_combinations <- counts_combinations |> 
  mutate(
    combination = fct_reorder(combination, n, .desc = TRUE)
  )


schools_dat |> 
  pivot_longer(
    cols = -school_id,
    names_to = 'subject',
    values_to = 'below',
    names_prefix = 'below_in_'
  )


counts_of_subjects <- schools_dat |> 
  pivot_longer(
    cols = -school_id,
    names_to = 'subject',
    values_to = 'below',
    names_prefix = 'below_in_'
  ) |> 
  summarize(
    counts = sum(below),
    .by = subject,
  ) 

counts_of_subjects


bar_chart <- counts_combinations  |> 
  ggplot(aes(x = combination, y = n)) +
  geom_col(width = 0.6, fill = 'dodgerblue4') +
  theme_minimal(
    base_size = 16, 
    base_family = 'Source Sans Pro'
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())


counts_combinations |> 
  mutate(
    subjects = map(
      combination, 
      ~str_split_1(as.character(.), ',')
    )
  ) |> 
  unnest(subjects)

points_data <- counts_combinations |> 
  mutate(
    subjects = map(
      combination, 
      ~str_split_1(as.character(.), ',')
    )
  ) |> 
  unnest(subjects) |> 
  mutate(
    subjects = factor(
      subjects, 
      levels = counts_of_subjects |> 
        arrange(counts) |> 
        pull(subject)
    ),
    combination = factor(
      combination,
      levels = levels(counts_combinations$combination)
    )
  ) |> 
  filter(!is.na(subjects)) # Remove the missing values


point_chart <- points_data |> 
  ggplot(aes(x = combination, y = subjects)) +
  geom_line(aes(group = combination), col = 'dodgerblue4') +
  geom_point(size = 10, col = 'dodgerblue4') +
  theme_minimal(base_size = 16, base_family = 'Source Sans Pro') +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(hjust = 0.5)
  ) +
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(drop = FALSE)


subject_bars <- counts_of_subjects |> 
  mutate(subject = fct_reorder(subject, counts)) |>
  ggplot(aes(x = -counts, y = subject)) +
  geom_col(width = 0.6, fill = 'dodgerblue4') +
  theme_minimal(base_size = 16, base_family = 'Source Sans Pro') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())


library(patchwork)

layout <- '
##AAAA
BBCCCC'

bar_chart + 
  subject_bars + 
  point_chart +
  plot_layout(design = layout) +
  plot_annotation(
    title = 'Some schools are below average in maths, reading\nor writing. Here\'s how that looks by subject.',
    caption = 'This uses fake data.',
  )
