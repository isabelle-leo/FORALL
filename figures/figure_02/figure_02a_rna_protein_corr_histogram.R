library(tidyverse)
library(ggrepel)
library(ggthemes)
library(svglite)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_02")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
data <- read_tsv(file = file.path("data", "mrna_protein", "ALL_corresponding_RNA_protein_correlation.txt")) %>%
  filter(!is.na(spearman_corr)) %>%
  mutate(positive = factor(as.numeric(spearman_corr >= 0), levels = c(1, 0)))

meta <- data.frame(mean = data$spearman_corr %>% mean(),
                   median = data$spearman_corr %>% median())

# plot
binWidth <- 0.05
p <- ggplot(data = data, aes(x = spearman_corr)) +
  geom_vline(xintercept = 0, size = 0.3) +
  geom_histogram(aes(y = stat(count / sum(count) / binWidth), fill = positive), binwidth = binWidth, alpha = 0.6, color = "white") +
  scale_fill_manual(values = c("0" = "tomato", "1" = "deepskyblue2")) +
  geom_density(color = NA, alpha = 0.3, fill = "steelblue") +
  geom_line(stat = "density", alpha = 0.9, color = "steelblue", size = 1) +
  geom_vline(xintercept = meta$mean, linetype = "dotted", color = "gray40", size = 0.5) +
  geom_text(data = meta, aes(x = mean * 0.65 , y = 1.4, label = paste0("Mean: ", round(mean, 2))), color = "#444444") +
  geom_vline(xintercept = meta$median, linetype = "dotted", color = "gray40", size = 0.5) +
  geom_text(data = meta, aes(x = median * 1.38, y = 1.4, label = paste0("Median: ", round(median, 2))), color = "#444444") +
  xlab("Spearman mRNA-protein correlation") +
  ylab("Density") +
  theme(legend.position = "none")
p

# export to svg
ggsave(filename = "protein_corr_histogram.svg",
       plot = p,
       device = NULL,
       path = outputFolder,
       scale = 1,
       width = 6,
       height = 4,
       units = c("cm"),
       dpi = 600,
       limitsize = FALSE)
  