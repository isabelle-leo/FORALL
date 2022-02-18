library(tidyverse)
library(ggrepel)
library(ggthemes)
library(svglite)
library(gtools)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_05")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# aesthetic constants
lineColor <- "#999999"
labelColor <- "#555555"
segmentColor <- "#AAAAAA"
xLabel <- expression(paste("delta mean sDSS of B-lineage vs. T-lineage"))
yLabel <- expression(paste(-log[10]~italic(adjp), "-value"))
densityColor <- "Red"

# load drug sensitivity for drugs
sDSS_data <- readRDS(file = file.path("data", "dsrt", "drugsens_2021-03-22.RDS"))

# labeling manual
#drugsToLabel <- c("FIMM115484") "FIMM115484...Venetoclax"

#if set to TRUE, it will overwrite drugsToLabel
cutoff_based_label  = TRUE # set to FALSE to disable 
abs_delta_sDSS_cutoff = 5
neglog10_pvalue_cutoff = -log10(0.02)

# calculate B versus T result ===============

#Set your comparison groups variables!
# based on types
groupA <- sDSS_data@phenoData@data$Cell.Line.R[grep("preB|B-ALL", sDSS_data@phenoData@data$Type)]
groupB <- sDSS_data@phenoData@data$Cell.Line.R[grep("T-ALL", sDSS_data@phenoData@data$Type)]
 
# make empty data frame for filling t tests, will put FIMM name
ttests <- data.frame("drug" = rownames(sDSS_data@assayData$exprs),
                     "t" = 1:528, 
                     "p.value" = 1:528,
                     "FDR" = 1:528,
                     "delta_mean_sDSS" = 1:528)
#will put drug name
ttests <- data.frame("drug" = sDSS_data@featureData@data$DRUG.NAME,
                     "t" = 1:528,
                     "p.value" = 1:528,
                     "FDR" = 1:528,
                     "delta_mean_sDSS" = 1:528)

#calculate statistics for each drug
for(i in 1:nrow(sDSS_data@assayData$exprs)){
  data_A = c(sDSS_data@assayData$exprs[i, groupA])
  data_B = c(sDSS_data@assayData$exprs[i, groupB])
  
  if(all(sd(data_A, na.rm = T) > 0, sd(data_B, na.rm = T) > 0)){
    tmp_ttest_output <- t.test(data_A, data_B)
    ttests[i, c("t", "p.value", "delta_mean_sDSS")] <- c(round(tmp_ttest_output$statistic, 3), signif(tmp_ttest_output$p.value, 3), round(mean(data_A) - mean(data_B), 3))
    }else{
      ttests[i, c("t", "p.value", "delta_mean_sDSS")] <- c(NA, NA,  round(mean(data_A) - mean(data_B), 3))
    }
}

#calculate FDR
ttests$FDR <- p.adjust(ttests$p.value,method = "fdr", n = nrow(ttests))

#filter results that can't be transformed or plotted
ttests.filt <- ttests[!is.na(ttests$p.value) , ]
            

#log transform the statistics for plotting
ttests.filt  <- mutate(ttests.filt, negativeLog10.p.value = -log10(as.numeric(p.value)))
ttests.filt  <- mutate(ttests.filt, negativeLog10.adjpvalue = -log10(as.numeric(FDR)))

if(cutoff_based_label){
  drugsToLabel <- ttests.filt$drug[(abs(ttests.filt$delta_mean_sDSS) >= abs_delta_sDSS_cutoff) & 
                                     (abs(ttests.filt$negativeLog10.adjpvalue) >= neglog10_pvalue_cutoff)]
}

# calculate density
# adapted from https://slowkow.com/notes/ggplot2-color-by-density/
getDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

ttests.filt$density <- getDensity(x = ttests.filt$delta_mean_sDSS,
                           y = ttests.filt$negativeLog10.p.value,
                           h = c(2, 2),
                           n = 200)

# plot
p <- ggplot(mapping = aes(x = delta_mean_sDSS, y = negativeLog10.adjpvalue, color = density)) +
  geom_hline(yintercept = -log10(0.05), color = lineColor, linetype = "dotted") +
  geom_vline(xintercept = 5, color = lineColor, linetype = "dotted") +
  geom_vline(xintercept = -5, color = lineColor, linetype = "dotted") +
  geom_point(data = ttests.filt %>% filter(!drug %in% drugsToLabel),
             size = 1.5,
             stroke = 0,
             alpha = 0.9) +
  scale_color_gradient_tableau(densityColor) +
  geom_label_repel(data = ttests.filt %>% filter(drug %in% drugsToLabel),
                   aes(label = drug),
                   color = labelColor,
                   size = 2,
                   direction = "both",
                   box.padding = 0.05,
                   point.padding = 0.2,
                   label.size = 0.15,
                   segment.size = 0.6,
                   segment.color = segmentColor,
                   max.overlaps = 30,
                   min.segment.length = 0.2) +
  geom_point(data = ttests.filt %>% filter(drug %in% drugsToLabel),
             size = 3,
             stroke = 0,
             alpha = 0.7,
             color = colorScheme$blue) +
  xlab(xLabel) +
  ylab(yLabel) +
  theme(legend.position = "none")


#view plot
p

outputFolder <- file.path("figures", "output")

# export to svg
ggsave(
  "figure_05_B_versus_T_sDSS_volcano_rev2.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 6,
  height = 6,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
