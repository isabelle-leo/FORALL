# dependencies
library(Biobase)
library(tidyverse)
library(ggpubr)
library(ggthemes)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_SI05")
# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# load data
drugsens <- readRDS(file = file.path("data", "dsrt", "drugsens_2021-03-22.RDS")) 

# filter example to remove B-ALL
#drugsens <- drugsens[, !pData(drugsens)$Cell.Line.R %in% c("TANOUE", "MN.60")] #If you want to remove certain cell lines

# parameters
drugId <- c("FIMM136453") # FIMM136453	Sabutoclax

# FIMM003707	Navitoclax
# FIMM136514	A-1155463
# FIMM115484	Venetoclax
# FIMM136513	A-1331852
# FIMM001823	AT-101
# FIMM136391	WEHI-539

# highlight a cell line
#cellLineHighlight <- c("X.697", "RCH.ACV", "KASUMI.2", "COG.319", "MHH.CALL.3", "COG.402")
cellLineHighlight <- c(NULL)

# tidy
df <- drugsens %>%
  exprs() %>%
  as.data.frame() %>%
  cbind(drugsens %>% fData()) %>%
  filter(FIMM.ID == drugId) %>%
  gather(key = "Cell.Line.R", value = "sdss", 1:ncol(drugsens)) %>%
  left_join(drugsens %>% pData(), by = "Cell.Line.R") %>%
  arrange(desc(sdss)) %>%
  mutate(Cell.Line.R = factor(Cell.Line.R, levels = .$Cell.Line.R))

# plot
# with subtype coloring
p <- df %>%
  ggplot(aes(x = Cell.Line.R, y = sdss, color = Subtype_Paper)) + #Change for coloring by #Type_Paper #Subtype_Paper need to change the scalecolor below 
  geom_hline(yintercept = 0,
             color = "gray50") +
  geom_segment(aes(xend = Cell.Line.R,
                   y = 0,
                   yend = sdss),
               size = 1,
               alpha = 0.6) + # lollipop plot
  geom_point(data = df %>% filter(Cell.Line.R %in% cellLineHighlight),
             size = 5,
             color = colorScheme$red,
             alpha = 0.6) + # scatter and lollopop plot
  geom_point(size = 3) + # scatter and lollopop plot #change lolli-size
  # geom_bar(stat = "identity") + # bar plot
  scale_color_manual(values = colorScheme$subtype.alt) +      #Change this if you change coloring legend "colorScheme$type" "colorScheme$subtype.alt"
    xlab("Cell line") +
  ylab("sDSS") +
  ggtitle(df$DRUG.NAME %>% unique() %>% paste(collapse = ", ")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  #geom_hline(data = summary_statistics,aes(yintercept=value,linetype=Statistics,colour=NULL),size=0.5)
  geom_hline(yintercept=8, linetype="dashed", 
            color = "firebrick", size=0.3)

# View the plot
p

# export to svg
library(svglite)
ggsave(
  "figure_SI04_drugSensPerSample_Sabutoclax_rev.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 12,
  height = 12,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
