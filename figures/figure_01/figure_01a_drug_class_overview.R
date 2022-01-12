library(tidyverse)
library(ggthemes)
library(svglite)

source(file.path("functions", "getColorScheme.R"))
source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

# constants
useVersion   <- "publication"
colorScheme  <- getColorScheme()

# Set outputfolder
outputFolder <- file.path("figures", "output", "figure_01")

# create maijor output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)



df <- as.data.frame(read.table(file = file.path("meta", "drug_classes.txt"), sep = "\t", header = T))

drugColors <- colorScheme$Class.explained %>%
  set_names(names(.) %>%
              gsub("^[A-Z]\\.", "", .) %>%
              gsub("\\.", " ", .) %>%
              gsub("Conv Chemo", "Conv. Chemotherapy", .) %>%
              gsub("Combo", "combo", .) %>%
              gsub("Protease proteasome", "Protease/proteasome", .))


#Version 1: Y-axis = Drug Type
df$Drug.Type <- factor(df$Drug.Type, levels = df$Drug.Type[order(df$Number.of.Drugs)])
# pdf(file="Figure_1A_drug_class_v1.pdf", width = 8, height = 4)
p <- ggplot(df, aes(y = Drug.Type, x = Number.of.Drugs, fill = Drug.Type)) +
  geom_bar(stat = "identity") +
  geom_text(label= df$Number.of.Drugs, color = "black", hjust = -0.2, size = 2) +
  scale_fill_manual(guide = FALSE, values = drugColors) +
  xlab("Number of drugs") +
  ylab("Drug class") +
  theme(axis.text = element_text(size = 6, family = "Arial"))
# dev.off()
p

#Version 2: Drug Type as legend
# df$Drug.Type <- factor(df$Drug.Type, levels = c("Kinase inhibitor", "Differentiating epigenetic modifier", "Conv. Chemotherapy", "Other", "Apoptotic modulator", "Hormone therapy", "Immunomodulatory", "Metabolic modifier", "Protease/proteasome inhibitor", "HSP inhibitor", 	"Rapalog","Hedgehog inhibitor","Kinesin inhibitor", "NSAID", "Conv. Chemotherapy combo"))
# pdf(file="Figure_1A_drug_class_v2.pdf", width = 8, height = 4)
# ggplot(df, aes(x = Drug.Type, y = Number.of.Drugs, fill = Drug.Type))+
#   geom_bar(stat = "identity")+
#   geom_text(label= df$Number.of.Drugs, color = "Black", vjust= -0.3)+
#   scale_fill_manual(name = "Drug Class",values=c("#ffafbb","#fe968b","#8dc993","#b5bdda","#4abdde","#4790db","#fe7e75","#49dcb5","#fefd90","#fded9c","#969b96","#c9c9fb","#49bd93","#fde190","#d7b6e5"))+
#   ylab("Number of Drugs")+
#   xlab("Drug Class")+
#   theme_classic()+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# dev.off()

#c("Conv. Chemotherapy combo" = "#d7b6e5", "NSAID" = "#fde190", "Immunomodulatory" = "#fe7e75","Metabolic modifier"  = "#49dcb5","Differentiating.epigenetic modifier"  = 	"#fe968b",	"Kinase inhibitor" = "#ffafbb",	"Other" = "#b5bdda",	"Apoptotic modulator" = "#4abdde",	"Conv.Chemo" = 	"#8dc993",	"Protease.proteasome inhibitor" = "#fefd90",	"Hormone therapy"  = 	"#4790db",	"HSP inhibitor"  = "#fded9c",	"Hedgehog inhibitor"  = "#c9c9fb",	"Kinesin inhibitor"  = "#49bd93",	"Rapalog"  = "#969b96")

# export to svg
ggsave(
  "drug_classes_overview.svg",
  plot = p,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = 7,
  height = 4,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
