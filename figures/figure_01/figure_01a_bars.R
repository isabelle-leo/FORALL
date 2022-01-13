library(ggplot2)
library(scales)
library(ggthemes)
library(svglite)

source(file.path("figures", "all_theme.R"))
theme_set(all_theme())

width <- 5
height <- 4.5

# set output folder
outputFolder <- file.path("figures", "output", "figure_01")

# create major output folder, if it doesn't exist
if (!is.null(outputFolder) && !dir.exists(outputFolder)) dir.create(outputFolder)

# peptide table
df_a <- as.data.frame(c(279351-95644, 95644))
colnames(df_a) <- "Peptides"
df_a$Type <- c("Peptides", "Peptides")
df_a$Peptide <- c(">1 sets", "All sets")

#protein table
df_b <- as.data.frame(c(12446, 9100, 8981))
colnames(df_b) <- "Gene Symbols"
df_b$Type <- c("Total Proteins","Proteins \n in all sets","Overlap \n Proteins-mRNA")

#mrna table
df_c <- as.data.frame(c(55239, 19583, 16213, 13140, 1307))
colnames(df_c)<- "mRNA"
df_c$Type <- c("Transcripts","Protein Coding","lncRNA","Pseudogenes","miRNA")
df_c$label <- c("White","White","White","White","Black")

#figures
p_peps <- ggplot(df_a,aes(factor(Type),Peptides,fill=factor(Peptide,c(">1 sets","All sets"))))+
  geom_bar(stat="identity")+
  geom_text(label= df_a$Peptides, color = "White", vjust=2)+
  scale_fill_manual(values = c("#742c55","#443158"),name="Quantification")+
  xlab(element_blank())+
  scale_y_continuous(labels = comma)

ggsave(
  filename = "bars_peptides.svg",
  plot = p_peps,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width + 0.7,
  height = height,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)


p_prots <- ggplot(df_b,aes(factor(Type,levels = c("Total Proteins","Proteins \n in all sets","Overlap \n Proteins-mRNA")),`Gene Symbols`,fill=factor(Type,levels = c("Overlap \n Proteins-mRNA","Proteins \n in all sets","Total Proteins"))))+
  geom_bar(stat="identity",show.legend = F)+
  geom_text(label= df_b$`Gene Symbols`, color = "White", vjust=1.5)+
  scale_fill_manual(values = c("#474c64","#443158","#742c55"))+
  xlab(element_blank())

ggsave(
  filename = "bars_proteins.svg",
  plot = p_prots,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width - 1.3,
  height = height,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)


p_mrna <- ggplot(df_c,aes(factor(Type,levels = c("Transcripts","Protein Coding","lncRNA","Pseudogenes","miRNA")),mRNA,fill=factor(Type,levels =  c("miRNA","Pseudogenes","lncRNA","Protein Coding","Transcripts"))))+
  geom_bar(stat="identity",show.legend = F)+
  geom_text(data = df_c[df_c$mRNA != "1307",], label = df_c$mRNA[1:4],color = "white", vjust=1.5)+
  geom_text(data= df_c[df_c$mRNA == "1307",], label = df_c$mRNA[5],color = "black", vjust=-0.5)+
  scale_fill_manual(values = c("#f3e5ab","#cad7ae","#a2cab2","#7abcb5","#28a1bc"))+
  xlab(element_blank())

ggsave(
  filename = "bars_mrna.svg",
  plot = p_mrna,
  device = NULL,
  path = outputFolder,
  scale = 1,
  width = width,
  height = height,
  units = c("cm"),
  dpi = 600,
  limitsize = FALSE)
