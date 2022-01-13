getColorScheme <- function () {
  colorScheme <- list()
  
  colorPanel1                                <- c("#8dd3c7", "#b3b300", "#bebada", "#fb8072", "#80b1d3",
                                                  "#fdb462", "#b3de69", "#fccde5")
  names(colorPanel1)                         <- c(1:8)
  
  colorPanel2                                <- c("#75b4c0", "#438644", "#8bc2a1", "#f0c0b5", "#ebe5e6",
                                                  "#216100", "#8ab22e", "#bc8e46", "#e3d9ad", "#633e5a",
                                                  "#cb95aa", "#fdae55", "#658cc8", "#007d92", "#62a399",
                                                  "#fd7e52", "#f3166b", "#ab4252", "#fea617", "#f9c099")
  names(colorPanel2)                         <- c(1:20)
  
  colorScheme$blue                           <- "deepskyblue3"
  colorScheme$white                          <- "white"
  colorScheme$red                            <- "firebrick3"
  
  colorScheme$green                          <- "#b7dd75"
  colorScheme$brick                          <- "#e07b39"
  colorScheme$lagoon                         <- "#69bdd2"
  
  colorScheme$samplesClusters$kmeans         <- colorPanel1
  colorScheme$samplesClusters$hierar         <- colorPanel1
  
  colorScheme$featureClusters$hierar         <- colorPanel2
  
  colorScheme$type                           <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#7c827c")
  names(colorScheme$type)                    <- c("T-ALL", "B-ALL", "BCP-ALL", "EBV")
  
  colorScheme$subtype                        <- c("#71BC78", "#FF9BAA", "#7c827c", "#FD5E53", "#FD7C6E",
                                                  "#1974D2", "#1DACD6", "#1CD3A2", "#1CAC78", "#FCD975",
                                                  "#FCE883", "#FDFC74", "#bbbbfa", "#CDA4DE", "#A2ADD0",
                                                  "#D68A59", "#ECEABE", "#9D81BA", "#FFA089", "#FC6C85",
                                                  "#EFCDB8", "#F664AF", "#E6A8D7", "#FC6C85", "gray",
                                                  "#5C80BC", "#CDD1C4", "#E0AFA0", "#1B264F", "#B5CA8D", "#885053")
  names(colorScheme$subtype)                 <- c("BCL11B.TLX3", "BCR.ABL1", "EBV", "ETV6.PDGFRB", "ETV6.RUNX1",      
                                                  "Hypodiploid", "Hypodiploid.LMO2", "IGH.CRLF2", "IGH.MYC", "KMT2A.AFF1",
                                                  "KMT2A.FOXO4", "KMT2A.MLLT1", "LMO1.TCRD", "LMO2", "LMO2.STAG2",
                                                  "MEF2D.HNRNPUL1", "Near.haploid", "NUP214.ABL1", "PAX5.ETV6", "SFPQ.ABL",
                                                  "SIL.SCL", "TCF3.HLF", "TCF3.PBX1", "TRA.MYC", "Unknown",
                                                  "ABL1.ZMIZ1", "BUB1B+", "PTMA.TMSB4X", "TRB.LCK", "SFPQ.ABL1", "CERS2+/IL32+")
  
  colorScheme$subtype.alt                    <- c("#71BC78", "#FF9BAA", "#7c827c", "#FD5E53", "#FD7C6E",
                                                  "#1974D2", "#1DACD6", "#1CD3A2", "#1CAC78", "#FCD975",
                                                  "#FCE883", "#FDFC74", "#bbbbfa", "#CDA4DE")
  names(colorScheme$subtype.alt)             <- c("ABL1-ZMIZ1", "B-Other", "BCL11B-TLX3", "BCR-ABL1", "EBV",  
                                                  "ETV6-RUNX1", "IGH-MYC",  "KMT2A-AFF1", "KMT2A-MLLT1", "MEF2D-HNRNPUL1", 
                                                  "NUP214-ABL1", "PAX5-ETV6", "T-Other", "TCF3-PBX1")
  
  colorScheme$Class.explained                <- c("#71BC78", "#FF9BAA", "#7c827c", "#FD5E53", 
                                                  "#FD7C6E", "#1974D2", "#1DACD6", "#1CD3A2", 
                                                  "#1CAC78", "#FCD975", "#FCE883", "#FDFC74", 
                                                  "#bbbbfa", "#CDA4DE", "#A2ADD0")
  names(colorScheme$Class.explained)         <- c("A.Conv.Chemo", "B.Kinase inhibitor", "C.Rapalog", "D.Immunomodulatory",
                                                  "E.Differentiating.epigenetic modifier", "F.Hormone therapy", "G.Apoptotic modulator", "H.Metabolic modifier",
                                                  "I.Kinesin inhibitor", "J.NSAID", "K.HSP inhibitor", "L.Protease.proteasome inhibitor",
                                                  "M.Hedgehog inhibitor", "N.Conv.Chemo Combo", "X.Other")
  
  colorScheme$tissue                         <- c("#ff9999", "#e6ccb3", "#ff9999")
  names(colorScheme$tissue)                  <- c("PB", "BM", "PE")
  
  colorScheme$gender                         <- c("#ff99c2", "#99ccff")
  names(colorScheme$gender)                  <- c("F", "M")
  
  colorScheme$age                            <- c("#99ebff", "#bbbb77") # for gradient
  
  colorScheme$AssignedStages                 <- c("pre-B other" = "#fb8072", "pre-pro-B" = "#80b1d3",
                                                  "pro-B" = "#fdb462", "early pre-B" = "#fccde5", "late pre-B" = "#b3b300", 
                                                  "DNTT+ B-like" = "gray90", "immature B" = "gray20")
  
  return(colorScheme)
}