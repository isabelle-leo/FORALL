makeExpressionSet <- function (expressionsFileName,
                               type = "protein",
                               mode = "gene_symbol_centric",
                               quantRegex,
                               geneSymbolColumn = "Protein accession",
                               metaFileName,
                               metaSheetName,
                               tmtString = "_tmt10plex_") {
  
  # load raw expression dataset
  expressions.df <- getData(fileName = expressionsFileName)
  
  # load meta data
  metaRaw <- read_excel(path = metaFileName,
                        sheet = metaSheetName) %>%
    dplyr::rename(proteomics_id = `Proteomics sample`) %>%
    dplyr::rename(rnaseq_id = `RNA-seq sample`) %>%
    dplyr::rename(cell_line = `Cell Line Name`) %>%
    dplyr::rename(Cell.Line.R = `Cell Line.R`) %>%
    dplyr::mutate(Age = gsub(",", "\\.", Age)) %>%
    dplyr::mutate(Age = gsub("-", NA, Age)) %>%
    dplyr::mutate(isReplicate = as.logical(isReplicate)) %>%
    dplyr::select(proteomics_id, everything())
  
  # expand rows with multiple proteomics_ids
  multiRows <- grep(";", metaRaw$proteomics_id)
  allNewRows <- lapply(X = multiRows,
                       FUN = function (row, metaRaw) {
                         newRowIds <- strsplit(x = metaRaw$proteomics_id[[row]],
                                               split = ";",
                                               fixed = TRUE)[[1]]
                         
                         content <- metaRaw[row, ] %>%
                           dplyr::select(-proteomics_id)
                         
                         newRows <- lapply(X = newRowIds,
                                           FUN = function (id, content) {
                                             data.frame(proteomics_id = id,
                                                        content,
                                                        stringsAsFactors = FALSE,
                                                        check.names = FALSE)
                                           },
                                           content = content)
                       },
                       metaRaw = metaRaw)
  
  df.intermediate <- do.call(what = rbind,
                             args = unlist(allNewRows, recursive = FALSE))
  
  meta <- rbind(metaRaw %>% filter(!grepl(";", proteomics_id)), df.intermediate) %>%
    dplyr::mutate(position = proteomics_id) %>%
    as.data.frame()
  
  row.names(meta) <- meta$position
  
  # filter data for accession, description and tmt columns
  if (type == "protein" & mode == "gene_symbol_centric") {
    expressions.df.filt <- expressions.df %>%
      dplyr::select(c("gene_symbol" = !!geneSymbolColumn),
                    c("description" = Description),
                    c("protein_id" = `Protein ID(s)`),
                    matches(quantRegex),
                    ends_with("_quanted_psm_count"))
    
    # set the rownames
    row.names(expressions.df.filt) <- expressions.df.filt$gene_symbol
  } else if (type == "protein" & mode == "gene_id_centric") {
    expressions.df.filt <- expressions.df %>%
      dplyr::select(c("gene_id" = !!geneSymbolColumn),
                    c("description" = Description),
                    matches(quantRegex),
                    ends_with("_quanted_psm_count"))
    
    # get gene symbols via biomaRt
    # start biomaRt
    ensembl <- useEnsembl(biomart = "ensembl",
                          dataset = "hsapiens_gene_ensembl")
    
    # get translations to gene names
    lookup.raw <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = expressions.df.filt$gene_id,
                        mart = ensembl)
    
    # make a lookup
    lookup <- lookup.raw$external_gene_name
    names(lookup) <- lookup.raw$ensembl_gene_id
    
    # get gene names in tpm table
    expressions.df.filt$gene_symbol <- lookup[expressions.df.filt$gene_id]
    
    # set the rownames
    row.names(expressions.df.filt) <- expressions.df.filt$gene_id
  } else if (type == "peptide" & mode == "gene_symbol_centric") {
    expressions.df.filt <- expressions.df %>%
      dplyr::select(c("peptide" = `Peptide sequence`),
                    c("gene_id" = `Gene(s)`),
                    c("gene_symbol" = `Associated gene ID(s)`),
                    c("description" = `Description(s)`),
                    c("protein_id" = `Protein(s)`),
                    matches(quantRegex),
                    ends_with("_quanted_psm_count")) %>%
      filter(!grepl(";", gene_symbol)) %>% # no multiple gene_symbols allowed
      dplyr::mutate(gene_symbol = ifelse(gene_symbol == "", "nogenesymbol", gene_symbol))
    
    # set the rownames
    row.names(expressions.df.filt) <- expressions.df.filt$peptide
  } else {
    stop("Type and/or mode not supported!")
  }
  
  # extract minimum PSM counts
  expressions.df.filt$min_quanted_psms <- expressions.df.filt %>%
    dplyr::select(ends_with("_quanted_psm_count")) %>%
    dplyr::mutate_all(~ ifelse(. == 0, NA, .)) %>%
    as.matrix() %>%
    rowMins(na.rm = TRUE)
  
  # change column names to exclude sample names
  colnames(expressions.df.filt) <- gsub("^.*?_Set", "Set", colnames(expressions.df.filt))
  
  # get quantification colnames and numbers
  quantColsNum <- grep(quantRegex, colnames(expressions.df.filt))
  quantCols <- colnames(expressions.df.filt)[quantColsNum]
  
  # order meta rows according to protein columns
  meta.ord <- meta[quantCols, ]
  
  # check for consistency between meta file and proteins file
  if (!identical(row.names(meta.ord), quantCols)) stop("Inconsistency between meta and protein data!")
  
  # build up eSet
  expressions.raw <- ExpressionSet(assayData = as.matrix(expressions.df.filt[, quantColsNum]),
                                   phenoData = AnnotatedDataFrame(data = meta.ord,
                                                                  varMetadata = data.frame(labelDescription = colnames(meta.ord))),
                                   featureData = AnnotatedDataFrame(data = expressions.df.filt[, -quantColsNum],
                                                                    varMetadata = data.frame(labelDescription = colnames(expressions.df.filt[, -quantColsNum]))))
  
  return(expressions.raw)
}
