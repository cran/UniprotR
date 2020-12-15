#' Connect and parse UniProt information
#'
#' This function is used for Enrichment analysis of given list of genes or proteins 
#'
#' @usage Enrichment.gprofiler(Accs, sources=NULL, OS="hsapiens", p_value=0.05, directorypath=NULL)
#'
#' @param Accs Vector of UniProt Accession/s or genes 
#' 
#' @param sources a vector of data sources to use, these include GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC
#' 
#' @param OS  organism name Example: human - 'hsapiens', mouse - 'mmusculus'
#' 
#' @param p_value custom p-value threshold for significance, default = 0.05
#'
#' @export
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@57357.com} and Ali Mostafa \email{ali.mo.anwar@std.agr.cu.edu.eg}
Enrichment.gprofiler <- function(Accs, sources = NULL, OS = "hsapiens", p_value = 0.05)
{
  AccList <- as.character(unique(Accs))
  Enrich.object <- gost(query = Accs, sources = sources, organism = OS, user_threshold = p_value)
  
  gostplot(Enrich.object, capped = TRUE, interactive = T)
}
#' Connect and parse UniProt information
#'
#' This function is used for Enrichment analysis of given list of genes or proteins based on Reactome API 
#'
#' @usage Enrichment.React(AccList, organism = "human", pvalueCutoff = 0.05, directorypath = NULL)
#'
#' @param AccList Vector of UniProt Accession/s or genes 
#' 
#' @param organism one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
#' 
#' @param pvalueCutoff Cutoff value of pvalue, default = 0.05 
#'
#' @param directorypath path to save excel file and plot returened by the function
#'
#' @export
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@57357.com} and Ali Mostafa \email{ali.mo.anwar@std.agr.cu.edu.eg}
Enrichment.React <- function(AccList, organism = "human", pvalueCutoff = 0.05, directorypath = NULL)
{
  AccList <- as.character(unique(AccList))
  
  Entrez <- ConvertID(ProteinAccList = AccList, ID_to = "P_ENTREZGENEID")
  Entrez <- as.character(na.omit(as.numeric(as.character(Entrez$`To P_ENTREZGENEID`))))
  Enrichment <- enrichPathway(gene=Entrez, pvalueCutoff = pvalueCutoff, organism = organism)
  Pathways <- Enrichment@result
  
  TopPathways <- Pathways[order(Pathways$Count , decreasing = T),]
  TopPathways <- TopPathways[1:10,]
  
  
  p <- ggplot(TopPathways, aes(y= reorder(TopPathways$Description , TopPathways$Count), x= TopPathways$Count , fill = `p.adjust`)) +
    geom_bar(stat="identity") + theme_classic()  + ylab("Description") + xlab("Count") 
  plot(p)
  if(!is.null(directorypath))
  {
    ggsave(plot = p, filename = paste0(directorypath,"/Top 10 pathways.jpeg"), device = "jpeg", width = 8, height = 6, dpi = 300)
    write.csv(x = Pathways, file = paste0(directorypath,"/Pathways.csv"))
  }
}