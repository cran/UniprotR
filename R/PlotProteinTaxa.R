#' Connect and parse UniProt information.
#'
#' This Function is used to plot different taxas found of the accessions.
#'
#' @usage PlotProteinTaxa(ProteinsData , directorypath = NULL)
#'
#' @param ProteinsData input a Dataframe of proteins as rownames.
#'
#' @param directorypath path to save excel file containig results returened by the function.
#'
#'
#' @author Mohmed Soudy and Ali Mostafa
#'
#' @export
#'
PlotProteinTaxa <- function(ProteinsData , directorypath = NULL)
{
  TaxaPlot <- ggplot(ProteinsData , aes(x = rownames(ProteinsData) , fill = ProteinsData$Organism)) +
    geom_bar() + theme(axis.text.x  = element_text(angle = 90 , hjust = 1 , vjust = 0.5) , axis.title.y = element_blank() , axis.text.y = element_blank()) + xlab("Protein Accession") +
    guides(fill=guide_legend(title="Orgainsims"))
  plot(TaxaPlot)
  if (!is.null(directorypath)){
  ggsave(paste0(directorypath, "//" ,"Proteins Taxonomy.png") ,plot = TaxaPlot , device = "png")
  }
}

#' Connect and parse UniProt information.
#'
#' This Function is used to plot different locations where accession/s can be found in.
#'
#' @usage PlotProteinsLoc(ProteinsData, directorypath = NULL)
#'
#' @param ProteinsData input a Dataframe of proteins as rownames.
#'
#' @param directorypath path to save excel file containig results returened by the function.
#'
#' @author Mohmed Soudy and Ali Mostafa
#'
#' @export
#'
PlotProteinsLoc <- function(ProteinsData, directorypath = NULL)
{
  TaxaPlot <- ggplot(ProteinsData , aes(x = rownames(ProteinsData) , fill = ProteinsData$Proteomes)) +
    geom_bar() + theme(axis.text.x  = element_text(angle = 90 , hjust = 1 , vjust = 0.5) , axis.title.y = element_blank() , axis.text.y = element_blank()) + xlab("Protein Accession") +
    guides(fill=guide_legend(title="Chromosomes"))
  TaxaPlot
  if (!is.null(directorypath)){
  ggsave(paste0(directorypath,"/","Proteins Locations.png") ,plot = TaxaPlot , device = "png")
  }
}
