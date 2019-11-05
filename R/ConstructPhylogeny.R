#' Connect and parse UniProt information.
#'
#' This Function is used to construct phylogenetic tree from protein sequence applying msa on the sequence list.
#'
#' @usage ConstructPhylogeny(ProteinDataObject , directorypath = NULL)
#'
#' @param ProteinDataObject  Dataframe returned from GetSequence function.
#'
#' @param directorypath path to save files containig results returened by the function ( default = NA ).
#'
#' @note To use this package you MUst INSTALL msa package from Bioconductor ,if no dir_path was given ( default = NA ) the function will only view the plot and will not save it
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@57357.com} and Ali Mostafa \email{ali.mo.anwar@std.agr.cu.edu.eg}
#'
#' @export
#'
ConstructPhylogeny <- function(ProteinDataObject , directorypath = NULL)
{
  #Get sequences
  Seqlist <- as.character(ProteinDataObject$Sequence)

  #Apply multiple sequence alignment
  MSAresult <- msa(Seqlist , type = "protein")
  #Convert algniment results to tree
  MSAtree <- msaConvert(MSAresult, type="seqinr::alignment")
  # generate a distance matrix using seqinr package
  MSAdistance <- dist.alignment(MSAtree, "identity")
  #neighbor-joining tree estimation
  mTree <- nj(MSAdistance)
  mTree$tip.label <- as.character(rownames(ProteinDataObject))

  # pile up the functions to make a new tree
  NTree <- nj(dist.alignment(MSAtree, "identity"))
  NTree$tip.label <- mTree$tip.label
  if (!is.null(directorypath))
  {
    pdf(paste0(directorypath , "/" ,"Phylogenymod.pdf") , width = 10 , height = 13)
  }
  #Plot new tree
  ModTree <- plot(NTree, "f", FALSE, cex = 0.7 , main="Phylogenetic Tree of proteins")

  RegularTree <- plot(mTree)
  dev.off()
}
