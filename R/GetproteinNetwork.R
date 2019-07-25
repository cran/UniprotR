#' Connect and parse stringdb information.
#'
#' This function is connecting to stringdb and retrieve all possible interactions
#' for the searched protein/s.
#'
#' @usage GetproteinNetwork(ProteinAccList , directorypath = NULL)
#'
#' @param ProteinAccList input a vector of UniProt Accession/s
#'
#' @param directorypath path to save excel file containig results returened by the function.
#'
#' @usage GetproteinNetwork(ProteinAccList)
#'
#' @author Mohmed Soudy and Ali Mostafa
#'
#' @export
GetproteinNetwork <- function(ProteinAccList , directorypath = NULL)
{
  pdf(paste0(directorypath , "/","Protin Network.pdf"))
  baseUrl <- "https://string-db.org/api/image/network?identifiers="
  for (identifier in ProteinAccList)
  {
    ProteinString <- paste0(baseUrl , identifier)
    Request <- GET(ProteinString)
    if (Request$status_code == 200)
    {
      Network <- image_read(ProteinString)
      plot(Network)
    }else{
      dev.off()
      HandleBadRequests(Request$status_code)
      pdf(paste0(directorypath , "/","Protin Network.pdf"))
    }
  }
  ProteinList <- paste0(baseUrl,ProteinAccList , collapse = "%0d" , "&add_color_nodes=20&network_flavor=actions&block_structure_pics_in_bubbles=1")
  WholeNetwork <- image_read(ProteinList)
  plot(WholeNetwork)
  dev.off()
}
