#' Connect and parse UniProt information.
#'
#' The function is work to retrieving GetpdbStructure and download it to user directory.
#'
#' @param ProteinAccList input a vector of UniProt Accession/s
#'
#' @param directorypath path to save excel file containig results returened by the function.
#'
#' @usage GetpdbStructure(ProteinAccList , directorypath = NULL)
#'
#'
#' @export
#'
#' @author Mohmed Soudy and Ali Mostafa
#'
GetpdbStructure <- function(ProteinAccList , directorypath = NULL)
{
  baseUrl <- "https://swissmodel.expasy.org/repository/uniprot/"
  for (ProteinAcc in ProteinAccList)
  {
    RequiredAcc <- paste0(baseUrl , ProteinAcc ,".pdb")
    Request_Response <- GET(RequiredAcc)
    if(Request_Response == 200)
    {
      download.file(RequiredAcc , paste0(ProteinAcc,".pdb"))
    }
    else{
      HandleBadRequests(Request_Response)
    }
  }
}
