#' Connect and parse UniProt General Information.
#'
#' The function is work to retrieve General Information data from UniProt
#' for a list of proteins accessions.For more information about what included in the 
#' General Information data see https://www.uniprot.org/help/uniprotkb_column_names.
#'
#' @usage GetGeneral_Information(ProteinAccList , directorypath = NULL)
#'
#' @param ProteinAccList Vector of UniProt Accession/s
#'
#' @param directorypath path to save excel file containig results returened by the function.
#'
#' @return DataFrame where rows names are the accession
#'      and columns contains the General Information of protein from the UniProt
#'
#' @note The function also, Creates a csv file with the retrieved information.
#'
#' @examples Obj <- GetGeneral_Information("O14520")
#'
#' @export
#'
#' @author Mohmed Soudy and Ali Mostafa

GetGeneral_Information <- function(ProteinAccList , directorypath = NULL)
{
  ProteinInfoParsed_total = data.frame()
  Colnames = "created,last-modified,sequence-modified,version(entry)"
  baseUrl <- "http://www.uniprot.org/uniprot/"
  for (ProteinAcc in ProteinAccList)
  {
    #to see if Request == 200 or not
    Request <- GET(paste0(baseUrl , ProteinAcc,".xml"))

    #this link return information in tab formate (format = tab)
    ProteinName_url <- paste0("?query=accession:",ProteinAcc,"&format=tab&columns=",Colnames)
    RequestUrl <- paste0(baseUrl , ProteinName_url)
    if (Request$status_code == 200){
      # parse the information in DataFrame
      ProteinDataTable <- tryCatch(read.table(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
      if (!is.null(ProteinDataTable))
      {
        ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
        # add Dataframes together if more than one accession
        ProteinInfoParsed_total <- rbind(ProteinInfoParsed_total, ProteinInfoParsed)
      }
    }else {
      HandleBadRequests(Request$status_code)
    }
  }
  if(!is.null(directorypath))
  {
    write.csv(ProteinInfoParsed_total , paste0(directorypath , "/" , "General Information.csv"))
  }
  return(ProteinInfoParsed_total)
}