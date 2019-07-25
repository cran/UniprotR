#' Connect and parse UniProt Sequences information.
#'
#' The function is work to retrieve Sequences data from UniProt for a list of proteins accessions.
#' For more information about what included in the Sequences data
#' see https://www.uniprot.org/help/uniprotkb_column_names.
#'
#' @usage GetSequences(ProteinAccList, directorypath = NULL)
#'
#' @param ProteinAccList Vector of UniProt Accession/s
#'
#' @param directorypath path to save excel file containig results returened by the function
#'
#' @return DataFrame where rows names are the accession
#'      and columns contains the information retrieved from the UniProt
#'
#' @examples Obj <- GetSequences("O14520")
#'
#' @note The function also, Creates a csv file with the retrieved information.
#'
#' @author Mohmed Soudy and Ali Mostafa
#'
#' @export

GetSequences <- function(ProteinAccList, directorypath = NULL){

  # Sequences information to be collected
  columns <- c("fragment","encodedon","comment(ALTERNATIVE PRODUCTS)",
               "comment(ERRONEOUS GENE MODEL PREDICTION)",
               "comment(ERRONEOUS INITIATION)","comment(ERRONEOUS TERMINATION)",
               "comment(ERRONEOUS TRANSLATION)","comment(FRAMESHIFT)",
               "comment(MASS SPECTROMETRY)","comment(POLYMORPHISM)",
               "comment(RNA EDITING)","comment(SEQUENCE CAUTION)",
               "length","mass","sequence","feature(ALTERNATIVE SEQUENCE)",
               "feature(NATURAL VARIANT)","feature(NON ADJACENT RESIDUES)",
               "feature(NON STANDARD RESIDUE)","feature(NON TERMINAL RESIDUE)",
               "feature(SEQUENCE CONFLICT)","feature(SEQUENCE UNCERTAINTY)",
               "version(sequence)")


  baseUrl <- "http://www.uniprot.org/uniprot/"
  ProteinInfoParsed_total = data.frame()
  for (ProteinAcc in ProteinAccList)
  {
    #to see if Request == 200 or not
    Request <- GET(paste0(baseUrl , ProteinAcc,".xml"))

    #this link return information in tab formate (format = tab)
    #columns = what to return from all of the information (see: https://www.uniprot.org/help/uniprotkb_column_names)
    ProteinName_url <- paste0("?query=accession:",ProteinAcc,"&format=tab&columns=",columns)

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
    write.csv(ProteinInfoParsed_total , paste0(directorypath, "/" ,"Sequences Information.csv"))
  }
  return(ProteinInfoParsed_total)
}