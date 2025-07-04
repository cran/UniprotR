#' Connect and parse UniProt PTM_Processsing information.
#'
#' The function is work to retrieve PTM_Processsing data from UniProt for a list of proteins accessions.
#' For more information about what included in the PTM_Processsing data
#' see https://www.uniprot.org/help/return_fields.
#'
#' @usage GetPTM_Processing(ProteinAccList, directorypath = NULL)
#'
#' @param ProteinAccList Vector of UniProt Accession/s
#'
#' @param directorypath path to save excel file containig results returened by the function
#'
#' @return DataFrame where rows names are the accession
#'      and columns contains the information retrieved from the UniProt
#'
#' @examples Obj <- GetPTM_Processing("O14520" )
#'
#' @note The function also, Creates a csv file with the retrieved information.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@57357.com} and Ali Mostafa \email{ali.mo.anwar@std.agr.cu.edu.eg}
#'
#' @export


GetPTM_Processing<- function(ProteinAccList, directorypath = NULL ){
  
  if(!has_internet())
  {
    message("Please connect to the internet as the package requires internect connection.")
    return()
  }
  options(timeout=10000)
  message("Please wait we are processing your accessions ...")
  pb <- progress::progress_bar$new(total = length(ProteinAccList))
  # PTM_Processsing information to be collected
  columns <- c("ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,cc_ptm,ft_propep,ft_signal,ft_transit")
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  ProteinInfoParsed_total = data.frame()
  for (ProteinAcc in ProteinAccList)
  {
    #to see if Request == 200 or not
    Request <- tryCatch(
      {
        GET(paste0(baseUrl , ProteinAcc,"&format=tsv") , timeout(7))
      },error = function(cond)
      {
        message("Internet connection problem occurs and the function will return the original error")
        message(cond)
      }
    )
    #this link return information in tab formate (format = tab)
    #columns = what to return from all of the information (see: https://www.uniprot.org/help/uniprotkb_column_names)
    ProteinName_url <- paste0(ProteinAcc,"&format=tsv&fields=",columns)
    
    RequestUrl <- paste0(baseUrl , ProteinName_url)
    RequestUrl <- URLencode(RequestUrl)
    if (length(Request) == 0)
    {
      message("Internet connection problem occurs")
      return()
    }
    if (Request$status_code == 200){
      # parse the information in DataFrame
      ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
      if (!is.null(ProteinDataTable))
      {
        ProteinDataTable <- ProteinDataTable[1,]
        ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
        # add Dataframes together if more than one accession
        ProteinInfoParsed_total <- rbind(ProteinInfoParsed_total, ProteinInfoParsed)
      }
    }else {
      HandleBadRequests(Request$status_code)
    }
    pb$tick()
    
  }
  if(!is.null(directorypath))
  {
    write.csv(ProteinInfoParsed_total ,paste0(directorypath,"/","PTM_Processsing Information.csv"))
  }
  return(ProteinInfoParsed_total)
}
