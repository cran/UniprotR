% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetStructureInfo.R
\name{GetStructureInfo}
\alias{GetStructureInfo}
\title{Connect and parse UniProt protein Structure information.}
\usage{
GetStructureInfo(ProteinAccList, directorypath = NULL)
}
\arguments{
\item{ProteinAccList}{Vector of UniProt Accession/s.}

\item{directorypath}{path to save excel file containig results returened by the function.}
}
\value{
DataFrame where rows names are the accession
     and columns contains the Structural information of protein from the UniProt
}
\description{
The function is work to retrieve Structral data from UniProt for a list of proteins accessions.
For more information about what included in the structral 
data see https://www.uniprot.org/help/return_fields.
}
\note{
The function also, Creates a csv file with the retrieved information.
}
