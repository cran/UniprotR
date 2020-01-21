## -----------------------------------------------------------------------------
library(UniprotR)
# User can Open a browser to select CVS or txt contains accession list or Enter accessions manually
#Select accession from CSV or txt
#Accessions <- GetAccessionList(file.choose()) 
#Enter Accessions manually 
Accessions <- c("Q8TAP6", "Q9P209", "P40121", "Q8NEF3", "Q5VTM2", "P24864", "Q9UHD1", "P09564",
                "P09603", "Q9BS18", "O15078","Q9BZP3", "Q5M9N0", "Q96LT6", "Q9HCU4", "O43852",  "Q9BWT7", "Q53HC0", "O15234", "Q96KN2", "P63098", "Q9WVC3","Q8BGK2", "Q70KF4", "Q6IR41", "Q9QXJ4", "Q8C3W1", "Q8K2W9", "O88870", "Q3V3N7", "Q8BGD0", "Q3V037", "Q8CG72","Q8CDN1", "Q9QY84", "Q9CWF6", "Q61271", "Q8C142", "P04186", "Q6EBV9", "Q9WUE3", "Q9D1D6", "Q9WUB6", "Q8C1Z7","P09470", "Q99M07", "Q8BHN7", "Q9QY83", "Q91ZV8", "Q4ZJN1", "Q80TI0", "Q4V7F0", "Q5XIS7", "P30120", "Q9QZI7","Q2TGJ1", "Q6QMY6", "P08542", "Q9JJW1", "Q6AY86", "Q6AYE4", "Q9Z1F2", "A2VCW5", "Q6GV27", "Q62941", "Q9JKY3","Q63665", "Q5J3E5", "P60841", "Q9Z288", "Q6P6U0", "P17425", "Q62661")

## ----fig.width=10 , fig.height=8----------------------------------------------
Taxa <- GetNamesTaxa(Accessions , getwd()) #getwd() Could be replaced with user's path where csv file contains Names & Taxonomy will be  saved
#Plot Names and taxanomic information
PlotProteinTaxa(Taxa , getwd()) #User should specify path to save plots 

## ----fig.width=12 , fig.height=9----------------------------------------------
Miscellaneous <- GetMiscellaneous(Accessions , getwd())
#Plot protein scoring
PlotproteinScore(Miscellaneous , getwd())

## ----fig.width=12 , fig.height=9----------------------------------------------
SeqInfo <- GetSequences(Accessions , getwd())
#Plot summary of GO Information
PlotPhysicochemical(SeqInfo , getwd())

## -----------------------------------------------------------------------------
GetproteinNetwork(Accessions , getwd())

