#' Connect and parse UniProt information.
#'
#' This Function is used to plot Gene ontolgy summary in the data of the accession/s.
#'
#' @usage plotGoannotation(ProteinDataObject,directorypath = NULL)
#'
#' @param ProteinDataObject input a Dataframe returned from GetProteinGOInfo function
#'
#' @param directorypath path to save files returened by the function.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@57357.com} and Ali Mostafa \email{ali.mo.anwar@std.agr.cu.edu.eg}
#'
#' @export

plotGoannotation <- function(ProteinDataObject, directorypath = NULL) {
  
  Goparse <- function(GOObj , index = 3)
  {
    GO_df_obj_bio <- toString(na.omit(GOObj[,index]))
    
    GO_df_obj_bio <- strsplit(GO_df_obj_bio,";|,")
    
    GO_df_obj_bio_df <- data.frame(GO_df_obj_bio)
    
    colnames(GO_df_obj_bio_df) <- "bio"
    
    
    trim <- function (x) gsub("^\\s+|\\s+$", "", x)
    
    GO_df_obj_bio_df$bio <- trim(GO_df_obj_bio_df$bio)
    
    test1 <- strsplit(as.character( GO_df_obj_bio_df$bio ), ".\\[")
    
    test2 <- lapply(test1, function(x) x[[1]][1])
    
    occurences<-table(unlist(test2))
    
    occurences<- as.data.frame(occurences)
    
    occurences <-  occurences[order(-occurences$Freq),]
    
    colnames(occurences) <- c("Goterm","Frequences")
    
    occurences %>%
      mutate(freq = percent(occurences$Freq / length(rownames(GOObj)))) -> occurences
    return(occurences)
  }
  # Get GO annotation data for each category
  BiologicalDF <- Goparse(ProteinDataObject, 3)[1:10,]
  MolecularDF <- Goparse(ProteinDataObject, 4)[1:10,]
  CellularDF  <- Goparse(ProteinDataObject, 5)[1:10,]
  
  # Add 'group' labels
  BiologicalDF$group <- "Biological process"
  MolecularDF$group  <- "Molecular function"
  CellularDF$group   <- "Subcelluar localization"
  
  # Combine all and ensure correct columns
  GoOntology <- bind_rows(BiologicalDF, MolecularDF, CellularDF)
  GoOntology <- na.omit(GoOntology)
  GoOntology <- GoOntology[,c("Goterm", "Frequences", "group")]
  
  # Ensure data types are correct
  GoOntology$Goterm    <- as.character(GoOntology$Goterm)
  GoOntology$group     <- factor(GoOntology$group, levels = c("Biological process", "Molecular function", "Subcelluar localization"))
  GoOntology$Frequences <- as.numeric(GoOntology$Frequences)
  
  # Arrange data and create 'empty bars' for visual separation
  GoOntology <- GoOntology %>% arrange(group, desc(Frequences))
  empty_bar  <- 3
  to_add     <- data.frame(Goterm=rep(NA, empty_bar * nlevels(GoOntology$group)),
                           Frequences=NA,
                           group=rep(levels(GoOntology$group), each=empty_bar))
  data <- bind_rows(GoOntology, to_add) %>% arrange(group)
  data$id <- seq_len(nrow(data))
  # Calculate angles and hjust for circular text
  number_of_bar <- nrow(data)
  angle <- 90 - 360 * (data$id - 0.5) / number_of_bar
  data$hjust <- ifelse(angle < -90, 1, 0)
  data$angle <- ifelse(angle < -90, angle + 180, angle)
  
  # Prepare data for group base lines and grid
  base_data <- data %>%
    group_by(group) %>%
    summarize(start = min(id), end = max(id) - empty_bar, .groups = "drop") %>%
    mutate(title = mean(c(start, end)))
  
  grid_data <- base_data
  grid_data$end <- c(grid_data$end[nrow(grid_data)], head(grid_data$end, -1)) + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  # Only show labels on real (non-empty) bars
  real_label_data <- data[!is.na(data$Frequences), ]
  
  # Set offsets for label position
  offset <- max(real_label_data$Frequences) * 0.04
  bar_label <- NULL
  hjust <- NULL 
  group <- NULL
  Frequences <- NULL
  id <- NULL
  # -----
  # Build the label as "GO term (Count)"
  real_label_data$bar_label <- paste0(real_label_data$Goterm, " (", real_label_data$Frequences, ")")
  
  # Plot
  p <- ggplot(data, aes(x = as.factor(id), y = Frequences, fill = group)) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.85) +
    # ... (grid lines etc)
    geom_text(
      data = real_label_data,
      aes(x = id, y = Frequences + offset, label = bar_label, angle = angle, hjust = hjust),
      color = "black", fontface = "bold", alpha = 0.85, size = 3.2, inherit.aes = FALSE
    ) +
    ylim(-max(real_label_data$Frequences, na.rm=TRUE)*0.2, max(real_label_data$Frequences, na.rm=TRUE)*1.15) +
    theme_minimal() +
    theme(
      legend.position = c(0.9, 0.45),  # Move legend closer to the plot
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 13, face = "bold"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar() +
    geom_segment(
      data = base_data,
      aes(x = start, y = -5, xend = end, yend = -5),
      colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE
    )
  
  print(p)
  
  # Optionally save plot
  if (!is.null(directorypath)) {
    ggsave(filename = file.path(directorypath, "Go_annotation.jpeg"), plot = p, device = "jpeg", width = 20, height = 15)
    ggsave(filename = file.path(directorypath, "Go_annotation.tiff"), plot = p, device = "tiff", width = 20, height = 15)
  }
  return(p)
}
