library(ggplot2)
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}
plot_total <- function (anchors_file, query_fai, ref_fai, query_label, ref_label, output_file){
  # read fai file. (first and second column)
  query_fai <- read.table(query_fai, header=F)
  query_length <- query_fai[, 1:2]
  colnames(query_length) <- c("queryChr", "length")
  query_length$queryChr <- as.character(query_length$queryChr)
  query_length$length <- as.integer(query_length$length)
  query_string_vector <- query_length$queryChr
  
  ref_fai <- read.table(ref_fai, header=F)
  ref_length <- ref_fai[, 1:2]
  colnames(ref_length) <- c("refChr", "length")
  ref_length$refChr <- as.character(ref_length$refChr)
  ref_length$length <- as.integer(ref_length$length)
  ref_string_vector <- ref_length$refChr
  
  # get blank df
  blank_df <- data.frame(matrix(nrow = 0, ncol = 4))
  for (i in 1:nrow(ref_length)) {
    refRowData <- ref_length[i, ]
    refChr <- refRowData[1, "refChr"]
    refChrLength <- refRowData[1, "length"]
    for (j in 1:nrow(query_length)) {
      queryRowData <- query_length[j, ]
      queryChr <- queryRowData[1, "queryChr"]
      queryChrLength <- queryRowData[1, "length"]
      new_row <- c(refChr, refChrLength, queryChr, queryChrLength)
      new_row_zero <- c(refChr, 0, queryChr, 0)
      blank_df <- rbind(blank_df, new_row)
      blank_df <- rbind(blank_df, new_row_zero)
    }
  }
  
  blank_colnames <- c("refChr", "refLength", "queryChr", "queryLength")
  colnames(blank_df) <- blank_colnames
  
  # convert column's type(factor and integer)
  blank_df$refLength <- as.integer(blank_df$refLength)
  blank_df$queryLength <- as.integer(blank_df$queryLength)
  blank_df$refChr <- factor(blank_df$refChr, levels = ref_string_vector)
  blank_df$queryChr <- factor(blank_df$queryChr, levels = query_string_vector)
  
  # read anchors file generated from AnchorWave.
  data = read.table(anchors_file, header=T)
  print(data)
  data$refChr <- as.character(data$refChr)
  data$queryChr <- as.character(data$queryChr)
  data = data[which(data$refChr %in% ref_string_vector),]
  data = data[which(data$queryChr %in% query_string_vector),]
  
  data$refChr = factor(data$refChr, levels = ref_string_vector)
  data$queryChr = factor(data$queryChr, levels = query_string_vector)
  
  plot = ggplot(data=data, aes(x=queryStart, y=referenceStart))+
    facet_grid(refChr~queryChr, scales = "free", space = "free")+
    geom_point(size=0.5, aes(color=strand)) + 
    geom_blank(data=blank_df, aes(x=queryLength, y=refLength)) +
    theme_grey(base_size = 30) +
    labs(x=query_label, y=ref_label)+scale_x_continuous(labels=changetoM, expand=c(0, 0)) + scale_y_continuous(labels=changetoM, expand=c(0, 0)) +
    theme(axis.line = element_blank(),
          panel.spacing = unit(0, "mm"),
          strip.background = element_rect(color = "white"),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          axis.text.y = element_text( colour = "black"),
          legend.position='none',
          axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )
  png(output_file , width=2000, height=1500)
  print(plot)
  dev.off()
}