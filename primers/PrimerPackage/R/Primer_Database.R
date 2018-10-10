#' primer_database_create
#'
#' @param primer_file 
#' @param colors 
#'
#' @return
#' @export
#'
#' @examples
primer_database_create <- function(primer_file,colors=c("purple","orange")){
  primers<-read.csv(primer_file,header = TRUE, na.strings = c("","NA"),stringsAsFactors = FALSE)
  
  primer_done<-primer_scheme_ends(primers)
  primer_names<-primer_scheme_names(primers)
  
  rownames(primer_done)<- primer_names
  
  pools <- c()
  primer_color <- c()
  height <- c(200,1000)
  y_height<-c()
  count = 1
  index = 1
  for (i in 1:nrow(primer_done)){
    if (primers[count, "Pool"] == 1){
      pools[index] = "Pool1"
      primer_color[index] <- colors[1]
      y_height[index] <- height[1]
    }else{
      pools[index] = "Pool2"
      primer_color[index] <- colors[2]
      y_height[index] <- height[2]
    }
    count = count + 2
    index = index + 1
    
  }
  
  primer_done <- cbind(primer_done, pools,primer_color)
  
  
  return(primer_done)
  
}


#' primer_scheme_ends
#'
#' @param primer_file 
#'
#' @return
#' @export
#'
#' @examples
primer_scheme_ends <- function(primer_file){
  
  primer_range_begin <- primer_file$End[grepl(".*LEFT",primer_file$Primer.Name,perl=TRUE)]
  primer_range_end <- primer_file$End[grepl(".*RIGHT",primer_file$Primer.Name,perl=TRUE)]
  
  return(cbind(primer_range_begin,primer_range_end))
}


#' primer_scheme_names
#'
#' @param primer_file 
#'
#' @return
#' @export
#'
#' @examples
primer_scheme_names<- function(primer_file){
  
  full_names <- primer_file$Primer.Name[grepl(".*RIGHT",primer_file$Primer.Name,perl=TRUE)]
  split_primer_names <- as.vector(strsplit(x = as.character(full_names), split = "\\_R",perl = TRUE))
  primer_names <- apply(as.array(split_primer_names),1,function(x) x[[1]][1])
  
  return(primer_names)
}
