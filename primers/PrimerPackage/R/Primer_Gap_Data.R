#' file_names
#'
#' @param file_directory 
#'
#' @return
#' @export
#'
#' @examples
file_names <- function(file_directory){
  #list.files(path ="/Users/Curly/Dropbox/SLE/SLE_Genomes/Seq_Test",pattern = "RT.*._Hits",full.names = TRUE )
  
  return(list.files(path =file_directory,pattern = ".*._[Hh]its.*",full.names = TRUE ))
}


#' RTnum_Contain
#'
#' @param file_array 
#'
#' @return
#' @export
#'
#' @examples
RTnum_Contain<-function (file_array){ 
  for( i in file_array[[1]]){
    c <- grepl("RT...",i,perl=TRUE)
    if(isTRUE(c)){
      return(i)
    }
  }
}


#' read_files_and_sum
#'
#' @param File_Names 
#' @param RT_Names 
#'
#' @return
#' @export
#'
#' @examples
read_files_and_sum <- function(File_Names,RT_Names){
 
  count = 1
  RT_Counts <- list()
  RT_sums <- list()
  for(i in File_Names ){
      
    RT_Counts[[RT_Names[count]]] <- as.data.frame(read.table(textConnection(gsub("\t",",", readLines(i))), header = TRUE,sep = ",",stringsAsFactors = FALSE))
    count = count + 1
  }

  #sum across all columns
  for ( j in RT_Names){

    RT_sums[[paste(j,"_sum",sep ="")]] <- apply(RT_Counts[[j]],1,sum)
    
  }
  return(RT_sums)
  
}

#' Find_Sequence_Gaps
#'
#' @param less_than_indexes 
#'
#' @return
#' @export
#'
#' @examples
Find_Sequence_Gaps <- function(less_than_indexes){
  total = length(less_than_indexes)
  start_end = c(0,0)
  final = c()
  count = 0
  
  for (i in 1:total){
    
    
    if(less_than_indexes[i]+1 == less_than_indexes[i+1] && count == 0){
      count = count + 1
      start_end[1] <- less_than_indexes[i]
    }
    
    else if(i == total){
      
      start_end[2] <- less_than_indexes[i]
      final = rbind(final, start_end)
      
    }  
    else if (less_than_indexes[i]+1 != less_than_indexes[i+1]){
      
      start_end[2] <- less_than_indexes[i]
      final = rbind(final, start_end)
      start_end[1] <- less_than_indexes[i+1]
      
    }
    
  }
  
  # colnames(final)<-c("Begin","End")
  # rownames(final)<-c(as.character(seq(1,nrow(final))))
  rownames(final) <-c()
  return(final)
}



#Finds regions that have low coverage.
#' Primer_Gap_Data
#'
#' @param sum_vector 
#' @param primers 
#' @param max_count 
#'
#' @return
#' @export
#'
#' @examples
Low_Coverage <- function(sum_vector,primers,max_count){
  RT_LC <- list()
  
  
  for ( i in names(sum_vector)){
    
    less_than_indexes <- which(sum_vector[[i]] < max_count)
    Gap_vector<-Find_Sequence_Gaps(less_than_indexes)
    LC_name <-paste(i,"_LC",sep ="")
    
    
    Gap_vector_diffs <- Gap_vector[,2] - Gap_vector[,1] 
    
    RT_LC[[LC_name]] <- Gap_vector[Gap_vector_diffs > 20,]
    
    primer_row_names <- rownames(primers)
    
    
    primer_gap <- c()
    count <- 1
    for(j in 1:nrow(RT_LC[[LC_name]])){
      
      
      if( as.numeric(RT_LC[[LC_name]][j,][1])<=as.numeric(primers[,1]) && count == 1 ){
        primer_gap[count] <-primer_row_names[count]
        count <- count + 1
        
      }else if (as.numeric(RT_LC[[LC_name]][j,][1])>tail(as.numeric(primers[,2]),n=1)){
        primer_gap[count] <-c("Gap out of range of primers")
        count <- count + 1
        
      }else if(as.numeric(RT_LC[[LC_name]][j,][1])<=as.numeric(primers[,2]) && as.numeric(RT_LC[[LC_name]][j,][1]) > as.numeric(primers[,2])){
        primer_gap[count] <- tail(primer_row_names,n=1)
        count <- count + 1
        
      }else{     
        
        Gap_primers<-primer_row_names[(as.numeric(primers[,1]) <= as.numeric(RT_LC[[LC_name]][j,][1]) & as.numeric(primers[,2]) >= as.numeric(RT_LC[[LC_name]][j,][1]))] #| ( as.numeric(primers[,1])<=as.numeric(RT_LC[[LC_name]][j,][1])& as.numeric(primers[,2])>=as.numeric(RT_LC[[LC_name]][j,][2]))]
        
        if(length(Gap_primers) > 1){
          primer_gap[count] <- paste(Gap_primers,collapse = "_")
          count = count + 1
        }else{
          
          
          primer_gap[count] <- Gap_primers
          
          count = count + 1
          
        }
      }     
      
    }
    Gap_Size  <- RT_LC[[LC_name]][,2] - RT_LC[[LC_name]][,1]
    RT_LC[[LC_name]] <- cbind(t(t(primer_gap)),RT_LC[[LC_name]],Gap_Size)
    colnames(RT_LC[[LC_name]]) <- c("Primer","Gap_Begins", "Gap_Ends", "Gap_Size")
    
  }
  
  return(RT_LC)
  
}


#' Primer_Gap_Data
#'
#' @param Hits_Dir 
#' @param primer_file 
#' @param max_count 
#'
#' @return
#' @export
#'
#' @examples
Primer_Gap_Data <- function(Hits_Dir,primer_file,max_count){
  
  print("cool")
  primers<-PrimerPackage::primer_database_create(primer_file = primer_file,max_count)
  RT_Files<-PrimerPackage::file_names(Hits_Dir)

  split_filePath <- as.array(strsplit(RT_Files, split = "\\/|\\/|\\_"))
  Rt_names<-apply(X = split_filePath, MARGIN = 1, RTnum_Contain)


 
  Rt_Sum <- PrimerPackage::read_files_and_sum(RT_Files,Rt_names)

  Primer_Gaps <- PrimerPackage::Low_Coverage(Rt_Sum,primers = primers,max_count = max_count)
 
  return(Primer_Gaps)
}

