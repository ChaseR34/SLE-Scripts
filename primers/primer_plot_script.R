new_primer <- primer_database_create("~/Downloads/primers (5).csv")
old_primer <- primer_database_create("~/Dropbox/SLE/Primers/SLEV_NA_TiledAmpliconScheme.csv",colors = c("red","blue"))
order_primer<- primer_database_create("~/Dropbox/SLE/Primers/Primer_Amp_Seq_10918/Primers_For_Amp_100918.csv",colors = c("green","gold"))
head(primer_db)

mar(par)

plot_primers <- function(...){
  
  primer_db <- list(...)
  names <- strsplit(deparse(substitute(...[1:2])), split = "\\[|\\,")[[1]][1:length(primer_db)]
  plot(NULL,xlim = c(1,11000), ylim = c(0,2))

  for(j in 1:length(primer_db) ){

      for ( i in 1:nrow(primer_db[[j]])){
    
        if (primer_db[[j]][i,3] == "Pool1"){
          lines(x = c(as.numeric(primer_db[[j]][i,1]),as.numeric(primer_db[[j]][i,2])), y = rep(0.3 + ((j-1)/10),length(c(70,100))),col = primer_db[[j]][i,4],lwd = 6)
          text(x = (as.numeric(primer_db[[j]][i,1]) +as.numeric(primer_db[[j]][i,2]))/2, y = rep(0.35 + ((j-1)/10),length(c(70,100))),labels = row.names(primer_db[[j]])[i],cex = 0.5)
          #abline(v = as.numeric(primer_db[[j]][i,1]),col = primer_db[[j]][i,4])
          #abline(v = as.numeric(primer_db[[j]][i,2]),col = primer_db[[j]][i,4])
         }
        if(primer_db[[j]][i,3] == "Pool2"){
          lines(x = c(as.numeric(primer_db[[j]][i,1]),as.numeric(primer_db[[j]][i,2])), y = rep(0 + ((j-1)/10),length(c(70,100))),col = primer_db[[j]][i,4],lwd = 6)
          text(x = (as.numeric(primer_db[[j]][i,1]) +as.numeric(primer_db[[j]][i,2]))/2, y = rep(0.05 + ((j-1)/10),length(c(70,100))),labels = row.names(primer_db[[j]])[i],cex = 0.5)
          #abline(v = as.numeric(primer_db[[j]][i,1]),col = primer_db[[j]][i,4])
          #abline(v = as.numeric(primer_db[[j]][i,2]),col = primer_db[[j]][i,4])
          
          }
      }
    legend(4000 + ((j-1)*3000), 1, c("Pool1","Pool2"), col = unique(primer_db[[j]][,4]),lty = c(1,1),title = names[j],cex = 0.6)
  }
 
}
plot_primers(old_primer,new_primer,order_primer)
