#https://www.biostars.org/p/335187/
#the function need dataframe as input
#row gene
#column sample
c2tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


# c2rpkm <- function(counts, lengths) {
#   rate <- counts / lengths 
#   rate / sum(counts) * 1e6
# }


#remotes::install_github("plger/RNAontheBENCH")
#function 
#fpkm2tpm
