getProfilesMetadata <- function(colNames){
  l1000_names <- as.vector(colNames)
  l1000_compounds <- unlist(lapply(strsplit(l1000_names, '_'), function(X){paste0(X[1:(length(X)-3)], collapse = '_')}))
  l1000_cellLine <- unlist(lapply(strsplit(l1000_names, '_'), function(X){X[(length(X)-2)]}))
  l1000_concentration <- unlist(lapply(strsplit(l1000_names, '_'), function(X){X[(length(X)-1)]}))
  l1000_time <- unlist(lapply(strsplit(l1000_names, '_'), function(X){X[length(X)]}))
  l1000_metadata <- data.frame(compounds = l1000_compounds,
                               cellLine = l1000_cellLine,
                               concentration = l1000_concentration,
                               time = l1000_time,
                               names = l1000_names)
  return(l1000_metadata)
}
