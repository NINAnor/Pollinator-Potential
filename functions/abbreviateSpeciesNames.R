abbreviateSpeciesNames <- function(names){
  result <- sapply(names, function(name){
    x <- stringr::str_split(name, " ")[[1]]
    ret <- paste(abbreviate(x[1], 1), x[2], sep = ". ")
  })
  return(result)
}