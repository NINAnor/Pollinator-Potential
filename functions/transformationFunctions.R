mycloglog <- function(x){
  return(log(-log(1-x)))
}

myinvcloglog <- function(x){
  return(1- exp(-exp(x)))
}