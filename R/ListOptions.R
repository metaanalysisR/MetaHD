# 

sumList <- function (list){
  res <- 0
  for (i in seq(list)){
    res <- res + list[[i]]
  }
  res
}


getList <- function (object){
  if (is.list(object))
    object
  else list(object)
}


dropList <- function (object){
  if (is.list(object) && length(object) == 1L)
    object[[1L]]
  else object
}
