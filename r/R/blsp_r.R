blsp_test <- function(a, b, c) {
  result <- .Call(C_generate_data, a, b, c)#, PACKAGE = "blspR")
  return(result)
}
