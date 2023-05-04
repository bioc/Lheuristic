#' read2
#'
#' \code{read2} Wrapper to read 2 files with the same format, dimensions and row and column names.
#'
#' @param expresFName Name of first file, expected to contain expression values.
#' @param metFName Name of second file, expected to contain methylation values.
#' @param dataDirectory Name of directory where the files are staore. Defaults to ".".
#' @param sepChar Name of character used to separate matrix columns. Defaults to ";".
#' @param decChar Name of character used as decimal point. Defaults to ".".
#' 
#' @return a list with a methylation and an expression matrix
#' 
#' @importFrom utils read.table
#' @export read2
#' @examples
#'
##'
##' # Example usage of the read2 function
##' # This example assumes that the data files 'expression_data.txt' and 'methylation_data.txt' are
##' # located in the working directory.
##' 
##' #data <- read2(expresFName = "expression_data.txt", 
##' #             metFName = "methylation_data.txt")
##' 

read2 <- function (expresFName, metFName,
                   dataDirectory=".", sepChar=";", decChar= "."){
  expres<- utils::read.table(file=file.path(dataDirectory, expresFName), header=TRUE,
                      sep=sepChar,dec=decChar, row.names = 1)
  mets <-utils::read.table(file=file.path(dataDirectory, metFName), header=TRUE,
                    sep=sepChar,dec=decChar, row.names = 1)
  if (checkPairing(expres, mets))
    result <- list(expres,mets)
  else
    result <- NULL
  return(result)
}
