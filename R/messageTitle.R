#' messageTitle
#' \code{messageTitle} A  wrapper function for sending messages to console.
#' @param aMessage text of the message.
#' @param underChar character to use for underlining the message.
#'
#' @export messageTitle
#'
#' @examples
#' \dontrun{
#' messageTitle("Hello world", "$")
#'}
#'
messageTitle <- function(aMessage, underChar="-"){
  cat (aMessage,"\n")
  cat(paste(rep(underChar, nchar(aMessage)),collapse=""),"\n")
}


