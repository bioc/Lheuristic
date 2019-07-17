#' A transcript coordinate annotation tool
#'
#' \code{getGenesLocations} Given a vector with Gene Symbols, it will produce an object with ENTREZ ID, gene name, chromosome number,
#' start position and end position for Homo sapiens.
#'
#' @param geneSymbolsSEL vector containing Gene Symbols (unique abbreviation for the gene name) to be annotated
#' @param sortByChrom logical; indicating if the results have to be ssorted ascending by chromosome number. Default to TRUE
#' @param csvFileName default is NULL; if a name is given, a csv file will be written as output
#' @keywords annotation gene
#' @import Homo.sapiens TxDb.Hsapiens.UCSC.hg19.knownGene utils ensembldb
#' @export getGenesLocations
#' @examples
#'\dontrun{
#'transcriptCoords <- getGenesLocations(rownames(selectedGenes),
#'csvFileName="results/coordsselectedGenes.csv")
#'}
#'

getGenesLocations <- function (geneSymbolsSEL, sortByChrom=TRUE, csvFileName=NULL){
  if (length(geneSymbolsSEL)>0){
    anotacs<- select(Homo.sapiens, keys=geneSymbolsSEL, columns="ENTREZID", keytype="SYMBOL")
    transcriptCoords2 <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                keys = anotacs$ENTREZID,
                                columns=c('GENEID', 'TXNAME', 'TXID', 'TXCHROM', 'TXSTART', 'TXEND' ),
                                keytype="GENEID")
    names(transcriptCoords2) <- c('ENTREZID', 'TXNAME', 'TXID', 'TXCHROM', 'TXSTART', 'TXEND')
    anotacs2<- merge(anotacs, transcriptCoords2, 'ENTREZID')
    anotacs3a <- anotacs2[!duplicated(anotacs2$ENTREZID),]
    anotacs3b <- na.omit(anotacs3a)
    if (sortByChrom){
      anotacs3 <- anotacs3b[order(anotacs3b$TXCHROM), ]
    }else{
      anotacs3 <- anotacs3b
    }
    if(!is.null(csvFileName)) write.csv(anotacs3, file=csvFileName, row.names = FALSE)
  }else{
    anotacs3=NULL
  }
  return(anotacs3)
}

