
#' ploGenesInChromosomes: A function to plot specific genes onto the corresponding chromosomes
#'
#' \code{ploGenesInChromosomes} Given a list of genes with their transcript coordenates, the function will plot the genes
#' on their specific locations on each corresponding chromosme. It can also collocate the genes with CpG islands information and
#' and DNAseI hypersesitive sites.
#'
#' @param transcriptCoords object of class containing a list of genes annotated with the getGenesLocations
#' @param plotsFilename object with the name of the file that will be used to save the .pdf with the graphs.
#' @param minbase number for the smallest basepair position of the chromosome. First position on the chromosome from the 5' side.
#' @param maxbase number for the largest basepair position of the chromosome. Last position on the chomosome from the 5' side.
#' @param islandData object of class GRanges containing CpG islands position data.
#' @param dnaseData object of class GRanges containing DNAseI hypersensitive sites position data.
#' @keywords genes chromosomes plot
#' @import Gviz gtools GenomicRanges utils grDevices
#' @export plotGenesInChroms
#' @details The data for the cpG islands was obtained
#' from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz and the
#' data for the DNAseI hypersensitive sites from
#' http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/
#' For both datasets, positions were annotated with GenomicRanges::GRanges().
#' @examples
#'\dontrun{
#'transcriptCoords <- getGenesLocations(rownames(selectedGenes),
#'csvFileName="results/coordsselectedGenes.csv")
#'starts <- min(lapply(transcriptCoords, function(tc){as.numeric(min(tc$TXSTART))}))
#'ends <- max(lapply(transcriptCoords, function(tc){as.numeric(max(tc$TXEND))}))
#'
#'minbase <- starts - (0.25*(ends-starts)) #allow 25% extra space on the plot
#'maxbase <- ends + (0.25*(ends-starts))

#'islandHMM = read.csv(paste("dades", "model-based-cpg-islands-hg19.txt",sep="/"),
#'stringsAsFactors=FALSE, header=TRUE)
#'islandData <- GRanges(seqnames=Rle(islandHMM$chr),
#'                    ranges=IRanges(start=islandHMM$start, end=islandHMM$end),
#'                    strand=Rle(strand(rep("*",nrow(islandHMM)))))
#'
#'dnase <- read.csv(paste("dades","wgEncodeRegDnaseClusteredV3.bed",sep="/"),
#'stringsAsFactors=FALSE,header=FALSE)
#'dnaseData <- GRanges(seqnames=dnase[,1],
#'                   ranges=IRanges(start=dnase[,2], end=dnase[,3]),
#'                   strand=Rle(rep("*",nrow(dnase))),
#'                   data=dnase[,5])
#'
#'plotGenesInChromosomes(transcCoords, fileName, minbase, maxbase, islanData, dnaseData)
#'}
#'

plotGenesInChroms <- function (transcriptCoords, plotsFilename, minbase, maxbase, islandData, dnaseData){
  if(!is.null(transcriptCoords)){
    anotacs4<-transcriptCoords[,c('TXCHROM', 'TXSTART', 'TXEND')]
    anotacs4<-anotacs4[complete.cases(anotacs4),]

    #change column names to read it into GenomicRanges
    colnames(anotacs4)<-c("chromosome","start","end")
    genRangList<-makeGRangesFromDataFrame(anotacs4) #if we keep geneid we add keep.extra.columns=TRUE inside function

    axisT<-GenomeAxisTrack()

    data <- read.table(paste("dades","cytoBandIdeo.txt", sep="/"), header=F, sep="\t")
    colnames(data) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')


    pdf(file=plotsFilename, width= 8, height = 12)
    #draw chromosomes (it has to be done 1 by 1)
    for (i in 1:length(unique(anotacs4$chromosome))){
      chr<-mixedsort(unique(anotacs4$chromosome))[i]
      #draw axis list from our genomic ranges object containing the gene positions
      #read each gene per chromosome
      genList<-AnnotationTrack(anotacs4, name = "Genes", genome ="hg19", chromosome = chr,  stacking ="dense", col= "#5E2366", fill= "#5E2366")
      ideoT<-IdeogramTrack(chromosome = chr,genome="hg19", bands=data,  showId=FALSE)
      # CpG island track
      islandData2 <- islandData[seqnames(islandData) == chr &  (start(islandData) >= minbase & end(islandData) <= maxbase)]
      islandTrack <- AnnotationTrack(range=islandData2, genome="hg19", name="CpG Islands",
                                     chromosome=chr)
      # DNaseI hypersensitive site data track
      dnaseTrack <- DataTrack(range=dnaseData, genome="hg19", name="DNAseI",
                              type="gradient", chromosome=chr)
      plotTracks(list(ideoT,axisT,genList, dnaseTrack, islandTrack), sizes=c(1,2,2,1,20), main=chr, cex.main=1,littleTicks = TRUE, showTitle=TRUE)
    }
    dev.off()
  }
}

