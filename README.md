# Lheuristic
## Introduction

This vignette gives information on how to use the 'Lheuristic' package functions to select genes with an expression vs methylation scatterplot that have an L-shaped pattern. 

### Input data
The data provided is a subset of the methylation array and the expression array from the colon adenocarcinoma (COAD) from The Cancer Genome Atlas (TCGA). The full dataset is available from the TCGA website (https://tcga-data.nci.nih.gov/docs/publications/tcga/?).

The methylation and expression data are in 2 different files,  the columns are the samples (patients) and the rows are the genes. We have to check that the matrices  have matching rows and columns, and in the same order.

```{r}
data("TCGAexp")
data("TCGAmet")

checkPairing(TCGAexp, TCGAmet)
```

## Data analysis

There are two methods implemented in the package set for the selection of L-shaped genes: correlation and heuristic. Each method has various parameters which are optimized to select for L-shaped scatterplot representations from methylation and expression data, and with a negative correlation.
The various functions allow for parameter tunning and flexibility of scatterplot distributions; even though the methods have been optimized to select for L-shaped scatterplot distributions, the selection of any  distribution of choice is possible. 

### Correlation method

The correlation method is based on the selection of genes by significant negative correlation coeficient.

```{r}
correlation <- correlationSelection (TCGAexp, TCGAmet, pValCutoff=0.25,  rCutoff=-0.5, type="Spearman",adj=TRUE)
correlationL <-correlation[correlation$SigNegCorr,] 
correlationNoL <-correlation[!correlation$SigNegCorr,] 

cat("The number of genes selected with the correlation method is: ", sum(correlationL$SigNegCorr),"\n")
```

We can observe the distribution of the correlation coeficients:

```{r}
hist(correlationL[,1], xlim=c(-1,0), main="Significant correlations in TCGA dataset")
```

Depending on the *pvalue* and the *r* cutoff selected, a different number of genes will be retained.
The resulting output is a list with the genes lassified as selected and not selected, in table format with the following columns: *r coefficient* and the *p-value* for the chosen correlation, the *adjusted p value*, the *distance correlation* and a logical wheter the gene was classified as L-shaped or not according to the set parameters.

```{r}
kable(correlationL, caption = "Genes selected with L-shape with the correlation method")
```

The resulting genes can also be plotted and saved in a PDF file.

```{r}
plotGenesMat (mets=TCGAmet[rownames(correlationL),], 
              expres=TCGAexp[rownames(correlationL),], 
              fileName ="correlationLgenes.pdf",
              text4Title = correlationL[rownames(correlationL),""]) 
```



### Heuristic method

The heuristic method intends to select L-shaped scatterplots by superimposing a grid on the graph and defining cells which have to (or do not have to) contain a minimum (or maximum) percentage of points if the scatterplot is to be called L-shaped.

The method also computes a score in such a way that scores in selected regions (L region) score positively and points off the L region score negatively. An appropriate setting of scores and weights should yield positive scores for L-shaped scatterplots and negative scores for those that are not L-shaped. 

```{r}
sampleSize <- dim(TCGAmet)[2]
numGenes <-   dim(TCGAmet)[1]

 
reqPercentages <- matrix (c(2, 20, 5, 1, 40, 20, 0, 1, 2), nrow=3, byrow=TRUE)
(maxminCounts <- toReqMat(sampleSize, reqPercentages))

(theWeightMifL=matrix (c(2,-2,-sampleSize/5,1,0,-2,1,1,2), nrow=3, byrow=TRUE))
(theWeightMifNonL=matrix (c(0,-2,-sampleSize/5,0,0,-2,0,0,0), nrow=3, byrow=TRUE))
   
   

heur <- scoreGenesMat (mets=TCGAmet,
							              expres=TCGAexp,
                            aReqPercentsMat=reqPercentages,
                            aWeightMifL=theWeightMifL,
                            aWeightMifNonL=theWeightMifNonL )
  cat("Number of scatterplots scored  : ", dim(heur)[1],"\n")
  cat("Number of L-shape scatterplots : ", sum(heur[,1]),"\n")
  
heurL <- heur[heur$logicSc,]
heurNoL <- heur[!heur$logicSc,]
```

Results can be checked in the following table, were there is a logical value describing if the gene has or not an L-shape based on our criteria and the *numerSc* score:

```{r}
kable((heurL), caption = "Genes selected with L-shape with the heuristic method")
```

Next, the scatterplots of the selected genes can be visualized and saved  on a PDF file.

```{r}
plotGenesMat (mets=TCGAmet[rownames(heurL),], 
              expres=TCGAexp[rownames(heurL),], 
              fileName ="heurLgenes.pdf",
              text4Title = heurL[rownames(heurL),"numeriSc"]) 
```

## Comparison of selected L-shaped genes

We have for lists of genes that we have identified with the 4 different methods. We can create the intersection of all lists and then visualize the results with a Venn Diagram.

```{r}
myVenn<- venn.diagram(x=list(correlation=rownames(correlationL), 
                                heuristic = rownames(heurL),
                                filename=NULL, lty = "blank",  
                                fill=c("pink1", "skyblue"),
                       main="Genes in common between the 4 methods")
grid.newpage()
grid.draw(myVenn)
```

We can decide to choose the genes that have been selected by 2 or more selection criteria or methods, to have higher consistancy in the selection. 

```{r}
inCommonL <- intersect(rownames(correlationL), rownames(heurL))

commonL <- unique(inCommonL)
allL <- inCommonL[duplicated(inCommonL)]
```

We can also plot selected genes.

```{r}
par(mfrow=c(2,2))
myGene1 <-allL[1]
xVec<- as.numeric(TCGAmet[myGene1,])
yVec<-as.numeric(TCGAexp[myGene1,])
titleT <- paste (myGene1, "(May be GRM)")
plotGeneSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)


myGene2 <-allL[2]
xVec<- as.numeric(TCGAmet[myGene2,])
yVec<-as.numeric(TCGAexp[myGene2,])
titleT <- paste (myGene2, "(May be GRM)")
plotGeneSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)


myGene3 <-allL[3]
xVec<- as.numeric(TCGAmet[myGene3,])
yVec<-as.numeric(TCGAexp[myGene3,])
titleT <- paste (myGene3, "(May be GRM)")
plotGeneSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)


myGene5 <-allL[5]
xVec<- as.numeric(TCGAmet[myGene5,])
yVec<-as.numeric(TCGAexp[myGene5,])
titleT <- paste (myGene5, "(May be GRM)")
plotGeneSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)

```

### Annotate selected genes on the chromosome

To annotate genes on the corresponding chromosomes, we will get the transcript coordinates for each gene. We only need to annotate the genes once, since we can store the information in a .csv file.

```{r}
recalc<- TRUE
if (recalc) {
tccorrelation <- getGenesLocations(unique(rownames(correlationL)),
                              csvFileName="coordsLcorrelation.csv")
tcHeur <- getGenesLocations(unique(rownames(heurL)),
                              csvFileName="coordsLheur.csv")

  transcriptCoordsList <- 
     list(tccorrelation = tccorrelation, tcHeuristic=tcHeur)
  save(transcriptCoordsList, file="transcriptCoordsLgenes.Rda")
}else{
  load(file="transcriptCoordsLgenes.Rda")
}
```

For example the table of annotations corresponding to the \texttt{selectedcorrelationDA} gene list has the following aspect:

```{r}

kable(head(transcriptCoordsList[["tccorrelation"]]))

```

