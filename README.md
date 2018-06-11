# lpattern
# lpattern
## Introduction

This vignette gives information on how to use the 'lpattern' package functions to select genes with an expression vs methylation scatterplot that have an L-shaped pattern. 

### Input data
The data provided is a subset of the methylation array and the expression array from the colon adenocarcinoma (COAD) from The Cancer Genome Atlas (TCGA). The full dataset is available from the TCGA website (https://tcga-data.nci.nih.gov/docs/publications/tcga/?).

The methylation and expression data are in 2 different files,  the columns are the samples (patients) and the rows are the genes. We have to check that the matrices  have matching rows and columns, and in the same order.

```{r}
data("TCGAexp")
data("TCGAmet")

checkPairing(TCGAexp, TCGAmet)
```

## Data analysis

There are four methods implemented in the package set for the selection of L-shaped genes: naive, CMI, heuristic and scagnostics. Each methods has various parameters which are optimized to select for the L-shape in scatterplot data, and with a negative correlation.
The functions allow for parameter tunning, in which case, any other selection is possible.

### Naive method

The naive method is based on the selection of genes by significant negative correlation.

```{r}
naive <- naiveSelection (TCGAexp, TCGAmet, pValCutoff=0.25,  rCutoff=-0.5, type="Spearman",adj=TRUE)
naiveL <-naive[naive$SigNegCorr,] 
naiveNoL <-naive[!naive$SigNegCorr,] 

cat("The number of genes selected with the naive method is: ", sum(naiveL$SigNegCorr),"\n")
```

We can observe the distribution of the correlation coeficients:
```{r}
hist(naiveL[,1], xlim=c(-1,0), main="Significant correlations in TCGA dataset")
```

Depending on the *pvalue* and the *r* cutoff more or less genes will be retained.
The result is a list with the genes slected and not selected in a table format with columns *r coefficient* and the *p-value* for the chosen correlation, the *adjusted p value*, the *distance correlation* and a logical wheter the gene was classified as L-shaped or not according to the parameters set.

```{r}
kable(naiveL, caption = "Genes selected with L-shape with the naive method")
```

The resulting genes can also be plotted and saved in a PDF file.

```{r}
plotGenesMat (mets=TCGAmet[rownames(naiveL),], 
              expres=TCGAexp[rownames(naiveL),], 
              fileName ="naiveLgenes.pdf",
              text4Title = naiveL[rownames(naiveL),""]) 
```

### CMI method

The Conditional Mutual Information method was based on the expression and methylation values computed at different points between 0 and 1 reached a minimum. This minimum should be small enough according predefined thresholds. This minimum was considered to be the cutoff point for methylation. 

The cMI function computes cMI values for different t values, from 0 to 1 and a step of 0.01. The output is stored in a data frame. For each gene, the optimal threshold is the t-value that results the minimum CMI. We can run the CMI function with the default parameters, were $h=0.2$, $smallR=0.25$ and $minCMI=0.1$.

```{r}
cm <- cmiSelection (methData = TCGAmet, exprData = TCGAexp )
cmiL <- cm[cm[,"meth_regulated"],]
cmiNotL <- cm[!cm[,"meth_regulated"],]
cat("The number of genes selected with the CMI method is: ", sum(cmiL$meth_regulated),"\n")
```

The results are presented in a table format as follows, were there is the *minimum CMI*, the optimal *t*, and a logical describing if the gene has or not an L-shape based on our criteria:

```{r}
kable(cmiL, caption = "Genes selected with L-shape with the CMI method")
```

Next, we can visualize the scatterplots of the selected genes with the \text{plotGenesMat} function and save them on a PDF file.
```{r}
plotGenesMat (mets=TCGAmet[rownames(cmiL),], 
              expres=TCGAexp[rownames(cmiL),], 
              fileName ="cmiLgenes.pdf",
              text4Title = cmiL[rownames(cmiL),""]) 
```

### Heuristic method

The heuristic method intends to select L-shaped scatterplots by superimposing a grid on the graph and defining cells which have to (or don't have to) contain a minimum (or maximum) percentage of points if the scatterplot is to be called L-shaped.

The method also computes a score in such a way that scores in selected regions (L region) score positively and points off the region of interest score negatively. An appropriate setting of scores and weights should yield positive scores for L-shaped scatterplots and negative scores for those that are not. 

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

We can check the results in the following table, were there is a logical value describing if the gene has or not an L-shape based on our criteria and the *numerSc* score:

```{r}
kable((heurL), caption = "Genes selected with L-shape with the heuristic method")
```

Next, we can visualize the scatterplots of the selected genes and save them on a PDF file.
```{r}
plotGenesMat (mets=TCGAmet[rownames(heurL),], 
              expres=TCGAexp[rownames(heurL),], 
              fileName ="heurLgenes.pdf",
              text4Title = heurL[rownames(heurL),"numeriSc"]) 
```


### Scagnostics method

The scagnostics provides seven coeficients for each gene that describe the scatterplot. These are: outlying, skewed, clumpy, sparse, striated, convex, skinny, stringy and monotonic. After a fine-tunning exercise, the parameters that we select from the DA dataset are: Monotonic, Convex, Striated and Clumpy. The parameters that we select from the TCGA dataset are: Monotonic, Convex, Skinnyand Clumpy. For the GEO dataset, none the parameters seems to be able to discriminate between TRUE and FALSE. That is also why we did a visual selection of L-shaped genes, however the results did not improve.

```{r}
###need to change that for the functions

si <- c()
for (i in 1:nrow(TCGAexp)){
  si <- cbind(si, scagnostics(TCGAmet[i,], TCGAexp[i,] ))

}

colnames(si) <- rownames(TCGAexp)


si_t <- as.data.frame(t(si))
dim(si_t)

scagL <- si_t[which( (si_t$Outlying < 0.06 ) & (si_t$Convex <  0.09) & (si_t$Sparse > 0.20) & (si_t$Clumpy > 0.17)), ]

plotGenesMat (mets=TCGAmet[rownames(scagL),], 
              expres=TCGAexp[rownames(scagL),], 
              fileName ="scagLgenes.pdf",
              text4Title = rownames(scagL)) 
```

## Comparison of selected L-shaped genes

We have for lists of genes that we have identified with the 4 different methods. We can create the intersection of all lists and then visualize the results with a Venn Diagram.

```{r}
myVenn<- venn.diagram(x=list(naive=rownames(naiveL), 
                                CMI = rownames(cmiL), 
                                heuristic = rownames(heurL),
                                scagnostics = rownames(scagL)), 
                                filename=NULL, lty = "blank",  
                                fill=c("pink1", "skyblue", "mediumorchid", "lightgoldenrod"),
                       main="Genes in common between the 4 methods")
grid.newpage()
grid.draw(myVenn)
```

We can decide to choose the genes that have been selected by 2 or more methods, for example, to have higher consistancy in the selection. 

```{r}
inCommonL1 <- intersect(rownames(heurL), rownames(cmiL))
inCommonL2 <- intersect(rownames(naiveL), rownames(heurL))
inCommonL3 <- intersect(rownames(heurL), rownames(scagL))

commonAll <- c(inCommonL1, inCommonL2, inCommonL3)
commonL <- unique(commonAll)
allL <- commonAll[duplicated(commonAll)]
```

We can also plot selected genes.
```{r}
par(mfrow=c(2,2))
myGene1 <-allL[1]
xVec<- as.numeric(TCGAmet[myGene1,])
yVec<-as.numeric(TCGAexp[myGene1,])
titleT <- paste (myGene1, "(May be GRM)")
plotGenSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)


myGene2 <-allL[2]
xVec<- as.numeric(TCGAmet[myGene2,])
yVec<-as.numeric(TCGAexp[myGene2,])
titleT <- paste (myGene2, "(May be GRM)")
plotGenSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)


myGene3 <-allL[3]
xVec<- as.numeric(TCGAmet[myGene3,])
yVec<-as.numeric(TCGAexp[myGene3,])
titleT <- paste (myGene3, "(May be GRM)")
plotGenSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)


myGene5 <-allL[5]
xVec<- as.numeric(TCGAmet[myGene5,])
yVec<-as.numeric(TCGAexp[myGene5,])
titleT <- paste (myGene5, "(May be GRM)")
plotGenSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)

```

### Annotate selected genes on the chromosome

To annotate genes on the corresponding chromosomes, we will get the transcript coordinates for each gene. We only need to annotate the genes once, since we can store the information in a .csv file.

```{r}
recalc<- TRUE
if (recalc) {
tcNaive <- getGenesLocations(unique(rownames(naiveL)),
                              csvFileName="coordsLnaive.csv")
tcCMI <- getGenesLocations(unique(rownames(cmiL)),
                              csvFileName="coordsLcmi.csv")
tcHeur <- getGenesLocations(unique(rownames(heurL)),
                              csvFileName="coordsLheur.csv")
tcScag <- getGenesLocations(rownames(rownames(scagL)),
                              csvFileName="coordsLscag.csv")

  transcriptCoordsList <- 
     list(tcNaive = tcNaive, tcCMI=tcCMI, 
          tcHeuristic=tcHeur,tcScagnostics=tcScag)
  save(transcriptCoordsList, file="transcriptCoordsLgenes.Rda")
}else{
  load(file="transcriptCoordsLgenes.Rda")
}
```

For example the table of annotations corresponding to the \texttt{selectedNaiveDA} gene list has the following aspect:

```{r}

kable(head(transcriptCoordsList[["tcNaive"]]))

```

