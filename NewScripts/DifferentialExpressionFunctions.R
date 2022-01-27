
#____________Download/Load packages________#

CRAN_packages<- c("pheatmap", "RColorBrewer", "ggrepel", "gplots", "reshape2", "ggplot2", "factoextra", "cluster", "NbClust", "tidyverse", "reticulate", "dplyr", "tidyr")
new.packages <- CRAN_packages[!(CRAN_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages)
}



biocManagerPackages = c("DESeq2", "vsn", "DEGreport", "pathview")

new.packages <- biocManagerPackages[!(biocManagerPackages %in% installed.packages()[,"Package"])]

downloadBiocManager = !requireNamespace("BiocManager", quietly = TRUE)
if(downloadBiocManager){
  install.packages("BiocManager")
}

for (package in new.packages){
  BiocManager::install(package)
}

for(package in c(CRAN_packages, biocManagerPackages)){
  library(package, character.only = TRUE)
}



#__________________________________________#

colorVector <- c("#ffff00", "#ff0000", "#00ffff", "#0000ff", "#ff00ff", "#ff8000", "#00ff00",
                    "#ffc0cb", "#cbffc0", "#ff5e7a", "#ff8c00", "#c0ebff", "#764100")


checkColumnsMatch <- function(coldata, cts){
  if(!all(rownames(coldata) %in% colnames(cts)))
  {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

setColumnsInOrder <- function(coldata, cts){
  if(all(rownames(coldata) != colnames(cts))){
    cts <- cts[, rownames(coldata)]
  }
}


# input:
# dds -> deseq2 object
# gene -> string of gene name
# condition1 -> the x axis that counts will be plotted across
# condition2 -> the grouping the data will be grouped in
# xlabel -> the x-axis label default to condition1
# grouplabel -> the label for the group legend, default conditoin2
#
# ouput:
# ggplot object
countPlotOver2Variables <- function(dds, gene, condition1, condition2, xlabel = condition1, grouplabel = condition2){

  intgroup = c(condition1, condition2)
  LNvRMtime <- DESeq2::plotCounts(dds, gene,
                          intgroup = c(condition1, condition2), returnData = TRUE)

  return(ggplot(LNvRMtime,
                aes(x = LNvRMtime[,condition1], y = LNvRMtime[,"count"], color = LNvRMtime[,condition2], group = LNvRMtime[,condition2])) +
           geom_point() + stat_summary(fun.y=mean, geom="line") +
           scale_y_log10() + ggtitle(gene) +
          labs(x = xlabel, y = "count", color = grouplabel))

}

# input:
# pcaData -> data.frame with columns PC1 and PC2 minimum.
# color -> categorical column for coloring
# shape -> categorical column for shape
# scale -> makes figure larger by a factor
#
# output:
# returns ggplot object
pcaPlotFormated <- function(pcaData, color = "", shape = "", legendTitleC = color, legendTitleS = shape, scale=1){

  percentVar <- round(100 * attr(pcaData, "percentVar"))

  if(color != "" & shape != ""){

    return(ggplot(pcaData, aes(x=PC1, y=PC2, color=pcaData[,color], shape=pcaData[,shape])) +
             geom_point(size=3*scale) +
             xlab(paste0("PC1: ",percentVar[1],"% variance")) +
             ylab(paste0("PC2: ",percentVar[2],"% variance")) +
             labs(color = legendTitleC, shape = legendTitleS) +
             coord_fixed() +
             theme(axis.title=element_text(size=16*scale),
                   legend.text=element_text(size=11*scale),
                   legend.key.height = unit(1*scale, units = "cm"),
                   legend.key.width = unit(1*scale, units = "cm"),
                   legend.title = element_text(size=9*scale)))

  } else if(color == "" & shape == "") {

    return(ggplot(pcaData, aes(x=PC1, y=PC2)) +
                    geom_point(size=3*scale) +
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed() +
                    theme(axis.title=element_text(size=16*scale),
                          legend.text=element_text(size=11*scale),
                          legend.key.height = unit(1*scale, units = "cm"),
                          legend.key.width = unit(1*scale, units = "cm"),
                          legend.title = element_text(size=9*scale)))
  } else if(color == ""){

    return(ggplot(pcaData, aes(x=PC1, y=PC2, shape=pcaData[,shape])) +
             geom_point(size=3*scale) +
             xlab(paste0("PC1: ",percentVar[1],"% variance")) +
             ylab(paste0("PC2: ",percentVar[2],"% variance")) +
             labs(shape = legendTitleS) +
             coord_fixed() +
             theme(axis.title=element_text(size=16*scale),
                   legend.text=element_text(size=11*scale),
                   legend.key.height = unit(1*scale, units = "cm"),
                   legend.key.width = unit(1*scale, units = "cm"),
                   legend.title = element_text(size=9*scale)))
  } else {

    return(ggplot(pcaData, aes(x=PC1, y=PC2, color=pcaData[,color])) +
                    geom_point(size=3*scale) +
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    labs(color = legendTitleC) +
                    coord_fixed() +
                    theme(axis.title=element_text(size=16*scale),
                          legend.text=element_text(size=11*scale),
                          legend.key.height = unit(1*scale, units = "cm"),
                          legend.key.width = unit(1*scale, units = "cm"),
                          legend.title = element_text(size=9*scale)))
  }
}


# input:
# rld -> DESeq Transform object
# topGeneNumber -> subsets genes for "n" most variable.
#
# output:
# dataframe of positons on each principal component.
generateGenePCAValues <- function(rld, topGeneNumber = 500, groups=NULL){

    assay_table <- as.data.frame(assay(rld))
    if (!is.null(groups)){
      assay_table <- assay_table[, (names(assay_table) %in% groups)]
    }

    assay_table <- as.matrix(assay_table)
    rv <- rowVars(assay_table) #Get variation
    select <- order(rv, decreasing = TRUE)[seq_len(min(topGeneNumber, length(rv)))]

    test.PCA <- prcomp(t(assay_table[select, ]))

    return(test.PCA$rotation)

}



# input:
# rld -> DESeq TRansform object
# coldata -> annotation data, used for row an column labels of sampleDistMatrix
# column -> the column from anotation to label the rows and columns of the sampleDistMatrix
# scale -> default 2, just scales the figure.
#
# output:
# returns pheatmap object
# heatmap is a distance matrix with a dendrogram showing how similar each sample is.
generateDistMatrix <- function(rld, coldata, column, scale=2, groups=NULL) {

    assay_table <- as.data.frame(assay(rld))
    if(!is.null(groups)){
      assay_table <- assay_table[, (names(assay_table) %in% groups)]
      coldata[,"holder"] <- 1
      coldata <- coldata[(rownames(coldata) %in% groups), ]
      coldata[, "holder"] <- NULL
    }
    assay_table <- as.matrix(assay_table)

    sampleDists <- dist(t(assay_table))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(coldata[,column])
    colnames(sampleDistMatrix) <- paste(coldata[,column])

    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

    (pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colors) +
                 theme(axis.title=element_text(size=20),
                       legend.text=element_text(size=15),
                       legend.key.height = unit(1, units = "cm"),
                       legend.key.width = unit(1, units = "cm"),
                       legend.title = element_text(size=13)))

}


## this is for bar plot generation. Cleans up column names
detectGroups <- function (x){  # x are col names
  tem <- gsub("Counts$", "", x)
  tem <- gsub("counts$", "", tem)
  # tem <-gsub(".", "", tem)
  tem <- gsub("[0-9]*$","",tem) # Remove all numbers from end
  #tem = gsub("_Rep|_rep|_REP","",tem)
  tem <- gsub("_$","",tem); # remove "_" from end
  tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
  tem <- gsub("_rep$","",tem); # remove "_rep" from end
  tem <- gsub("_REP$","",tem)  # remove "_REP" from end
  return( tem )
}

# input:
# cts -> count matrix
#
# output:
# ggplot object, barplot for each condition type total amount of reads.
generateBarplotForTotalReadCounts <- function(cts, groups=NULL){

    if(!is.null(groups)){
      cts <- cts[, (names(cts) %in% groups)]

    }

    # print(colnames(cts))
    # print(detectGroups(colnames(cts)))
    conditions <- as.factor(detectGroups(colnames(cts)))
    col = rainbow(nlevels(conditions))[ conditions ]
    par(mar=c(9,5,3,3))
    return(barplot( colSums(cts)/1e6,
         col=col,las=3, main="Total read counts (millions)", ylab="CPM", ylim=c(0,30)))

}

# input:
# cts -> count matrix
#
# output:
# ggplot object, boxplot for each conditions normalized counts.
generateNormalizedBoxplot <- function(cts, groups=NULL){

    if(!is.null(groups)){
      cts <- cts[, (names(cts) %in% groups)]

    }
    pseudoCount = rlog(data.matrix(cts))
    pseudoCount = as.data.frame(pseudoCount)
    df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix

    tem <- gsub("Counts$", "", df$Samples)
    tem <- gsub("counts$", "", tem)
    tem <- gsub("[0-9]*$","",tem) # Remove all numbers from end
    tem <- gsub("_$","",tem); # remove "_" from end
    tem <- gsub(".$","",tem)


    df = data.frame(df, Condition = tem)

    return(ggplot(df, aes( y = count, fill = Condition)) + geom_boxplot() + xlab("") +
        ylab(expression(rlog[](count))))

}


# input:
# cts -> count count matrix
#
# output:
# ggplot object, density plot of the normalized counts
generateNormalizedDensityPlot <- function(cts, groups=NULL){

    if(!is.null(groups)){
      cts <- cts[, (names(cts) %in% groups)]

    }

    pseudoCount = rlog(data.matrix(cts))
    pseudoCount = as.data.frame(pseudoCount)
    df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix
    # tem <- gsub("Counts$", "", df$Samples)
    # tem <- gsub("counts$", "", tem)
    # tem <- gsub("[0-9]*$","",tem) # Remove all numbers from end
    # tem <- gsub("_$","",tem); # remove "_" from end
    # tem <- gsub(".$","",tem)
    df = data.frame(df, Condition = "")

    return(ggplot(df, aes(x = count, colour = Samples, fill = Samples)) + ylim(c(0, .4)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(rlog[](count))))
}


# input:
# cts -> count matrix
#
# output:
# dataframe the the variance of each gene in each condition
generateGeneCountVarianceVaraince <- function(cts){

    cts_Matrix <- data.matrix(cts)
    logcounts = rlog(cts_Matrix)
    rownames(logcounts) <- rownames(cts_Matrix)
    var_genes <- apply(logcounts, 1, var)

    return(var_genes)

}


# input:
# cts -> count matrix
#
# output:
# dataframe with the rlog counts
generaterlogNormalizedCountDataFrame <- function(cts, groups=NULL, method="rlog"){

    if(!is.null(groups)){
      cts <- cts[, (names(cts) %in% groups)]

    }

    cts_Matrix <- data.matrix(cts)
    if(method == "vst"){
      rlogcounts = vst(cts_Matrix)
    } else {
      rlogcounts = rlog(cts_Matrix)

    }
    rownames(rlogcounts) <- rownames(cts_Matrix)

    return(rlogcounts)
}


# input:
# rlogcounts -> data.frame reutned from generaterlogNormalizedCountDataFrame
# coldata -> annotation data.frame
# column -> column that samples will be labeled with.
# geneNumber -> top n most variable genes.
#
# output:
# heatmap object for variable genes.
generateGeneCountVarianceHeatmap <- function(rlogcounts, coldata, column, geneNumber, margins=c(10,6), labRow = NULL, labCol = NULL, cexRow=1){
    # cts_Matrix <- data.matrix(cts)
    # logcounts = rlog(cts_Matrix)
    # rownames(logcounts) <- rownames(cts_Matrix)

    # print(coldata[,column])
    numberOfGroups <- length(coldata[,column][!duplicated(coldata[,column])])
    # print(numberOfGroups)

    var_genes <- apply(rlogcounts, 1, var)

    select_var <- names(sort(var_genes, decreasing=TRUE))[1:geneNumber]

    highly_variable_lcpm <- rlogcounts[select_var,]


    mypalette <- brewer.pal(11,"RdYlBu")
    morecols <- colorRampPalette(mypalette)
    # print(morecols(50))



    colorGrouptoColor <- vector(mode="list", length=numberOfGroups)
    groups <- c()
    for (group in coldata[,column]){
      if (!(group %in% groups)){
        groups <- append(groups, group)
      }
    }
    
    names(colorGrouptoColor) <- groups
    
    
    for(i in 1:numberOfGroups){
      # colors <- append(colors, colorVector[i])
      colorGrouptoColor[[i]] <- colorVector[i]
    }
    
    
    colors <- c()
    for (group in coldata[,column]){
      colors <- append(colors, colorGrouptoColor[[group]])
    }
    
    # print(colors[1:numberOfGroups])
    # col.cell <- colors[coldata[,column]]
    col.cell <- colors
    # print(class(coldata[,column]))
    # print(coldata[,column])
    # print(colors)
    # print(col.cell)
    # print(rownames(highly_variable_lcpm))
    return(gplots::heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", ColSideColors=col.cell,scale="row", margins=margins, labRow = labRow, labCol = labCol, cexRow=cexRow))
    # return(gplots::heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",scale="row"))

}

getPerThreshold <- function(res, per, variable){
  return(quantile(res[complete.cases(res),variable], c(per)))


}

#Dont use
getPercentileIndex <- function(res, per){

    return(ceiling(nrow(res[complete.cases(res),])*per))

}

#dont use
getPercentileValue <- function(index, df, column){
    holder <- df[order(df[,column]),]
    holder <- holder[complete.cases(holder),]

    return(holder[index,column])
}


differentialExpression <- function(control, treated, variable, dds, alpha=0.01){
    resp.01 <- results(dds, contrast=c(variable,treated,control), alpha = alpha)
    resp.01$score <- abs(resp.01$log2FoldChange * -log(resp.01$padj))

    return(resp.01)

}

combinedResultDataFrames <- function(formatedCSVThereshold0, res, index, control, treated){
    # this code adds the LFC and padj for each tested condition to one master file
    names(formatedCSVThereshold0)[index] <- paste("LFC_", treated, "v", control, sep="")
    names(formatedCSVThereshold0)[(index+1)] <- paste("padj_", treated, "v", control, sep="")
    names(formatedCSVThereshold0)[(index+2)] <- paste("Score_", treated, "v", control, sep="")
    formatedCSVThereshold0[,index] <- res[,"log2FoldChange"]
    formatedCSVThereshold0[,(index+1)] <- res[,"padj"]
    formatedCSVThereshold0[,(index+2)] <- res[,"score"]

    return(formatedCSVThereshold0)
}


generateDataframeForVolcanoPlot <- function(res, threshold, label, topNumberOfGenes){


    testingRES_Vol <- res
    threshold_OE <- testingRES_Vol$score <= threshold
    testingRES_Vol$threshold <- threshold_OE

    # making dataframe
    # here we store the 3 columns we care about as vectors.
    log2FoldChange <- testingRES_Vol$log2FoldChange
    padj <- testingRES_Vol$padj
    threshold <- testingRES_Vol$threshold

    # we then make a dataframe from those 3 columns.
    testDataframe.data <- data.frame(log2FoldChange,padj,threshold)
    # give the dataframe the gene rownames.
    row.names(testDataframe.data) <- rownames(testingRES_Vol)

    ## Sort dataframe by ordered Score
    if(label == "score"){
        res_tableOE_ordered <- testDataframe.data[rev(order(abs(testDataframe.data$log2FoldChange * -log(testDataframe.data$padj)))), ]
    } else if(label == "LFC"){
        res_tableOE_ordered <- testDataframe.data[rev(order(abs(testDataframe.data$log2FoldChange))), ]
    } else if(label == "padj"){
        res_tableOE_ordered <- testDataframe.data[rev(order(abs(-log(testDataframe.data$padj)))), ]
    }

    # return(res_tableOE_ordered)
    #
    res_tableOE_ordered <- res_tableOE_ordered[complete.cases(res_tableOE_ordered),]
    ## Create a column to indicate which genes to label (Top 10 genes by LFC)
    res_tableOE_ordered$genelabels <- ""
    res_tableOE_ordered$genelabels[1:topNumberOfGenes] <- rownames(res_tableOE_ordered)[1:topNumberOfGenes]
    return(res_tableOE_ordered)

}


generateVolcanoPlot <- function(res_tableOE_ordered, y_limits, x_limits, x_breaks){

    return(ggplot(res_tableOE_ordered) +
        geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
        geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = genelabels)) +
        ggtitle(paste(treated, "v", control," Volcano", sep="")) +
        xlab("log2 fold change") +
        ylab("-log10 adjusted p-value") +
        #This is y axis scale
        scale_y_continuous(limits = y_limits) +
        #This x axis scall.
        scale_x_continuous(limits = x_limits, breaks = x_breaks) +
        theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))))

}


getAllComparisonsVector <- function(formatedCSVThereshold0){
    splitColumns <- strsplit(names(formatedCSVThereshold0), "_")
    holder <- c()
    for(i in 1:length(splitColumns)){
        # print(i)
        column <- splitColumns[[i]][[2]]
        if(!(column %in% holder)){
            holder <- append(holder, splitColumns[[i]][[2]])
        }
    }
    return(holder)
}


keggEnrichmentOnPercentiles <- function(formatedCSVThereshold0, percentiles, metric, outputDir){
    holder <- py_getKEGGDataframe(formatedCSVThereshold0)

    comparisons <- getAllComparisonsVector(holder)
    filenames <- c()
    for (comp in comparisons){
        for (per in percentiles){
            #can change this so that it is dependent on score, padj or LFC
            if(metric == "score"){
                test = holder[order(holder[[paste("Score_",comp, sep="")]]),]
            } else if(metric == "padj"){
                test = holder[order(-log(holder[[paste("padj_",comp, sep="")]])),]
            } else if(metric == "LFC"){
                test = holder[order(abs(holder[[paste("LFC_",comp, sep="")]])),]

            }
            test = test[complete.cases(test),]
            index = ceiling(per*nrow(test))
            genes = rownames(test)
            upreg= c()
            downreg = c()
            for (i in index:nrow(test)){
                print(genes[i])
                if(test[genes[i],paste("LFC_",comp, sep="")] > 0){
                    upreg = append(upreg, genes[i])
                } else {
                    downreg = append(downreg, genes[i])
                }
            }
            
            if(!is.null(upreg)){
              filenames <- c(filenames, pythonKeggEnrichmentOnRVector(upreg, paste(outputDir, "/", metric, "_per_", per, "_", comp, "_UP", sep="")))
            }
            if(!is.null(downreg)){
              filenames <- c(filenames, pythonKeggEnrichmentOnRVector(downreg, paste(outputDir, "/", metric, "_per_", per, "_", comp, "_DOWN", sep="")))
            }
        }
    }

    return(filenames)
}


enrichmentOnPercentiles <- function(formatedCSVThereshold0, percentiles, metric, outputDir){

    comparisons <- getAllComparisonsVector(formatedCSVThereshold0)
    filenames <- c()
    for (comp in comparisons){
        for (per in percentiles){
            #can change this so that it is dependent on score, padj or LFC
            if(metric == "score"){
                test = formatedCSVThereshold0[order(formatedCSVThereshold0[[paste("Score_",comp, sep="")]]),]
            } else if(metric == "padj"){
                test = formatedCSVThereshold0[order(-log(formatedCSVThereshold0[[paste("padj_",comp, sep="")]])),]
            } else if(metric == "LFC"){
                test = formatedCSVThereshold0[order(abs(formatedCSVThereshold0[[paste("LFC_",comp, sep="")]])),]

            }
            test = test[complete.cases(test),]
            index = ceiling(per*nrow(test))
            genes = rownames(test)
            upreg= c()
            downreg = c()
            for (i in index:nrow(test)){
                print(genes[i])
                if(test[genes[i],paste("LFC_",comp, sep="")] > 0){
                    upreg = append(upreg, genes[i])
                } else {
                    downreg = append(downreg, genes[i])
                }
            }

            filenames <- c(filenames, pythonEnrichmentOnRVector(upreg, paste(outputDir, "/", metric, "_per_", per, "_", comp, "_UP", sep="")))
            filenames <- c(filenames, pythonEnrichmentOnRVector(downreg, paste(outputDir, "/", metric, "_per_", per, "_", comp, "_DOWN", sep="")))
        }
    }

    return(filenames)
}


#Something is wrong with the union here.
# this function is messy should change. Maybe have a funcitno that returns one list and then
# on the front end combined the lists manually?
generateEnrichmentObject <- function(enrichmentResultFileNames, comparisons, alpha=0.01){
    vectorOfEnrichmentData = list()
    for (i in seq_along(enrichmentResultFileNames)){
        # print(i)
        # print(enrichmentResultFileNames[i])
        vectorOfEnrichmentData[[i]] <- list(enrichmentResultFileNames[i], read.csv(file = enrichmentResultFileNames[i], row.names = 2))

    }
    # return(vectorOfEnrichmentData)

    for (i in seq_along(vectorOfEnrichmentData)){
        # print(vectorOfEnrichmentData[[i]][[1]])
        # print(i)
        df <- vectorOfEnrichmentData[[i]][[2]]
        holderlist <- c()
        for(row in 1:nrow(df)){
            # print(df[row, "Pvalue"])
            if(df[row, "Pvalue"] < alpha){
                print(df[row, "Pvalue"])
                print(rownames(df)[row])
                holderlist <- append(holderlist, rownames(df)[row])
            }
        }
        vectorOfEnrichmentData[[i]] <- list(vectorOfEnrichmentData[[i]][[1]], vectorOfEnrichmentData[[i]][[2]], holderlist)
    }

    unionOfGOTerms <- NULL
    for(i in seq_along(vectorOfEnrichmentData)){
        # print(i)
        # print(vectorOfEnrichDataUp[[i]][[3]])
        unionOfGOTerms = union(unionOfGOTerms, vectorOfEnrichmentData[[i]][[3]])
    }



    listOfVectorspvalues <- list()
    for(i in 1:length(vectorOfEnrichmentData)){
        listOfVectorspvalues[[i]] <- list()
    }

    for(term in unionOfGOTerms){
        for(i in 1:length(vectorOfEnrichmentData)){
            if(vectorOfEnrichmentData[[i]][[2]][term, "Pvalue"] < 0.01){
                listOfVectorspvalues[[i]] <- append(listOfVectorspvalues[[i]], TRUE)
            } else {
                listOfVectorspvalues[[i]] <- append(listOfVectorspvalues[[i]], FALSE)

            }
        }
    }

    for(i in 1:length(vectorOfEnrichmentData)){
        vectorOfEnrichmentData[[i]] <- append(vectorOfEnrichmentData[[i]],
                    list(listOfVectorspvalues[[i]]))
    }

    columns <- c()

    for(i in seq_along(vectorOfEnrichmentData)){
        filename <- vectorOfEnrichmentData[[i]][[1]]
        splitList <- strsplit(filename, "_")

        for(comp in comparisons){
            if(comp %in% splitList[[1]]){
                index <- match(comp, splitList[[1]])
                label <- paste(splitList[[1]][[index-1]], splitList[[1]][[index]], splitList[[1]][[index+1]], sep="_")
            }
        }

        vectorOfEnrichmentData[[i]] <- append(vectorOfEnrichmentData[[i]], list(label))
    }

    return(list(vectorOfEnrichmentData, unionOfGOTerms))
}

rader_plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, groups=NULL)
{
  assay_table <- as.data.frame(assay(object))
  coldata_table <- as.data.frame(colData(object))

  # print(head(assay_table))
  # print(head(coldata_table))
  # print(names(assay_table))
  # print(groups)
  # print((names(assay_table) %in% groups))
  if(!is.null(groups)){
    assay_table <- assay_table[, (names(assay_table) %in% groups)]
    coldata_table[,"holder"] <- 1
    coldata_table <- coldata_table[(rownames(coldata_table) %in% groups), ]
    coldata_table[, "holder"] <- NULL
  }

  # print(head(assay_table))
  # print(head(coldata_table))

  rv <- rowVars(as.matrix(assay_table))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay_table[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(coldata_table))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(coldata_table[, intgroup,
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    coldata_table[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
                  intgroup.df, name = colnames(assay_table))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *
                                                                                                            100), "% variance")) + coord_fixed()
}

rowFilter <- function(df, groups){

  if(is.null(groups)){
    df <- droplevels(df)
    return(df)
  }
  df[,"holder"] <- 1
  df <- df[(rownames(df) %in% groups), ]
  df[, "holder"] <- NULL
  df <- droplevels(df)
  return(df)
}
