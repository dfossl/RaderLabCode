

# Creating DESeq2Dataset
# Must pick design, the variable you want actually contrast should be second.
# Following link provides a comprehensive description of how different designs can work.
# See: https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
# I have preset work flows for different analysis for looking if there is some variable specific effect, comparing at specific timepoints, and controling for time.


#________________Time Course Analysis___________________#
#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments

# Preforms a likelhood ratio test where we remove the condition-specific
# differences over Time. Genes with small p values
# from this test are those which in one or more timepoints after
# 1 hour showed a condition-specific effect.
# this means genes that moved up or down in both samples over
# time will not get small p-values

# Creates a deseq object.
# design factors should be columns in the annotation file

# uncomment following lines for data with time points

#ddsTimeSeries<-DESeqDataSetFromMatrix( countData = cts,
#                             colData = coldata,
#                             design = ~ timepoint + condition + timepoint:condition)


# ddsTimeSeries <- DESeq(ddsTimeSeries, test="LRT", reduced = ~ timepoint + condition)

# Gets results with signifgant alpha being 0.01
#resTimeSeriesp.01 <- results(ddsTimeSeries, alpha = 0.01)
#resTimeSeriesp.01

#how many gense differ due to chosen variable "time"
#sum(resTimeSeriesp.01$padj < 0.01, na.rm=TRUE)

# orders most signifigant to least
#resTSOrdered <- resTimeSeriesp.01[order(resTimeSeriesp.01$padj),]

# stores the ordered data in CSV
#filename <- "TimeSeries.csv"
#write.csv(resTSOrdered, file = paste(OutputFileDirectory, "/", gsub(":", "-", gsub(" ", "_", Sys.time())), filename, sep=""))



# Goes through top 10 genes of most variation over time.
# change n to change number of genes you want to save
#filename <- "OverTime.tiff"
#n <- 1
#for(i in 1:n){
#  index = as.numeric(i)

  #Note you need to have the ordered data "resTSOrdered"
#  geneName <- as.character(rownames(resTSOrdered)[index])

#  p <- countPlotOver2Variables(ddsTimeSeries, geneName, "timepoint", "condition")

#  tiff(paste(OutputFileDirectory, "/", geneName, filename, sep=""), width=12, height=12, res=300)
#  print(p)
#  dev.off()

#}

#____________________END of Time Course analysis_________________#



#____________________Analysis at specific timepoints_____________#
# If time course shows genes are variable over time (or some other variable)
# then good idea to do analysis at each time
source("NewFullRNAseqPipeline_metadata.R")



dds<-DESeqDataSetFromMatrix( countData = cts,
                             colData = coldata,
                             design = configDesign)

# sorts using mincount set in first section
{
  keep<-rowSums(DESeq2::counts(dds))>=minimumCount
  dds<-dds[keep,]
}

dds <- DESeq(dds)

#_________PreProcessing___________#


counts <- DESeq2::counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))

rld <- rlog(dds, blind=T)
colData(rld)
#PCA
# pcaData <- plotPCA(rld, intgroup=c("Groups"), returnData=TRUE)

# rader_plotPCA is the same as plotPCA but allows you to give a paraemeter "groups" That is the list of columns of count data you want to use.
# pcaData <- rader_plotPCA(rld, intgroup=c("Groups"), returnData=TRUE)

# groupsToLookAt = c("CmLS_24.1Counts","CmLS_24.2Counts", "CmLS_24.3Counts",
#                    "CmLS_1.1Counts","CmLS_1.2Counts","CmLS_1.3Counts",
#                    "RM_1h_1","RM_1h_2","RM_1h_3",
#                    "RM_24h_1", "RM_24h_2", "RM_24h_3",
#                    "CmRM_60.1Counts", "CmRM_60.2Counts", "CmRM_60.3Counts")

# This is a variable for the plotPCA function that lets you subset what you want to look at
#NULL means all columns in this case.
groupsToLookAt = NULL

#Note default this is top 500 most variable genes.
pcaData <- rader_plotPCA(object=rld,
#                         intgroup=c("conditionAndtimepoint"),
                         intgroup=c(configDesignStr),
                         returnData=TRUE,
                         groups=groupsToLookAt,
                         ntop = 500)
pcaData_FileName = paste(experimentIdentifier, "_PCA.csv", sep="")
write.csv(pcaData, file.path(OutputFileDirectory, pcaData_FileName))


p <- pcaPlotFormated(pcaData,
  color=configDesignStr,
  legendTitleC = str_to_title(configDesignStr))

tiff(file.path(OutputFileDirectory,paste(experimentIdentifier, "_PCA_Plot_rlog.tiff"), sep=""), width=6, height=5,
     units = 'in', res = 150)
print(p)
dev.off()

#This code will get the weights for each gene to determine biggset contributers
#To a principal component
pcaResultsGenes <- generateGenePCAValues(rld, groups=groupsToLookAt)
holderDF <- as.data.frame(pcaResultsGenes)
rownames(holderDF) <- rownames(pcaResultsGenes)
genePCAData_name = paste(experimentIdentifier,"_Gene_PCA_Data.csv", sep="")
write.csv(holderDF, file.path(OutputFileDirectory,genePCAData_name))

tiff(file.path(OutputFileDirectory, paste(experimentIdentifier,"_Distance_Plot_rlog_ALL_timepointAsControl.tiff", sep="")), width=575*2, height=575*2, res=200)
#generateDistMatrix(rld, coldata, "conditionAndtimepoint", groups=groupsToLookAt)
generateDistMatrix(rld, coldata, configDesignStr, groups=groupsToLookAt)
dev.off()

tiff(file.path(OutputFileDirectory, paste(experimentIdentifier,"_exploratoryBarPlottemp.tiff", sep="")), width=575*2, height=575*2, res=200)
generateBarplotForTotalReadCounts(cts, groups=groupsToLookAt)
dev.off()


p <- generateNormalizedBoxplot(cts, groups=groupsToLookAt)
tiff(file.path(OutputFileDirectory,paste(experimentIdentifier,"_exploratoryBoxPlottemp.tiff", sep="")), width=575*2, height=575*2, res=200)
print(p)
dev.off()

p <- generateNormalizedDensityPlot(cts, groups=groupsToLookAt)
tiff(file.path(OutputFileDirectory,paste(experimentIdentifier,"_exploratoryDensityPlottemp.tiff", sep="")), width=575*3, height=575*2, res=200)
print(p)
dev.off()

colnames(cts)

rlogcounts <- generaterlogNormalizedCountDataFrame(cts, groups=groupsToLookAt)
coldata_filtered <- rowFilter(coldata, groups=groupsToLookAt)
coldata_filtered

numberOfGenes <- c(10,100,500,1000)

for(geneNumber in numberOfGenes){
    tiff(file.path(OutputFileDirectory,paste(experimentIdentifier,"_Z-Score.HeatMap", geneNumber, "MostVariableGenes_ALL.tiff", sep="")),
         width = 5,
         height = 6,
         units='in',
         res=250)
  
  if(geneNumber < 100){
    generateGeneCountVarianceHeatmap(rlogcounts, coldata_filtered, configDesignStr, geneNumber)
  }else{
    generateGeneCountVarianceHeatmap(rlogcounts, coldata_filtered, configDesignStr, geneNumber, labRow = FALSE)
  }
    
    dev.off()
}

#________End_PreProcessing__________#

#________Differential Expression Analysis_______#

# control vector and treatedVector should be the same size and be in the order you
# are comparing the to conditions.
# not can do any number of variable combinations.
# Note: The items of this vector should be entries in the column of coldata you are analysizing

#controlVector <- c("RM1h", "RM2h", "RM24h", "RM49h")
#treatedVector <- c("LP1h", "LP2h", "LP24h", "LP49h")

#Moved to metadata file. maybe keep here?
# controlVector <- c("RM")
# treatedVector <- c("Au")

# The variable should be the name of the column in the annotation file that holds
# the terms found in controlVector and treatedVector
#variable <- "conditionAndtimepoint"
variable <- configDesignStr


firstEntry <- TRUE
index <- 1

for(i in 1:length(controlVector)){

    control <- controlVector[i]
    print(control)
    treated <- treatedVector[i]
    print(treated)

    res <- differentialExpression(control, treated, variable, dds, alpha=0.01)

    write.csv(py_addProteinColumnToDataframe(as.data.frame(res)), file = file.path(OutputFileDirectory, paste(experimentIdentifier,"_",treated, "v", control,"_DESeq2ResultsThreshold0.csv", sep = "")))
    write.csv(py_addProteinColumnToDataframe(as.data.frame(subset(res, padj < 0.01))), file = file.path(OutputFileDirectory, paste(experimentIdentifier,"_",treated, "v", control,"_DESeq2ResultsThreshold0p0.01.csv", sep = "")))

    #Generate MA plot
    py_generateMAPlot(as.data.frame(res), OutputFileDirectory, "score", percentile=.9, comparison = paste(experimentIdentifier,"_", treated, "v", control, sep=""))
    py_pvalueVSbasemean(as.data.frame(res), OutputFileDirectory, comparison = paste(experimentIdentifier,"_",treated, "v", control, sep=""), lfc_transform="max")
    

    # this is code that is constructing the a file of all the results.
    print(nrow(res))
    if(firstEntry){
        formatedCSVThereshold0 <- data.frame(matrix(ncol=length(controlVector)*3, nrow =nrow(res) ))
        rownames(formatedCSVThereshold0) <- rownames(res)

        # formatedCSVThereshold1.99 <- data.frame(matrix(ncol=length(controlVector)*2, nrow =nrow(resLFCThreshold1.99p.01) ))
        # rownames(formatedCSVThereshold1.99) <- rownames(resLFCThreshold1.99p.01)

        comparison <- c()
        thresholdVector <- c()


        firstEntry <- FALSE
    }

    formatedCSVThereshold0 <- combinedResultDataFrames(formatedCSVThereshold0, res, index, control, treated)

    index <- index + 3

    # per_Index <- getPercentileIndex(res, 0.90)
    # print(per_Index)
    per_threshold <- getPerThreshold(res, .90, "score")
    print(per_threshold)
    
    percentiles = c(.90, .95, .99)
    for(per in percentiles){
      py_makeHistPlot(as.data.frame(res), OutputFileDirectory, "score", percentile=per, comparison=paste(experimentIdentifier,"_",treated, "v", control, sep=""))
    }
    # holder <- res[order(res[,"score"]),]
    # holder <- holder[complete.cases(holder),]
    # head(holder)
    # holder[per_Index,"score"]
    

    comparison <- append(comparison, paste(treated, "v", control, sep=""))
    thresholdVector <- append(thresholdVector, per_threshold)
    thresholdValueHolder <- data.frame(comparison, thresholdVector)

    res_ordered <- generateDataframeForVolcanoPlot(res, per_threshold, "LFC", 10)
    # print(res_ordered)
    
    scale <- 2.5
    p <- generateVolcanoPlot(res_ordered, y_limits = c(0,300), x_limits = c(-7.5,12.5), x_breaks = c(-5,0,5,10,15))
    tiff(file.path(OutputFileDirectory,paste(experimentIdentifier,"_",treated, "v", control, "_VolcanoPlot.tiff", sep="")),
        width = scale*550,
        height=scale*500,
        res=scale*100)
    print(p)
    dev.off()
    
    if(i == length(controlVector)){
      write.csv(thresholdValueHolder, file = file.path(OutputFileDirectory,paste(experimentIdentifier,"_Comparison_Thresholds.csv", sep = "")))
      write.csv(py_addProteinColumnToDataframe(formatedCSVThereshold0), file = file.path(OutputFileDirectory,paste(experimentIdentifier,"_AllLFCandPadj_Thereshold0.csv", sep = "")))
      
    }
}


#columns = c("Score_LP1hVRM1h", "Score_LP2hVRM2h", "Score_LP24hVRM24h", "Score_LP49hVRM49h")


columns = c()
for(e in comparison){
  print(e)
  columns = append(columns, paste("Score_",e,sep=""))
}


percentiles <- c(.99, .95, .90)
for(per in percentiles){
  dist_data = py_densityplotCompareColumns(formatedCSVThereshold0, columns, percentile=per, comparison=str_to_title(configDesignStr), xlabel="score", outputfiledir=OutputFileDirectory, legend=TRUE, expLabel=experimentIdentifier)
  
}


##____dot plots___Not Finished_dontRUN##

#res will hold the last group comparison in this case 49 hour timepoint
genes = py_getlistofGenes_Padj_range(as.data.frame(res), min_=600)

plotCounts(dds, gene=genes[1], intgroup=c(configDesignStr))

cnts <- counts(dds)[genes[1],]
group <- colData(dds)[[c(configDesignStr)]]
class(group)
data <- data.frame(count=cnts, group=as.integer(group))
data
plot(data$group + runif(ncol(dds),-.05,.05), data$count, xlim=c(.5,max(data$group)+.5), log="", xaxt="n", xlab="group", ylab="count", main=genes[1])
axis(1, at=seq_along(levels(group)), levels(group))

#__________Python Enrichment_________#
percentiles <- c(.99, .95, .90)
test <- keggEnrichmentOnPercentiles(formatedCSVThereshold0, percentiles, metric="score", outputDir=OutputFileDirectory)
test[0:2]

pathways <- py_GetAllPathwaysOfSig(test[0:6], p_cutoff=0.05)

keggDataFrame <- py_getKEGGDataframe(formatedCSVThereshold0)

colnames(keggDataFrame)

#kegg_LFC <- keggDataFrame[,c("LFC_LP1hVRM1h", "LFC_LP2hVRM2h", "LFC_LP24hVRM24h", "LFC_LP49hVRM49h")]
kegg_LFC <- keggDataFrame[,c("LFC_AuVRM")]

#So that legend shows green as positive and red negative.
kegg_LFC <- -1*kegg_LFC



#pathview(gene.data = kegg_LFC, pathway.id = "cme03030", species = "cme", gene.idtype = "KEGG",plot.col.key=FALSE)

#kegg_LFC["CYME_CMK133C",]

pv.out.list <- sapply(pathways, function(pid) pathview(gene.data = kegg_LFC, pathway.id = pid, species = "cme", gene.idtype = "KEGG",plot.col.key=FALSE))

# Get all pathways that are of interest. visualize eachone. Need LFCs for each timepoint, not sure if I need to rescale

#__GO Enrichment__#
# enrichmentResultFileNames <- c()
enrichmentResultFileNames <- enrichmentOnPercentiles(formatedCSVThereshold0, percentiles, metric="score", outputDir=OutputFileDirectory)
print(enrichmentResultFileNames[[1]][[1]])

comparisons <- getAllComparisonsVector(formatedCSVThereshold0)


enrichmentResultsObjectHolder <- generateEnrichmentObject(enrichmentResultFileNames, comparisons, alpha=0.01)
enrichmentResultsObject <- enrichmentResultsObjectHolder[[1]]
goTerms <- enrichmentResultsObjectHolder[[2]]
print((enrichmentResultsObject[[6]][[3]]))

mm <- data.frame(matrix(NA, nrow = length(enrichmentResultsObject[[2]][[4]]), ncol=0))

for(i in seq_along(enrichmentResultsObject)){
    mm[,enrichmentResultsObject[[i]][[5]]] <- unlist(enrichmentResultsObject[[i]][[4]], use.names = FALSE)
}
rownames(mm) <- goTerms
print(colnames(mm))
write.csv(mm, paste(OutputFileDirectory,"/", "unionofenrich_results.csv", sep=""))

# percentile and regulation are used to find the appropriate columns in mm
# you can give only percentile or only regulation if you want.
# you can also just provide a list to the "columns" variable to directly select columns.
# columns = c("col1", "col2")
# py_subsetEnrichmentResults
# you give this function mm and then filtering options
# 1) percentile ex: " percentile=0.9 "
# OR
# 2) regulation ex: " regulation="up" "
# OR
# 3) percentile and regulation ex: " percentile=0.9, regulation="up" "
# OR
# 4) A list of columns. ex: " columns=c("col1","col2") "
# The goal is to have atleast 2 columns with a MAX of 6 columns in the return.
# if you have more then that you will throw an error.
colnames(mm)
UP <- py_subsetEnrichmentResults(mm, regulation="up")

if (length(UP) < 2 | length(UP) > 6){
  print("You'll Throw an Error on Visualization. Column number is less then 2 or greater then 6.")
}

# layout = tight will try and auto scale. 
py_enrichmentVisualization(UP, OutputFileDirectory, "AuvRM_UP_enrichmentVisualized", layout="tight")

# If you get a warning that the plot can't fit you may need to scale to you liking
# bottom and left give space to labels. scale is 0-1, bigger the number more space.
py_enrichmentVisualization(UP, OutputFileDirectory, "90_AuvRM_Up_enrichmentVisualized", bottom=.15, left=.95)


# give Venn diagram using the provided dataframe. Cannot have more then 6 columns
py_venDiagramGoTerms(UP, OutputFileDirectory, "vennGO_UP")


DOWN <- py_subsetEnrichmentResults(mm, percentile=0.90, regulation="down")


py_enrichmentVisualization(DOWN, OutputFileDirectory, "90_AuvRM_Down_enrichmentVisualized", layout="tight")

py_enrichmentVisualization(DOWN, OutputFileDirectory, "90_AuvRM_Down_enrichmentVisualized", bottom=.15, left=.95)


# give Venn diagram using the provided dataframe. Cannot have more then 6 columns
py_venDiagramGoTerms(DOWN, OutputFileDirectory, "vennGO_DOWN")


# given the percentile and/or regulation OR a list of columns, intersect returns the intersect of those columns
py_intersectOfColumns_GoTerms(mm, percentile=0.90, regulation="down")

colnames(mm)
py_intersectOfColumns_GoTerms(mm, columns=c("0.9_AuVRM_UP", "0.9_AuVRM_DOWN"))

# returns a set of goterms that is EXCLUSIVELY shared among the provided columns
# careful need to use a subset of columns to start.
py_exclusiveSet(test, c("0.9_LP49hVRM49h_UP", "0.9_LP24hVRM24h_UP"))
#___Generate Ven Diagrams?__#

#__________Exploring Overlapping Genes________#
# py_subsetDifferentialExpressionResults(formatedCSVThereshold0, columns=c("Score_LP1hVRM1h", "Score_LP2hVRM2h"), percentile=.90, regulation="down")

py_geneVennDiagram(formatedCSVThereshold0, columns=c("Score_LP1hVRM1h", "Score_LP2hVRM2h", "Score_LP24hVRM24h", "Score_LP49hVRM49h"), percentile=.90, regulation="down", outputdir=OutputFileDirectory, filename="gene_venn_90_down_score")

py_geneVennDiagram(formatedCSVThereshold0, columns=c("Score_LP1hVRM1h", "Score_LP2hVRM2h", "Score_LP24hVRM24h", "Score_LP49hVRM49h"), percentile=.90, regulation="up", outputdir=OutputFileDirectory, filename="gene_venn_90_up_score")

shared_1_2_hour <- py_intersectOfColumns_SigGenes(formatedCSVThereshold0, columns=c("Score_LP1hVRM1h", "Score_LP2hVRM2h"), percentile=.90, regulation="down")

pythonEnrichmentOnRVector(geneList=shared_1_2_hour, filename=paste(OutputFileDirectory,"/enrichment_1_2h_down",sep=""))
#write funciton to do set operations and then enrichment
#_____________________________________________#


#_________k means and enrichment______#
#Most of this is in python!

#Subsets genes to the top percentiles in each group for a given variable
genesToKeep = (py_sortcolumns(formatedCSVThereshold0, .90, "score"))

sig_gene_subset_formatedCSVThereshold0 <-py_subsetDataOnGenes(formatedCSVThereshold0, genesToKeep)

#elbow plot can be used to visually pick  a cluster, optimal K is a algorthum that pick it
py_elbowPlot(sig_gene_subset_formatedCSVThereshold0, "Score", OutputFileDirectory, label="LP", max_k=50)

# if you dont want to use elbow plot this will give a suggested k value.
optimal_k <- (py_optimalK(sig_gene_subset_formatedCSVThereshold0, "lfc", nrefs=5, maxClusters=20)[[2]])
print(optimal_k)


kmeans_result <- py_performKmeans(sig_gene_subset_formatedCSVThereshold0, "lfc", 7)
write.csv(kmeans_result, paste(OutputFileDirectory, "/", "Kmeans_results.csv", sep=""))


# py_visualizeKmeans(kmeans_result, "lfc", OutputFileDirectory)

# refernece column ins a column in the provided dataframe
normalized_LFC <- py_generateNormalized_LFC(kmeans_result, "lfc", reference_column="LFC_LP1hVRM1h")
write.csv(normalized_LFC, "Normalized_LFC_kmeans_results.csv")


# py_visualizeKmeans_Normalized(normalized_LFC, "lfc", OutputFileDirectory)

# the returned format is not ideal. But is in such a form to make it easier to group and graph in python.
# the functions needs alot of inputs, reference_column is a column in coldata, reference entry is and entry in that reference column
# color_column is a column in coldata that the samples will be grouped and colored on.
# x_axis_group_column is a column in coldata that will group samples on the x axis of the plot
normalizedCountData <- py_generateNormalized_Count(df_kmeans = kmeans_result,
                                                   df_counts = cts,
                                                   coldata = coldata,
                                                   reference_column = "conditionAndtimepoint",
                                                   reference_entry = "RM1h",
                                                   color_column = "condition",
                                                   x_axis_group_column = "timepoint",
                                                   outputdir = OutputFileDirectory)
write.csv(normalizedCountData, "Normalized_count_kmeans_results.csv")
max(normalizedCountData[,"count"])

# py_countPlotKmeans(normalizedCountData, OutputFileDirectory)

py_enrichmentOnKmeans(kmeans_result)


py_kmeansGraphsCombined(kmeans_result, kmeans_result, normalizedCountData, outputdir=OutputFileDirectory, label="_NewTestMean_divide_scale", estimator="median")



#TEstING SOME THINGS


# testing_df <- py_exploringTheNormalization(df_kmeans = kmeans_result,
#                              df_counts = cts,
#                              coldata = coldata,
#                              reference_column = "TimeAndCondition",
#                              reference_entry = "RM1h")
# head(testing_df)
# describe_test <- py_lookingAtNormalizaiton(testing_df, cluster=2)
# 
# py_histogramForColumnsAndCluster(testing_df, cluster=2)
# 
# 
# x <- py_getMaxRange(df_kmeans = kmeans_result)
# 
# y <- py_getMaxRangeOfCluster(df_kmeans = kmeans_result, cluster=0)
# 
# testing <- (x - y[[1]])/2
# 
# y[[2]] + testing
# y[[3]] - testing

#____________Out of Regular Pipeline____________#

#____________Explore genes that dont make signfigant cut but still highly expressed________#
# Rader was intersted in genes that had High LFCs but didn't make the top 10%, you can ingnore this.

lfc_set_pos_1h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP1hVRM1h", threshold=3, direction="above")
percent_score_set_1h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP1hVRM1h", percentile=.90)

leftOut_HighLFC_pos_1h <- py_Difference(lfc_set_pos_1h, percent_score_set_1h)

lfc_set_neg_1h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP1hVRM1h", threshold=-3, direction="below")
percent_score_set_1h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP1hVRM1h", percentile=.90)

leftOut_HighLFC_neg_1h <- py_Difference(lfc_set_neg_1h, percent_score_set_1h)




lfc_set_pos_2h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP2hVRM2h", threshold=3, direction="above")
percent_score_set_2h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP2hVRM2h", percentile=.90)

leftOut_HighLFC_pos_2h <- py_Difference(lfc_set_pos_2h, percent_score_set_2h)

lfc_set_neg_2h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP2hVRM2h", threshold=-3, direction="below")
percent_score_set_2h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP2hVRM2h", percentile=.90)

leftOut_HighLFC_neg_2h <- py_Difference(lfc_set_neg_2h, percent_score_set_2h)




lfc_set_pos_24h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP24hVRM24h", threshold=3, direction="above")
percent_score_set_24h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP24hVRM24h", percentile=.90)

leftOut_HighLFC_pos_24h <- py_Difference(lfc_set_pos_24h, percent_score_set_24h)

lfc_set_neg_24h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP24hVRM24h", threshold=-3, direction="below")
percent_score_set_24h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP24hVRM24h", percentile=.90)

leftOut_HighLFC_neg_24h <- py_Difference(lfc_set_neg_24h, percent_score_set_24h)



lfc_set_pos_49h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP49hVRM49h", threshold=3, direction="above")
percent_score_set_49h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP49hVRM49h", percentile=.90)

leftOut_HighLFC_pos_49h <- py_Difference(lfc_set_pos_49h, percent_score_set_49h )

lfc_set_neg_49h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="LFC_LP49hVRM49h", threshold=-3, direction="below")
percent_score_set_49h <- py_subsetColumnOnThreshold(formatedCSVThereshold0, column="Score_LP49hVRM49h", percentile=.90)

leftOut_HighLFC_neg_49h <- py_Difference(lfc_set_neg_49h, percent_score_set_49h )

py_venn(c(leftOut_HighLFC_pos_1h, leftOut_HighLFC_pos_2h, leftOut_HighLFC_pos_24h, leftOut_HighLFC_pos_49h), c("LN1h", "LN2h", "LN24h", "LN49h"), OutputFileDirectory, "positive_subgroup.png")
py_venn(c(leftOut_HighLFC_neg_1h, leftOut_HighLFC_neg_2h, leftOut_HighLFC_neg_24h, leftOut_HighLFC_neg_49h), c("LN1h", "LN2h", "LN24h", "LN49h"), OutputFileDirectory, "negitive_subgroup.png")


#____________Exploring genes that are in that subgroup?___________________________#
LN1h_exclusive <- lfc_set_pos_1h
for(genes in c(lfc_set_pos_2h, lfc_set_pos_24h, lfc_set_pos_49h)){
  LN1h_exclusive <- py_Difference(LN1h_exclusive, genes)
}
LN1h_exclusive
pythonEnrichmentOnRVector(LN1h_exclusive, paste(OutputFileDirectory, "/subset_LN1h_pos_exclusive", sep=""))

intersectAll <- lfc_set_pos_1h
for(genes in c(lfc_set_pos_2h, lfc_set_pos_24h, lfc_set_pos_49h)){
  intersectAll <- py_Intersect(intersectAll, genes)
}
intersectAll
pythonEnrichmentOnRVector(intersectAll, paste(OutputFileDirectory, "/subset_pos_allConditionIntersect", sep=""))



#________generic code for exploring a gene subset__________#

#____________________Enrichment on subset genes____________________#
pythonEnrichmentOnRVector(leftOut_HighLFC_pos_1h, paste(OutputFileDirectory, "/leftOUthigh_pos_1h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN1hvRM1h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_pos_1h, "testing_1h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_pos_1h, "testing_1h")

pythonEnrichmentOnRVector(leftOut_HighLFC_neg_1h, paste(OutputFileDirectory, "/leftOUthigh_neg_1h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN1hvRM1h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_neg_1h, "testing_1h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_neg_1h, "testing_1h")




pythonEnrichmentOnRVector(leftOut_HighLFC_pos_2h, paste(OutputFileDirectory, "/leftOUthigh_pos_2h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN2hvRM2h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_pos_2h, "testing_2h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_pos_2h, "testing_2h")

pythonEnrichmentOnRVector(leftOut_HighLFC_neg_2h, paste(OutputFileDirectory, "/leftOUthigh_neg_2h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN2hvRM2h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_neg_2h, "testing_2h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_neg_2h, "testing_2h")





pythonEnrichmentOnRVector(leftOut_HighLFC_pos_24h, paste(OutputFileDirectory, "/leftOUthigh_pos_24h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN24hvRM24h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_pos_24h, "testing_24h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_pos_24h, "testing_24h")

pythonEnrichmentOnRVector(leftOut_HighLFC_neg_24h, paste(OutputFileDirectory, "/leftOUthigh_neg_24h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN24hvRM24h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_neg_24h, "testing_24h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_neg_24h, "testing_24h")




pythonEnrichmentOnRVector(leftOut_HighLFC_pos_49h, paste(OutputFileDirectory, "/leftOUthigh_pos_49h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN49hvRM49h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_pos_49h, "testing_49h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_pos_49h, "testing_49h")

pythonEnrichmentOnRVector(leftOut_HighLFC_neg_49h, paste(OutputFileDirectory, "/leftOUthigh_neg_49h", sep=""))

df <- read.csv("/Users/dfossl/Documents/Work_And_Projects/Research/RaderLab/ManuscriptNewFiguresEdits/NewFunctions_LN/LN49hvRM49h_DESeq2ResultsThreshold0.csv", row.names = 1)
py_makeHistPlotFromList(df, OutputFileDirectory, "baseMean", leftOut_HighLFC_neg_49h, "testing_49h")
py_makeHistPlotFromList(df, OutputFileDirectory, "score", leftOut_HighLFC_neg_49h, "testing_49h")


#__________________Exploring subsets and doing enrichment on those____________________#



#____________________________________________________________________________________#

geneName <- "CMH019C"
filename <- "test.png"
py_plotCounts(cts, coldata, "TimeAndCondition", gene=geneName, groups=c("RM1h", "LN1h"), outputdir=OutputFileDirectory, label="newtest")

p <- countPlotOver2Variables(dds, geneName, "timepoint", "condition")
tiff(paste(OutputFileDirectory, "/", geneName, filename, sep=""),width=1500, height=1500, res=300)
print(p)
dev.off()
