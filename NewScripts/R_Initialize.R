# Working Directory (Default is current working directory.)
workingDirectory = "C:\\Users\\dfoss\\Documents\\Projects\\RaderLab\\RaderLabCode\\NewScripts"
setwd(workingDirectory)

# This loads all the function from the DifferentialExpressionFunctions.R script. Make sure its the right directory for you
source("R_initialize_metadata.R")
source("DifferentialExpressionFunctions.R")


set.seed(123)

#Setting up python instance
#cmd- which python3 gets path
# python working directory is not R working directory but is the directory RStudio is opened in
# lame
{
use_python(pythonInstanceDir, required = T)
py_config()
source_python(pythonModuleDir)
}


#Input Data


# Setting Working Directory
# Shows current working directory.
# getwd()

# If working Directory should change place path here.
# setwd(workingDirectory)

# This location is where ouput Files will Start. Should be name of folder in working directory
OutputFileDirectory <- file.path(workingDirectory, outputDirName)

# Location of data
# This is usually in the working directory but doesn't have to be.
#Keep all data in one directory and padting the specific directory here
# can be useful for not having duplicates
CountFileDirectory <- CountFileDirectory
# AnnotationFileDirectory <- "/Users/dfossl/OneDrive/Documents/Dylan_School_Cloud/Rader Lab/Analysis-LP/Deseq2Analysis_LowPhosphorous/LP_annotation.csv"
AnnotationFileDirectory <- AnnotationFileDirectory

cts <- read.csv(CountFileDirectory, row.names=1)
coldata <- read.csv(AnnotationFileDirectory, row.names=1)
coldata <- droplevels(coldata)

# WARNING: A common error is having numbers in the annotation makes those columns not be considered Factors
# You can force a column to be considered a factor with the following code.
# coldata$column <- factor(coldata$column)


#Makes sure the columns in data match rows in annotation.
checkColumnsMatch(coldata, cts)
rownames(coldata) %in% colnames(cts)
if(!checkColumnsMatch(coldata, cts)){
  print("ERROR, you have rows and columns with different names or exrta rows or columns.
          Therefore reformat Data till TRUE. ")
}

#expVariables holds list of variable conditions
expvariables <- colnames(coldata)
expvariables

# Set what the minimum count you wish each row to sum too.
# Deseq2 Documentation claims more robust filtering in making the Deseq object
# So I am trusting them on this,
minimumCount <- 10

