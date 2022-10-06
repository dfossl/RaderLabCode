# RaderLabCode
General Code Storage for Rader Lab

# Current up to date Work

Folders NewScripts and Annotations are the most current relivent scripts.

## Description of Running Process.
Okay so everything should work, which are some famous last words.

Important Files:

I split the main code into two scripts "NewScripts/R_Initialize.R" and "NewScripts/NewFullRNAseqPipeline.R". The idea is that one is for initializing the directories and loading python and all that and the other actually has analysis. "R_Initialize.R" should be run first with "NewFullRNAseqPipeline.R" being second. There are three metadata files in NewScripts "global_python_meta.py", "R_initialize_metadata.R", and "NewFullRNAseqPipeline_metadata.R". 

"global_python_meta.py" only has one variable which is the annotation directory.

"R_initialize_metadata.R" holds other directory information as well as the directory information for the count and annotation files.

"NewFullRNAseqPipeline_metadata.R" has the DESeq design, an experiment identifier used for tagging outputs and the control and treated vectors. 

"DifferentialExpressionFunctions.R" This holds all the R functions.

"Deseq2EnrichmentWorkFlowFunctions.py" This hold all the python functions.

"NewFullRNAseqPipeline.R" is the DEseq pipeline.

## Data

There is a data folder. It has LP data and a file with all experiments I have seen. Dr. Rader has all the data this is hear to have something to play and test with primarily

## Running Scripts

Metadata should be set first. "global_python_meta.py" and "R_initialize_metadata.R" are essentially just setting directories and pointing a file locations. In "global_python_meta.py" the AnnotationFileDirectory points to the DESeq annotation file that labels conditions. Rest should be self explanitory.

"NewFullRNAseqPipeline_metadata.R" is used with the goal of seperating experiment specific things from the pipeline. So here you set the DESeq design and some other experiment specific information.

"R_Initialize.R" can be ran first. This is often the first hiccup, it will prompt to install the necessary packages BUT I find R quite bad with dependencies and will often error out saying it failed to install. I find that simply reading the error and manually installing the dependency listed in the error with install.package() or bioconductor gets the job done.

I don't have a script that installs the python packages or the dependency list but "Deseq2EnrichmentWorkFlowFunctions.py" has all the packages at the top.

You should only have to run "NewFullRNAseqPipeline.R" and you will run it line by line.

"##____dot plots___Not Finished_dontRUN##" This marks the end where everything is basicly automated. Enrichment hasn't been automated yet so there is still some manually things that need to be done and filled in.

Notably, KEGG related stuff is between 318-343. Essentially that pathview function on line 340 is used to construct Kegg pathways. However, as far as I could look into at the time the visualizations were very hard to customize. That package with the pathview function can be looked into more. OR other visualization packages can be explored.

"#____________Out of Regular Pipeline____________#" marks the end of the pipeline. 