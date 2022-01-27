# These values should be the same name and should be the name
# of a column in annotation
configDesign = ~condition
configDesignStr = "condition"


# This will be tagged at the front of all generated files
# as a way to identify what analysis is being done.
experimentIdentifier = "Au"
# Dylan -> experimentIdentifier = "LPvRm"


# Differential expression comparisons
# Look at the annotation file column that mathches design and pick
# values to compare.
# Dylan -> controlVector <- c("RM")
# Dylan -> treatedVector <- c("LP")
controlVector <- c("RM")
treatedVector <- c("Au")
 