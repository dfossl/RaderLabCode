# These values should be the same name and should be the name
# of a column in annotation
configDesign = ~conditionAndtimepoint
configDesignStr = "conditionAndtimepoint"


# This will be tagged at the front of all generated files
# as a way to identify what analysis is being done.
# experimentIdentifier = "Au"
experimentIdentifier = "LPvRm"


# Differential expression comparisons
# Look at the annotation file column that mathches design and pick
# values to compare.
controlVector <- c("RM1h", "RM2h", "RM24h", "RM49h")
treatedVector <- c("LP1h", "LP2h", "LP24h", "LP49h")
# controlVector <- c("RM")
# treatedVector <- c("Au")
 