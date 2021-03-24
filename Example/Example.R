library(ISFrag)
MS1directory <- "X:/Users/Sam_Shen/test_20201221/fullscan"
MS2directory <- "X:/Users/Sam_Shen/test_20201221/DDA"
type <- "multi"

featureTable <- generate.featuretable(MS1directory = MS1directory, MS2directory = MS2directory, type = type)
featureTable <- ms2.tofeaturetable(MS2data = MS2data, featureTable = featureTable, type = type)
level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1.files, featureTable = featureTable, type = type)
level2 <- find.level2(ISFtable = level3)
level1 <- find.level1(ISF_putative = level2)
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)

plot.tree.all(ISFresult = results, directory = MS2directory)
export.ISFrag.results(ISFresult = ISFresult, directory = MS2directory)
