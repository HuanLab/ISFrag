library(ISFrag)

MS1directory <- "F:/Jian_Guo/ISFrag_NISTPlasma_Parameters_20210205/CollisionEnergy/9/3/f"
MS2directory <- "F:/Jian_Guo/ISFrag_NISTPlasma_Parameters_20210205/CollisionEnergy/9/3/d"
type <- "single"
lib_directory <- "X:/Users/Sam_Shen/Library"
lib_name <- "convertedLibraryPos.msp"
ft_directory <- "F:/Jian_Guo/ISFrag_NISTPlasma_Parameters_20210205/CollisionEnergy/9/3/f"
ft_name <- "collision9_3.csv"

xcmsFT <- XCMS.featuretable(MS1directory = MS1directory, type = type, peakwidth = c(5,20))
customFT <- custom.featuretable(ft_directory = ft_directory, ft_name = ft_name)

featureTable <- ms2.assignment(MS2directory = MS2directory, XCMSFT = xcmsFT, customFT = customFT)
# featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)

level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1.files, featureTable = featureTable, type = type)
level2 <- find.level2(ISFtable = level3)
level1 <- find.level1(ISF_putative = level2)
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)
output_dir <- "F:/Jian_Guo/ISFrag_NISTPlasma_Parameters_20210205/CollisionEnergy/9/3/d"
plot.tree.single(ISFresult = results, featureID = "F7539", directory = output_dir)

