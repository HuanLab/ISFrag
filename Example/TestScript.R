library(ISFrag)

MS1directory <- "X:/Users/Jian_Guo/DaDIA_20200710/Leukemiaapplication_20200819/DaDIA/DIA"
MS2directory <- "X:/Users/Jian_Guo/DaDIA_20200710/Leukemiaapplication_20200819/DaDIA/DDA"
# type <- "multi"
# lib_directory <- "X:/Users/Sam_Shen/Library"
# lib_name <- "convertedLibraryPos.msp"
# ft_directory <- "F:/Jian_Guo/ISFrag_NISTPlasma_Parameters_20210205/CollisionEnergy/9/3/f"
# ft_name <- "collision9_3.csv"

# xcmsFT <- XCMS.featuretable(MS1directory = MS1directory, type = type, peakwidth = c(5,20))
# customFT <- custom.featuretable(ft_directory = ft_directory, ft_name = ft_name)
#
# featureTable <- ms2.assignment(MS2directory = MS2directory, XCMSFT = xcmsFT, customFT = customFT)
# # featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)

xcmsFT <- XCMS.featuretable(MS1directory = MS1directory, type = "multi", peakwidth = c(5,80))
featureTable <- ms2.assignment(MS2directory = MS2directory, XCMSFT = xcmsFT)
# Identify level 3 in-source fragments.
level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1.files, featureTable = featureTable, type = "multi")

# Identify level 2 in-source fragments.
level2 <- find.level2(ISFtable = level3)

# Identify level 1 in-source fragments.
level1 <- find.level1(ISF_putative = level2)
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)
resultFT <- export.ISFrag.results(ISFresult = results)
setwd(MS1directory)
write.csv(resultFT, file = "hi.csv", row.names = F, col.names = T)
