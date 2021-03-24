library(ISFrag)
MS1directory <- "X:/Users/Sam_Shen/ISFtest20210127/HILIC(+)/HILIC(+)3/fullscan"
MS2directory <- "X:/Users/Sam_Shen/ISFtest20210127/HILIC(+)/HILIC(+)3/DDA"
type <- "single"
lib_directory <- "E:/SAM"
lib_name <- "convertedLibraryPos.msp"
ft_directory <- "X:/Users/Sam_Shen/ISFtest20210127/HILIC(+)/HILIC(+)3/fullscan"
ft_name <- "HILIC(+)featuretable.csv"

featureTable <- generate.featuretable(MS1directory = MS1directory, type = type, peakwidth = c(10,60))
featureTable <- add.features(featureTable = featureTable, ft_directory = ft_directory, ft_name = ft_name)
featureTable <- ms2.tofeaturetable(MS2directory = MS2directory, featureTable = featureTable)
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)
level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1.files, featureTable = featureTable, type = type)
level2 <- find.level2(ISFtable = level3)
level1 <- find.level1(ISF_putative = level2)
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)

plot.tree.all(ISFresult = results, directory = MS2directory)
export.ISFrag.results(ISFresult = results, directory = MS2directory)
plot.tree.single(ISFresult = results, featureID = "F4260", directory = MS2directory)






ISF_List <- level1
for(a in 1:length(ISF_List)){
  names(ISF_List)[a] <- rownames(ISF_List[[a]][1,])
  ISF_List[[a]] <- ISF_List[[a]][ISF_List[[a]]$ISF_level != "Level_3", ]
  if(nrow(ISF_List[[a]]) == 1){
    ISF_List[a] <- NA
  }
}
ISF_List <- ISF_List[!is.na(ISF_List)]


summary <- data.frame(matrix(nrow = length(results[["TreesList"]]), ncol = 5))
summary[is.na(summary)] <- ""
for(i in 1:length(results[["TreesList"]])){
  loc <- which(names(ISF_List) == names(results[["TreesList"]])[i])
  currTable <- ISF_List[[loc]]
  summary[i, 1] <- rownames(currTable)[1]
  summary[i, 2] <- currTable$mz[1]
  summary[i, 3] <- currTable$rt[1]
  for(j in 2:nrow(currTable)){
    if(currTable$ISF_level[j] == "Level_1"){
      summary[i, 4] <- paste0(summary[i, 4], round(currTable$mz[j], digits = 4), " (",round(currTable$ppcor[j], digits = 2), "); ")
    }else{
      summary[i, 5] <- paste0(summary[i, 5], round(currTable$mz[j], digits = 4), " (",round(currTable$ppcor[j], digits = 2), "); ")
    }
  }
}
colnames(summary) <- c("Feature ID",  "mz", "rt", "Level 2 ISF", " Level 1 ISF")

setwd("X:/Users/Sam_Shen/IsFrag_StandardsTest_20210106")
write.csv(summary, file = "output(-)20200114.csv", col.names = T, row.names = F)



#
#
# png(file = "tmp.png", width = 480, height = 480)
# plot(peak_smooth(tmpEIC))
# lines(peak_smooth(currEIC))
# dev.off()
#
# png(file = paste0(featureTable$mz[988],"_", featureTable$RT[988],".png"), width = 480, height = 480)
# eic <- plotEIC(xraw[[1]], mzrange = c((featureTable$mz[988] - mz.tol), (featureTable$mz[988] + mz.tol)),
#                rtrange = c((featureTable$RT[988] - rt.tol), (featureTable$RT[988] + rt.tol)))
# dev.off()
#
# png(file = paste0(featureTable$mz[2015],"_", featureTable$RT[2015],".png"), width = 480, height = 480)
# eic <- plotEIC(xraw[[3]], mzrange = c((featureTable$mz[1617] - mz.tol), (featureTable$mz[1617] + mz.tol)),
#                rtrange = c((featureTable$RT[1617] - rt.tol), (featureTable$RT[1617] + rt.tol)))
# eic <- plotEIC(xraw[[3]], mzrange = c((featureTable$mz[2015] - mz.tol), (featureTable$mz[2015] + mz.tol)),
#                rtrange = c((featureTable$RT[2015] - rt.tol), (featureTable$RT[2015] + rt.tol)))
# dev.off()
