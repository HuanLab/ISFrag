#ISFrag Functions--------------------------------------------------------------------------------------------

#Helper function: match ms2 to features
matchMS2 <- function(x,featuretable, expandRt = 0, expandMz = 0, ppm = 0) {
  pks <- featuretable
  pks <- cbind(pks, 0)
  colnames(pks)[ncol(pks)] <- "mzmin"
  pks <- cbind(pks, 0)
  colnames(pks)[ncol(pks)] <- "mzmax"
  if (ppm != 0){
    mz.diff <- pks[, "mz"] * ppm / 1e6
    expandMz <- rep(expandMz, nrow(pks))
  }else{
    mz.diff <- rep(0, nrow(pks))
    expandMz <- rep(expandMz, nrow(pks))
  }
  pks[, "mzmin"] <- pks[, "mz"] - expandMz - mz.diff
  pks[, "mzmax"] <- pks[, "mz"] + expandMz + mz.diff
  pks[, "rtmin"] <- pks[, "rt"] - expandRt
  pks[, "rtmax"] <- pks[, "rt"] + expandRt

  peak_ids <- rownames(pks)
  sps <- spectra(x)
  pmz <- precursorMz(x)
  rtm <- rtime(x)
  res <- vector(mode = "list", nrow(pks))
  for (i in 1:nrow(pks)) {
    if (is.na(pks[i, "mz"])) next()
    idx <- which(pmz >= pks[i, "mzmin"] & pmz <= pks[i, "mzmax"] &
                   rtm >= pks[i, "rtmin"] & rtm <= pks[i, "rtmax"])
    if (length(idx)) {
      res[[i]] <- lapply(sps[idx], function(z) {
        z
      })
    }
  }
  names(res) <- peak_ids
  return(res)
}

#Feature detection using XCMS
generate.featuretable <- function(MS1directory, type, ppm=10, peakwidth=c(10,120), mzdiff = 0.01,
                              snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100, bw = 5,
                              mzwid = 0.015, max = 100, CAMERA = F){

  #XCMS feature detection
  cwp <- CentWaveParam(ppm = ppm,
                       peakwidth = peakwidth,
                       mzdiff = mzdiff,
                       snthresh = snthresh,
                       integrate = integrate,
                       prefilter = prefilter,
                       noise = noise)

  if(type == "single"){
    setwd(MS1directory)
    MS1.files <- list.files(pattern = ".mzXML")
    MS1data <- readMSData(MS1.files, msLevel. = 1, mode = "onDisk")
    MS1data <- findChromPeaks(MS1data, param = cwp)
    data_filtered <- filterMsLevel(MS1data, msLevel = 1L)
    xset <- as(data_filtered, 'xcmsSet')

    featureTable <- as.data.frame(xset@peaks)
    featureTable <- featureTable[order(featureTable$mz, decreasing = TRUE ),]
    rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")
    featureTable <- featureTable[, c("mz", "rt", "rtmin", "rtmax", "maxo")]

    if(CAMERA){
      xsa<-xsAnnotate(xset)
      anF <- groupFWHM(xsa, perfwhm = 0.6)
      anI <- findIsotopes(anF, mzabs = 0.01)
      anIC <- groupCorr(anI, cor_eic_th = 0.75)
      anFA <- findAdducts(anIC, polarity="positive")
      peaklist <- getPeaklist(anFA)
      peaklist <- peaklist[order(peaklist$mz),]
      featureTable <- cbind(featureTable, peaklist$isotopes)
      colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
      featureTable <- cbind(featureTable, peaklist$adduct)
      colnames(featureTable)[ncol(featureTable)] <- "Adduct"
      featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
      colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
    }

  } else{
    setwd(MS1directory)
    MS1.files <- list.files(pattern = ".mzXML")
    MS1data <- readMSData(MS1.files, msLevel. = 1, mode = "onDisk")
    MS1data <- findChromPeaks(MS1data, param = cwp, SnowParam())
    data_filtered <- filterMsLevel(MS1data, msLevel = 1L)
    xset <- as(data_filtered, 'xcmsSet')

    xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid, minsamp = 1, max = max)
    xset <- retcor(xset, method = "obiwarp", profStep = 1)
    xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid, minsamp = 1, max = max)
    xset <- fillPeaks(xset)
    XCMt <- data.frame(xset@groups)
    xcmI <- groupval(xset, value = "maxo")
    featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI)
    colnames(featureTable)[1:4] <- c("mz", "rt", "rtmin", "rtmax")
    featureTable <- as.data.frame(featureTable)
    featureTable <- featureTable[order(featureTable$mz, decreasing = TRUE ),]
    rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")

    if(CAMERA){
      xsa<-xsAnnotate(xset)
      anF <- groupFWHM(xsa, perfwhm = 0.6)
      anI <- findIsotopes(anF, mzabs = 0.01)
      anIC <- groupCorr(anI, cor_eic_th = 0.75)
      anFA <- findAdducts(anIC, polarity="positive")
      peaklist <- getPeaklist(anFA)
      peaklist <- peaklist[order(peaklist$mz),]
      featureTable <- cbind(featureTable, peaklist$isotopes)
      colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
      featureTable <- cbind(featureTable, peaklist$adduct)
      colnames(featureTable)[ncol(featureTable)] <- "Adduct"
      featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
      colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
    }
  }
  MS1.files <<- MS1.files
  MS1data <<- MS1data
  return(featureTable)
}

#Add additional features (CSV) to XCMS feature table
add.features <- function(ft_directory, ft_name, featureTable = NA){
  setwd(ft_directory)
  addFT <- read.csv(ft_name, header = T, stringsAsFactors = F)

  if(is.na(featureTable)){
    colnames(addFT)[1:4] <- c("mz", "rt", "rtmin", "rtmax")
    featureTable <- addFT
  }else{
    colnames(addFT) <- colnames(featureTable)
    featureTable <- rbind(featureTable, addFT)
  }

  dereplicatedFT <- data.frame(matrix(ncol = ncol(featureTable), nrow = 0)) #generate data frame with dereplicated features
  colnames(dereplicatedFT) <- colnames(featureTable)
  for(m in (1:nrow(featureTable))) {
    mass.lower.limit <- featureTable$mz[m] - 0.01
    mass.upper.limit <- featureTable$mz[m] + 0.01
    rt.lower.limit <- featureTable$rt[m] - 30
    rt.upper.limit <- featureTable$rt[m] + 30
    temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & dereplicatedFT$mz <= mass.upper.limit,]
    temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
    if(nrow(temp) == 0) {
      dereplicatedFT[nrow(dereplicatedFT) + 1,] = featureTable[m,]
    }else{
      index <- which(dereplicatedFT$mz == temp$mz[1])[1]
      if(sum(dereplicatedFT[index, (5:ncol(featureTable))]) < temp[1, (5:ncol(featureTable))]){
        dereplicatedFT[index, ] <- temp[1, ]
      }
    }
  }

  featureTable <- dereplicatedFT
  featureTable <- featureTable[order(featureTable$mz, decreasing = TRUE ),]
  rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")
  return(featureTable)
}

#MS2 Assignment
ms2.tofeaturetable <- function(MS2directory, featureTable){
  setwd(MS2directory)
  MS2.files <- list.files(pattern = ".mzXML")
  MS2data <- readMSData(MS2.files, mode = "onDisk")

  MS2spectra <- matchMS2(MS2data, featureTable, expandRt = 10, expandMz = 0.01, ppm = 0)
  featureTable <- cbind(featureTable,F,0,0,0,0)
  colnames(featureTable)[(ncol(featureTable)-4):ncol(featureTable)] <- c("MS2_match", "MS2mz", "MS2int",
                                                                         "PeaksCount", "fromFile")
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "ISF_level"

  for (i in 1:nrow(featureTable)) {
    if(!is.null(MS2spectra[[i]])){
      tmpSpectra <- MS2spectra[[i]]
      for (j in 1:length(tmpSpectra)){
        if(tmpSpectra[[j]]@peaksCount == 0){
          tmpSpectra[[j]] <- NA
        }
      }
      tmpSpectra <- tmpSpectra[is.na(tmpSpectra)==FALSE]
      if(length(tmpSpectra) > 0){
        currInt = tmpSpectra[[1]]@precursorIntensity
        currIdx = 1
        for(k in 1:length(tmpSpectra)){
          if(tmpSpectra[[k]]@precursorIntensity > currInt){
            currIdx = k
            currInt = tmpSpectra[[k]]@precursorIntensity
          }
        }
        finalSpectra = tmpSpectra[[currIdx]]
        indices <- finalSpectra@intensity >= 1000;
        finalSpectra@intensity <- finalSpectra@intensity[indices];
        finalSpectra@mz <- finalSpectra@mz[indices];
        featureTable$MS2_match[i] <- TRUE
        featureTable$MS2mz[i] <- paste(round(finalSpectra@mz,4),collapse = ";")
        featureTable$MS2int[i] <- paste(finalSpectra@intensity, collapse = ";")
        featureTable$PeaksCount[i] <- finalSpectra@peaksCount
        featureTable$fromFile[i] <- finalSpectra@fromFile
      }
    }
  }
  featureTable[featureTable == ""] <- 0
  return(featureTable)
}

#Feature Annotation
feature.annotation <- function(featureTable, lib_directory, lib_name, dp = 0.7, ms1.tol = 0.01, ms2.tol = 0.02){
  # Dot product function
  dp.score <- function(x,y){
    if(nrow(x)==0 | nrow(y)==0){return(0)}
    x[,2] <- 100*x[,2]/max(x[,2])
    y[,2] <- 100*y[,2]/max(y[,2])
    alignment <- data.frame(matrix(nrow=nrow(x), ncol=3))
    alignment[,1:2] <- x[,1:2]
    y1 <- y  ##in case one row in y can be selected multiple times
    for(i in 1:nrow(x)){
      mass.diff <- abs(y1[,1] - x[i,1])
      if(min(mass.diff) <= ms2.tol){
        alignment[i,3] <- y1[mass.diff==min(mass.diff),2][1]
        y1[mass.diff==min(mass.diff),1][1] <- NA   # after matched, NA assigned
        y1 <- y1[complete.cases(y1),]
        if(is.null(nrow(y1)) ==TRUE) break
        if(nrow(y1)==0) break
      }
    }
    alignment <- alignment[complete.cases(alignment),]
    if(nrow(alignment)==0){score <- 0}
    if(nrow(alignment)>0){
      #dot product calculation
      AB <- sum(alignment[,2]*alignment[,3])
      A <- sum(x[,2]^2)
      B <- sum(y[,2]^2)
      dp.score <- AB/sqrt(A*B)
      score <- as.numeric(dp.score)
    }
    match_No <- nrow(alignment)
    return  <- c(score,match_No)
    return(return)
  }

  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  setwd(lib_directory)
  database <- read.msp(lib_name, only.org = FALSE)
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "Annotation"
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "DPscore"

  rez <- foreach(x = 1:nrow(featureTable)) %dopar% {
    premass.Q <- featureTable$mz[x]     ###query precursor ion mass

    if(featureTable$MS2mz[x] == 0){
      # featureTable$Annotation[x] <- "unknown"
      return(c("unknown", 0))
    }

    ms2.Q <- data.frame(m.z = strsplit(featureTable$MS2mz[x], ";")[[1]],
                        int = strsplit(featureTable$MS2int[x], ";")[[1]])  ###query ms2 input, ncol = 2, m.z & int
    ms2.Q$m.z <- as.numeric(as.character(ms2.Q$m.z))
    ms2.Q$int <- as.numeric(as.character(ms2.Q$int))

    output <- data.frame(matrix(ncol=3))
    colnames(output) <- c('std.name','DP.score','match_No')
    h <- 1
    for(i in 1:length(database)){
      if(is.null(database[[i]]$PrecursorMZ)==TRUE) next # no precursor mass

      premass.L <- database[[i]]$PrecursorMZ # database precursor
      if(abs(premass.L-premass.Q) > ms1.tol) next # precursor filter

      ms2.L <- as.data.frame(database[[i]]$pspectrum) # database spectrum
      name.L <- database[[i]]$Name

      output[h,1] <- name.L
      output[h,2] <- dp.score(ms2.Q,ms2.L)[1]
      output[h,3] <- dp.score(ms2.Q,ms2.L)[2]

      h <- h + 1
    }
    output <- output[complete.cases(output),]

    # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
    output <- output[output[,2] >= dp,]
    # output <- output[output[,3] >= match.number.threshold,]

    if(nrow(output)==0) {
      # featureTable$Annotation[x] <- "unknown"
      return(c("unknown", 0))
    } else {
      output <- output[order(-output[,2]),] # sort by scores
      # featureTable$Annotation[x] <- output[1,1]
      # featureTable$DPscore[x] <- output[1,2]
      return(c(output[1,1], output[1,2]))
    }
  }

  for(x in 1:nrow(featureTable)){
    featureTable[x, c("Annotation", "DPscore")] <- rez[[x]]
  }
  featureTable$DPscore <- as.numeric(featureTable$DPscore)
  #clean up the cluster
  stopImplicitCluster()

  return(featureTable)
}

#Find Level3 ISF
find.level3 <- function(MS1directory, MS1.files, featureTable, type, peakCOR = 0.8, loss = 10,
                     mz.tol = 0.01, rt.tol = 30){

  #EIC peak smoothing
  peak_smooth <- function(x,level=2){
    n <- level
    if(length(x) < 2*n){
      return(x)
    } else if(length(unique(x))==1){
      return(x)
    } else{
      y <- vector(length=length(x))
      for(i in 1:n){
        y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
      }
      for(i in (n+1):(length(y)-n)){
        y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
      }
      for(i in (length(y)-n+1):length(y)){
        y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
      }
      return(y)
    }
  }

  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  if(type == "single"){
    setwd(MS1directory)
    xraw <- xcmsRaw(MS1.files, profstep=0)
    ISFtable <- foreach(i = (1:nrow(featureTable)), .packages = c("xcms", "MSnbase", "dplyr")) %dopar% {
      currFeature <- featureTable[i,]
      minRT <- currFeature$rt - 10
      maxRT <- currFeature$rt + 10
      similarFeatures <- featureTable[featureTable$rt > minRT & featureTable$rt < maxRT,]
      if(nrow(similarFeatures) == 0){
        return(NA)
      }
      similarFeatures <- similarFeatures[similarFeatures$mz <= currFeature$mz[1] - loss,]
      if(nrow(similarFeatures) == 0){
        return(NA)
      }
      for(u in 1:nrow(similarFeatures)){
        if((sum(similarFeatures[u, 5:(5+length(MS1.files)-1)])/(length(MS1.files))) >
           (sum(currFeature[1, 5:(5+length(MS1.files)-1)])/(length(MS1.files)))){
          similarFeatures[u,] <- NA
        }
      }
      similarFeatures <- similarFeatures[complete.cases(similarFeatures),]

      if(nrow(similarFeatures) == 0){
        return(NA)
      }

      curr.mass.lower.limit <- currFeature$mz[1] - mz.tol
      curr.mass.upper.limit <- currFeature$mz[1] + mz.tol
      curr.rt.lower.limit <- currFeature$rt[1] - rt.tol
      curr.rt.upper.limit <- currFeature$rt[1] + rt.tol
      # filter the features out of the retention time range
      if(curr.rt.lower.limit > tail(xraw@scantime, n=1) | curr.rt.upper.limit > tail(xraw@scantime, n=1)){
        return (NA)
      }
      if(curr.rt.lower.limit < xraw@scantime[1]+1){
        curr.rt.lower.limit <- xraw@scantime[1]+1
      }
      if(curr.rt.lower.limit < 1){
        curr.rt.lower.limit <- 1
      }
      if(curr.rt.upper.limit > tail(xraw@scantime, n=1)){
        curr.rt.upper.limit <- tail(xraw@scantime, n=1) -1
      }
      # filter the features out of the m/z range
      if(curr.mass.lower.limit < xraw@mzrange[1]){
        return (NA)
      }
      if(curr.mass.upper.limit > xraw@mzrange[2]){
        return (NA)
      }
      mzRange <- as.double(cbind(curr.mass.lower.limit, curr.mass.upper.limit))
      RTRange <- as.integer(cbind(curr.rt.lower.limit, curr.rt.upper.limit))
      eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange) #extracted EIC object
      currEIC <-eeic[["intensity"]]
      # png(file = paste0("TMP_",".png"), width = 480, height = 480)
      # eic <- plotEIC(xraw, mzrange = mzRange, rtrange = RTRange)
      # dev.off()

      putativeISF <- data.frame(matrix(ncol = ncol(featureTable)+1, nrow = 0))
      ppcor <- 0
      putativeISF <- rbind(putativeISF, cbind(currFeature, ppcor))

      for(j in 1:nrow(similarFeatures)){
        mass.lower.limit <- similarFeatures$mz[j] - mz.tol
        mass.upper.limit <- similarFeatures$mz[j] + mz.tol
        rt.lower.limit <- curr.rt.lower.limit
        rt.upper.limit <- curr.rt.upper.limit

        # filter the features out of the m/z range
        if(mass.lower.limit < xraw@mzrange[1]) next()
        if(mass.upper.limit > xraw@mzrange[2]) next()
        mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
        RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
        eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange) #extracted EIC object
        tmpEIC <- eeic[["intensity"]]

        #normalize length
        currEIC <- currEIC[1:min(length(currEIC), length(tmpEIC))]
        tmpEIC <- tmpEIC[1:min(length(currEIC), length(tmpEIC))]

        if(is.na(cor(peak_smooth(currEIC), peak_smooth(tmpEIC)))) next()

        #peak peak correlation
        ppcor <- cor(peak_smooth(currEIC), peak_smooth(tmpEIC))
        if(ppcor >= peakCOR){
          putativeISF <- rbind(putativeISF, cbind(similarFeatures[j,], ppcor))
        }
      }
      if(nrow(putativeISF) == 1){
        return(NA)
      }else{
        putativeISF$ISF_level[2:nrow(putativeISF)] <- "Level_3"
        putativeISF$ISF_level[1] == "Parent"

        dereplicatedFT <- data.frame(matrix(ncol = ncol(putativeISF), nrow = 0)) #generate data frame with dereplicated features
        colnames(dereplicatedFT) <- colnames(putativeISF)
        for(m in (1:nrow(putativeISF))) {
          mass.lower.limit <- putativeISF$mz[m] - 0.01
          mass.upper.limit <- putativeISF$mz[m] + 0.01
          rt.lower.limit <- putativeISF$rt[m] - 30
          rt.upper.limit <- putativeISF$rt[m] + 30
          temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & dereplicatedFT$mz <= mass.upper.limit,]
          temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
          if(nrow(temp) == 0) {
            dereplicatedFT[nrow(dereplicatedFT) + 1,] = putativeISF[m,]
          }else{
            index <- which(dereplicatedFT$mz == temp$mz[1])[1]
            if(sum(dereplicatedFT[index, (5:5+length(MS1.files)-1)]) < temp[1, (5:5+length(MS1.files)-1)]){
              dereplicatedFT[index, ] <- temp[1, ]
            }
          }
        }
        putativeISF <- dereplicatedFT
        putativeISF <- putativeISF[order(putativeISF$mz, decreasing = TRUE ),]
        putativeISF[putativeISF == ""] <- 0
        return(putativeISF)
      }
    }

  } else{
    setwd(MS1directory)
    xraw <- list()
    for(w in 1:length(MS1.files)){
      xraw[[w]] <- xcmsRaw(MS1.files[w], profstep=0)
    }
    ISFtable <- foreach(i = (1:nrow(featureTable)), .packages = c("xcms", "MSnbase")) %dopar% {
      currFeature <- featureTable[i,]
      minRT <- currFeature$rt - 10
      maxRT <- currFeature$rt + 10
      similarFeatures <- featureTable[featureTable$rt > minRT & featureTable$rt < maxRT,]
      if(nrow(similarFeatures) == 0){
        return(NA)
      }
      similarFeatures <- similarFeatures[similarFeatures$mz <= currFeature$mz[1] - loss,]
      if(nrow(similarFeatures) == 0){
        return(NA)
      }
      for(u in 1:nrow(similarFeatures)){
        if((sum(similarFeatures[u, 5:(5+length(MS1.files)-1)])/(length(MS1.files))) >
           (sum(currFeature[1, 5:(5+length(MS1.files)-1)])/(length(MS1.files)))){
          similarFeatures[u,] <- NA
        }
      }
      similarFeatures <- similarFeatures[complete.cases(similarFeatures),]

      if(nrow(similarFeatures) == 0){
        return(NA)
      }
      # similarFeatures[1:nrow(similarFeatures),12:16] <- 0

      putativeISF <- data.frame(matrix(ncol = ncol(featureTable)+1, nrow = 0))
      ppcor <- 0
      putativeISF <- rbind(putativeISF, cbind(currFeature, ppcor))

      for(j in 1:nrow(similarFeatures)){
        ppcor <- 0
        count <- 0
        for(k in 1:length(MS1.files)){
          if(currFeature[1, 4+k] == 0) next()
          if(similarFeatures[j, 4+k] == 0) next()

          #TARGET
          curr.mass.lower.limit <- currFeature$mz[1] - mz.tol
          curr.mass.upper.limit <- currFeature$mz[1] + mz.tol
          curr.rt.lower.limit <- currFeature$rt[1] - rt.tol
          curr.rt.upper.limit <- currFeature$rt[1] + rt.tol

          # filter the features out of the retention time range
          if(curr.rt.lower.limit > tail(xraw[[k]]@scantime, n=1) |
             curr.rt.upper.limit > tail(xraw[[k]]@scantime, n=1)) next()
          if(curr.rt.lower.limit < xraw[[k]]@scantime[1]+1){
            curr.rt.lower.limit <- xraw[[k]]@scantime[1]+1
          }
          if(curr.rt.lower.limit < 1){
            curr.rt.lower.limit <- 1
          }
          if(curr.rt.upper.limit > tail(xraw[[k]]@scantime, n=1)){
            curr.rt.upper.limit <- tail(xraw[[k]]@scantime, n=1) -1
          }
          # filter the features out of the m/z range
          if(curr.mass.lower.limit < xraw[[k]]@mzrange[1]) next()
          if(curr.mass.upper.limit > xraw[[k]]@mzrange[2]) next()
          mzRange <- as.double(cbind(curr.mass.lower.limit, curr.mass.upper.limit))
          RTRange <- as.integer(cbind(curr.rt.lower.limit, curr.rt.upper.limit))
          eeic <- rawEIC(xraw[[k]], mzrange=mzRange, rtrange=RTRange) #extracted EIC object
          currEIC <-eeic[["intensity"]]

          #ISF
          mass.lower.limit <- similarFeatures$mz[j] - mz.tol
          mass.upper.limit <- similarFeatures$mz[j] + mz.tol
          rt.lower.limit <- curr.rt.lower.limit
          rt.upper.limit <- curr.rt.upper.limit

          # filter the features out of the m/z range
          if(mass.lower.limit < xraw[[k]]@mzrange[1]) next()
          if(mass.upper.limit > xraw[[k]]@mzrange[2]) next()
          mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
          RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
          eeic <- rawEIC(xraw[[k]], mzrange=mzRange, rtrange=RTRange) #extracted EIC object
          tmpEIC <- eeic[["intensity"]]

          #normalize length
          currEIC <- currEIC[1:min(length(currEIC), length(tmpEIC))]
          tmpEIC <- tmpEIC[1:min(length(currEIC), length(tmpEIC))]

          if(is.na(cor(peak_smooth(currEIC), peak_smooth(tmpEIC)))) next()

          #peak peak correlation
          ppcor <- ppcor + cor(peak_smooth(currEIC), peak_smooth(tmpEIC))
          count <- count + 1
        }
        if(count == 0) next()
        ppcor <- ppcor/count
        if(ppcor >= peakCOR){
          putativeISF <- rbind(putativeISF, cbind(similarFeatures[j,], ppcor))
        }
      }
      if(nrow(putativeISF) == 1){
        return(NA)
      }else{
        putativeISF$ISF_level[2:nrow(putativeISF)] <- "Level_3"
        putativeISF$ISF_level[1] <- "Parent"

        dereplicatedFT <- data.frame(matrix(ncol = ncol(putativeISF), nrow = 0)) #generate data frame with dereplicated features
        colnames(dereplicatedFT) <- colnames(putativeISF)
        for(m in (1:nrow(putativeISF))) {
          mass.lower.limit <- putativeISF$mz[m] - 0.01
          mass.upper.limit <- putativeISF$mz[m] + 0.01
          rt.lower.limit <- putativeISF$rt[m] - 30
          rt.upper.limit <- putativeISF$rt[m] + 30
          temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & dereplicatedFT$mz <= mass.upper.limit,]
          temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
          if(nrow(temp) == 0) {
            dereplicatedFT[nrow(dereplicatedFT) + 1,] = putativeISF[m,]
          }else{
            index <- which(dereplicatedFT$mz == temp$mz[1])[1]
            if(sum(dereplicatedFT[index, (5:5+length(MS1.files)-1)]) < temp[1, (5:5+length(MS1.files)-1)]){
              dereplicatedFT[index, ] <- temp[1, ]
            }
          }
        }
        putativeISF <- dereplicatedFT
        putativeISF <- putativeISF[order(putativeISF$mz, decreasing = TRUE ),]
        putativeISF[putativeISF == ""] <- 0
        return(putativeISF)
      }
    }
  }
  #clean up the cluster
  stopImplicitCluster()

  ISFtable <- ISFtable[!is.na(ISFtable)]
  for(a in 1:length(ISFtable)){
    names(ISFtable)[a] <- paste0(rownames(ISFtable[[a]][1,]), "_", round(ISFtable[[a]]$mz[1], digits = 2), "_",
                                  round(ISFtable[[a]]$rt[1], digits = 0))
  }
  return(ISFtable)
}

#Find Level2 ISF
find.level2 <- function(ISFtable){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)


  ISF_list <- foreach(k = (1:length(ISFtable)), .packages = c("xcms", "MSnbase", "dplyr")) %dopar% {
    currTable <- ISFtable[[k]]
    if(currTable$MS2mz[1] == 0){
      return(currTable)
    } else{
      fragmentMZ <- as.numeric(strsplit(currTable$MS2mz[1], ";")[[1]])
      for(t in 2:nrow(currTable)){
        diff <- abs(fragmentMZ - currTable$mz[t])
        if(min(diff) <= 0.01){
          currTable$ISF_level[t] <- "Level_2"
        }
      }
      return(currTable)
    }
  }
  ISF_putative <- ISF_list[!is.na(ISF_list)]

  #clean up the cluster
  stopImplicitCluster()

  for(a in 1:length(ISF_putative)){
    names(ISF_putative)[a] <- paste0(rownames(ISF_putative[[a]][1,]), "_", round(ISF_putative[[a]]$mz[1], digits = 2), "_",
                                 round(ISF_putative[[a]]$rt[1], digits = 0))
  }
  return(ISF_putative)
}

#Find Level1 ISF
find.level1 <- function(ISF_putative, ms2.tol = 0.02){
  #Reverse dot product (x -> y, y as a template)
  reverse_dp <- function(x,y){ # x: experimental MS2, y: reference MS2
    x <- x[x[,2]>0,]
    y <- y[y[,2]>0,]
    if(nrow(x)==0 | nrow(y)==0){return(0)}
    x[,2] <- 100*x[,2]/max(x[,2])
    y[,2] <- 100*y[,2]/max(y[,2])
    alignment <- data.frame(matrix(nrow=nrow(x), ncol=3))
    alignment[,1:2] <- x[,1:2]
    y1 <- y  # in case one row in y can be selected multiple times
    for(i in 1:nrow(x)){
      mass.diff <- abs(y1[,1] - x[i,1])
      if(min(mass.diff) <= ms2.tol){
        alignment[i,3] <- y1[mass.diff==min(mass.diff),2][1]
        y1[mass.diff==min(mass.diff),1][1] <- NA   # after matched, NA assigned
        y1 <- y1[complete.cases(y1),]
        if(is.null(nrow(y1)) ==TRUE) break
        if(nrow(y1)==0) break
      }
    }
    alignment <- alignment[complete.cases(alignment),]
    if(nrow(alignment)==0){score <- 0}
    if(nrow(alignment)>0){
      #dot product calculation
      AB <- sum(alignment[,2]*alignment[,3])
      A <- sum(alignment[,2]^2)
      B <- sum(y[,2]^2)
      score <- as.numeric(AB/sqrt(A*B))
    }
    match_No <- nrow(alignment)
    matched_ratio <- as.numeric(nrow(alignment)/nrow(y))
    return  <- c(score,match_No,matched_ratio)
    return(return)
  }

  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  ISF_confirmed <- foreach(a = (1:length(ISF_putative)), .packages = c("xcms", "MSnbase", "dplyr")) %dopar% {
    currTable <- ISF_putative[[a]]
    if(currTable$MS2mz[1] == 0){
      return(currTable)
    }else{
      X <- data.frame(m.z = strsplit(currTable$MS2mz[1], ";")[[1]],
                      int = strsplit(currTable$MS2int[1], ";")[[1]]) #query MS2
      X$m.z <- as.numeric(as.character(X$m.z))
      X$int <- as.numeric(as.character(X$int))
      for(b in 2:nrow(currTable)){
        if(currTable$ISF_level[b] != "Level_2") next()
        if(currTable$MS2mz[b] != 0){
          Y <- data.frame(m.z = strsplit(currTable$MS2mz[b], ";")[[1]],
                          int = strsplit(currTable$MS2int[b], ";")[[1]]) #query MS2 input
          Y$m.z <- as.numeric(as.character(Y$m.z))
          Y$int <- as.numeric(as.character(Y$int))

          # currTable$DPscore[b] <- reverse_dp(X,Y)[1]
          # currTable$matched_ratio[b] <- reverse_dp(X,Y)[3]
          if(reverse_dp(X,Y)[1] > 0.5 |reverse_dp(X,Y)[3] > 0.7){
            currTable$ISF_level[b] <- "Level_1"
          }
        }
      }
      return(currTable)
    }
  }
  for(a in 1:length(ISF_confirmed)){
    names(ISF_confirmed)[a] <- paste0(rownames(ISF_confirmed[[a]][1,]), "_", round(ISF_confirmed[[a]]$mz[1], digits = 2), "_",
                                     round(ISF_confirmed[[a]]$rt[1], digits = 0))
  }

  #clean up the cluster
  stopImplicitCluster()
  return(ISF_confirmed)
}

#Produce ISFrag Results
get.ISFrag.results <- function(ISF_List, featureTable){
  results <- list()
  for(a in 1:length(ISF_List)){
    names(ISF_List)[a] <- rownames(ISF_List[[a]][1,])
    ISF_List[[a]] <- ISF_List[[a]][ISF_List[[a]]$ISF_level != "Level_3", ]
    if(nrow(ISF_List[[a]]) == 1){
      ISF_List[a] <- NA
    }
  }
  ISF_List <- ISF_List[!is.na(ISF_List)]
  treeList <- list()
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "Num_Level2"
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "Num_Level1"
  for(m in 1:length(ISF_List)){
    index <- which(rownames(featureTable) == names(ISF_List)[[m]])
    if(featureTable$ISF_level[index] == 0){
      featureTable$ISF_level[index] <- paste(names(ISF_List)[[m]],"Parent",sep = "<-")
      assign(paste0(names(ISF_List)[[m]], "_", names(ISF_List)[[m]]), Node$new(paste0(rownames(featureTable)[index], ": mz=",
                                                   round(featureTable$mz[index], digits = 2), ", rt=",
                                                   round(featureTable$rt[index], digits = 0), ", ParentFeature")))
      searchSpace <- ISF_List[[m]][2:nrow(ISF_List[[m]]),]
      searchSpace <- cbind(searchSpace, rownames(ISF_List[[m]][1,]))
      colnames(searchSpace)[ncol(searchSpace)] <- "Parent"
      searchSpace <- cbind(searchSpace, rownames(searchSpace))
      colnames(searchSpace)[ncol(searchSpace)] <- "ID"
      while (nrow(searchSpace) != 0) {
        tmpFeature <- searchSpace[1,]
        searchSpace[1,] <- NA
        searchSpace <- searchSpace[complete.cases(searchSpace),]

        tmpindex <- which(rownames(featureTable) == tmpFeature$ID)
        if(featureTable$ISF_level[tmpindex] == 0){
          featureTable$ISF_level[tmpindex] <- paste0(names(ISF_List)[[m]], "<-",
                                                     tmpFeature$ISF_level)
        }else{
          featureTable$ISF_level[tmpindex] <- paste0(featureTable$ISF_level[tmpindex], ";",
                                                     names(ISF_List)[[m]], "<-",
                                                     tmpFeature$ISF_level)
        }
        if(tmpFeature$ISF_level == "Level_2"){
          featureTable$Num_Level2[tmpindex] <- featureTable$Num_Level2[tmpindex] + 1
        }else{
          featureTable$Num_Level1[tmpindex] <- featureTable$Num_Level1[tmpindex] + 1
        }
        assign(paste0(tmpFeature$ID, "_", names(ISF_List)[[m]]),
               eval(as.name(paste0(tmpFeature$Parent, "_", names(ISF_List)[[m]])))$AddChild(paste0(tmpFeature$ID, ": mz=",
                                                                round(tmpFeature$mz, digits = 2), ", rt=",
                                                                round(tmpFeature$rt, digits = 0), ", PPcor=",
                                                                round(tmpFeature$ppcor, digits = 2), ", level=",
                                                                tmpFeature$ISF_level)))

        newindex <- which(names(ISF_List) == tmpFeature$ID)
        if(length(newindex) == 1){
          toBind <- ISF_List[[newindex]]
          toBind <- toBind[2:nrow(toBind),]
          toBind <- cbind(toBind, tmpFeature$ID)
          colnames(toBind)[ncol(toBind)] <- "Parent"
          toBind <- cbind(toBind, rownames(toBind))
          colnames(toBind)[ncol(toBind)] <- "ID"
          uniqueToBind <- toBind
          # uniqueToBind <- anti_join(toBind, searchSpace, by = c("mz" = "mz"))
          searchSpace <- rbind(searchSpace, uniqueToBind)
        }
      }
      treeitem <- eval(as.name(paste0(names(ISF_List)[[m]], "_", names(ISF_List)[[m]])))
      treeList[[length(treeList)+1]] <- treeitem
      names(treeList)[length(treeList)] <- names(ISF_List)[[m]]
    }
  }

  results[[1]] <- treeList
  results[[2]] <- featureTable
  names(results) <- c("TreesList", "FeatureTable")
  return(results)
}

#Plot ISF Tree by FeatureID
plot.tree.single <- function(ISFresult, featureID, directory){
  setwd(directory)
  ID <- paste0(featureID, "_", featureID)
  treeIndex <- which(names(ISFresult[["TreesList"]]) == ID)
  tree <- ISFresult[["TreesList"]][[treeIndex]]
  SetGraphStyle(tree, rankdir = "TB")
  SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2)
  SetNodeStyle(tree, style = "filled,rounded", shape = "box", fillcolor = "steelblue1",
               tooltip = GetDefaultTooltip)
  Do(tree$children, function(node) SetNodeStyle(node, style = "filled,rounded", shape = "box", fillcolor = "tan1",
                                                tooltip = GetDefaultTooltip))
  treeplot <- ToDiagrammeRGraph(tree)
  export_graph(treeplot,
               file_name = paste0(featureID,"_",
                                  round(ISFresult[["FeatureTable"]][featureID, "mz"], digits = 2), "_",
                                  round(ISFresult[["FeatureTable"]][featureID, "rt"], digits = 0), ".png"),
               file_type = "png",
               title = paste0("Tree Diagram of Feature ", featureID))
  render_graph(treeplot)
}

#Plot all ISF Trees
plot.tree.all <- function(ISFresult, directory){
  setwd(directory)
  for(i in 1:length(ISFresult[["TreesList"]])){
    tree <- ISFresult[["TreesList"]][[i]]
    featureID <- names(ISFresult[["TreesList"]])[i]
    SetGraphStyle(tree, rankdir = "TB")
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2)
    SetNodeStyle(tree, style = "filled,rounded", shape = "box", fillcolor = "steelblue1",
                 tooltip = GetDefaultTooltip)
    Do(tree$children, function(node) SetNodeStyle(node, style = "filled,rounded", shape = "box", fillcolor = "tan1",
                                                  tooltip = GetDefaultTooltip))
    treeplot <- ToDiagrammeRGraph(tree)
    export_graph(treeplot,
                 file_name = paste0(featureID,"_",
                                    round(ISFresult[["FeatureTable"]][featureID, "mz"], digits = 2), "_",
                                    round(ISFresult[["FeatureTable"]][featureID, "rt"], digits = 0), ".png"),
                 file_type = "png",
                 title = paste0("Tree Diagram of Feature ", featureID))
    # render_graph(treeplot)
  }
}

#Export ISFrag Results
export.ISFrag.results <- function(ISFresult, directory){
  setwd(directory)
  write.csv(ISFresult[["FeatureTable"]], file = "ISFrag_Results.csv", row.names = T, col.names = T)
}

#Export Detailed ISFrag Feature Info
export.ISFrag.detailed <- function(ISF_List, featureID, directory){
  setwd(directory)
  write.csv(ISF_List[[which(grepl(featureID, names(ISF_List), fixed=TRUE)==T)]],
            file = paste0(featureID, "_Detailed.csv"), row.names = T, col.names = T)
}


