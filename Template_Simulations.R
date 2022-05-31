###################################################
##### --------- GERMAN UNGULATES ------------ #####
##### --- SCR-VIOLATIONS SIMULATION STUDY --- #####
###################################################

## This script presents the template for the simulation and analysis of SCR data
## based on different ungulate populations characteristics.
## This particular simulation mimicks a population with a specific age- and sex-
## composition and corresponding space-use patterns (represented by variable sigma/p0 
## parameters of the half-normal detection function.
rm(list=ls())
gc()


## ------ LIBRARIES ------
library(rgdal)
library(rgeos)
library(raster)
library(nimble)
library(nimbleSCR)
library(abind)
library(coda)
library(sf)
library(fasterize)
library(R.utils)
library(plotrix)
library(extraDistr)



## ------ WORKING DIRECTORIES ------
WD <- file.path("C:/Users/pidu/OneDrive - Norwegian University of Life Sciences/PROJECTS/SIMULATION_STUDY")
dir.create(WD)
path.IN <- file.path(WD, "InFiles")
dir.create(path.IN)
path.OUT <- file.path(WD, "OutFiles")
dir.create(path.OUT)


## ------ FUNCTIONS ------
sourceDirectory(file.path(WD, "functions"),
                modifiedOnly = FALSE)


## -----------------------------------------------------------------------------
## ------   1. SET SIMULATION CHARACTERISTICS ------
## Load parameter space
parm.df <- read.csv2( file = file.path(WD, "PARAM.SPACE/param.space.csv"),
                      header = T, sep = ",", dec = ".")

temp <- MakeParameterSpace(list.param = list("Grouping" = c(FALSE, TRUE),
                                             "Heterogeneity" = c(FALSE, TRUE),
                                             "Adaptive" = c(FALSE, TRUE),
                                             "searchEffort" = c("high","low"),
                                             "species" = c("chamois","red_deer","wild_boar")),
                           n.rep = 1)
parm.df <- merge(temp, parm.df, by = "species")

## Assign scenario number to each simulation
ref <- expand.grid(list(G = c(F,T), H = c(F,T), A = c(F,T)))
for(s in 1:nrow(parm.df)){
  parm.df$scenario[s] <- which(ref$G == parm.df$Grouping[s] &
                                 ref$H == parm.df$Heterogeneity[s] &
                                 ref$A == parm.df$Adaptive[s])
}#s

## Identify when the simulation changes species
changeSpecies <- unlist(lapply(unique(parm.df$species),function(x)min(which(parm.df$species == x))))


for(sim in 1:nrow(parm.df)){
  ## Set study area characteristics
  species <- parm.df$species[sim]
  site <- parm.df$site[sim]
  habitatFile <- parm.df$habitatFile[sim]
  detectorFile <- parm.df$detectorFile[sim]
  D <- parm.df$D[sim]
  habitatRes <- 500
  detectorRes <- 200
  scenario <- parm.df$scenario[sim]
  
  ## Set population characteristics
  pop.space <- read.csv2( file = file.path(WD, "PARAM.SPACE/",
                                           paste0(species, ".csv")),
                          header = T, sep = ",", dec = ".")
  numGroups <- nrow(pop.space)
  propGroups <- pop.space$groupProp
  prop <- pop.space[ ,grep("prop_", names(pop.space))]
  p0 <- pop.space[ ,grep("p0", names(pop.space))]
  if(parm.df$searchEffort[sim] == "high") p0 <- p0*1.5
  sigma <- pop.space[ ,grep("sigma", names(pop.space))]
  sex <- pop.space[ ,grep("sex", names(pop.space))]
  GS <- pop.space[ ,grep("GS", names(pop.space))]
  betaP0 <- 0.2
  
  ## SIMULATE A NEW HABITAT ONLY ONCE PER SPECIES TO SAVE TIME
  if(sim %in% changeSpecies){
    ## ------     1.1. STUDY AREA ------
    studyArea <- readOGR(file.path(WD, "PARAM.SPACE/StudyArea.shp"))
    plot(studyArea, col = "gray60")
    
    
    ## ------     1.2. DETECTORS ------
    searchTracks <- readOGR(file.path(WD, "PARAM.SPACE/Transects.shp"))
    searchTracks <- spTransform(searchTracks, CRSobj = proj4string(studyArea))
    searchTracks <- gSimplify(searchTracks, tol = 50)
    plot(searchTracks, add = T, col = "firebrick3")
    
    detectors <- MakeSearchGrid( data = searchTracks,
                                 resolution = detectorRes,
                                 fasterize = T,
                                 plot = F)
    
    detCoords <- coordinates(detectors$detector.sp)
    detCovs <- detectors$detector.sp$count
    numDetectors <- dim(detCoords)[1]
    
    
    ## ------     1.3. HABITAT ------
    ## CREATE POLYGON OF SUITABLE HABITAT
    ## by considering a buffer of 3*max(sigma) around the search tracks
    ## and joining with the study area polygon
    habPoly <- gBuffer(spgeom = searchTracks, width = 3*max(sigma, na.rm = T))
    habPoly <- union(habPoly, studyArea)
    habPoly <- aggregate(habPoly)
    plot(habPoly, add = T, col = adjustcolor("forestgreen", alpha.f = 0.3))
    
    habArea <- gArea(habPoly)/1e6
    
    ## CREATE A RASTER OF SUITABLE HABITAT
    habRaster <- raster(extent(habPoly), resolution = habitatRes)
    habRaster <- rasterize(habPoly, habRaster)
    plot(habRaster, add = T, col = adjustcolor("forestgreen", alpha.f = 0.3))
    
    ## CREATE A MATRIX OF SUITABLE HABITAT
    habMatrix <- as.matrix(habRaster)
    habMatrix[is.na(habMatrix)] <- 0
    
    ## IDENTIFY WHICH RASTER CELLS ARE SUITABLE HABITAT
    isHabitat <- which(habRaster[] %in% 1)
    
    
    
    ## ------     1.4. SCALE DETECTORS & HABITAT ------
    ## SCALE HABITAT GRID CELLS
    scaledHabCoords <- scaleCoordsToHabitatGrid(
      coordsData = coordinates(habRaster)[isHabitat, ],
      coordsHabitatGridCenter = coordinates(habRaster)[isHabitat, ])$coordsDataScaled
    
    scaledUpperHabCoords <- scaledHabCoords + 0.5
    scaledLowerHabCoords <- scaledHabCoords - 0.5
    
    numHabCells <- dim(scaledHabCoords)[1]
    
    ## CREATE AN HABITAT GRID
    habGrid <- matrix(0, dim(habRaster)[1], dim(habRaster)[2])
    for(i in 1:numHabCells){
      habGrid[trunc(scaledHabCoords[i,2])+1,
              trunc(scaledHabCoords[i,1])+1] <- i
    }#i
    
    ## SCALE DETECTOR COORDINATES
    scaledDetCoords <- scaleCoordsToHabitatGrid(
      coordsData = coordinates(detectors$detector.sp),
      coordsHabitatGridCenter = coordinates(habRaster)[isHabitat, ])$coordsDataScaled
    
    
    
    ## ------     1.5. LOCAL EVALUATION APPROACH ------
    localDetectors <- getLocalObjects( habitatMask = habMatrix,
                                       coords = scaledDetCoords,
                                       dmax = 4*max(sigma, na.rm = T)/habitatRes,
                                       resizeFactor = 1,
                                       plot.check = F)
    
    
    
  }#if
  
  
  
  ## -------------------------------------------------------------------------
  ## ------   2. SET POPULATION CHARACTERISTICS ------
  ## If no grouping ; set all group sizes to 1
  if(!parm.df$Grouping[sim]){
    unGroup <- TRUE
  }
  
  ## If no heterogeneity ; set all detection parameters to average pop values.
  if(!parm.df$Heterogeneity[sim]){
    p0[] <- sum(p0 * (propGroups * prop), na.rm = T)
    sigma[] <- sum(sigma * (propGroups * prop), na.rm = T)
  }
  
  ## If no adaptive sampling; set betaP0 = 0
  if(!parm.df$Adaptive[sim]){
    betaP0 <- 0
  }
  
  
  
  ## ------     2.1. SIMULATE GROUP SIZES ------
  ## Get expected population size
  exp.N <- D * habArea
  
  ## Sample group sizes
  groupSize <- groupType <- NULL
  index <- 1
  for(g in 1:numGroups){
    temp.N <- exp.N * propGroups[g]
    while(temp.N > 0){
      groupSize[index] <- rtpois(1, lambda = GS[g,2], a = GS[g,1], b = GS[g,3])
      groupType[index] <- g
      temp.N <- temp.N - groupSize[index]
      index <- index + 1
    }
  }#g
  
  if(unGroup){
    groupType <- rep(groupType, groupSize)
    groupSize <- rep(1, sum(groupSize))
  }
  
  
  ## Get realized population size
  sim.N <- sum(groupSize)
  
  ## Get realized density
  sim.D <- sim.N/habArea
  
  
  ## ------     2.2. SIMULATE GROUP ACs ------
  habIntensities <- rep(1, numHabCells)
  logIntensities <- log(habIntensities)
  
  ## Simulate AC locations
  groupACs <- do.call(rbind, lapply(1:length(groupSize),
                                    function(x){
                                      rbernppAC(
                                        n = 1,
                                        lowerCoords = scaledLowerHabCoords,
                                        upperCoords = scaledUpperHabCoords,
                                        logIntensities = logIntensities,
                                        logSumIntensity = log(sum(habIntensities)),
                                        habitatGrid = habMatrix,
                                        numGridRows = nrow(habMatrix),
                                        numGridCols = ncol(habMatrix))
                                    }))
  
  ## Rescale AC coordinates for plots
  groupACs.rescaled <- scaleCoordsToHabitatGrid(coordsData = groupACs,
                                                coordsHabitatGridCenter = coordinates(habRaster)[isHabitat, ],
                                                scaleToGrid = F)$coordsDataScaled
  groupACs.rescaled <- SpatialPointsDataFrame(coords = groupACs.rescaled,
                                              data = as.data.frame(groupACs.rescaled))
  
  
  
  ## ------     2.3. SIMULATE INDIVIDUAL STATUS ------
  ## Assign individuals to groups
  idGroupID <- rep(1:length(groupSize), groupSize)
  idGroupType <- rep(groupType, groupSize)
  
  ## Assign individual activity centers
  idACs <- groupACs[idGroupID, ]
  
  idACs.rescaled <- groupACs.rescaled[idGroupID, ]
  
  ## Sample individual status
  idStatus <- unlist(lapply( X = 1:sim.N,
                             FUN = function(x)sample(x = 1:dim(prop)[2],
                                                     size = 1,
                                                     replace = T,
                                                     prob = prop[idGroupType[x], ])))
  idSex <- idSigma <- idP0 <- NULL
  for(i in 1:sim.N){
    ## Assign individual sex
    idSex[i] <- sex[idGroupType[i],idStatus[i]]
    ## Assign individual sigma
    idSigma[i] <- sigma[idGroupType[i],idStatus[i]]
    idSigma[i] <- idSigma[i]/habitatRes
    ## Assign individual p0
    idP0[i] <- p0[idGroupType[i],idStatus[i]]
  }#i
  
  ## Check population composition
  exp.compo <- propGroups * prop
  
  sim.compo <- table(idGroupType,idStatus)
  sim.compo <- round(sim.compo/sum(sim.compo),2)
  
  ## Check population sex ratio
  exp.SR <- sum(unlist(lapply(1:nrow(exp.compo),
                              function(x){
                                y <- which(sex[x, ] %in% c("m"))
                                return(sum(exp.compo[x,y]))
                              })), na.rm = T)
  
  sim.SR <- sum(idSex == "m")/sim.N
  
  
  
  ## ------     2.4. SIMULATE ADAPTIVE SAMPLING ------
  ## Calculate individual space-use
  UD <- do.call(rbind, lapply(1:sim.N,
                              function(i){
                                d2 <- (scaledDetCoords[ ,1] - idACs[i,1])^2 +
                                  (scaledDetCoords[ ,2] - idACs[i,2])^2
                                ud <- exp(-d2/(2*idSigma[i]*idSigma[i]))
                                return(ud)
                              }))
  
  ## Calculate collective space use (i.e. density)
  density <- colSums(UD)
  adSampInt <- scale(density)
  
  
  
  ## ------     2.5. GENERATE SCR DATA ------
  y <- do.call(rbind, lapply(1:sim.N,
                             function(i){
                               p0i <- ilogit(logit(idP0[i]) + betaP0*adSampInt)
                               pi <- p0i*UD[i, ]
                               yi <- rbinom(n = length(pi), size = 1, prob = pi)
                               return(yi)
                             }))
  numDetTot <- sum(y)
  
  
  
  ## ------     2.6. PLOT CHECK ------  
  if(sim %in% changeSpecies){
    my.col <- matrix(terrain.colors(11)[c(1,2,3,5,6,7,9,10,11)],
                     nrow = 3, ncol = 3, byrow = T)
    rgb.col <- col2rgb(my.col)/255
    my.colors <- matrix(rgb(rgb.col[1,], rgb.col[2,], rgb.col[3,], 0.3),
                     nrow = 3, ncol = 3)
  ## Grouping patterns
  plot(habPoly, col = "gray80", main = "", border = F)
  plot(studyArea, col = "gray60", add = T, border = F)
  for(i in 1:sim.N){
    idLoc <- jitter(coordinates(idACs.rescaled)[i, ], amount = 100)
    points(x = idLoc[1], y = idLoc[2],
           pch = 21, cex = 1,
           bg = my.col[idGroupType[i]+idStatus[i]],
           col = my.colors[idGroupType[i],idStatus[i]])
  }#i

  levels <- unique(paste0(idGroupType,".",idStatus))
  legend( x = 592000, y = 5656000, bty = "n",
          legend = levels,
          title = "id status",
          y.intersp = 1.2,
          x.intersp = 2,
          pch = 21,
          col = as.vector(t(my.colors))[c(1,2,3,4,7)],
          pt.bg = as.vector(t(my.colors))[c(1,2,3,4,7)],
          pt.cex = 1)


  ## Adaptive sampling
  p0.test <- unique(idP0)
  test <- matrix(NA, length(p0.test), length(adSampInt))
  for(i in 1:length(p0.test)){
    test[i, ] <-  ilogit(logit(p0.test[i]) + betaP0*adSampInt)
    if(i == 1){
      plot(adSampInt, test[i, ], col = i, ylim = c(0,1))
    }else{points(adSampInt, test[i, ], col = i)}
  }

  ## Use shades of colors
  plot(habPoly, col = "gray80", main = "Adaptive sampling intensity")
  plot(studyArea, col = "gray60", add = T)
  for(d in 1:dim(detectors$detector.sp)[1]){
    plot( detectors$detector.sp[d, ], add = TRUE, pch = 21, cex = 1,
          bg = adjustcolor("firebrick3", alpha.f = adSampInt[d]),
          col = adjustcolor("firebrick3", alpha.f = adSampInt[d]))
  }#d

  ## Use jittered points for detections
  plot(habPoly, col = "gray80", main = "Detections")
  plot(studyArea, col = "gray60", add = T)
  for(i in 1:dim(y)[1]){
    dets <- which(y[i, ] >= 1)
    if(length(dets) > 0){
      detLocs <- jitter(coordinates(detectors$detector.sp)[dets, ], amount = 100)
      points(detLocs,
             pch = 21, cex = 0.5,
             bg = adjustcolor(hcl.colors(dim(y)[1])[i], alpha.f = 0.5),
             col = adjustcolor(hcl.colors(dim(y)[1])[i], alpha.f = 0.5))
    }
  }#d
  }#if
  
  
  ## Number of ids detected per detector
  numIdsPerDet <- mean(colSums(y))
  
  ## Number of detections per individual
  numDetsPerId <- mean(rowSums(y))
  
  ## % of individuals detected
  propDet <- mean(rowSums(y) > 0)
  
  
  
  ## ------     2.7. AUGMENT & SPARSE DATA ------
  ## Augment y by a factor 2
  y <- rbind(y, matrix(0, dim(y)[1], dim(y)[2]))
  ## Transform into sparse matrix
  ySparse <- getSparseY(y)
  
  ## Keep sex for detected individuals only
  sex.data <- c(idSex, rep(NA,dim(y)[1]))
  sex.data[rowSums(y) == 0] <- NA
  sex.data[sex.data == "f"] <- 0
  sex.data[sex.data == "m"] <- 1
  sex.data <- as.numeric(sex.data)
  
  
  
  ## ---------------------------------------------------------------------------
  ## ------   3. SCR NIMBLE MODEL ------
  ## ------     3.1. DEFINE NIMBLE MODEL ------
  nimModel <- nimbleCode({
    
    ##-----------------------------##
    ##------ SPATIAL PROCESS ------##
    ##-----------------------------##
    for(i in 1:n.individuals){
      s[i,1:2] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:n.cells,1:2],
        upperCoords = upperHabCoords[1:n.cells,1:2],
        logIntensities = logIntensities[1:n.cells],
        logSumIntensity = logSumIntensity,
        habitatGrid = habitatGrid[1:y.max, 1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)
    }#i
    
    
    ##-------------------------------##
    ##----- DEMOGRAPHIC PROCESS -----##
    ##-------------------------------##
    sexRatio ~ dunif(0.2,0.8)
    psi ~ dunif(0.1,0.9)
    for(i in 1:n.individuals){
      sex[i] ~ dbern(sexRatio)
      sex1[i] <- sex[i]+1
      z[i] ~ dbern(psi)
    }#i
    N.females <- sum(sex[1:n.individuals]*z[1:n.individuals])
    N.males <- N - N.females
    N <- sum(z[1:n.individuals])
    
    
    ##-----------------------------##
    ##----- DETECTION PROCESS -----##
    ##-----------------------------##
    sigma[1] ~ dunif(0,10)
    sigma[2] ~ dunif(0,10)
    
    p0[1] ~ dunif(0,0.7)
    p0[2] ~ dunif(0,0.7)
    
    for(i in 1:n.individuals){
      y[i,1:maxNumDets] ~ dbinomLocal_normal(
        detNums = detNums[i],
        detIndices = detIndices[i,1:maxNumDets],
        size = trials[1:n.detectors],
        p0 = p0[sex1[i]],
        sigma = sigma[sex1[i]],
        s = s[i,1:2],
        trapCoords = trapCoords[1:n.detectors,1:2],
        localTrapsIndices = localTrapsIndices[1:numResizeCells, 1:numLocalDetsMax],
        localTrapsNum = localTrapsNum[1:numResizeCells],
        resizeFactor = resizeFactor,
        habitatGrid = resizeHabGrid[1:y.resize, 1:x.resize],
        indicator = z[i])
    }#i
  })
  
  
  
  ## ------     3.2. DEFINE NIMBLE CONSTANTS ------
  nimConstants <- list(
    x.max = dim(habGrid)[2],
    y.max = dim(habGrid)[1],
    x.resize = dim(localDetectors$habitatGrid)[2],
    y.resize = dim(localDetectors$habitatGrid)[1],
    n.individuals = dim(ySparse$y)[1],
    n.cells = numHabCells,
    n.detectors = numDetectors,
    maxNumDets = ySparse$maxDetNums,
    resizeFactor = localDetectors$resizeFactor,
    numLocalDetsMax = dim(localDetectors$localIndices)[2],
    numResizeCells = dim(localDetectors$localIndices)[1])
  
  
  
  ## ------     3.3. DEFINE NIMBLE DATA ------
  nimData <- list(
    y = ySparse$y[,,1],
    sex = sex.data,
    sex1 = sex.data + 1,
    detIndices = ySparse$detIndices[,,1],
    detNums = ySparse$detNums[,1],
    trials = rep(1,numDetectors),
    lowerHabCoords = scaledLowerHabCoords,
    upperHabCoords = scaledUpperHabCoords,
    trapCoords = scaledDetCoords,
    localTrapsIndices = localDetectors$localIndices,
    localTrapsNum = localDetectors$numLocalIndices,
    logIntensities = rep(0, numHabCells),
    logSumIntensity = log(numHabCells),
    habitatGrid = habGrid,
    resizeHabGrid = localDetectors$habitatGrid)
  
  
  
  ## ------     3.4. DEFINE NIMBLE INITS -------
  s.init <- matrix(NA, nimConstants$n.individuals, 2)
  for(i in 1:nimConstants$n.individuals){
    if(i <= nimConstants$n.individuals/2){
      s.init[i, ] <- idACs[i, ]
    }else{
      s.init[i, ] <- rbernppAC( n = 1,
                                lowerCoords = nimData$lowerHabCoords,
                                upperCoords = nimData$upperHabCoords,
                                logIntensities = nimData$logIntensities,
                                logSumIntensity = nimData$logSumIntensity,
                                habitatGrid = nimData$habitatGrid,
                                numGridRows = nrow(nimData$habitatGrid),
                                numGridCols = ncol(nimData$habitatGrid))
    }
  }#i
  
  sex.init <- rbinom(n = nimConstants$n.individuals, 1, prob = 0.5)
  sex.init[!is.na(nimData$sex)] <- NA
  sex1.init <- sex.init + 1
  
  nimInits <- list(
    sexRatio = 0.5,
    sex = sex.init,
    sex1 = sex1.init,
    s = s.init,
    z = c(rep(1,nimConstants$n.individuals/2),
          rep(0,nimConstants$n.individuals/2)),
    sigma = c(3,3),
    p0 = c(0.4,0.4),
    psi = 0.5)
  
  
  
  ## ------     3.5. DEFINE NIMBLE PARAMETERS ------
  nimParams <- c( "N", "sigma", "p0", "psi", "sexRatio","s", "z")
  
  
  
  ## ------     3.6. SAVE THE INPUT ------
  ## Store useful info
  parm.df$sim.D[sim] <- sim.D
  parm.df$exp.N[sim] <- exp.N
  parm.df$sim.N[sim] <- sim.N
  
  parm.df$exp.SR[sim] <- exp.SR
  parm.df$sim.SR[sim] <- sim.SR
  
  parm.df$propDet[sim] <- propDet
  parm.df$numDetsPerId[sim] <- numDetsPerId
  parm.df$numIdsPerDet[sim] <- numIdsPerDet
  parm.df$numDetTot[sim] <- numDetTot
  
  paramspace <- parm.df[sim, ]
  
  thisFileName <- paste( "NimbleInFile_Set", parm.df$set_ID[sim],
                         "_Rep", parm.df$rep_ID[sim],
                         ".RData", sep = "")
  save( nimModel,
        nimData,
        nimConstants,
        nimInits,
        nimParams,
        paramspace,
        exp.compo,
        sim.compo,
        file = file.path(path.IN, thisFileName))
  print(sim)
}#sim


## -----------------------------------------------------------------------------
## ------   4. RUN THE MODEL ------
inFiles <- list.files(path.IN)
for(m in 1:length(inFiles)){
  load(file.path(path.IN, inFiles[m]))
  
  Rmodel <- nimbleModel( code = nimModel,
                         constants = nimConstants,
                         data = nimData,
                         inits = nimInits,
                         check = F,
                         calculate = F)
  cmodel <- compileNimble(Rmodel)
  cmodel$calculate()
  MCMCconf <- configureMCMC(model = Rmodel,
                            monitors = nimParams,
                            control = list(reflective = TRUE),
                            thin = 1)
  MCMC <- buildMCMC(MCMCconf)
  cMCMC <- compileNimble(MCMC, project = Rmodel, resetFunctions = TRUE)
  MCMCRuntime <- system.time(nimOutput <- runMCMC(
    mcmc = cMCMC,
    nburnin = 0,
    niter = 1000,
    nchains = 1,
    samplesAsCodaMCMC = TRUE))
  MCMCRuntime

  save(nimOutput, MCMCRuntime,
       file = file.path(path.OUT, paste0("OutFor", inFiles[m])))
}#m



## -----------------------------------------------------------------------------
## ------   5. PLOT RESULTS -----
load(file.path(WD, "parm.df.RData"))
load(file.path(WD, "results.RData"))

myCols <- terrain.colors(8)
rgb.col <- col2rgb(myCols)/255
colors <- rgb(rgb.col[1,], rgb.col[2,], rgb.col[3,], 0.3)


## PLOTS N
{ 
  graphics.off()
  pdf(file = file.path(WD, "results_N.pdf"),
      width = 15, height = 8)
  
  species <- unique(parm.df$species)
  scenarios <- 1:8
  valNoise <- seq(-0.35,0.35, length.out = 8)
  par(mfrow = c(1,2))
  
  ## RELATIVE BIAS
  for(se in c("high", "low")){
    ylim <- 0.5
    plot(1,1, xlim = c(0.5,3.5), ylim = c(-ylim,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "RB(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(-ylim,ylim,length.out = 5),
         labels = seq(-ylim,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(h = 0, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    
    p <- 1
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc  &
                                 myResults[[p]]$searchEffort == se, ]
        try(plot.violins2( dat.list = list(na.omit((temp$mean - temp$sim.N)/temp$sim.N)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  
  ## COEFFICIENT OF VARIATION
  for(se in c("high", "low")){
    ylim <- 0.15
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "CV(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se,  ]
        try(plot.violins2( dat.list = list(na.omit(temp$CV)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  ## COVERAGE
  for(se in c("high", "low")){
    ylim <- 1
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "95%CI coverage", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.2))
    abline(h = 0.95, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se, ]
        
        coverage <- mean(temp$lci <= temp$sim.N & temp$uci >= temp$sim.N, na.rm =T)
        try(points( y = coverage,
                    x = sp + valNoise[sc],
                    pch = 21, cex = 2,
                    bg = colors[sc],
                    col = myCols[sc]), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 0.5, y = 0.35,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  graphics.off()
}

## PLOTS Sex Ratio
{ 
  graphics.off()
  pdf(file = file.path(WD, "results_SR.pdf"), width = 15, height = 8)
  
  species <- unique(parm.df$species)
  scenarios <- 1:8
  valNoise <- seq(-0.35,0.35, length.out = 8)
  par(mfrow = c(1,2))
  
  ## RELATIVE BIAS
  for(se in c("high", "low")){
    ylim <- 0.5
    plot(1,1, xlim = c(0.5,3.5), ylim = c(-ylim,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "RB(Sex Ratio)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(-ylim,ylim,length.out = 5),
         labels = seq(-ylim,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.1))
    abline(h = 0, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    p <- 3
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc  &
                                 myResults[[p]]$searchEffort == se, ]
        try(plot.violins2( dat.list = list(na.omit((temp$mean - temp$sim.SR)/temp$sim.SR)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 2.8, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  
  ## COEFFICIENT OF VARIATION
  for(se in c("high", "low")){
    ylim <- 0.15
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "CV(N)", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.1))
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se,  ]
        try(plot.violins2( dat.list = list(na.omit(temp$CV)),
                           x = sp + valNoise[sc],
                           at = sp + valNoise[sc],
                           add = T,
                           col = myCols[sc],
                           violin.width = 0.05,
                           alpha = 0.9,
                           border.col = myCols[sc],
                           scale.width = F), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 3, y = ylim,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  ## COVERAGE
  for(se in c("high", "low")){
    ylim <- 1
    plot(1,1, xlim = c(0.5,3.5), ylim = c(0,ylim),
         type = "n", xaxt = "n", main = se,
         ylab = "95%CI coverage", xlab = "", axes = F)
    axis(1, at = 1:3, labels = c("Chamois","Red Deer", "Wild boar"))
    axis(2,
         at = seq(0,ylim,length.out = 5),
         labels = seq(0,ylim,length.out = 5))
    polygon(x = c(1.5,1.5,2.5,2.5),
            y = c(-ylim,ylim,ylim,-ylim),
            border = F, col = adjustcolor("gray80", alpha.f = 0.1))
    abline(h = 0.95, lty = 2, lwd = 2)
    abline(v = 1.5, lty = 2, lwd = 1)
    abline(v = 2.5, lty = 2, lwd = 1)
    for(sp in 1:length(species)){
      for(sc in scenarios){
        temp <- myResults[[p]][myResults[[p]]$species == species[sp] &
                                 myResults[[p]]$scenario == sc &
                                 myResults[[p]]$searchEffort == se, ]
        
        coverage <- mean(temp$lci <= temp$sim.SR & temp$uci >= temp$sim.SR, na.rm = T)
        try(points( y = coverage,
                    x = sp + valNoise[sc],
                    pch = 21, cex = 2,
                    bg = colors[sc],
                    col = myCols[sc]), silent = F)
      }#sc
    }#sp
  }#se
  legend(x = 0.5, y = 0.35,
         legend = c("m0", "G", "H", "GH", "A", "GA", "HA", "GHA"),
         fill = colors)
  
  graphics.off()
}


## -----------------------------------------------------------------------------