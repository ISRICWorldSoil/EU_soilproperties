## Fit models for As, Cd, Cu, Pb, Zn, SOC, ph CaCl2, pH H2O, clay, sand and silt
## Tom.Hengl@isric.org & Gerard.Heuvelink@wur.nl
## predictions only: "Arable land", "Grassland"

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "R.utils", "dplyr", "parallel", "dismo", "snowfall", "lattice", "ranger", "mda", "psych", "stringr", "caret", "plotKML", "maptools", "maps", "stringr", "R.utils", "grDevices")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/Integrator")
load(".RData")
library(plyr)
library(dplyr)
library(stringr)
library(sp)
library(rgdal)
#library(e1071)
#library(randomForest)
library(devtools)
#devtools::install_github('dmlc/xgboost')
library(xgboost)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger", ref="forest_memory") ## version to deal with Memory problems
library(ranger)
library(raster)
library(caret)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
library(snowfall)
library(utils)
library(plotKML)
library(R.utils)
library(GSIF)
library(parallel)
library(doParallel)

plotKML.env(convert="convert", show.env=FALSE)
gdalwarp = "gdalwarp"
gdalbuildvrt = "gdalbuildvrt"
system("gdal-config --version")

## load from github: https://github.com/ISRICWorldSoil/SoilGrids250m/tree/master/grids/models
source_https <- function(url, ...) {
  require(RCurl)
  cat(getURL(url, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), file = basename(url))
  source(basename(url))
}
source_https("https://raw.githubusercontent.com/ISRICWorldSoil/SoilGrids250m/master/grids/models/saveRDS_functions.R")
source_https("https://raw.githubusercontent.com/ISRICWorldSoil/SoilGrids250m/master/grids/models/wrapper.predict_cs.R")

## training points:
SPROPS.GEMAS = readRDS("/data/soilstorage/SoilData/GEMAS/pnts_GEMAS.rds") ## GEMAS
load("./points/SPROPS.LUCAS_LC.rda") ## LUCAS
pnts = dplyr::bind_rows(list(SPROPS.GEMAS,SPROPS.LUCAS_LC))
pnts = pnts[!is.na(pnts$LONWGS84)&!is.na(pnts$LATWGS84),]
str(pnts)
coordinates(pnts) = ~ LONWGS84 + LATWGS84
proj4string(pnts) = "+proj=longlat +datum=WGS84"

## Fix missing values for LC1:
cor.r = raster("/mnt/cartman/CORINE/g100_clc12_V18_5.tif")
proj4string(cor.r) = "+init=epsg:3035"
cor.pnts = raster::extract(cor.r, spTransform(pnts, CRS("+init=epsg:3035")))
## land cover codes 
leg.LC1 = read.csv("./points/legend_LC1.csv", na.strings = c("","#NA","#N/A","NA"))
leg.cor = read.csv("./points/CLC_2012_legend.csv", na.strings = c("","#NA","#N/A","NA"))
LC1.cor = join(leg.cor, leg.LC1)
x = join(data.frame(GRID_CODE=cor.pnts), LC1.cor, match="first")$AGGR_LC1
summary(leg.LC1$AGGR_LC1)
pnts$CLC = join(pnts@data, leg.LC1)$AGGR_LC1
pnts$CLC = factor(ifelse(is.na(pnts$CLC), paste(x), paste(pnts$CLC)))
summary(pnts$CLC)
#Arable land   Grassland          NA      Nature       Urban 
#391       10543        5161         317        7330         131 
#Wetlands 
#48 

## spatia overlay (50 mins):
covs250m.lst = list.files(path="/data/stacked250m", glob2rx("*.tif$"), full.names=TRUE)
covs250m.lst <- covs250m.lst[-unlist(sapply(c("QUAUEA3","BICUSG5.tif","LCEE10","B08CHE3","B09CHE3","S01ESA4","S02ESA4","S11ESA4","S12ESA4","CSCMCF5"), function(x){grep(x, covs250m.lst)}))]
extract.tif = function(x, y){
  r = raster(x)
  if(!is.na(proj4string(r))){
    y = spTransform(y, proj4string(r))
  }  
  out = raster::extract(y=y, x=r)
  return(out)
}
## overlay in parallel TAKES 20 MINS:
sfInit(parallel=TRUE, cpus=56)
sfExport("pnts", "covs250m.lst", "extract.tif")
sfLibrary(rgdal)
sfLibrary(raster)
ov <- data.frame(sfClusterApplyLB(covs250m.lst, function(i){try( extract.tif(i, pnts) )}))
sfStop()
names(ov) = basename(covs250m.lst)

ovA <- cbind(as.data.frame(pnts), ov)
## 23921 obs. of  219 variables

## Data inspection (final checks)
hist(log1p(ovA$CECSUM))
hist(ovA$CLYPPT) ## probably many missing values are 0?
hist(log1p(ovA$ORCDRC))
hist(ovA$PHIHOX)
hist(ovA$PHICAL) ## also many missing values
summary(ovA$DEPTH)
summary(ovA$ORCDRC > 80) ## 5% of profiles >8% ORC

## Convert land cover to indicators:
clc.mat = data.frame(model.matrix(~CLC-1, ovA))
clc.mat$CLCNA = NULL
ovA = cbind(ovA, clc.mat)
ovA$LOC_ID <- as.factor(paste("ID", ovA$LONWGS84, ovA$LATWGS84, sep="_"))

write.csv(ovA, file="ov.SPROPS_EU_SoilGrids250m.csv")
unlink("ov.SPROPS_EU_SoilGrids250m.csv.gz")
gzip("ov.SPROPS_EU_SoilGrids250m.csv")
saveRDS(ovA, "ov.SPROPS_EU_SoilGrids250m.rds")

## ------------- MODEL FITTING -----------

## use log-transformed values to improve prediction of lower values
for(i in c("ORCDRC","As", "Cd", "Cu", "Pb", "Zn")){   ovA[,paste0("log_",i)] = log1p(ovA[,i]) }
t.vars <- c("ORCDRC", "PHIHOX", "PHICAL", "SNDPPT", "SLTPPT", "CLYPPT", "As", "Cd", "Cu", "Pb", "Zn", "log_As", "log_Cd", "log_Cu", "log_Pb", "log_Zn", "log_ORCDRC")
col.m <- lapply(ovA[,t.vars], quantile, probs=c(0.01,0.5,0.99), na.rm=TRUE)
col.m <- as.data.frame(col.m)
write.csv(col.m, "vars_quantiles.csv")

z.min <- as.list(c(0,20,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
names(z.min) = t.vars
z.max <- as.list(c(700,950,950,100,100,100,150,10,200,200,500,800,800,800,800,800,800))
names(z.max) = t.vars
## FIT MODELS:
pr.lst <- basename(covs250m.lst)
## remove some predictors that might lead to artifacts (buffer maps and land cover):
formulaString.lst = lapply(t.vars, function(x){as.formula(paste(x, ' ~ ', paste(pr.lst, collapse="+"), "+", paste(names(clc.mat), collapse="+")))})
#all.vars(formulaString.lst[[1]])
save.image()

## cleanup:
# unlink(list.files("./models/",pattern=glob2rx("^mrf.*.rds"), full.names = TRUE))
# unlink(list.files("./models/",pattern=glob2rx("^t.mrf.*.rds"), full.names = TRUE))
# unlink(list.files("./models/",pattern=glob2rx("^Xgb.*"), full.names = TRUE))
# unlink(list.files("./models/",pattern=glob2rx("^mgb.*.rds"), full.names = TRUE))
# unlink(list.files("./models/",pattern=glob2rx("^mrf.*.rds"), full.names = TRUE))
# unlink(list.files("./models/",pattern=glob2rx("^mrfX_*_*.rds"), full.names = TRUE))
# unlink(list.files("./models/",pattern=glob2rx("^RF_fit_*.csv.gz"), full.names = TRUE))

## sub-sample to speed up model fitting:
## TAKES >2 hrs to fit all models
Nsub <- 2.5e3 
## Caret training settings (reduce number of combinations to speed up):
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5), nrounds = c(50,100,150), max_depth = 2:4, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
rf.tuneGrid <- expand.grid(mtry = seq(5,50,by=5))
## Initiate cluster
cl <- makeCluster(56)
doParallel::registerDoParallel(cl)
## Takes 1 hour to fit all models:
for(j in 1:length(t.vars)){
  out.file = paste0("./models/", t.vars[j],"_resultsFit.txt")
  if(!file.exists(out.file)){
    cat("Results of model fitting 'randomForest / XGBoost':\n\n", file=out.file)
    cat("\n", file=out.file, append=TRUE)
    cat(paste("Variable:", all.vars(formulaString.lst[[j]])[1]), file=out.file, append=TRUE)
    cat("\n", file=out.file, append=TRUE)
    out.rf <- paste0("./models/mrf.",t.vars[j],".rds")
    if(!file.exists(out.rf)){
      LOC_ID <- ovA$LOC_ID
      dfs <- ovA[,all.vars(formulaString.lst[[j]])]
      sel <- complete.cases(dfs)
      dfs <- dfs[sel,]
      ## optimize mtry parameter:
      if(!file.exists(gsub("./models/mrf","./models/t.mrf",out.rf))){
        t.mrfX <- caret::train(formulaString.lst[[j]], data=dfs[sample.int(nrow(dfs), Nsub),], method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
        saveRDS.gz(t.mrfX, file=gsub("./models/mrf", "./models/t.mrf", out.rf))
      } else {
        t.mrfX <- readRDS.gz(gsub("./models/mrf","./models/t.mrf",out.rf))
      }
      ## fit RF model using 'ranger' (fully parallelized)
      ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
      mrfX <- ranger(formulaString.lst[[j]], data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)  
      saveRDS.gz(mrfX, file=paste0("./models/mrf.",t.vars[j],".rds"))
      ## Top 15 covariates:
      sink(file=out.file, append=TRUE, type="output")
      print(mrfX)
      cat("\n Variable importance:\n", file=out.file, append=TRUE)
      xl <- as.list(ranger::importance(mrfX))
      print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:35]])))
      ## save fitting success vectors:
      fit.df <- data.frame(LOC_ID=LOC_ID[sel], observed=dfs[,1], predicted=predictions(mrfX))
      unlink(paste0("./models/RF_fit_", t.vars[j], ".csv.gz"))
      write.csv(fit.df, paste0("./models/RF_fit_", t.vars[j], ".csv"))
      gzip(paste0("./models/RF_fit_", t.vars[j], ".csv"))
      if(!file.exists(paste0("./models/mgb.",t.vars[j],".rds"))){
        ## fit XGBoost model (uses all points):
        mgbX <- caret::train(formulaString.lst[[j]], data=dfs, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid) 
        saveRDS.gz(mgbX, file=paste0("./models/mgb.",t.vars[j],".rds"))
        ## save also binary model for prediction purposes:
        xgb.save(mgbX$finalModel, paste0("./models/Xgb.",t.vars[j]))
      } else {
        mgbX <- readRDS.gz(paste0("./models/mgb.",t.vars[j],".rds"))
      }
      importance_matrix <- xgb.importance(mgbX$coefnames, model = mgbX$finalModel)
      cat("\n", file=out.file, append=TRUE)
      print(mgbX)
      cat("\n XGBoost variable importance:\n", file=out.file, append=TRUE)
      print(importance_matrix[1:25,])
      cat("--------------------------------------\n", file=out.file, append=TRUE)
      sink()
    }
  }
}
rm(mrfX); rm(mgbX)
stopCluster(cl); closeAllConnections()
save.image()
#x = list.files(pattern=glob2rx(".*.rds$"), all.files = TRUE, full.names=TRUE)
#file.rename(from=x, to=gsub("./.t.","./models/t.",x))

## ------------- PREDICTIONS -----------

system(paste0('gdalwarp /mnt/cartman/CORINE/g100_clc12_V18_5.tif ./stacked1km/g100_clc12_V18_5.tif -tr 1000 1000 -r \"near\" -te 857000 855000 7407000 5438000'))

## resample to 1 km resolution:
sfInit(parallel=TRUE, type="SOCK", cpus=54)
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("gdalwarp", "covs250m.lst", "cor.r")
out <- sfClusterApplyLB(covs250m.lst, function(x){ if(!file.exists(paste0('./stacked1km/', gsub(".tif", "_1km.tif", basename(x))))){ system(paste0(gdalwarp, ' ', x, ' ./stacked1km/', gsub(".tif", "_1km.tif", basename(x)), ' -co \"COMPRESS=DEFLATE\" -r \"average\" -t_srs \"', proj4string(cor.r), '\"  -tr 1000 1000 -te 857000 855000 7407000 5438000')) }})
sfStop()
## test:
x = stack(list.files("./stacked1km/", glob2rx("*.tif$"), full.names=TRUE))

## Tiling system:
obj <- GDALinfo("./stacked1km/g100_clc12_V18_5.tif")
tile.lst <- getSpatialTiles(obj, block.x=200000, return.SpatialPolygons=TRUE)
proj4string(tile.lst)
tile.tbl <- getSpatialTiles(obj, block.x=200000, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
## ?? tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
unlink("tiles200km.shp")
writeOGR(tile.pol, "tiles200km.shp", "tiles200km", "ESRI Shapefile")
plot(spTransform(pnts, CRS(proj4string(tile.lst))), pch="+")
lines(tile.lst, col="red")
## Overlay tiles and admin units (fully parallelized):
system(paste('gdal_translate ./stacked1km/g100_clc12_V18_5.tif g100_clc12_V18_5.sdat -ot \"Byte\" -of \"SAGA\"'))
system(paste0('saga_cmd -c=56 shapes_grid 2 -GRIDS=\"g100_clc12_V18_5.sgrd\" -POLYGONS=\"tiles200km.shp\" -PARALLELIZED=1 -RESULT=\"ov_tiles200km.shp\"'))
ov_tiles = readOGR("ov_tiles200km.shp", "ov_tiles200km")
str(ov_tiles@data)
summary(sel.t <- !is.na(ov_tiles$g100_clc12_))
## 293 tiles with values
ov_tiles = ov_tiles[sel.t,]
plot(ov_tiles)
unlink("ov_tiles200km.shp")
writeOGR(ov_tiles, "ov_tiles200km.shp", "ov_tiles200km", "ESRI Shapefile")

## Prepare covariates as tiles:
t.sel = as.character(ov_tiles$ID)
new.dirs <- paste0("./tiled/T", t.sel)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

mean.val = lapply(ovA[,all.vars(formulaString.lst[[1]])[-1]], function(x){quantile(x, probs=.5, na.rm=TRUE)})
names(mean.val) = all.vars(formulaString.lst[[1]])[-1]
## Function to make predictions locations
make_newdata <- function(i, tile.tbl, in.covs, out.path="./tiled", mask.l="./stacked1km/g100_clc12_V18_5.tif", LC1.cor, mean.val){
  out.rds <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/T", tile.tbl[i,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(fname=mask.l, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    m$band1 = ifelse(m$band1>39, NA, m$band1)
    sel.p = !is.na(m$band1)
    if(sum(sel.p)>1){
      m$CLC = join(data.frame(GRID_CODE=m$band1), LC1.cor, match="first")$AGGR_LC1
      m = as(m["CLC"], "SpatialPixelsDataFrame")
      clc.t = data.frame(model.matrix(~CLC-1, m@data))
      #clc.t$CLCNA = NULL
      m@data = cbind(m@data, clc.t)
      for(j in 1:length(in.covs)){
        cov.name = gsub("_1km", "", basename(in.covs[j]))
        m@data[,cov.name] <- round(readGDAL(in.covs[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index])
      }
      ## Fill-in missing values (if necessary):
      sel.mis = sapply(m@data, function(x){sum(is.na(x))>0})
      if(sum(sel.mis)>0){
        x = which(sel.mis)
        xn = attr(x, "names")
        for(j in 1:length(x)){
          if(!is.factor(m@data[,x[j]])){
            if(length(grep(pattern="USG5", names(m)[x[j]]))>0 | length(grep(pattern="GLC5", names(m)[x[j]]))>0){ 
              repn = 0 
            } else {
              if(sum(is.na(m@data[,x[j]]))<length(m)){
                repn = quantile(m@data[,x[j]], probs=.5, na.rm=TRUE)
              } else {
                repn = mean.val[[xn[j]]]
              }
            }
            m@data[,x[j]] = ifelse(is.na(m@data[,x[j]]), repn, m@data[,x[j]])
          }
        }
      }
      ## limit to complete data sets only
      m = m[complete.cases(m@data),]
      saveRDS(m, out.rds)
      gc(); gc()
    }
  }
}
## test it:
#make_newdata(i=249, tile.tbl=tile.tbl, in.covs=paste0("./stacked1km/", gsub(".tif", "_1km.tif", pr.lst)), LC1.cor=LC1.cor, mean.val=mean.val)
#x = readRDS("./tiled/T249/T249.rds")
#which(!all.vars(formulaString.lst[[1]])[-1] %in% names(x))

#del.lst <- list.files(path="./tiled", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)
#unlink(del.lst)
## TAKES 15 mins
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus=54)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(plyr)
snowfall::sfExport("make_newdata", "t.sel", "tile.tbl", "pr.lst", "LC1.cor", "mean.val")
out <- snowfall::sfClusterApplyLB(as.numeric(t.sel), function(i){ make_newdata(i, tile.tbl=tile.tbl, in.covs=paste0("./stacked1km/", gsub(".tif", "_1km.tif", pr.lst)), LC1.cor=LC1.cor, mean.val=mean.val) })
snowfall::sfStop()

## Predict per tile:
pr.dirs <- basename(dirname(list.files(path="./tiled", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
str(pr.dirs)

#save(pr.dirs, file="pr.dirs.rda")
## 292 dirs
sel.vars <- c("ORCDRC", "PHIHOX", "PHICAL", "SNDPPT", "SLTPPT", "CLYPPT", "log_As", "log_Cd", "log_Cu", "log_Zn")
type.lst <- c("Int16", "Byte", "Byte", "Byte", "Byte", "Byte", "Int16", "Int16", "Int16", "Int16", "Int16")
names(type.lst) = sel.vars
mvFlag.lst <- c(-32768, 255, 255, 255, 255, 255, -32768, -32768, -32768, -32768, -32768)
names(mvFlag.lst) = sel.vars

#del.lst <- list.files(path="./tiled", pattern=glob2rx("*.tif$"), recursive = TRUE, full.names = TRUE)
#unlink(del.lst)
#del.lst <- list.files(path="./tiled", pattern=glob2rx("log_*.tif$"), recursive = TRUE, full.names = TRUE)
#unlink(del.lst)

## Run per property (TAKES ABOUT 30 mins in total)
library(ranger)
library(xgboost)
library(tools)
library(parallel)
library(doParallel)
library(rgdal)
library(plyr)
## Run per property (RF = 20 tiles per minute)
for(j in sel.vars){
  try( detach("package:snowfall", unload=TRUE), silent=TRUE)
  try( detach("package:snow", unload=TRUE), silent=TRUE)
  if(j %in% c("PHIHOX","PHICAL")){ multiplier = 10 }
  if(j %in% c("log_As","log_Cd","log_Cu","log_Zn")){ multiplier = 100 }
  if(j %in% c("ORCDRC","SNDPPT","SLTPPT","CLYPPT")){ multiplier = 1 }
  ## Random forest predictions:
  gm = readRDS.gz(paste0("./models/mrf.", j,".rds"))
  gm1.w = 1/gm$prediction.error
  ## Estimate amount of RAM needed per core
  cpus = unclass(round((500-50)/(3.5*(object.size(gm)/1e9))))
  require(parallel)
  cl <- makeCluster(ifelse(cpus>54, 54, cpus), type="FORK")
  x = parLapply(cl, pr.dirs, fun=function(x){ if(any(!file.exists(paste0("./tiled/", x, "/", j, "_M_", x, ".tif")))){ try( split_predict_n(x, gm, in.path="./tiled", out.path="./tiled", split_no=NULL, varn=j, method="ranger", depths=FALSE, multiplier=multiplier, rds.file=paste0("./tiled/", x, "/", x,".rds")) ) } } )
  stopCluster(cl)
  gc(); gc()
  ## XGBoost:
  gm = readRDS.gz(paste0("./models/mgb.", j,".rds"))
  gm2.w = 1/(min(gm$results$RMSE, na.rm=TRUE)^2)
  cpus = unclass(round((500-30)/(3.5*(object.size(gm)/1e9))))
  require(parallel)
  cl <- makeCluster(ifelse(cpus>54, 54, cpus), type="FORK")
  x = parLapply(cl, pr.dirs, fun=function(x){ if(any(!file.exists(paste0("./tiled/", x, "/", j, "_M_", x, ".tif")))){ try( split_predict_n(x, gm, in.path="./tiled", out.path="./tiled", split_no=NULL, varn=j, method="xgboost", depths=FALSE, multiplier=multiplier, rds.file=paste0("./tiled/", x, "/", x,".rds")) ) } } )
  stopCluster(cl)
  gc(); gc()
  ## sum up predictions:
  library(snowfall)
  sfInit(parallel=TRUE, cpus=56)
  sfExport("pr.dirs", "sum_predict_ensemble", "j", "z.min", "z.max", "gm1.w", "gm2.w", "type.lst", "mvFlag.lst")
  sfLibrary(rgdal)
  sfLibrary(plyr)
  x <- sfClusterApplyLB(pr.dirs, fun=function(x){ try( sum_predict_ensemble(x, in.path="./tiled", out.path="./tiled", varn=j, num_splits=NULL, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]], depths=FALSE, rds.file=paste0("./tiled/", x, "/", x,".rds")) ) } )
  sfStop()
  gc(); gc()
}

## Mosaic all maps:
mosaic_tif = function(j, in.path="./tiled", out.path="./output/"){
  out.tif <- paste0(out.path, j, "_topsoil_1km.tif")
  if(!file.exists(out.tif)){
    tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(j, "_M_*.tif$")), full.names=TRUE, recursive=TRUE)
    out.tmp <- tempfile(fileext = ".txt")
    vrt.tmp <- tempfile(fileext = ".vrt")
    cat(tmp.lst, sep="\n", file=out.tmp)
    system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
    system(paste0('gdalwarp ', vrt.tmp, ' ', out.tif, ' -r \"near\" -co \"COMPRESS=DEFLATE\" -te 857000 855000 7407000 5438000'))
    unlink(vrt.tmp)
    unlink(out.tmp)
  }
}

sfInit(parallel=TRUE, cpus=ifelse(length(sel.vars)>56, 56, length(sel.vars)))
sfExport("sel.vars", "mosaic_tif")
out <- sfClusterApplyLB(sel.vars, function(x){try( mosaic_tif(x) )})
sfStop()

## summary statistics:
m = readGDAL("./output/PHIHOX_topsoil_1km.tif")
m@data[,"x"] = readGDAL("./stacked1km/g100_clc12_V18_5.tif")$band1
m$CLC = join(data.frame(GRID_CODE=m$x), LC1.cor, match="first")$AGGR_LC1
sum.pH = plyr::ddply(m@data, .(CLC), summarize, pH_mean=mean(band1, na.rm=TRUE)/10, pH_sd=sqrt(var(band1/10, na.rm=TRUE)))

