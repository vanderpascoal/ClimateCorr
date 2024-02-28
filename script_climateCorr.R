# Algorithm developed to generate indicators for correlating climate anomalies within a 0.5
# degree grid. The generated products include minimum temperature, average temperature,
# median temperature, maximum temperature, standard deviation, temperature variation, and
# the frequency of temperature surpassing the maximum climatological normal.
# DOI: 10.5281/zenodo.10592772
# author: "Vanderlei Matos"
# date: "2024-01-30"

#configure your workbook
rm(list = ls())
setwd("###")
library("raster")
library(sf)
library(parallel)
library(doParallel)
library(tidyverse)
library("raster")
library(geojsonio)
library(compiler)
#Load climatological normal data
nclima01 <- raster("normal/wc2.1_30s_tavg_01.tif")
nclima02 <- raster("normal/wc2.1_30s_tavg_02.tif")
nclima03 <- raster("normal/wc2.1_30s_tavg_03.tif")
nclima04 <- raster("normal/wc2.1_30s_tavg_04.tif")
nclima05 <- raster("normal/wc2.1_30s_tavg_05.tif")
nclima06 <- raster("normal/wc2.1_30s_tavg_06.tif")
nclima07 <- raster("normal/wc2.1_30s_tavg_07.tif")
nclima08 <- raster("normal/wc2.1_30s_tavg_08.tif")
nclima09 <- raster("normal/wc2.1_30s_tavg_09.tif")
nclima10 <- raster("normal/wc2.1_30s_tavg_10.tif")
nclima11 <- raster("normal/wc2.1_30s_tavg_11.tif")
nclima12 <- raster("normal/wc2.1_30s_tavg_12.tif")
#Configure the number of cores used
n.cores <- parallel::detectCores() - 5
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
#check cluster definition
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
#Loading geometry
shp <- sf::read_sf("br_2010.shp")
#Configuring the grid
g = st_make_grid(shp, cellsize = 0.5, square = TRUE)
g<-g[shp]
g<-as_Spatial(g)
#Clip the climatological normal for each grid value
nclima01<-crop(nclima01, g)
nclima02<-crop(nclima02, g)
nclima03<-crop(nclima03, g)
nclima04<-crop(nclima04, g)
nclima05<-crop(nclima05, g)
nclima06<-crop(nclima06, g)
nclima07<-crop(nclima07, g)
nclima08<-crop(nclima08, g)
nclima09<-crop(nclima09, g)
nclima10<-crop(nclima10, g)
nclima11<-crop(nclima11, g)
nclima12<-crop(nclima12, g)
#
temp.grid<-function(g1=g[9], ano = 0){
  files.temperaute <- list.files(path = "../update_obs/ECMWF/diario/AL/", pattern =
                                   "2m_temperature")
  pos <- stringr::str_detect(string = files.temperaute,pattern = as.character(ano))
  files.temperaute <- files.temperaute[pos]
  dt <- list()
  value1 <-list()
  i<-1
  for(i in 1:length(files.temperaute)){
    print(files.temperaute[i])
    p1 <- raster(paste0("../update_obs/ECMWF/diario/AL/",files.temperaute[i]))
    data1 <- crop(x = p1,
                  y = g1)
    bandas_stackeadas <- stack(c(data1)) # Creating the stack from list objects
    temperatura <- overlay(bandas_stackeadas,
                           fun=function(x){(x-273.15)})
    Stack.crop = as.matrix(temperatura)
    #Treating inconsistent values with leaps and bounds
    if(mean(x = Stack.crop) < -60){next()}
    value1 <- append(value1, mean(x = Stack.crop))
    dt <- append(dt, p1@z[[1]])
  }
  tmp <- data.frame(temp = unlist(value1), date = as.Date(unlist(dt)))
  rm(p1, value1, dt)
  tmp$anomes<-paste0(format(tmp$date,"%Y"),format(tmp$date,"%m"))
  ##### Cutting by the normal with respect to average
  nclima01.1<-crop(nclima01, g1)
  nclima01.1<-mean(as.matrix(nclima01.1))
  nclima02.1<-crop(nclima02, g1)
  nclima02.1<-mean(as.matrix(nclima02.1))
  nclima03.1<-crop(nclima03, g1)
  nclima03.1<-mean(as.matrix(nclima03.1))
  nclima04.1<-crop(nclima04, g1)
  nclima04.1<-mean(as.matrix(nclima04.1))
  nclima05.1<-crop(nclima05, g1)
  nclima05.1<-mean(as.matrix(nclima05.1))
  nclima06.1<-crop(nclima06, g1)
  nclima06.1<-mean(as.matrix(nclima06.1))
  nclima07.1<-crop(nclima07, g1)
  nclima07.1<-mean(as.matrix(nclima07.1))
  nclima08.1<-crop(nclima08, g1)
  nclima08.1<-mean(as.matrix(nclima08.1))
  nclima09.1<-crop(nclima09, g1)
  nclima09.1<-mean(as.matrix(nclima09.1))
  nclima10.1<-crop(nclima10, g1)
  nclima10.1<-mean(as.matrix(nclima10.1))
  nclima11.1<-crop(nclima11, g1)
  nclima11.1<-mean(as.matrix(nclima11.1))
  nclima12.1<-crop(nclima12, g1)
  nclima12.1<-mean(as.matrix(nclima12.1))
  tmp$nclima<-NA
  tmp$nclima[which(format(tmp$date,"%m")=="01")] <- nclima01.1
  tmp$nclima[which(format(tmp$date,"%m")=="02")] <- nclima02.1
  tmp$nclima[which(format(tmp$date,"%m")=="03")] <- nclima03.1
  tmp$nclima[which(format(tmp$date,"%m")=="04")] <- nclima04.1
  tmp$nclima[which(format(tmp$date,"%m")=="05")] <- nclima05.1
  tmp$nclima[which(format(tmp$date,"%m")=="06")] <- nclima06.1
  tmp$nclima[which(format(tmp$date,"%m")=="07")] <- nclima07.1
  tmp$nclima[which(format(tmp$date,"%m")=="08")] <- nclima08.1
  tmp$nclima[which(format(tmp$date,"%m")=="09")] <- nclima09.1
  tmp$nclima[which(format(tmp$date,"%m")=="10")] <- nclima10.1
  tmp$nclima[which(format(tmp$date,"%m")=="11")] <- nclima11.1
  tmp$nclima[which(format(tmp$date,"%m")=="12")] <- nclima12.1
  ##### Cutting by the normal with respect to max
  nclima01.2<-crop(nclima01, g1)
  nclima01.2<-max(as.matrix(nclima01.2))
  nclima02.2<-crop(nclima02, g1)
  nclima02.2<-max(as.matrix(nclima02.2))
  nclima03.2<-crop(nclima03, g1)
  nclima03.2<-max(as.matrix(nclima03.2))
  nclima04.2<-crop(nclima04, g1)
  nclima04.2<-max(as.matrix(nclima04.2))
  nclima05.2<-crop(nclima05, g1)
  nclima05.2<-max(as.matrix(nclima05.2))
  nclima06.2<-crop(nclima06, g1)
  nclima06.2<-max(as.matrix(nclima06.2))
  nclima07.2<-crop(nclima07, g1)
  nclima07.2<-max(as.matrix(nclima07.2))
  nclima08.2<-crop(nclima08, g1)
  nclima08.2<-max(as.matrix(nclima08.2))
  nclima09.2<-crop(nclima09, g1)
  nclima09.2<-max(as.matrix(nclima09.2))
  nclima10.2<-crop(nclima10, g1)
  nclima10.2<-max(as.matrix(nclima10.2))
  nclima11.2<-crop(nclima11, g1)
  nclima11.2<-max(as.matrix(nclima11.2))
  nclima12.2<-crop(nclima12, g1)
  nclima12.2<-max(as.matrix(nclima12.2))
  tmp$nclima_max<-NA
  tmp$nclima_max[which(format(tmp$date,"%m")=="01")] <- nclima01.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="02")] <- nclima02.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="03")] <- nclima03.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="04")] <- nclima04.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="05")] <- nclima05.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="06")] <- nclima06.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="07")] <- nclima07.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="08")] <- nclima08.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="09")] <- nclima09.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="10")] <- nclima10.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="11")] <- nclima11.2
  tmp$nclima_max[which(format(tmp$date,"%m")=="12")] <- nclima12.2
  ##### Cutting by the normal with respect to min
  nclima01.3<-crop(nclima01, g1)
  nclima01.3<-min(as.matrix(nclima01.3))
  nclima02.3<-crop(nclima02, g1)
  nclima02.3<-min(as.matrix(nclima02.3))
  nclima03.3<-crop(nclima03, g1)
  nclima03.3<-min(as.matrix(nclima03.3))
  nclima04.3<-crop(nclima04, g1)
  nclima04.3<-min(as.matrix(nclima04.3))
  nclima05.3<-crop(nclima05, g1)
  nclima05.3<-min(as.matrix(nclima05.3))
  nclima06.3<-crop(nclima06, g1)
  nclima06.3<-min(as.matrix(nclima06.3))
  nclima07.3<-crop(nclima07, g1)
  nclima07.3<-min(as.matrix(nclima07.3))
  nclima08.3<-crop(nclima08, g1)
  nclima08.3<-min(as.matrix(nclima08.3))
  nclima09.3<-crop(nclima09, g1)
  nclima09.3<-min(as.matrix(nclima09.3))
  nclima10.3<-crop(nclima10, g1)
  nclima10.3<-min(as.matrix(nclima10.3))
  nclima11.3<-crop(nclima11, g1)
  nclima11.3<-min(as.matrix(nclima11.3))
  nclima12.3<-crop(nclima12, g1)
  nclima12.3<-min(as.matrix(nclima12.3))
  tmp$nclima_min<-NA
  tmp$nclima_min[which(format(tmp$data,"%m")=="01")] <- nclima01.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="02")] <- nclima02.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="03")] <- nclima03.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="04")] <- nclima04.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="05")] <- nclima05.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="06")] <- nclima06.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="07")] <- nclima07.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="08")] <- nclima08.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="09")] <- nclima09.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="10")] <- nclima10.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="11")] <- nclima11.3
  tmp$nclima_min[which(format(tmp$data,"%m")=="12")] <- nclima12.3
  rm(nclima01.1, nclima02.1, nclima03.1, nclima04.1, nclima05.1,
     nclima06.1, nclima07.1, nclima08.1, nclima09.1, nclima10.1,
     nclima11.1, nclima12.1, nclima01.2, nclima02.2, nclima03.2, nclima04.2, nclima05.2,
     nclima06.2, nclima07.2, nclima08.2, nclima09.2, nclima10.2,
     nclima11.2, nclima12.2,
     nclima01.3, nclima02.3, nclima03.3, nclima04.3, nclima05.3,
     nclima06.3, nclima07.3, nclima08.3, nclima09.3, nclima10.3,
     nclima11.3, nclima12.3)
  tmp <-
    tmp %>%
    group_by(anomes) %>%
    mutate(anoma_f = (temp<nclima), anoma_q = (temp>nclima),
           d35 = (temp>=35), anoma_min = (temp < nclima_min),
           anoma_max = (temp > nclima_max), d20 = (temp >= 20)) %>%
    summarise(min=min(temp),
              avg = mean(temp),
              median = median(temp),
              max=max(temp),
              sd = sd(temp),
              var = var(temp),
              nclimax = unique(nclima_max),
              ncli = unique(nclima),
              nclimin = unique(nclima_min),
              anomq=sum(anoma_q==TRUE),
              aonmf=sum(anoma_f==TRUE),
              dia20 = sum(d20==TRUE),
              dia35=sum(d35==TRUE),
              anom_max=sum(anoma_max == TRUE),
              anom_min=sum(anoma_min == TRUE)
    )
  tmp<-reshape2::melt(tmp)
  tmp$variable<-as.character(tmp$variable)
  tmp<-data.frame(variable = paste0(tmp$anomes,"",tmp$variable), value = tmp$value)
  field_name<-tmp$variable
  tmp<-as.data.frame(t(tmp))
  names(tmp)<-field_name
  tmp<-tmp[-1,]
  row.names(tmp)<-row.names(g1)
  tmp<-SpatialPolygonsDataFrame(Sr = g1,
                                data = tmp)
  return(tmp)
}
temp.grid <- compiler::cmpfun(temp.grid)
system.time(
  importance.scores <- foreach(
    i = 1:length(g),
    .combine = 'rbind'
  ) %dopar% {
    ano<-0
    tmp <-NA
    try(tmp <- temp.grid(g = g[i], ano = ano))
    try(names(tmp) <- make.names(names(tmp)))
    library(geojsonio)
    shp_json <- geojson_json(tmp)
    #folder with results in gson format
    geojson_write(shp_json, file = paste0("clima/temp_",ano,"_",i,".geojson"))
    #file for checking and controlling parallelism
    write.table(x = paste(Sys.time()," - ",i), file = "log.txt", append = T, row.names = T, col.names =
                  F)
    return(tmp)
  }
)