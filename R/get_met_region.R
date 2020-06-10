get_worldclim_rasters <- function(topath, clean=FALSE){
  
  download_worldclim <- function(basen, topath){
    
    wc_fn_full <- file.path(topath, basen)
    
    if(!file.exists(wc_fn_full)){
      message("Downloading WorldClim 10min layers ... ", appendLF=FALSE)
      download.file(file.path("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur",basen),
                    wc_fn_full, mode="wb")
      message("done.")
    }
    
    u <- unzip(wc_fn_full, exdir=topath)
    
    return(u)
  }
  
  download_worldclim("tmean_10m_esri.zip", topath)
  download_worldclim("prec_10m_esri.zip", topath) 
  
  if(clean){
    unlink(c(wc_fn_full,dir(file.path(topath,"tmean"),recursive=TRUE)))
    unlink(c(wc_fn_full,dir(file.path(topath,"prec"),recursive=TRUE)))
  }
  
  # Read the rasters into a list
  tmean_raster <- list()
  prec_raster <- list()
  #message("Reading Worldclim rasters ... ", appendLF = FALSE)
  for(i in 1:12){
    tmean_raster[[i]] <- raster(file.path(topath, sprintf("tmean/tmean_%s", i)))
    prec_raster[[i]] <- raster(file.path(topath, sprintf("prec/prec_%s", i)))
  }
  #message("done.")
  
  return(list(tmean_raster=tmean_raster, prec_raster=prec_raster))
}


get_worldclim_prectemp <- function(data, topath=tempdir(), return=c("all","summary"), worldclim=NULL){
  
  return <- match.arg(return)
  
  if(is.null(worldclim)){
    worldclim <- get_worldclim_rasters(topath)
  }
  tmean_raster <- worldclim$tmean_raster
  prec_raster <- worldclim$prec_raster
  
  #extract worldclim data; extract the gridCell ID for observations
  tmeanm <- precm <- matrix(ncol=12, nrow=nrow(data))
  for(i in 1:12){
    tmeanm[,i] <- 0.1 * extract(tmean_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
    precm[,i] <- extract(prec_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
  }
  colnames(tmeanm) <- paste0("tmean_",1:12)
  colnames(precm) <- paste0("prec_",1:12)
  
  pxy <- cbind(data, as.data.frame(tmeanm), as.data.frame(precm))
  names(pxy)[2:3] <- c("longitude","latitude")
  
  pxy$MAT <- apply(pxy[,grep("tmean_",names(pxy))],1,mean)
  pxy$MAP <- apply(pxy[,grep("prec_",names(pxy))],1,sum)
  
  # 
  if(return == "all")return(pxy)
  
  if(return == "summary"){
    
    dfr <- suppressWarnings(with(pxy, data.frame(species=unique(data$species),
                                                 n=nrow(data),
                                                 lat_mean=mean(latitude,na.rm=TRUE),
                                                 long_mean=mean(longitude,na.rm=TRUE),
                                                 MAT_mean=mean(MAT,na.rm=TRUE),
                                                 MAT_q05=quantile(MAT,0.05,na.rm=TRUE),
                                                 MAT_q95=quantile(MAT,0.95,na.rm=TRUE),
                                                 MAP_mean=mean(MAP,na.rm=TRUE),
                                                 MAP_q05=quantile(MAP,0.05,na.rm=TRUE),
                                                 MAP_q95=quantile(MAP,0.95,na.rm=TRUE))))
    rownames(dfr) <- NULL
    return(dfr)
  }
}

library(raster)
# 
coord.df = read.csv('data/part2 of Provenance _grass_seeds.csv',
                    skip=1)
coord.df$latitude = coord.df$Latidude
coord.df$longitude = coord.df$Longitude

out.df = get_worldclim_prectemp(coord.df,'download')

# write.csv(out.df,'grass_clim_envl_part1.csv',row.names = F)
write.csv(out.df,'grass_clim_envl_part2.csv',row.names = F)


