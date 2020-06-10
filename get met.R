# packages:
library(raster)

# Functions
# 
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
# 
get_worldclim_rasters <- function(topath){
  
  download_worldclim("tmean_10m_esri.zip", topath)
  download_worldclim("prec_10m_esri.zip", topath) 
  
  # Read the rasters into a list
  tmean_raster <- list()
  prec_raster <- list()
  #message("Reading Worldclim rasters ... ", appendLF = FALSE)
  for(i in 1:12){
    tmean_raster[[i]] <- raster(file.path(topath, sprintf("tmean/tmean_%s", i)))
    prec_raster[[i]] <- raster(file.path(topath, sprintf("prec/prec_%s", i)))
  }
  
  return(list(tmean_raster=tmean_raster, prec_raster=prec_raster))
  saveRDS(list(tmean_raster=tmean_raster, prec_raster=prec_raster),
          'p_t.rds')
}
# 
get_worldclim_prectemp <- function(data, topath='download'){
  if(!dir.exists(topath)){dir.create(topath)}
  
  if(!file.exists('p_t.rds')){
    worldclim <- get_worldclim_rasters(topath)
  }else{
    worldclim <- readRDS('p_t.rds')
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
  
  out.df <- data.frame(longitude = pxy$longitude,
                       latitude = pxy$latitude,
                       MAT = pxy$MAT,
                       MAP = pxy$MAP )
  
  return(out.df)
}

# the data has to have lon and lat
get_worldclim_prectemp(data.frame(longitude = 120,latitude = 30))
