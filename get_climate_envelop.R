library(raster)
## Loading required package: sp
library(rgbif)  # when using GBIF
library(ALA4R)  # when using ALA

# functions
get_occurrences_gbif <- function(species){
  spdat <- occ_search(scientificName=species,
                      limit=50000,
                      fields =c('name','decimalLatitude','decimalLongitude'),
                      hasGeospatialIssue=FALSE,
                      return="data")
  
  # remove obs with no lats and longs
  spdat <- as.data.frame(na.omit(spdat))
  
  # rename
  names(spdat) <- c("species","longitude","latitude")
  
  return(spdat)
}

get_occurrences_ala <- function(species){
  
  spdat <- occurrences(taxon=species, download_reason_id=7)
  
  spdat <- subset(spdat$data, 
                  !is.na(longitude),
                  select = c(species,longitude,latitude))
  
  # Returns other species as well, for some reason
  spdat <- spdat[spdat$species == species,]
  
  return(spdat)  
}

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

rasterize_occurrences <- function(spdat){
  
  # make a new raster same size as worlclim but each gridcellhas ID number
  gridcellID <- raster(nrow=900,ncol=2160,extent(c(-180,180,-60,90)))
  gridcellID[] <- 1:1944000
  
  # get centerpoint of gridcells where species occur, 1 observation for each gridcell
  spdat$GridID <- raster::extract(gridcellID,cbind(spdat$longitude,spdat$latitude),
                                  method='simple')
  
  # extract the gridCell ID for observations
  presence <- as.data.frame(cbind(1,unique(spdat$GridID)))
  
  # make dataframe of cells present in
  colnames(presence) <- c("val","gridID")
  
  # raster of presence
  presence_raster <- raster::subs(gridcellID,
                                  as.data.frame(table(presence$gridID)))
  
  # convert lat and longs of gridcells where species is present
  pxy <- as.data.frame(rasterToPoints(presence_raster))
  pxy <- pxy[,-3]
  names(pxy) <- c("longitude","latitude")
  
  pxy <- cbind(data.frame(species=unique(spdat$species)), pxy)
  
  return(pxy)
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

worldclim_presence <- function(species, database=c("ALA","GBIF"), topath=tempdir(), 
                               return=c("summary","all")){
  
  return <- match.arg(return)
  database <- match.arg(database)
  
  worldcl <- get_worldclim_rasters(topath)
  
  l <- list()
  
  for(i in seq_along(species)){
    
    if(database == "GBIF"){
      spocc <- get_occurrences_gbif(species[i])
    } else if(database == "ALA"){
      spocc <- get_occurrences_ala(species[i])
    }
    
    r <- rasterize_occurrences(spocc)
    
    l[[i]] <- get_worldclim_prectemp(r, topath=topath, return=return, worldclim=worldcl)
  }
  
  if(return == "all"){
    if(length(species) == 1)
      return(l[[1]])
    else{
      names(l) <- species
      return(l)
    }
  } else {
    return(do.call(rbind,l))
  }
  
}


# read 
species.df <- read.csv('Sorted_grass_species_list_25062019.csv')
# species.df <-data.frame(Binomial.name = c('Lachnagrostis filiformis'))

# library('Taxonstand')
# sp.nm <- as.character(species.df$Binomial.name)[1]
# TPL(sp.nm,infra = TRUE, 
#     corr = TRUE)
# 
# TPL(genus = strsplit(sp.nm,' ')[[1]][1], species = strsplit(sp.nm,' ')[[1]][2],
#     corr = TRUE)

condition.ls <- list()

for(i in seq_along(species.df$Binomial.name)){
  condition.ls[[i]] <- try(worldclim_presence(species.df$Binomial.name[i],return="summary",topath='download'))
  
  if(class(condition.ls[[i]])=="try-error"){condition.ls[[i]]<-NA}
}
out.df <- do.call(rbind,condition.ls)
# condition.df <- worldclim_presence(species.df$Binomial.name[1],return="summary",topath='download')
write.csv(out.df,'clim_enve.csv',row.names = F)
# condition.df <- worldclim_presence(species.df$Binomial.name[5],return="summary",topath='download')
# 
# see <- occurrences('species:Austrostipa bigeniculata ',download_reason_id = 7)
# see.df <- t(do.call(rbind,see[[1]]))
