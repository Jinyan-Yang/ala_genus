# Downloads the script into your working directory
url <- "https://bitbucket.org/!api/2.0/snippets/remkoduursma/pp6R4/263ee377b87bf487b3687c90755f7ef658f71151/files/ala_gbif_worldclim_species_envelope.r"
download.file(url, basename(url))

# and source
source("ala_gbif_worldclim_species_envelope.r")

library(raster)
## Loading required package: sp
library(rgbif)  # when using GBIF
library(ALA4R)  # when using ALA


ala.genus.func <- function(pft,gen.in){
  
  out.fn <- sprintf('data/%s/%s.df',pft,gen.in)
  if(!file.exists(sprintf('data/%s.df.rds',gen.in))){
    Agrostis <- occurrences(paste0("genus:",gen.in),download_reason_id=7,
                            fields=c("longitude","latitude",
                                     "genus","family",'species',
                                     'country','state'))
    
    
    Agrostis.df <- t(do.call(rbind,Agrostis[[1]]))
    
    Agrostis.df <- as.data.frame(Agrostis.df)
    Agrostis.df <- Agrostis.df[Agrostis.df$country == 'Australia',]
    Agrostis.df <- Agrostis.df[,c("longitude","latitude",
                                  "genus","family",'species',
                                  'country','state')]
    
    Agrostis.df <- Agrostis.df[complete.cases(Agrostis.df),]
    
    saveRDS(Agrostis.df,out.fn)
  }else{
    print(paste0(out.fn,' already exists'))
  }
}

C3.vec <- c('Agrostis', 'Avena', 'Bromus', 'Dactylis', 'Danthonia', 
            'Festuca', 'Holcus', 'Hordeum', 'Lolium', 'Phalaris', 
            'Phleum', 'Poa', 'Rytidosperma', 'Triticum', 'Zoysia') 

# c3
for(nm.in in C3.vec){
  ala.genus.func('c3',nm.in)
}

temp.ls <- list()

file.ls <- list.files('data/c3')

for(i in seq_along(file.ls)){
  temp.ls[[i]] <- readRDS(file.path('data','c3',file.ls[i]))
}

c3.df <- do.call(rbind,temp.ls)
c3.df$longitude <- as.numeric(as.character(c3.df$longitude))
c3.df$latitude <- as.numeric(as.character(c3.df$latitude))

saveRDS(c3.df,'data/c3.rds')

# c4
C4.vec <- c('Andropogon', 'Aristida', 'Astrebla', 'Cenchrus', 
            'Cynodon', 'Echinochloa', 'Enteropogon', 'Eragrostis', 
            'Panicum', 'Paspalum', 'Sorghum', 'Sporobolus', 'Themeda', 'Triodia')

for(nm.in in C4.vec){
  ala.genus.func('c4',nm.in)
}

temp.ls <- list()

file.ls <- list.files('data/c4')

for(i in seq_along(file.ls)){
  temp.ls[[i]] <- readRDS(file.path('data','c4',file.ls[i]))
}

c4.df <- do.call(rbind,temp.ls)

c4.df$longitude <- as.numeric(as.character(c4.df$longitude))
c4.df$latitude <- as.numeric(as.character(c4.df$latitude))

saveRDS(c4.df,'data/c4.rds')

c3.df <- readRDS('data/c3.rds')
c4.df <- readRDS('data/c4.rds')
# plot
tran.func <- function(cols,alpha=0.1){
  targe.vec <- col2rgb(cols)
  out.col <- apply(targe.vec,2,FUN = function(x){
    rgb(red = x["red"]/255, green=x["green"]/255, blue=x["blue"]/255,alpha=alpha)})
  return(out.col)
}


library(maps)
map('world','Australia',mar = c(1, 1, 1, 1))

points(latitude~longitude,data = c3.df,pch=15,col=tran.func('blue'),cex = 0.2)
points(latitude~longitude,data = c4.df,pch=15,col=tran.func('coral'),cex = 0.2)

par(mfrow=c(1,2))
smoothScatter(c3.df$longitude,c3.df$latitude,
              colramp = colorRampPalette(c("white", 'blue')),pch = " ",
              xlim=c(110,160),
              ylim=c(-50,-10),
              xlab='lon',ylab='lat')
# par(new=TRUE)
smoothScatter(c4.df$longitude,c4.df$latitude,
              colramp = (colorRampPalette(c("white", 'red'))),pch = " ",
              xlim=c(110,160),
              ylim=c(-50,-10),
              xlab='lon',ylab='lat')

