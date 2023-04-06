####################################################################################
########### Vegetation index extractions using RGB orthomosaic images ##############
####################################################################################

### Setting up the working directory 
rm(list=ls())
my.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(my.path)

if(!require("pacman")) install.packages("pacman")

pacman::p_load(FIELDimageR,
               raster,
               sp,
               dplyr,
               rgdal,
               lidR,
               parallel,
               foreach,
               doParallel,
               plyr,
               stringr,
               tidyverse,
               mapview
               
) 

# install.packages("mapview")
# install.packages("devtools")
# devtools::install_github("filipematias23/FIELDimageR", force = T)
# install.packages('raster', repos='https://rspatial.r-universe.dev')
# install.packages('terra', repos='https://rspatial.r-universe.dev')



# Basic pipeline

## Setting the work directories

### Folder directory containing the orthomosaics

dir_ortho <- "./MOSAIC"


### Folder directory containing the shapefiles

folder_shp <- "./Shapefile"


## Getting the files

imgFiles <-list.files(path = dir_ortho, pattern="*.tif$",full.names = T) #get the Orthosaics. 
imgFiles #Files that their name ends in group1.tif (Change all file names to otimization)

imgFiles_name <-list.files(path = dir_ortho, pattern="*.tif$") # Aux file names
imgFiles_name <- gsub(".tif", "", imgFiles_name, ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE) # Removing the extension
imgFiles_name



## Loading and Reading the plots Shapefile

# Loading the shape files
layer_prefix_shp <- "Export_Output" # name with the .shp extension
## Read the plots Shapefile
indPlots <- readOGR(dsn = folder_shp, layer = layer_prefix_shp)
## Polygon plot ID list from the shapefiles
names(indPlots)


## Getting the RGB bands

img.rgb.1 <-  stack(imgFiles[1])
img.rgb.1

##Ploting the orthomosaic
plotRGB(img.rgb.1)

plot(indPlots,add=T,col="red") # merge the external shp file and orthomosaic


## Reducing the image resoution if necessary

img.rgb.2<-aggregate(img.rgb.1, fact=4) ## reducing the resoluation. Using high numbers for fact is not desired. 


## Cropping the plot image

plot_crop_idx_1 <- crop(img.rgb.1,extent(indPlots[1,]))


## Plotting the plot image

plotRGB(plot_crop_idx_1)


## Building vegetation indices

img.index <- fieldIndex(mosaic = plot_crop_idx_1, index = c("BI", "GLI", "HI", "NGRDI", "BGI", "VARI", "SCI")) 
plot(img.index$NGRDI)


## Removing the soil

EX.Mask1 <- fieldMask(mosaic = plot_crop_idx_1, index = "NGRDI", cropValue = 0.05, cropAbove = F)

EX.Mask2 <- fieldMask(mosaic = plot_crop_idx_1, index = "GLI", cropValue = 0.05, cropAbove = F)

EX.Mask3 <- fieldMask(mosaic = plot_crop_idx_1, index = "HI", cropValue = 0.6, cropAbove = T) #all values above the cropValue will be accounted to make the mask.

EX.Mask4 <- fieldMask(mosaic = plot_crop_idx_1, index = "HI", cropValue = 0.8, cropAbove = T)
EX.Mask5 <- fieldMask(mosaic = plot_crop_idx_1, index = "HI", cropValue = 0.6, cropAbove = T) ## This model will be select to perform the analysis



##  Building vegetation indices


#Investigating the raw image plot
img.index.hi <- fieldIndex(mosaic = plot_crop_idx_1, index = c("HI") )
hist(img.index.hi$HI) # Image segmentation start from 0.7 (soil and plants)

#Investigating the clean image plot
img.index.hi.clean <- fieldIndex(mosaic = EX.Mask5$newMosaic, index = c("HI") )
hist(img.index.hi.clean$HI) # Image segmentation start from 0.7 (soil and plants)


## Orthomosaic image


ortho.index.hi.clean<- fieldIndex(mosaic = img.rgb.1, index = c("HI"))
#dev.off()

ortho.HI.RemSoil<- fieldMask(mosaic = img.rgb.1, Red = 1, Green = 2, Blue = 3, 
                             index = "HI", cropValue = 0.7, cropAbove = T) 



## Extracting the pixels values


##Selection the VI from the clean orthomosaic
Veg.Indices<-fieldIndex(mosaic = ortho.HI.RemSoil$newMosaic,
                        index=c("GLI"))

## Conf. the CRS projections
proj4string(indPlots) <- proj4string(Veg.Indices) ##Makes coordinate system of .tif with indices and without soil and .shp the same

ortho.Info1<-fieldInfo(mosaic=Veg.Indices,
                       fieldShape=indPlots,
                       n.core=10) ## Extracts all vegetation indices (layers) for each .shp file

VI_data<- as_tibble(ortho.Info1$fieldShape@data)

write.csv(VI_data, "VI_data_GLI.csv",  row.names=F) 

hist(ortho.Info1$fieldShape$GLI)


## Graphic visualization of trait values for each plot using the fieldShape file and the orthomosaic


#names(indPlots)
### Interpolating colors: c("white","black")
fieldPlot(fieldShape=ortho.Info1$fieldShape,fieldAttribute="MaturityDa", mosaic=Veg.Indices, color=c("white","black"), alpha = 0.5)

### Interpolating colors: c("red","blue")
fieldPlot(fieldShape=ortho.Info1$fieldShape,fieldAttribute="GLI", mosaic=Veg.Indices, color=c("red","blue"), alpha = 0.5)


---------------------------------------------
  
  ## Getting started - Complete data set
  
  ## Setting up the working directory 
  
  rm(list=ls())
my.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(my.path)


# Setting the work directories

dir_ortho <- "./MOSAIC"


### Folder directory containing the shapefiles

folder_shp <- "./Shapefile"


# Getting the files

imgFiles <-list.files(path = dir_ortho, pattern="*.tif$",full.names = T) #get the Orthosaics. 
imgFiles #Files that their name ends in group1.tif (Change all file names to otimization)

imgFiles_name <-list.files(path = dir_ortho, pattern="*.tif$") # Aux file names
imgFiles_name <- gsub(".tif", "", imgFiles_name, ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE) # Removing the extension
imgFiles_name



# Loading and Reading the plots Shapefile

# Loading the shape files
layer_prefix_shp <- "Export_Output"
## Read the plots Shapefile
indPlots <- readOGR(dsn = folder_shp, layer = layer_prefix_shp)
## Polygon plot ID list from the shapefiles
names(indPlots)

# Subsetting plots to speed up the loop
indPlots<- indPlots[1:5,"MNPlot"]
indPlots


# Setting the parameters from the loop

## Plots ID

## Data frame with plot names
Plot_ID<- as.data.frame(indPlots[,"MNPlot"])
indPlots # Plot identification propose


## VIs list

## Import the functions file
source("VIs_RGB-aux.R")

#Reading 38 VIs from RGB image
myIndex_list<- myIndex_list
myIndex_list_name<- myIndex_list_name


## Functions of extractions

func_list<- c('mean', 'median', 'sd')
# There are many other functions (SD, VAR, Quantiles, etc)



## Starting the parallel function
### Please, take your time and study how this function works and if this function is for you.


gc() #Cleaning unusual memmory

# Number of cores
n.core<-detectCores() # or detectCores()

# Starting parallel
cl <- makeCluster(n.core, output = "")
registerDoParallel(cl)
getDoParWorkers()



## Running the loops

## Using the code inside the system time to get the total time used.
system.time(
  for(i in 1:length(imgFiles)){ #loop through images
    
    message("Processing ortho: ",paste(imgFiles_name[i]))
    #i.h<-aggregate(stack(imgFiles[k]), fact=aggregateCells) to reduce the image size if the case
    i.h <-  stack(imgFiles[i])
    
    for(v in 1:length(myIndex_list)){ #loop through VIs
      
      message("Using VI: ",paste(myIndex_list_name[v]))
      
      for(f in 1:length(func_list)){  #loop through extractions functions
        
        message("Using function: ",paste(func_list[f]))
        
        results<- foreach(p = 1:length(indPlots),  #loop through plots numbers
                          .packages = c("raster", "FIELDimageR"), 
                          .combine = rbind) %dopar% {
                            
                            # Step 1: crop the plot from the orthomosaci
                            h.c <-  crop(i.h, extent(indPlots[p,]))
                            
                            # Step 2: remove the soil using the HUE VI with a threshold = 0.6
                            m.h <-  fieldMask(mosaic=h.c,
                                              Red=1,
                                              Green=2,
                                              Blue=3,
                                              index="HUE",
                                              cropValue=0.6, #or 0.6 | 0.8
                                              cropAbove=T, ## Removes any instance of soil from tif file
                                              plot = FALSE) 
                            
                            # Step 3: Obtain the VI from the vegetation only
                            Veg.Indices<-fieldIndex(mosaic = m.h$newMosaic, myIndex=myIndex_list[v],
                                                    plot = FALSE)
                            
                            projection(indPlots)<-projection(Veg.Indices)
                            
                            # Step 4: Extraction the pixel value from each individual plot using the function 
                            raster::extract(x = Veg.Indices$myIndex, y = indPlots[p,], fun = eval(parse(text = func_list[f])),  
                                    buffer = buffer, na.rm = T, df = T)
                            
                          }
        
        results$Func<-func_list[f]
        
        # Step 5: Saving the results into a data frame with the lenght of the Plots ID
        results$ID <- 1:length(indPlots)
        
        if(f==1){results.1<-results}else{results.1<-rbind(results.1, results)}
        
      }
      
      # Step 6: Saving the vegetation index name for each function
      results.1$VIs<-myIndex_list_name[v]
      
      if(v==1){results.2<-results.1}else{results.2<-rbind(results.2, results.1)}
      
    }
    
    # # Step 6: Saving the vegetation index name from each function for each orthomosaic
    results.2$imgFiles_name<-imgFiles_name[i]
    
    if(i==1){results.3<-results.2}else{results.3<-rbind(results.3, results.2)}
    
  })

parallel::stopCluster(cl) # Stopping the parallel function


###### the end ###########


## Adjusting the results and saving

as_tibble(results.3) # Saving into a table

# Adjusting the data for VIs and Func columns
results.final<-pivot_wider(results.3, names_from   = c("VIs", "Func"),
                           values_from = "myIndex")

head(results.final)

# Saving
write.csv(results.final, "test01_MRF_RGB.csv", quote = F, row.names = F)





























