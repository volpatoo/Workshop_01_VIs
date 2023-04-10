# Vegetation index extractions using MS reflectance images
> Recommended to large data set
> 
```r
####################################################################################
########### Vegetation index extractions using MS reflectance images ##############
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


####Change it###
dir_multispec <- "C:\\temp_vi\\Sample_MS"
folder_shp <- "C:\\temp_vi\\Sample_MS\\Shp"
layer_prefix_shp <- "Shapefile_MRC_02"
###Get all the files to process from the WD
imgFiles <-list.files(path = dir_multispec, pattern="*.tif$",full.names = T) #get the Orthosaics. Files that their name ends in group1.tif (Change all file names to otimization)
imgFiles

imgFiles_name <-list.files(path = dir_multispec, pattern="*.tif$")
imgFiles_name <- gsub(".tif", "", imgFiles_name, ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
imgFiles_name

## Read the plots Shapefile
indPlots <- readOGR(dsn = folder_shp, layer = layer_prefix_shp)
## Polygon plot ID list from the shapefiles
names(indPlots)
## Data frame with plot names
VIs.Table <- as.data.frame(indPlots[,"name"])
#VIs list to use in the loop
myIndex_list_MS<- c("PSRI","NDVI","GNDVI","RVI","NDRE","TVI","CVI","CIG","CIRE","DVI","EVI")

imageMulti_names_list <- imgFiles_name %>% substr( start = 1, stop = 7)

imageMulti_levels<- levels(as.factor(imageMulti_names_list))

func_list<- c('mean', 'median', 'sd')

gc()

#If the parallel functions does no work well
#cl <- parallel::makeCluster(n.core, output = "", setup_strategy = "sequential")

# Number of cores
n.core<-detectCores()

# Starting parallel
cl <- makeCluster(n.core, output = "")
registerDoParallel(cl)
getDoParWorkers()

j=0
k=5

system.time(
  
  for(i in 1:(length(imgFiles_name)/5)){
    
    if(i == 1) {
      j = j+0 }
    else {
      j = j +5
    }   
    
    message("Stacking relectance images from: ", imageMulti_levels[i],"_MRC_MS")
    imageMulti <- imgFiles[(j+1):(i*k)]
    
    #i.h<-aggregate(stack(imgFiles[i]), fact=aggregateCells) 
    imageMulti.blue <- imageMulti[[1]]
    imageMulti.gree <- imageMulti[[2]]
    imageMulti.nir <- imageMulti[[3]]
    imageMulti.edge <- imageMulti[[4]]
    imageMulti.red <- imageMulti[[5]]
    
    i.h <- stack(imageMulti.blue,
                 imageMulti.gree,
                 imageMulti.nir,
                 imageMulti.edge,
                 imageMulti.red)
    
    for(v in 1:length(myIndex_list_MS)){ 
      
      message("Using VI: ",paste(myIndex_list_MS[v]))
      
      for(f in 1:length(func_list)){ 
        
        message("Using function: ",paste(func_list[f]))
        
        results<- foreach(i = 1:length(indPlots), 
                          .packages = c("raster", "FIELDimageR", "plyr", "dplyr"), 
                          .combine = rbind) %dopar% {
                            h.c <-  crop(i.h, extent(indPlots[i,]))
                            m.h <-  fieldMask(mosaic=h.c, #i.h@layers
                                              Red=5,
                                              Green=2,
                                              Blue=1,
                                              RedEdge = 4, 
                                              NIR = 3,
                                              index="NDVI",
                                              cropValue=0.7, 
                                              cropAbove=F, ## Removes any instance of soil from tif file
                                              plot = F) 
                            # stackMS@layers
                            Veg.Indices<-fieldIndex(mosaic = m.h$newMosaic, 
                                                    Red=5,Green=2,Blue=1,RedEdge=4,NIR=3,
                                                    index = myIndex_list_MS[v],
                                                    plot = FALSE)
                            
                            # projection(indPlots)<-projection(Veg.Indices)
                            
                            
                            raster::extract(x = Veg.Indices[[6]], y = indPlots[i,], fun = eval(parse(text = func_list[f])),  
                                            buffer = buffer, na.rm = T, df = T)
                            
                            
                            
                          }
        
        results$Func<-func_list[f]
        
        results$Plot_ID <- 1:length(indPlots)
        colnames(results) <- c("ID","value","Func","Plot_ID")
        
        if(f==1){results.1<-results}else{results.1<-rbind(results.1, results)}

        
      }
      
      
      results.1$VIs<-myIndex_list_MS[v]
      
      if(v==1){results.2<-results.1}else{results.2<-rbind(results.2, results.1)}
      
    }
    results.2$imgFiles_name<-imageMulti_names_list[i]
    
    if(k==5){results.3<-results.2}else{results.3<-rbind(results.3, results.2)}
    
  })

parallel::stopCluster(cl)

as_tibble(results.3)

results.final<-pivot_wider(results.3, names_from   = c("VIs", "Func"),
                           values_from = "value")

head(results.final)

# saving
write.csv(results.final, "RGB_VIs_MRC_2021_test_R.csv", quote = F, row.names = F)


```


# Vegetation index extractions using Multispectral imagery 
> Recommended for small data set only

```r
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
               doParallel
               
) 

####Change it###
dir_multispec <- "C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\MS_reflectance"
folder_shp <- "C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Data_VIs\\Shapefile"
layer_prefix_shp <- "Shapefile_MRC_02"
###Get all the files to process from the WD
##Recommended to order the files
imgFiles <-list.files(path = dir_multispec, pattern="*.tif$",full.names = T) #get the Orthosaics. Files that their name ends in group1.tif (Change all file names to otimization)
imgFiles_name <-list.files(path = dir_multispec, pattern="*.tif$")
imgFiles
imgFiles_name
## Read the plots Shapefile
indPlots <- readOGR(dsn = folder_shp, layer = layer_prefix_shp)
## Polygon plot ID list from the shapefiles
names(indPlots)
## Data frame with plot names
VIs.Table <- as.data.frame(indPlots[,"name"])

myIndex_list_MS<- list("PSRI","NDVI","GNDVI","RVI","NDRE","TVI","CVI","CIG","CIRE","DVI","EVI")
#myIndex_list_name<- myIndex_list_name
# Number of cores

imageMulti_names_list <- imgFiles_name %>% substr( start = 1, stop = 7)

imageMulti_levels<- levels(as.factor(imageMulti_names_list))
#imageMulti_levels<-imageMulti_levels[2]
  
gc()
rasterOptions()
rasterOptions(chunksize = 1e+12)
rasterOptions(maxmemory = 1e+12)

#aggregateCells = 48
j=0
k=5

ptm <- proc.time() #6:00 PM 1-15-2022
for(i in 1:(length(imgFiles_name)/5)){
  
  if(i == 1) {
    j = j+0 }
  else {
    j = j +5
  }
  
  message("Stacking relectance images from: ", imageMulti_levels[i],"MRC_MS", sep="_")
  imageMulti <- imgFiles[(j+1):(i*k)] #function to select in a loop the 5 columns from each flight (data point)
                                      #image reflectance have to be in order 
  
  # imageMulti <- imgFiles %>%
  #   as_tibble() %>%
  #   dplyr::filter(stringr::str_detect(value, imageMulti_levels[1])) #image reflectance have to be in order 
  # 
  # imageMulti <- pull(imageMulti, value)
  
  
  imageMulti.blue <- imageMulti[[1]]
  imageMulti.gree <- imageMulti[[2]]
  imageMulti.nir <- imageMulti[[3]]
  imageMulti.edge <- imageMulti[[4]]
  imageMulti.red <- imageMulti[[5]]
  
  stackMS <- stack(imageMulti.blue,
                   imageMulti.gree,
                   imageMulti.nir,
                   imageMulti.edge,
                   imageMulti.red)
 
  # Removing the soil using index and mask from step 4:
  MS.RemSoil<-fieldMask(mosaic=stackMS,
                         Red=5,
                         Green=2,
                         Blue=1,
                         RedEdge = 4, 
                         NIR = 3,
                         index="NDVI",
                         cropValue=0.7, 
                         cropAbove=F, ## Removes any instance of soil from tif file
                         plot = F) 

 # stackMS@layers
  Veg.Indices<-fieldIndex(mosaic = MS.RemSoil$newMosaic, 
                          Red=5,Green=2,Blue=1,RedEdge=4,NIR=3,
                          #index = c("NDVI","NDRE"),
                          plot = FALSE)
  
  projection(indPlots)<-projection(Veg.Indices)
  
  EX1.Info1<-fieldInfo(mosaic=Veg.Indices,
                       fieldShape=indPlots,
                       fun = "std",
                       n.core=12,
                       projection = FALSE) ## Extracts all vegetation indices (layers) for each .shp file plot
  
  EX1.Info1_data<- EX1.Info1$fieldShape@data[,15:20] # Need to change here according to the number of columns from attribute shapefile
  colnames(EX1.Info1_data)<- paste(colnames(EX1.Info1_data),imageMulti_levels[i],sep="_")
  VIs.Table<- cbind(VIs.Table,EX1.Info1_data) #save all columns into a data frame
  
  rm(Veg.Indices,EX1.Info1,EX1.Info1_data)
  gc()
  
  ## 33 VIs will be used from RGB images!!!
  for( s in 1:(length(myIndex_list_MS))){
    
    message("Using VI: ",paste(myIndex_list_MS[s]))
    Veg.Indices2<-fieldIndex(mosaic = MS.RemSoil$newMosaic, 
                            Red=5,Green=2,Blue=1,RedEdge=4,NIR=3,
                            index=myIndex_list_MS[s],
                            plot = FALSE)
    
    projection(indPlots)<-projection(Veg.Indices2) ##Makes coordinate system of .tif with indices and without soil and .shp the same

    EX1.Info2<-fieldInfo(mosaic=Veg.Indices2[[6]],
                         fieldShape=indPlots,
                         fun = "std",
                         n.core=12,
                         projection = FALSE) ## Extracts all vegetation indices (layers) for each .shp file plot
    
    EX1.Info1_data2<- EX1.Info2$fieldShape@data[15]
    colnames(EX1.Info1_data2)<- paste(myIndex_list_MS[s],imageMulti_levels[i],"MRC_MS", sep="_")
    VIs.Table<- cbind(VIs.Table,EX1.Info1_data2) #save all columns into a data frame
    
    rm(EX1.Info2,EX1.Info1_data2,Veg.Indices2)
    gc()
    
  }
  
  rm(MS.RemSoil,imageMulti,stackMS)

}
    

proc.time() - ptm

data_beans_col <- as.data.frame(indPlots[,1:13])

VIs_MS_outp <- cbind(data_beans_col, VIs.Table)

write.csv(VIs_MS_outp, "MS_VIs_MRC_2021_mean.csv", quote = F, row.names = F)
```

# Adjusting the data set to have date of flight in factorial column
> Pivot data set by flight date
```r
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
               tidyr
) 

####Change it###
dir_multispec <- "C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\MS_reflectance"

VIs.out1<- read.csv("MS_VIs_MRC_2021_sd.csv")

###Get all the files to process from the WD
imgFiles <-list.files(path = dir_multispec, pattern="*.tif$",full.names = T) #get the Orthosaics. Files that their name ends in group1.tif (Change all file names to otimization)
imgFiles_name <-list.files(path = dir_multispec, pattern="*.tif$")
imgFiles_name
imgFiles

imageMulti_names_list <- imgFiles_name %>% substr( start = 1, stop = 7)

imageMulti_levels<- levels(as.factor(imageMulti_names_list))


# Taking the file names
names<- c("Blue", "Green", "NIR", "RedEdge", "Red", "HUE", "PSRI","NDVI",
          "GNDVI","RVI","NDRE","TVI","CVI","CIG","CIRE","DVI","EVI", 'Flight_Date')


for (i in 1:length(imageMulti_levels)) {

VIs.out1.sub1<-VIs.out1 %>%
  select(contains(imageMulti_levels[i]))
VIs.out1.sub1$Flight_Date<- imageMulti_levels[i]
# length(names)
# length(VIs.out1.sub1)

colnames(VIs.out1.sub1)<- c(names)

if(i==1){VIs.out1.sub2<-VIs.out1.sub1}else{VIs.out1.sub2<-rbind(VIs.out1.sub2, VIs.out1.sub1)}

}

VIs.out1.sub2 = VIs.out1.sub2 %>% select(Flight_Date, everything())
head(VIs.out1.sub2)

VIs.out1.adj<- cbind(VIs.out1[,1:13], VIs.out1.sub2)


write.csv(VIs.out1.adj, "MS_VIs_MRC_2021_sd_adj.csv", quote = F, row.names = F)
```
