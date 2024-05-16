library(tidyverse)
library(pliman) # requires Bioconductor package - EBImage
library(EBImage)
# upload file by using only one FOLDER!!!
#input1 <- "/Users/mackia/Desktop/Folder1"
#input2 <- "/Users/mackia/Desktop/Folder2"

input1 <- "./pic_prototype5/Folder1"
input2 <- "./pic_prototype5/Folder2"

img.files1 = list.files(path=input1,
                        pattern="*.jpg",
                        full.names = TRUE)
print(img.files1)

img1 = readImage(img.files1) # each images in an array 2040 1528 3 5
dim(img1) # each image is in an array of 2040 x1528 x 3 x 5

# resize images
img1 = resize(img1, w=1020, h=764)
dim(img1)

##
## save resized images in ./result folder (THIS part is NOT necessary!)
##
path_result <- paste(input1, "result", sep = "/")
print(path_result)
display(img1)

# naming files
files = paste(rep("Img",5),1:5,rep("resized.jpg",5), sep="")
paste(rep(path_result,5),files, sep="/")

writeImage(img1,
           files=paste(rep(path_result,5),files,sep="/")) 
# result/ folder should be created before saving...
##
##

##
## for calculating
##  

# adjust threshold and analyze picture
str(img1)
dim(img1)
colorMode(img1) = "Color"
colorMode(img1)
numberOfFrames(img1, type = "total") # 15 frames (RGB 3 frames x 5 images)
numberOfFrames(img1, type = "render") # 5 images
getFrame(img1,1, type="render") # "render" gets a single image
getFrames(img1, type="render")
##
##
micrometer <- function(a) { return(a*10000) }
Anabead <- function(a) { return(get_measures(a, dpi = 79016) %>% micrometer()) }  # 79016 dpi? 
Filterbead <- function(df) { 
  return(df %>% filter(circularity_norm > 6000 & circularity > 20000) %>% 
           select(id, diam_mean, area, circularity, circularity_norm) %>% 
           arrange(circularity_norm))
  } 
##

##
## each images in 'img1' object is called, analyzed...
##

npic = numberOfFrames(img1, type = "render") # total number of images
combine_sum = NULL # your final data?

for (i in 1:npic) {
   # 'i'th image in 'img1' object
   #i = 1
   x = analyze_objects(getFrame(img1,i,type="render"),
                       marker = "id", 
                       index = "I",          # <--- BI, I
                       object_size = "large", # <--- small, medium, large, elarge
                       lower_noise = 0.1,     # <--- change this?
                       lower_size = NULL,     # <--- change this?
                       upper_size = NULL,
                       lower_eccent = NULL,
                       lower_circ = NULL)
   # index https://tiagoolivoto.github.io/pliman/articles/indexes.html
   # 'BI' using brightness for detecting objects??
   # 'x' is the result of analyze_objects
   plot_measures(x, measure = "id")
   names(x)
   #[1] "results"          "statistics"       "object_rgb"       "object_index"     "efourier"        
   #[6] "efourier_norm"    "efourier_error"   "efourier_power"   "efourier_minharm" "veins"           
   #[11] "angles"           "width_at"         "mask"             "pcv"              "contours"        
   #[16] "parms
   x$results %>% 
     filter(circularity_norm > 0.6 & circularity > 2) %>%
     select(id, diam_mean, area, circularity, circularity_norm) %>% 
     arrange(circularity_norm)
   #
   y = get_measures(x, dpi=79016)
   names(y)
   # what does this number (79016) determine?
   # what for 'micrometer()' function?
   class(y)
   class(y) = "data.frame"
   y %>% 
     filter(circularity_norm > 0.6 & circularity > 2) %>%
     select(id, diam_mean, area, circularity, circularity_norm) %>% 
     arrange(circularity_norm)
   # collecting/combining results
   combine_sum = bind_rows(combine_sum, y, .id="img")
}

##
## histogram
##
  combinebead <- combine_sum %>% select(diam_mean)
  
  #==========================================================================
  #set function for creating HT
  BuildHT_SUM <- ggplot(data = combinebead)+
    aes(x = diam_mean, y = ..density..) +
    geom_histogram(fill = "azure3") +
    geom_density(alpha = 0.3, fill = "coral1") +
    labs(x="Bead size (um)", y="Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          legend.position = "bottom")
  
  
  BuildHT_DEN <- ggplot(data = combinebead)+
    aes(x = diam_mean, y = ..density..) +
    geom_density(alpha = 0.3, fill = "coral1") +
    labs(x="Bead size (um)", y="Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          legend.position = "bottom")
  
  BuildHT_HT <- ggplot(data = combinebead)+
    aes(x = diam_mean, y = ..density..) +
    geom_histogram(fill = "azure4") +
    labs(x="Bead size (um)", y="Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          legend.position = "bottom")
  
  
  #======================================================================
  
  #export Histogram_SUM
  png(paste(path_result, "Histogram summary.png", sep = "/"))
  plot(BuildHT_SUM)
  dev.off()
  dev.set()
  
  png(paste(path_result, "Histogram (Density).png", sep = "/"))
  plot(BuildHT_DEN)
  dev.off()
  dev.set()
  
  png(paste(path_result, "Histogram.png", sep = "/"))
  plot(BuildHT_HT)
  dev.off()
  dev.set()
  
  write_csv(combine_sum, paste(path_result, "summarize.csv", sep = "/"))
  
}



