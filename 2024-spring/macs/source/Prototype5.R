library(tidyverse)
library(pliman)

# upload file by using only one FOLDER!!!

input1 <- "/Users/mackia/Desktop/EEEEEE"
input2 <- "/Users/mackia/Desktop/RRR"


ip <-data.frame(input1, input2)


for(x in ip)
  {

  path <- x
  path_result <- paste(x, "result", sep = "/")
  folder <- image_import(image,
                         which = 5,
                         pattern = "Img",
                         path = path,
                         resize = F,
                         plot = F)
  npic <- length(folder)
  
  for(i in 1:npic) 
  {
    jpeg(paste0( 
      paste(path_result, "RAW_IMG", sep = "/"),  i, ".jpg"),
      unit = "px",
      width = 2040, 
      height = 1528)
    
    set.seed(i)
    
    if(i == 1){
      plot(folder$Img1.jpg)
      dev.off()
      dev.set()
    } else if(i==2){
      plot(folder$Img2.jpg)
      dev.off()
      dev.set()
    } else if(i==3){
      plot(folder$Img3.jpg)
      dev.off()
      dev.set()
    } else if(i ==4){
      plot(folder$Img4.jpg)
      dev.off()
      dev.set()
    } else if(i==5){
      plot(folder$Img5.jpg)
      dev.off()
      dev.set()
    }
    
  
  }
  #upload picture to R for calculating
  
  for(x in 1:npic){
    
    set.seed(x)
    assign(paste0("Img", x), image_import(paste(path_result, paste0("RAW_IMG", x, ".jpg"), sep = "/"),
                                          which = 1,
                                          pattern = NULL,
                                          path = NULL,
                                          resize = FALSE,
                                          plot = TRUE,
                                          nrow = NULL,
                                          ncol = NULL))}
  
  # adjust threshold and analyze picture
  
  for(y in 1:npic){
    set.seed(y)
    if(y ==1){assign("ImgR1", analyze_objects(Img1,
                                              marker = "id",
                                              index = "BI",
                                              lower_noise = 0.1))
    } else if(y==2){assign("ImgR2", analyze_objects(Img2,
                                                    marker = "id",
                                                    index = "BI",
                                                    lower_noise = 0.1))
    } else if(y==3){assign("ImgR3", analyze_objects(Img3,
                                                    marker = "id",
                                                    index = "BI",
                                                    lower_noise = 0.1))
    } else if(y==4){assign("ImgR4", analyze_objects(Img4,
                                                    marker = "id",
                                                    index = "BI",
                                                    lower_noise = 0.1))
    } else if(y==5){assign("ImgR5", analyze_objects(Img5,
                                                    marker = "id",
                                                    index = "BI",
                                                    lower_noise = 0.1))}
    
  }
  #==========================================================================
  
  micrometer <- function(a){
    return(a*10000)}
  
  Anabead <- function(a){
    return(get_measures(a, dpi = 79016) %>%
             micrometer())}
  
  for(i in 1: npic){
    set.seed(i)
    if(i ==1){assign("result1", Anabead(ImgR1))
    } else if(i ==2){assign("result2", Anabead(ImgR2))
    } else if(i ==3){assign("result3", Anabead(ImgR3))
    } else if(i ==4){assign("result4", Anabead(ImgR4))
    } else if(i ==5){assign("result5", Anabead(ImgR5))} }
  
  #==========================================================================
  library(tidyverse)
  Filterbead <- function(df){
    return(df %>% filter(circularity_norm > 6000 & circularity > 20000) %>% 
             select(id, diam_mean, area, circularity, circularity_norm) %>% 
             arrange(circularity_norm))} 
  
  for(i in 1:npic){
    if(i ==1){assign("resultF1", Filterbead(result1))
    } else if(i ==2){assign("resultF2", Filterbead(result2))
    } else if(i ==3){assign("resultF3", Filterbead(result3))
    } else if(i ==4){assign("resultF4", Filterbead(result4))
    } else if(i ==5){assign("resultF5", Filterbead(result5))}
  }
  
  #==========================================================================
  if(npic==5){
    assign("combine_sum", rbind(resultF1, resultF2, resultF3, resultF4, resultF5))
  } else if(npic ==4){
    assign("combine_sum", rbind(resultF1, resultF2, resultF3, resultF4))
  } else if(npic ==3){
    assign("combine_sum", rbind(resultF1, resultF2, resultF3))
  } else if(npic ==2){
    assign("combine_sum", rbind(resultF1, resultF2))
  } else if(npic ==1){
    assign("combine_sum", rbind(resultF1))
  } 
  
  #==========================================================================
  
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



