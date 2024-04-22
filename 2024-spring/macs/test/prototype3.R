library(tidyverse)
library(pliman)
#==============================================================

# upload file by using only one FOLDER!!!!

path <- "#INPUT FULL PATH OF FOLDER"
folder <- image_import(image,
                       which = 2,
                       pattern = "Img",
                       path = path,
                       resize = F,
                      plot = T
                    )


path2 <- paste(path, "result", sep = "/")

npic <- length(folder)

for(i in 1:npic) 
  {
  jpeg(paste( 
          paste(path2, "RAW_IMG", sep = "/"),  i, ".jpg"),
        unit = "px",
        width = 2040, 
        height = 1528)
  
  set.seed(i)
  
  if(i == 1){
    plot(folder$Prototype_Img1.jpg)
    dev.off()
    dev.set()
  } else if(i==2){
    plot(folder$Prototype_Img2.jpg)
    dev.off()
    dev.set()
  } else if(i==3){
    plot(folder$Prototype_Img3.jpg)
    dev.off()
    dev.set()
  } else if(i ==4){
    plot(folder$Prototype_Img4.jpg)
    dev.off()
    dev.set()
  } else if(i==5){
    plot(folder$Prototype_Img5.jpg)
    dev.off()
    dev.set()
  }

}

#==============================================================

#upload picture to R for calculating

      img_1 <- image_import(paste(path2, "RAW_IMG\ 1\ .jpg", sep = "/"),
                           which = 1,
                           pattern = NULL,
                           path = NULL,
                           resize = FALSE,
                           plot = TRUE,
                           nrow = NULL,
                           ncol = NULL)
      
      img_2 <- image_import(paste(path2, "RAW_IMG\ 2\ .jpg", sep = "/"),
                           which = 1,
                           pattern = NULL,
                           path = NULL,
                           resize = FALSE,
                           plot = TRUE,
                           nrow = NULL,
                           ncol = NULL)

      img_3 <- image_import(paste(path2, "RAW_IMG\ 3\ .jpg", sep = "/"),
                           which = 1,
                           pattern = NULL,
                           path = NULL,
                           resize = FALSE,
                           plot = TRUE,
                           nrow = NULL,
                           ncol = NULL)

      img_4 <- image_import(paste(path2, "RAW_IMG\ 4\ .jpg", sep = "/"),
                           which = 1,
                           pattern = NULL,
                           path = NULL,
                           resize = FALSE,
                           plot = TRUE,
                           nrow = NULL,
                           ncol = NULL)

      img_5 <- image_import(paste(path2, "RAW_IMG\ 5\ .jpg", sep = "/"),
                           which = 1,
                           pattern = NULL,
                           path = NULL,
                           resize = FALSE,
                           plot = TRUE,
                           nrow = NULL,
                           ncol = NULL)

#==========================================================================
# adjust threshold and analyze picture
  #1
      
      img_1bi <- image_binary(img_1,
      )
      
      
      img_1res <- analyze_objects(img_1,
                                  marker = "id",
                                  index = "BI",
                                  lower_noise = 0.1
      )
      
      #-------------------------------------
      #2
      img_2bi <- image_binary(img_2,
      )
      
      
      img_2res <- analyze_objects(img_2,
                                  marker = "id",
                                  index = "BI",
                                  lower_noise = 0.1
      )
      
      #-------------------------------------
      #3
      img_3bi <- image_binary(img_3,
      )
      
      
      img_3res <- analyze_objects(img_3,
                                  marker = "id",
                                  index = "BI",
                                  lower_noise = 0.1
      )
      
      #-------------------------------------
      #4
      img_4bi <- image_binary(img_4,
      )
      
      
      img_4res <- analyze_objects(img_4,
                                  marker = "id",
                                  index = "BI",
                                  lower_noise = 0.1
      )
      
      #-------------------------------------
      #5
      img_5bi <- image_binary(img_5,
      )
      
      
      img_5res <- analyze_objects(img_5,
                                  marker = "id",
                                  index = "BI",
                                  lower_noise = 0.1
      )
      
      #==========================================================================
      # analyze size and number of all counted objects.
      
      #Setting function
      
      #create function to change cm to um 
      micrometer <- function(a){
        return(a*10000)
      }
      
      # Final analyze
      Anabead <- function(a){
        return(get_measures(a, dpi = 157412.1) %>%
                 micrometer())
      }
      
      result1 <- Anabead(img_1res)
      result2 <- Anabead(img_2res)
      result3 <- Anabead(img_3res)
      result4 <- Anabead(img_4res)
      result5 <- Anabead(img_5res)
      
      #--------------------------------------------
      #filter result because I want only the circular beads
      
      library(tidyverse)
      Filterbead <- function(df){
        return(df %>% filter(circularity_norm > 6000 & circularity > 20000) %>% 
                 select(id, diam_mean, area, circularity, circularity_norm) %>% 
                 arrange(circularity_norm))
      } 
      
      result1_F <- Filterbead(result1)
      result2_F <- Filterbead(result2)
      result3_F <- Filterbead(result3)
      result4_F <- Filterbead(result4)
      result5_F <- Filterbead(result5)
      
      #===============================================================
      # combine picture
      
      combine_diameter <- c(result1_F$diam_mean, result2_F$diam_mean, result3_F$diam_mean,
                            result4_F$diam_mean, result5_F$diam_mean)
      
      combine_cirnorm <- c(result1_F$circularity_norm, result2_F$circularity_norm, result3_F$circularity_norm,
                           result4_F$circularity_norm, result5_F$circularity_norm)
      
      combine_cir <- c(result1_F$circularity, result2_F$circularity, result3_F$circularity,
                       result4_F$circularity, result5_F$circularity)
      
      combine_area <- c(result1_F$area, result2_F$area, result3_F$area,
                        result4_F$area, result5_F$area)
      
      
      combinebead <- data.frame(combine_diameter)
      
      combine_SUM <- data.frame(combine_diameter, combine_area, combine_cir, combine_cirnorm)
      
      #set function for creating HT
  BuildHT_SUM <- ggplot(data = combinebead)+
        aes(x = combine_diameter, y = ..density..) +
        geom_histogram(fill = "azure3") +
        geom_density(alpha = 0.3, fill = "coral1") +
        labs(x="Bead size (um)", y="Density") +
        theme_minimal() +
        theme(axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.position = "bottom")
  
  BuildHT_DEN <- ggplot(data = combinebead)+
    aes(x = combine_diameter, y = ..density..) +
    geom_density(alpha = 0.3, fill = "coral1") +
    labs(x="Bead size (um)", y="Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          legend.position = "bottom")
  
  BuildHT_HT <- ggplot(data = combinebead)+
    aes(x = combine_diameter, y = ..density..) +
    geom_histogram(fill = "azure4") +
    labs(x="Bead size (um)", y="Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          legend.position = "bottom")
      

#======================================================================
    
  #export Histogram_SUM
  png(paste(path2, "Histogram summary.png", sep = "/"))
  plot(BuildHT_SUM)
  dev.off()
  dev.set()
  
  png(paste(path2, "Histogram (Density).png", sep = "/"))
  plot(BuildHT_DEN)
  dev.off()
  dev.set()
  
  png(paste(path2, "Histogram.png", sep = "/"))
  plot(BuildHT_HT)
  dev.off()
  dev.set()
  
  write_csv(combine_SUM, paste(path2, "summarize.csv", sep = "/"))
      
  
      
      
      







