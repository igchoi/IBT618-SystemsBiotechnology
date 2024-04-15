# PROTOTYPE 2

* Right now, we can make analysis 5 microscopic pictures and combine to 1 histograms as picture and csv. file of raw analyzed data
* However, it have some limitation as below
  - we have to input 5 pictures manually and less than or more than 5 pictures is not acceptable.
    
---------
## 1. Setting dpi of picture with known distance. 
* In this case, I used scale in hemocytometer.

``` r
library(pliman)
# set scale from hemocytometer
scale_50um_20X <- image_import("# full path of hemocytometer picture",
                               resize = F,
                               plot = T)
dpi(scale_50um_20X)
# 157412.1

scale_50um_10X <- image_import("# full path of hemocytometer picture",
                               resize = F,
                               plot = T)
dpi(scale_50um_10X)
# 79016.64

```

## 2.Input picture (right now is manually)

```r

#1
img_1 <- image_import("#full path of prototype_Img1",
                      which = 1,
                      pattern = NULL,
                      path = NULL,
                      resize = FALSE,
                      plot = TRUE,
                      nrow = NULL,
                      ncol = NULL)

#2
img_2 <- image_import("#full path of prototype_Img2",
                      which = 1,
                      pattern = NULL,
                      path = NULL,
                      resize = FALSE,
                      plot = TRUE,
                      nrow = NULL,
                      ncol = NULL)

#3
img_3 <- image_import("#full path of prototype_Img3",
                      which = 1,
                      pattern = NULL,
                      path = NULL,
                      resize = FALSE,
                      plot = TRUE,
                      nrow = NULL,
                      ncol = NULL)

#4
img_4 <- image_import("#full path of prototype_Img4",
                      which = 1,
                      pattern = NULL,
                      path = NULL,
                      resize = FALSE,
                      plot = TRUE,
                      nrow = NULL,
                      ncol = NULL)

#5
img_5 <- image_import("#full path of prototype_Img5",
                      which = 1,
                      pattern = NULL,
                      path = NULL,
                      resize = FALSE,
                      plot = TRUE,
                      nrow = NULL,
                      ncol = NULL)

```

## 3. adjusted threshold and analyze picture

``` r
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

```

## 4. analyzed size and number of all counted objects.

```r 
#First, setting function is need

#create function to change cm to um 
micrometer <- function(a){
  return(a*10000)
}

#Final analyze
Anabead <- function(a){
  return(get_measures(a, dpi = 157412.1) %>%
           micrometer())
}

result1 <- Anabead(img_1res)
result2 <- Anabead(img_2res)
result3 <- Anabead(img_3res)
result4 <- Anabead(img_4res)
result5 <- Anabead(img_5res)

```

## 5. filtered results remaining only the circular beads

``` r
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

```

## 6. combined results of each picture
* because many microparticle numbers is needed for histogram creating

``` r

combine_diameter <- c(result1_F$diam_mean, result2_F$diam_mean, result3_F$diam_mean,
                      result4_F$diam_mean, result5_F$diam_mean)


combinebead <- data.frame(combine_diameter)

```

## 7. setting function and creating histogram
* histogram and raw data of analyzed file as .csv will be exported

``` r
BuildHT <- ggplot(data = combinebead)+
           aes(x = combinebead$combine_diameter, y = ..density..) +
           geom_histogram() +
           geom_density(alpha = 0.3, fill = "red") +
           labs(x="Bead size (um)", y="Density") +
           theme_minimal() +
           theme(axis.text.x = element_text(size=10),
                 axis.text.y = element_text(size=10),
                 legend.position = "bottom")


write_csv(combine_SUM, "SUMMARY OF DATA.csv")

getwd()

```


