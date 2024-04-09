## Microparticle analysis_prototype
-------
#### 1. Adding picture to R program

``` r
library(pliman)
img_test <- image_import("/Users/mackia/Desktop/TC_998.jpg",
  which = 1,
  pattern = NULL,
  path = NULL,
  resize = FALSE,
  plot = TRUE,
  nrow = NULL,
  ncol = NULL)

```

![TC_998](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/9317fd64-e9db-4a3e-bb66-e3398034a1b4)

-------
#### 2. Adjusting theshold of picture

``` r
img_test2 <- image_binary(img_test)
``` 

![binary_test](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/99af6bed-3211-4843-aa87-ba0fe9c79023)

-------
#### 3. Identifying microparticles by using bright intensity (BI)

``` r
img_res <- analyze_objects(img_test,
                           marker = "id",
                           index = "BI",
                           lower_noise = 1
                           )

``` 

![Resolution_test](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/01140df9-94e7-4bbd-b7ef-763eef0aab46)

-------
#### 4. Gaining results and creating result to histogram

``` r

result <- get_measures(img_res, dpi = 100)
library(ggplot2)
HT_size <- ggplot(data = result) +
  aes(x = result$perimeter, y = ..density..) +
  geom_histogram() +
  geom_density(alpha = 0.3, fill = "red") +
  labs(x="Bead size (cm)", y="Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.position = "bottom") +
  labs(title = "Size Distribution")
print(HT_size)

``` 

![Histogram_test](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/68a88cce-666d-4c06-8d1a-25719a1e93bb)

