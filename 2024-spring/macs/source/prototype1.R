
#### 1. Adding picture to R program


library(pliman)
img_test <- image_import("# put path of pic_prototype1",
  which = 1,
  pattern = NULL,
  path = NULL,
  resize = FALSE,
  plot = TRUE,
  nrow = NULL,
  ncol = NULL)

#### 2. Adjusting theshold of picture

img_test2 <- image_binary(img_test) 


#### 3. Identifying microparticles by using bright intensity (BI)

img_res <- analyze_objects(img_test,
                           marker = "id",
                           index = "BI",
                           lower_noise = 1
                           )

#### 4. Gaining results and creating result to histogram


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

print(HT_size)

