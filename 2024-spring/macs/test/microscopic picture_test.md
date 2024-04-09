* Microscopic picture acquiring from alginate microparticle production was added to process

![TC_998](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/9317fd64-e9db-4a3e-bb66-e3398034a1b4)

=====================================
# Adding picture to R program
library(pliman)
img_test <- image_import("/Users/mackia/Desktop/TC_998.jpg",
  which = 1,
  pattern = NULL,
  path = NULL,
  resize = FALSE,
  plot = TRUE,
  nrow = NULL,
  ncol = NULL)
  
=====================================
  # Adjusting theshold of picture

img_test2 <- image_binary(img_test)





==============================================
# Identifying microparticles by using bright intensity (BI)
img_res <- analyze_objects(img_test,
                           marker = "id",
                           index = "BI",
                           lower_noise = 1
                           )



=====================================
# gaining results and creating result to histogram

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











=====================================
