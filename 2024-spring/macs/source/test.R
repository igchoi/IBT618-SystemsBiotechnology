img_test <- image_import("picture_path",
                         which = 1,
                         pattern = NULL,
                         path = NULL,
                         resize = FALSE,
                         plot = TRUE,
                         nrow = NULL,
                         ncol = NULL
)

img_test2 <- image_binary(img_test)

img_res <- analyze_objects(img_test,
                           marker = "id",
                           index = "BI",
                           lower_noise = 0.01)

Is_px <- img_res$results$perimeter[1]
pixels_to_cm(px = Is_px, dpi = 100)

result <- get_measures(img_res, dpi = 100)


# Histogram of size

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

#====================================================


