library(hexSticker)
library(ggplot2)

curve_data <- data.frame(x = seq(0, 15, length.out = 100))
curve_data$y <- cos(curve_data$x) * curve_data$x + 10

p <- ggplot(curve_data, aes(x, y)) +
  geom_line(color = "white", size = 0.5) +
  theme_void() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(color = "white",arrow = arrow(length = unit(0.05, "inches"))))

library(showtext)
font_add_google("lobster", family = "Lobster")
showtext_auto()

library(magick)

# https://pixabay.com/vectors/dumbbell-gym-activity-fitness-body-7610150/
img <- image_read("scales.png")

info <- image_info(img)
new_height <- info$height - 635  # recorta 50 px desde abajo
img <- image_crop(img, geometry = paste0(info$width, "x", new_height, "+0+0"))
img_white <- image_colorize(img, opacity = 100, color = "#24693D")
image_write(img_white, path = "white_scale.png", format = "png")


sticker("white_scale.png",
        s_x = 1, s_y = 0.75, s_width = 0.5, s_height = 0.5,
        package = "WData",
        p_color = "#24693D", p_size = 30, p_family = "Lobster",
        h_fill = "#C5D4B9", h_color = "#24693D",  spotlight = FALSE,
        white_around_sticker = FALSE,
        url = "https://github.com/noeliasanchmrt/WData",
        u_color = "#24693D", u_size = 3.5, u_family = "sans",
        filename = "inst/extdata/WData_logo.png")

usethis::use_logo("inst/extdata/WData_logo.png") # To use the logo in the README
