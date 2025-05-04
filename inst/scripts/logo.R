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

sticker(p,
        s_x = 1, s_y = 0.75, s_width = 0.75, s_height = 0.75,
        package = "WData",
        p_color = "white", p_size = 30, p_family = "lobster",
        h_fill = "#CDCDC1", h_color = "#8B8B83",  spotlight = FALSE,
        white_around_sticker = FALSE,
        url = "https://github.com/noeliasanchmrt/WData",
        u_color = "white", u_size = 3.5, u_family = "sans",
        filename = "inst/extdata/WData_logo.png")

usethis::use_logo("inst/extdata/WData_logo.png") # To use the logo in the README
