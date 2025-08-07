# Creates shrub.data dataset from https://www.jstor.org/stable/2531448?seq=1&cid=pdf-reference#references_tab_contents
shrub.data <- read.csv("shrub.csv", sep = ";")

# output for the package
usethis::use_data(shrub.data, overwrite = TRUE)
