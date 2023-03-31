

# Document and create Package --------------------------------------------------
.rs.restartR()

# options(rmarkdown.html_vignette.check_title = FALSE)
Sys.setenv('PATH' = paste0('C:/Program Files/qpdf-10.3.1/bin;', Sys.getenv('PATH')))
library(here)
library(devtools)
library(usethis)
library(roxygen2)
library(RODBC)
devtools::document()
setwd("..")
install("sdmgamindex")
3
setwd(here::here())
# devtools::check()

# Create Package Readme.md -----------------------------------------------------
rmarkdown::render(paste0("./README.Rmd"),
                  output_dir = "./",
                  output_file = "README.md")

# Add stuff --------------------------------------------------------------------

# usethis::use_vignette("C-model-comparisons")

# Create Documentation GitHub-Pages --------------------------------------------

.rs.restartR()
# devtools::install_github("rstudio/fontawesome", force = T)
# library(fontawesome)
library(here)
library(usethis)
library(pkgdown)

# Run once to configure package to use pkgdown
# usethis::use_pkgdown()
# Run to build the website
pkgdown::build_site()


# If you’re using GitHub, we also recommend setting up GitHub actions to automatically build and publish your site:
# usethis::use_pkgdown_github_pages()

# Save Package tar.gz
version0 <- "0.1.0"
file.remove(paste0(dirname(here::here()), "/sdmgamindex_",version0,".tar.gz"))
file.remove(paste0((here::here()), "/sdmgamindex_",version0,".tar.gz"))
devtools::build()
file.copy(from = paste0(dirname(here::here()), "/sdmgamindex_",version0,".tar.gz"),
          to = paste0(here::here(), "/sdmgamindex_",version0,".tar.gz"),
          overwrite = TRUE)

