
########### Create Package Readme.md -------------------------------------------
rmarkdown::render(paste0("./README.Rmd"),
                  output_dir = "./",
                  output_file = paste0("README.md"))

########### Document and create Package ----------------------------------------
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
install("surveyIndex")
3
setwd(here::here())
# devtools::check()

########### Create Documentation GitHub-Pages ----------------------------------

.rs.restartR()
# devtools::install_github("rstudio/fontawesome", force = T)
# library(fontawesome)
library(here)
library(usethis)
library(pkgdown)

# devtools::install_github("r-lib/pkgdown")
# pkgdown::build_favicons()
# devtools::build_vignettes()
# usethis::use_pkgdown(config_file = "./pkgdown/_pkgdown.yml")

pkgdown::build_site(pkg = here::here())
# usethis::use_github_action("pkgdown")

# Save Package tar.gz
date0 <- "2022.10.01"
file.remove(paste0(dirname(here::here()), "/surveyIndex_",date0,".tar.gz"))
file.remove(paste0((here::here()), "/surveyIndex_",date0,".tar.gz"))
devtools::build()
file.copy(from = paste0(dirname(here::here()), "/surveyIndex_",date0,".tar.gz"),
          to = paste0(here::here(), "/surveyIndex_",date0,".tar.gz"),
          overwrite = TRUE)

