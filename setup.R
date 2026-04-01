library(devtools)
library(roxygen2)

#### NAME


desc::desc_set("Package", "uavRpheno")



#### AUTOR

usethis::use_author(given = "Marcelo",
                    family = "Bueno",
                    role = c("aut", "cre"),
                    email = "marcelobueno630@gmail.com")

#### DESRIPCION

desc::desc_set("Description", "Performs common high-throughput phenotyping (HTP) tasks using time series of multispectral imagery in the context of agronomic experimental trials, especially in randomized complete block designs (RCBD). The library enables derivation of canopy cover, canopy volume and vegetation indices, modeling of crop phenology using gaussian-based growth curves, and statistical modeling of derived features and field traits.")

### TITULO

# If you don't have it: install.packages("desc")
desc::desc_set("Title",
               "High-Throughput Phenotyping for Agronomic Trials using multispectral UAV images")


# DEPENDENCIES

libs <- c(
  "terra", "sf", "tidyverse", "exactextractr", "openxlsx", "signal",
  "zoo", "stringr", "lubridate", "gslnls", "minpack.lm", "skimr",
  "ggplot2", "lme4", "car", "emmeans", "multcomp", "multcompView",
  "broom.mixed", "MuMIn", "DHARMa", "mgcv"
)

sapply(libs, function(pkg) as.character(packageVersion(pkg)))



# --- Geospatial Dependencies ---
usethis::use_package("terra", min_version = "1.8.21")
usethis::use_package("sf", min_version = "1.0.19")
usethis::use_package("exactextractr", min_version = "0.10.0")

# --- Data Manipulation & Utilities ---


# These are the core engines of the tidyverse you are likely using
usethis::use_package("dplyr", min_version = "1.1.0")
usethis::use_package("tidyr", min_version = "1.3.0")
usethis::use_package("purrr", min_version = "1.0.0")
usethis::use_package("readr", min_version = "2.1.0")

usethis::use_package("openxlsx", min_version = "4.2.8")
usethis::use_package("stringr", min_version = "1.5.1")
usethis::use_package("lubridate", min_version = "1.9.4")
usethis::use_package("zoo", min_version = "1.8.12")
usethis::use_package("skimr", min_version = "2.2.2")
usethis::use_package("signal", min_version = "1.8.1")

# --- Modeling & Nonlinear Regression ---
usethis::use_package("gslnls", min_version = "1.4.1")
usethis::use_package("minpack.lm", min_version = "1.2.4")

# --- Statistical Analysis & Visualization ---
usethis::use_package("ggplot2", min_version = "3.5.1")
usethis::use_package("lme4", min_version = "1.1.36")
usethis::use_package("car", min_version = "3.1.3")
usethis::use_package("emmeans", min_version = "1.10.7")
usethis::use_package("multcomp", min_version = "1.4.28")
usethis::use_package("multcompView", min_version = "0.1.10")
usethis::use_package("broom.mixed", min_version = "0.2.9.6")
usethis::use_package("MuMIn", min_version = "1.48.4")
usethis::use_package("DHARMa", min_version = "0.4.7")

usethis::use_package("mgcv", min_version = "1.9.1")

### LICENCIA

usethis::use_gpl3_license()




#################
################# GIR IGNORE
#################



path <- "D:/DOCS/CIENCIA DE SUELOS/R_SUELOS/QUINOA/QUINOA_GITHUB/HTP_PROCESING/PhenoHTP"

all_files <- list.files(path, recursive = TRUE, full.names = TRUE)

print(all_files)



#################
################# TEST
#################

library(tidyverse)




devtools::document()

devtools::load_all()


test_zonal <- extract_htp_zonal(
  multi_path = 'PROCESING/multi/',
  dsm_path = 'PROCESING/dsm/',
  border_path = 'PROCESING/borders/trail_trigo.gpkg',
  dtm_file = 'PROCESING/dtm/dtm_14.tif',
  save_rasters = FALSE,
  treatment_var = "treatment",
  blocking_var = "blocking",
  daps = c("14", "35", "42", "56", "70","83","104","119","126","145"),
  band_indices = c(B = 1, G = 2, R = 3, RE = 4, NIR = 5, T = 6)
)

# write.csv(test_zonal, 'test_zonal.csv' )

htp_features <- extract_htp_pheno(
  data = test_zonal,
  indices = c("CC", "CV", "ExG", "NDVI"),
  time_var = "DAP"
)

##############


data("traits_trigo")

cor_table <- htp_correlations(traits = traits_trigo,
                              htp_features = htp_features,
                              trait = "yield_kg_plot")

reg_results <- htp_regression(traits = traits_trigo,
                              htp_features = htp_features,
                              trait = "yield_kg_plot")

anova_results <- run_glmer_anova(traits = traits_trigo,
                                 htp_features = htp_features,
                                 trait = 'yield_kg_plot')




###############
############### CHECKING RESULTS
###############

#
# library(ggplot2)
#
# ggplot(test_zonal, aes(x = DAP, y = CC)) +
#   geom_point(aes(color = treatment), alpha = 0.6) +
#   geom_line(aes(group = id), alpha = 0.3) + # Optional: connect points for the same ID
#   facet_grid(treatment ~ blocking, labeller = label_both) +
#   theme_bw() +
#   labs(
#     title = "Canopy Cover Development by Treatment and Block",
#     x = "Days After Planting (DAP)",
#     y = "Canopy Cover (CC)",
#     color = "Treatment"
#   ) +
#   theme(legend.position = "none") # Hide legend as facets already identify treatments
#
#
# sss <- sf::st_read('PROCESING/borders/trail_trigo.gpkg')
# sss %>% View()




#################
################# EXAMPLE DATA
#################

data <- sf::st_read('PROCESING/borders/trail_trigo.gpkg')
# data <- data %>% dplyr::rename(blocking = bloque )
# sf::st_write(data, 'PROCESING/borders/trail_trigo.gpkg',append= FALSE)

data <- data %>% sf::st_drop_geometry()
write.csv(data,'traits_trigo.csv')



traits_trigo <- data
usethis::use_data(traits_trigo, overwrite = TRUE)

library(uavRpheno)
data("traits_trigo")




