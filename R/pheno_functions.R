
#### MASK
#' Calculate Canopy Coverage (CC) Pixel Map
#'
#' Evaluates individual pixels to determine if they represent vegetation (canopy)
#' based on RGB ratios and optional intensity thresholds.
#'
#' @param r Numeric vector or matrix. Red channel values.
#' @param g Numeric vector or matrix. Green channel values.
#' @param b Numeric vector or matrix. Blue channel values.
#' @param p1 Numeric. Threshold for the Red/Green ratio. Default is 0.95.
#' @param p2 Numeric. Threshold for the Blue/Green ratio. Default is 0.95.
#' @param p3 Numeric. Minimum threshold for the Excess Green index (2*G - R - B). Default is 20.
#' @param p4 Numeric (Optional). Total RGB intensity threshold (R + G + B).
#'   Must be used in conjunction with \code{p5}.
#' @param p5 Numeric (Optional). Threshold for the Green-Blue difference (G - B).
#'   Must be used in conjunction with \code{p4}.
#'
#' @details
#' The function identifies canopy pixels using two main conditions:
#' \itemize{
#'   \item \strong{Condition 1 (Primary):} Based on the dominance of the green channel relative
#'   to red and blue, and a minimum "Excess Green" value (\code{p3}).
#'   \item \strong{Condition 2 (Secondary):} Only active if \code{p4} and \code{p5} are provided.
#'   This acts as an additional filter for specific lighting or intensity conditions.
#' }
#' Any \code{NA} values resulting from the calculation are coerced to \code{0} (non-canopy).
#'
#' @return A logical vector or matrix (same dimensions as input) where \code{TRUE}
#' indicates a canopy pixel and \code{FALSE} indicates background.
#'
#' @export
#'
#' @examples
#' # Basic usage with default parameters
#' r <- c(100, 50); g <- c(120, 150); b <- c(90, 40)
#' calculate_cc_pixel(r, g, b)
#'
#' # Using the secondary intensity condition
#' calculate_cc_pixel(r, g, b, p4 = 300, p5 = 50)
calculate_cc_pixel <- function(r, g, b, p1=0.95, p2=0.95, p3=20, p4=NULL, p5=NULL) {
  cond1 <- (r / g < p1) & (b / g < p2) & ((2 * g - r - b) > p3)
  cond2 <- FALSE  # <--- CHANGED FROM 0 TO FALSE
  if (!is.null(p4) && !is.null(p5)) {
    cond2 <- ((r + g + b) < p4) & ((g - b) > p5)
  }
  cc_map <- (cond1 | cond2)
  cc_map[is.na(cc_map)] <- 0
  return(cc_map)
}



#' Extract Canopy Cover Percentage from a Binary Mask
#'
#' This function calculates the percentage of canopy cover within given polygons
#' based on a binary canopy mask (where 1 represents canopy and 0 represents non-canopy).
#' It uses the \code{exactextractr} package for efficient extraction.
#'
#' @param cc_mask A SpatRaster or RasterLayer object. This should be a binary mask
#'   where cells with value 1 represent canopy cover.
#' @param polygons An sf or SpatVector object containing the boundaries for extraction.
#' @param id_col A character string specifying the column name in \code{polygons}
#'   to be used as a unique identifier. Default is "id".
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item \code{COD}: The identifier from the original polygon data.
#'   \item \code{CC_pct}: The calculated canopy cover percentage (0-100).
#' }
#'
#' @export
#'
#' @examples
#' # results <- extract_canopy_cover(my_mask, my_buffer_zones, id_col = "plot_id")
extract_canopy_cover <- function(cc_mask, polygons, id_col = "id") {
  if (!id_col %in% names(polygons)) {
    stop(paste("Column", id_col, "not found in polygon data!"))
  }

  stats <- exactextractr::exact_extract(cc_mask, polygons, function(values, coverage_fraction) {
    sum(values == 1, na.rm = TRUE) / length(values)
  }, progress = FALSE)

  results <- data.frame(
    COD = polygons[[id_col]],
    CC_pct = stats * 100
  )

  return(results)
}



#' Extract Canopy Volume from DSM and DTM
#'
#' This function calculates the total canopy volume (m³) within specified polygons.
#' It generates a Canopy Height Model (CHM) by subtracting the DTM from the DSM,
#' applies a height threshold, masks by canopy cover (optional), and integrates
#' pixel area to compute volume.
#'
#' @param dsm A SpatRaster representing the Digital Surface Model (top of canopy).
#' @param dtm A SpatRaster representing the Digital Terrain Model (ground surface).
#' @param polygons An sf or SpatVector object defining the areas of interest.
#' @param id_col A character string for the unique identifier column in \code{polygons}.
#'   Default is "id".
#' @param cc_mask An optional binary SpatRaster (1 = canopy, 0 = non-canopy).
#'   If provided, heights in non-canopy areas are set to 0. Default is NULL.
#'
#' @details
#' The function filters out noise by setting heights below 0.05m to zero.
#' Volume is calculated as: \eqn{Volume = \sum (Height \times Pixel Area)}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{COD}: The identifier from the polygon data.
#'   \item \code{Volume_m3}: The total estimated canopy volume in cubic meters.
#' }
#'
#' @export
extract_cv <- function(dsm, dtm, polygons, id_col = "id", cc_mask = NULL) {

  chm <- dsm - dtm
  chm[chm < 0.05] <- 0
  chm <- terra::mask(chm, cc_mask, maskvalues = 0, updatevalue = 0)
  pixel_area <- terra::res(chm)[1] * terra::res(chm)[2]
  vol_raster <- chm * pixel_area
  vol_sums <- exactextractr::exact_extract(vol_raster, polygons, 'sum')
  results <- data.frame(
    COD = polygons[[id_col]],
    Volume_m3 = vol_sums
  )
  return(results)
}


#' Extract Mean Excess Green Index (ExG)
#'
#' This function calculates the Excess Green Index (ExG) from RGB bands,
#' masks it to canopy areas, and extracts the mean value within specified polygons.
#'
#' @param red A SpatRaster or RasterLayer representing the Red channel.
#' @param green A SpatRaster or RasterLayer representing the Green channel.
#' @param blue A SpatRaster or RasterLayer representing the Blue channel.
#' @param cc_mask A binary mask SpatRaster (1 = canopy, 0 = non-canopy) used
#'   to restrict calculation to vegetation pixels.
#' @param polygons An sf or SpatVector object defining the areas for extraction.
#' @param id_col A character string for the unique identifier column in \code{polygons}.
#'   Default is "id".
#'
#' @details
#' The Excess Green Index is calculated using normalized RGB values:
#' \itemize{
#'   \item \eqn{r = R / (R + G + B)}
#'   \item \eqn{g = G / (R + G + B)}
#'   \item \eqn{b = B / (R + G + B)}
#'   \item \eqn{ExG = 2g - r - b}
#' }
#' This index is particularly effective for separating green vegetation from soil
#' and residue backgrounds.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{COD}: The identifier from the polygon data.
#'   \item \code{ExG}: The mean Excess Green Index value per polygon.
#' }
#'
#' @export
extract_exg_mean <- function(red, green, blue, cc_mask, polygons, id_col = "id") {

  total <- red + green + blue
  r_norm <- red / total
  g_norm <- green / total
  b_norm <- blue / total

  exg <- 2 * g_norm - r_norm - b_norm

  exg_masked <- terra::mask(exg, cc_mask, maskvalues = 0)

  mean_vals <- exactextractr::exact_extract(exg_masked, polygons, 'mean')

  results <- data.frame(
    COD = polygons[[id_col]],
    ExG = mean_vals
  )

  return(results)
}

#' Extract Mean Normalized Difference Vegetation Index (NDVI)
#'
#' This function calculates the NDVI from Red and Near-Infrared (NIR) bands,
#' masks the result to canopy areas, and extracts the mean value within specified polygons.
#'
#' @param red A SpatRaster or RasterLayer representing the Red channel (usually ~660nm).
#' @param nir A SpatRaster or RasterLayer representing the Near-Infrared channel (usually ~800nm).
#' @param cc_mask A binary mask SpatRaster (1 = canopy, 0 = non-canopy) used
#'   to restrict calculation to vegetation pixels.
#' @param polygons An sf or SpatVector object defining the areas (e.g., plots or tree crowns)
#'   for extraction.
#' @param id_col A character string for the unique identifier column in \code{polygons}.
#'   Default is "id".
#'
#' @details
#' The NDVI is calculated as:
#' \deqn{NDVI = \frac{NIR - Red}{NIR + Red}}
#' Values typically range from -1 to 1, where high positive values indicate
#' dense, healthy vegetation. By masking with \code{cc_mask}, the function
#' ensures the mean represents only the canopy pixels.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{COD}: The identifier from the polygon data.
#'   \item \code{NDVI}: The mean NDVI value per polygon.
#' }
#'
#' @export
extract_ndvi_mean <- function(red, nir, cc_mask, polygons, id_col = "id") {

  ndvi <- (nir - red) / (nir + red)

  ndvi_masked <- terra::mask(ndvi, cc_mask, maskvalues = 0)

  mean_vals <- exactextractr::exact_extract(ndvi_masked, polygons, 'mean')

  results <- data.frame(
    COD = polygons[[id_col]],
    NDVI = mean_vals
  )

  return(results)
}


########### WRAPPER

#' Full Remote Sensing Zonal Extraction Pipeline
#'
#' This function automates the processing of multi-temporal remote sensing data.
#' For each provided Day After Planting (DAP), it aligns DSMs and DTMs, computes
#' Canopy Cover (CC), Canopy Volume (CV), and spectral indices (ExG, NDVI),
#' and performs zonal extraction based on provided plot boundaries.
#'
#' @param multi_path Character string. Path to the directory containing multispectral .tif files.
#' @param dsm_path Character string. Path to the directory containing DSM .tif files.
#' @param border_path Optional; Character string. Path to a single shapefile/geojson
#'   used for all DAPs.
#' @param multiborder_path Optional; Character string. Path to a directory with
#'   DAP-specific border files. Used if \code{border_path} is NULL.
#' @param dtm_file Character string. Path to the static DTM file.
#' @param save_rasters Logical. If TRUE, saves intermediate rasters (CC, CV, ExG, NDVI) to disk.
#' @param output_path Character string. Root directory where results and rasters will be saved.
#' @param treatment_var Character string. The column name in the polygon data
#'   representing the treatment/unique ID. Default is "treatment".
#' @param blocking_var Character string. The column name for blocking/replication.
#'   Default is "blocking".
#' @param daps A vector (numeric or character) of Days After Planting to process.
#' @param band_indices A named vector defining the band order in the multispectral files.
#'   Default is \code{c(B=1, G=2, R=3, RE=4, NIR=5, T=6)}.
#'
#' @details
#' The function expects files to follow a specific naming convention ending in
#' \code{_DAP.tif} (e.g., \code{imagery_45.tif}). It automatically handles spatial
#' alignment by resampling the DSM and DTM to match the multispectral resolution
#' and extent using bilinear interpolation.
#'
#' @return A long-format data frame (tibble) containing extracted metrics (CV, ExG,
#'   NDVI, CC) for every plot and every DAP processed.
#' @export
extract_htp_zonal <- function(multi_path       = NULL,
                              dsm_path         = NULL,
                              border_path      = NULL,
                              multiborder_path = NULL,
                              dtm_file         = NULL,
                              save_rasters     = FALSE,
                              output_path      = NULL,
                              treatment_var    = "treatment",
                              blocking_var     = "blocking",
                              daps,
                              band_indices     = c(B=1, G=2, R=3, RE=4, NIR=5, T=6)) {

  message("--- Starting Full Remote Sensing Pipeline ---")

  dtm_raw <- terra::rast(fs::path(dtm_file))

  static_poly <- NULL
  if (!is.null(border_path)) {
    static_poly <- sf::st_read(as.character(border_path), quiet = TRUE)
    message("    > Using static border file for all DAPs")
  }

  if (save_rasters && !is.null(output_path)) {
    products <- c("CC", "CV", "ExG", "NDVI")
    for(prod in products) {
      fs::dir_create(fs::path(output_path, prod))
    }
  }

  all_results <- purrr::map_dfr(daps, function(dap) {

    message(paste("--- Processing DAP:", dap, "---"))

    pattern <- paste0("_", dap, "\\.tif$")
    m_file  <- fs::dir_ls(multi_path, regexp = pattern)[1]
    s_file  <- fs::dir_ls(dsm_path, regexp = pattern)[1]

    if (any(length(c(m_file, s_file)) == 0) || any(is.na(c(m_file, s_file)))) {
      warning(paste("Missing raster files for DAP", dap))
      return(NULL)
    }

    r_list <- terra::rast(m_file)
    dsm_raw <- terra::rast(s_file)

    if (!terra::compareGeom(r_list, dsm_raw, stopOnError = FALSE)) {
      dsm <- terra::resample(dsm_raw, r_list, method = "bilinear")
    } else {
      dsm <- dsm_raw
    }

    if (!terra::compareGeom(r_list, dtm_raw, stopOnError = FALSE)) {
      dtm_ref <- terra::resample(dtm_raw, r_list, method = "bilinear")
    } else {
      dtm_ref <- dtm_raw
    }

    if (!is.null(static_poly)) {
      poly_data <- static_poly
    } else {
      b_file <- fs::dir_ls(multiborder_path, regexp = pattern)[1]
      if (length(b_file) == 0 || is.na(b_file)) {
        warning(paste("No specific border file found for DAP", dap))
        return(NULL)
      }
      poly_data <- sf::st_read(as.character(b_file), quiet = TRUE)
    }

    b   <- r_list[[band_indices["B"]]]
    g   <- r_list[[band_indices["G"]]]
    r   <- r_list[[band_indices["R"]]]
    nir <- r_list[[band_indices["NIR"]]]

    cc_mask <- calculate_cc_pixel(r * 255, g * 255, b * 255, p1 = 1.05, p2 = 1.00, p3 = 5)

    chm <- dsm - dtm_ref
    chm[chm < 0.05] <- 0
    cv_raster <- terra::mask(chm, cc_mask, maskvalues = 0, updatevalue = NA)

    total <- r + g + b
    exg   <- 2 * (g/total) - (r/total) - (b/total)
    exg_masked <- terra::mask(exg, cc_mask, maskvalues = 0, updatevalue = NA)

    ndvi <- (nir - r) / (nir + r)
    ndvi_masked <- terra::mask(ndvi, cc_mask, maskvalues = 0, updatevalue = NA)


    if (save_rasters && !is.null(output_path)) {
      terra::writeRaster(cc_mask,   fs::path(output_path, "CC",   paste0("CC_DAP_", dap, ".tif")),   overwrite = TRUE)
      terra::writeRaster(cv_raster, fs::path(output_path, "CV",   paste0("CV_DAP_", dap, ".tif")),   overwrite = TRUE)
      terra::writeRaster(exg_masked,       fs::path(output_path, "ExG",  paste0("ExG_DAP_", dap, ".tif")),  overwrite = TRUE)
      terra::writeRaster(ndvi_masked,      fs::path(output_path, "NDVI", paste0("NDVI_DAP_", dap, ".tif")), overwrite = TRUE)
    }

    cv_res   <- extract_cv(dsm, dtm_ref, poly_data, treatment_var, cc_mask)
    exg_res  <- extract_exg_mean(r, g, b, cc_mask, poly_data, treatment_var)
    ndvi_res <- extract_ndvi_mean(r, nir, cc_mask, poly_data, treatment_var)
    cc_res   <- extract_canopy_cover(cc_mask, poly_data, treatment_var)

    dap_combined <- data.frame(
      id        = if("id" %in% colnames(poly_data)) poly_data$id else seq_len(nrow(poly_data)),
      treatment = poly_data[[treatment_var]],
      blocking  = if(blocking_var %in% colnames(poly_data)) poly_data[[blocking_var]] else NA,
      CV        = cv_res$Volume_m3,
      ExG       = exg_res$ExG,
      NDVI      = ndvi_res$NDVI,
      CC        = cc_res$CC_pct,
      DAP       = as.numeric(dap)
    )

    return(dap_combined)
  })

  message("--- Pipeline Complete ---")
  return(all_results)
}


##### PLOT RASTERS

#' Plot HTP Geospatial Products Grid
#'
#' @description
#' Generates a multi-panel grid visualization of RGB orthomosaics and phenotypic
#' index rasters (CC, CV, NDVI, ExG) for a specific experimental plot across
#' multiple time points (DAPs).
#'
#' @param output_path Path to the root directory containing processed index folders.
#' @param multi_path Path to the directory containing multispectral/RGB TIFF files.
#' @param treatment Character string identifying the treatment group to filter.
#' @param blocking Numeric or character identifying the block/rep to filter.
#' @param border_file Path to the geopackage (.gpkg) containing plot boundaries.
#' @param daps Character vector of Days After Planting to include in the columns.
#' @param plot_path Directory path where the final PNG image will be saved.
#'
#' @return A patchwork ggplot object containing the arranged raster visualizations.
#' @export
plot_htp<- function(output_path = "OUTPUT/",
                    multi_path  = "PROCESING/MULTI/",
                    treatment   = "N50",
                    blocking    = 1,
                    border_file = "PROCESING/borders/trail_trigo.gpkg",
                    daps        = c("14", "104", "126"),
                    plot_path   = "OUTPUT/PLOTS/") {

  boundaries <- st_read(border_file, quiet = TRUE)

  plot_poly <- boundaries %>%
    dplyr::filter(treatment == !!treatment, blocking == as.numeric(!!blocking))

  if(nrow(plot_poly) == 0) stop("No plot found for specified treatment/blocking.")

  v_poly <- vect(plot_poly)
  plot_list <- list()

  band_indices <- c(B = 1, G = 2, R = 3, RE = 4, NIR = 5, T = 6)
  rgb_idx <- c(band_indices["R"], band_indices["G"], band_indices["B"])

  for (dap in daps) {

    pattern <- paste0("_", dap, "\\.tif$")
    rgb_file <- fs::dir_ls(multi_path, regexp = pattern)[1]

    if (length(rgb_file) == 0 || is.na(rgb_file)) {
      warning(paste("RGB file missing for DAP", dap))
      next
    }

    rgb_raw  <- rast(rgb_file)[[rgb_idx]]
    rgb_crop <- crop(rgb_raw, v_poly) %>% mask(v_poly)
    rgb_bright <- stretch(rgb_crop, minv=0, maxv=255, minq=0.02, maxq=0.98)

    find_idx <- function(type, d) {
      dir_path <- fs::path(output_path, type)
      if(!dir.exists(dir_path)) return(NULL)
      fs::dir_ls(dir_path, regexp = paste0("_", d, "\\.tif$"))[1]
    }

    idx_files <- list(
      NDVI = find_idx("NDVI", dap),
      CC   = find_idx("CC", dap),
      CV   = find_idx("CV", dap),
      ExG  = find_idx("ExG", dap)
    )

    if (any(sapply(idx_files, is.null)) || any(sapply(idx_files, function(x) length(x) == 0 || is.na(x)))) {
      warning(paste("One or more index TIFFs missing for DAP", dap))
      next
    }

    ndvi_r <- rast(idx_files$NDVI) %>% crop(v_poly) %>% mask(v_poly)
    cc_r   <- rast(idx_files$CC)   %>% crop(v_poly) %>% mask(v_poly)
    cv_r   <- rast(idx_files$CV)   %>% crop(v_poly) %>% mask(v_poly)
    exg_r  <- rast(idx_files$ExG)  %>% crop(v_poly) %>% mask(v_poly)

    p_rgb  <- ggplot() + geom_spatraster_rgb(data = rgb_bright) + theme_void() + labs(title = paste("DAP", dap))
    p_cc   <- ggplot() + geom_spatraster(data = cc_r) + scale_fill_gradient(low="white", high="darkgreen", na.value=NA) + theme_void() + theme(legend.position="none")
    p_cv   <- ggplot() + geom_spatraster(data = cv_r) + scale_fill_viridis_c(option = "magma", na.value=NA) + theme_void() + theme(legend.position="none")
    p_ndvi <- ggplot() + geom_spatraster(data = ndvi_r) + scale_fill_gradientn(colors = rev(terrain.colors(10)), na.value=NA) + theme_void() + theme(legend.position="none")
    p_exg  <- ggplot() + geom_spatraster(data = exg_r) + scale_fill_gradient2(low="brown", mid="white", high="forestgreen", na.value=NA) + theme_void() + theme(legend.position="none")

    if (dap == daps[1]) {
      p_rgb  <- p_rgb  + ylab("RGB")  + theme(axis.title.y = element_text(angle=90, size=12))
      p_cc   <- p_cc   + ylab("CC")   + theme(axis.title.y = element_text(angle=90, size=12))
      p_cv   <- p_cv   + ylab("CV")   + theme(axis.title.y = element_text(angle=90, size=12))
      p_ndvi <- p_ndvi + ylab("NDVI") + theme(axis.title.y = element_text(angle=90, size=12))
      p_exg  <- p_exg  + ylab("ExG")  + theme(axis.title.y = element_text(angle=90, size=12))
    }

    plot_list[[paste0(dap, "_1")]] <- p_rgb
    plot_list[[paste0(dap, "_2")]] <- p_cc
    plot_list[[paste0(dap, "_3")]] <- p_cv
    plot_list[[paste0(dap, "_4")]] <- p_ndvi
    plot_list[[paste0(dap, "_5")]] <- p_exg
  }

  if(length(plot_list) == 0) stop("No data found to plot.")

  final_plot <- wrap_plots(plot_list, ncol = length(daps), byrow = FALSE) +
    plot_annotation(title = paste(treatment, "-", blocking)) &
    theme(plot.title = element_text(hjust = 0.5, size = 12))

  if(!is.null(plot_path)){
    dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(plot_path, paste0("geoProducts_", treatment, "_", blocking , ".png")),
           plot = final_plot, width = 10, height = 12, limitsize = FALSE)
  }

  return(final_plot)
}



###################### UTILS FITING


rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2, na.rm = TRUE))
}

r2 <- function(observed, predicted) {
  1 - sum((observed - predicted)^2, na.rm = TRUE) /
    sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
}



#' Extract High-Throughput Phenotyping (HTP) Metrics
#'
#' @description
#' Fits cubic spline functions (via GAMs) and calculates growth-rate derivatives to
#' extract 22 phenotypic features (F1 to F22) that characterize the full crop
#' life cycle across multiple vegetation indices.
#'
#' @details
#' The function extracts the following features:
#' \itemize{
#'   \item \bold{CC (Canopy Cover):} F1 (Max CC), F2 (Max growth rate), F3 (DAP at max growth rate), F4 (Duration over half maximum).
#'   \item \bold{CV (Canopy Volume):} F5 (Max CV), F6 (Max growth rate), F7 (DAP at max growth rate), F8 (Duration over half maximum).
#'   \item \bold{ExG (Excess Green):} F9 (Max ExG), F10 (DAP at max ExG), F11 (Max increasing slope), F12 (Value at peak), F13 (Duration of increasing period), F14 (Area of increasing period), F15 (Max decreasing slope), F16 (Value at peak), F17 (Duration of decreasing period), F18 (Area of decreasing period).
#'   \item \bold{NDVI:} F19 (Max NDVI), F20 (DAP at max NDVI), F21 (Increasing slope from first obs to peak), F22 (Decreasing slope from peak to last obs).
#' }
#'
#' @param data A data frame containing time-series vegetation indices.
#' @param indices Character vector of column names to process. Defaults to CC, CV, ExG, and NDVI.
#' @param group_cols Character vector of columns used to define unique experimental units (e.g., treatment and block).
#' @param time_var Character string specifying the time variable (e.g., "DAP").
#' @param k Integer specifying the basis dimension for the GAM spline. Defaults to 5.
#'
#' @return A data frame containing the original grouping columns and the calculated features F1 through F22.
#' @export
extract_htp_pheno <- function(data,
                              indices = c("CC", "CV", "ExG", "NDVI"),
                              group_cols = c("treatment", "blocking"),
                              time_var = "DAP",
                              k = 3) {

  groups <- unique(data[, group_cols])
  results_list <- list()

  for(i in 1:nrow(groups)) {
    # Dynamically filter for the current group
    sub_data <- data
    for(col in group_cols) {
      sub_data <- sub_data[sub_data[[col]] == groups[i, col], ]
    }

    group_info <- groups[i, , drop = FALSE]
    group_metrics <- group_info

    for(idx in indices) {
      current_sub <- sub_data[!is.na(sub_data[[idx]]), ]
      if(nrow(current_sub) < (k + 1)) next

      form <- as.formula(paste(idx, "~ s(", time_var, ", k =", k, ")"))

      model <- tryCatch({
        mgcv::gam(form, data = current_sub)
      }, error = function(e) return(NULL))

      if(is.null(model)) next

      t_seq <- seq(min(current_sub[[time_var]]), max(current_sub[[time_var]]), length.out = 1000)
      dt <- t_seq[2] - t_seq[1]

      pred_df <- data.frame(setNames(list(t_seq), time_var))
      y <- as.numeric(predict(model, newdata = pred_df))

      d1 <- diff(y) / dt
      d1 <- c(d1[1], d1)

      idx_max_y <- which.max(y)
      max_val <- y[idx_max_y]
      time_at_max <- t_seq[idx_max_y]

      inc_mask <- t_seq <= time_at_max
      dec_mask <- t_seq > time_at_max

      if(idx == "CC") {
        is_above_half <- y >= (max_val / 2)
        group_metrics$F1 <- max_val
        group_metrics$F2 <- max(d1, na.rm=TRUE)
        group_metrics$F3 <- t_seq[which.max(d1)]
        group_metrics$F4 <- if(any(is_above_half)) diff(range(t_seq[is_above_half])) else 0

      } else if(idx == "CV") {
        is_above_half <- y >= (max_val / 2)
        group_metrics$F5 <- max_val
        group_metrics$F6 <- max(d1, na.rm=TRUE)
        group_metrics$F7 <- t_seq[which.max(d1)]
        group_metrics$F8 <- if(any(is_above_half)) diff(range(t_seq[is_above_half])) else 0

      } else if(idx == "ExG") {
        group_metrics$F9  <- max_val
        group_metrics$F10 <- time_at_max
        group_metrics$F11 <- max(d1[inc_mask], na.rm=TRUE)
        group_metrics$F12 <- max_val
        group_metrics$F13 <- sum(inc_mask) * dt
        group_metrics$F14 <- sum(y[inc_mask]) * dt
        group_metrics$F15 <- min(d1[dec_mask], na.rm=TRUE)
        group_metrics$F16 <- max_val
        group_metrics$F17 <- sum(dec_mask) * dt
        group_metrics$F18 <- sum(y[dec_mask]) * dt

      } else if(idx == "NDVI") {
        y_first <- y[1]
        y_last  <- y[length(y)]
        t_first <- t_seq[1]
        t_last  <- t_seq[length(t_seq)]

        group_metrics$F19 <- max_val
        group_metrics$F20 <- time_at_max
        group_metrics$F21 <- (max_val - y_first) / (time_at_max - t_first)
        group_metrics$F22 <- (y_last - max_val) / (t_last - time_at_max)
      }
    }
    results_list[[i]] <- group_metrics
  }
  return(do.call(dplyr::bind_rows, results_list))
}



#' Plot HTP Fitted Functions with Performance Metrics
#'
#' @description
#' Visualizes the observed time-series data against fitted GAM splines for a
#' specific vegetation index. Calculates and displays R2 and RMSE for each
#' experimental unit within a faceted grid.
#'
#' @param data A data frame containing time-series vegetation indices.
#' @param indices Character string specifying the column name to plot (e.g., "NDVI").
#' @param treatment_var Character string for the column used as the row facet (e.g., "treatment").
#' @param blocking_var Character string for the column used as the column facet (e.g., "blocking").
#' @param time_var Character string specifying the time variable (e.g., "DAP").
#' @param k Integer specifying the basis dimension for the GAM spline. Defaults to 3.
#'
#' @return A ggplot object showing observed points, fitted lines, and accuracy metrics.
#' @export
plot_htp_fit <- function(data,
                         indices = "NDVI",
                         treatment_var = "treatment",
                         blocking_var = "blocking",
                         time_var = "DAP",
                         k = 3) {

  rmse_func <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
  r2_func <- function(obs, pred) {
    1 - sum((obs - pred)^2, na.rm = TRUE) /
      sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  }

  data$plot_group <- paste(data[[treatment_var]], data[[blocking_var]], sep = "_")

  all_preds <- list()
  all_stats <- list()
  groups <- unique(data$plot_group)

  for(g in groups) {
    sub_data <- data[data$plot_group == g & !is.na(data[[indices]]), ]
    if(nrow(sub_data) < (k + 1)) next

    form <- as.formula(paste(indices, "~ s(", time_var, ", k =", k, ")"))
    model <- tryCatch(gam(form, data = sub_data), error = function(e) NULL)

    if(!is.null(model)) {
      t_seq <- seq(min(sub_data[[time_var]]), max(sub_data[[time_var]]), length.out = 100)
      pred_line <- data.frame(time = t_seq)
      colnames(pred_line) <- time_var
      pred_line$fitted_val <- as.numeric(predict(model, newdata = pred_line))
      pred_line$plot_group <- g
      pred_line[[treatment_var]] <- sub_data[1, treatment_var]
      pred_line[[blocking_var]] <- sub_data[1, blocking_var]

      obs_pred <- as.numeric(predict(model, newdata = sub_data))
      m_rmse <- rmse_func(sub_data[[indices]], obs_pred)
      m_r2   <- r2_func(sub_data[[indices]], obs_pred)

      stat_df <- data.frame(
        treatment = sub_data[1, treatment_var],
        blocking = sub_data[1, blocking_var],
        label = paste0("R2: ", round(m_r2, 3), "\nRMSE: ", round(m_rmse, 4)),
        stringsAsFactors = FALSE
      )
      colnames(stat_df)[1:2] <- c(treatment_var, blocking_var)

      all_preds[[g]] <- pred_line
      all_stats[[g]] <- stat_df
    }
  }

  plot_data <- do.call(rbind, all_preds)
  stats_data <- do.call(rbind, all_stats)

  ggplot(data, aes_string(x = time_var, y = indices)) +
    geom_point(color = "grey40", alpha = 0.8) +
    geom_line(data = plot_data, aes(y = fitted_val), color = "blue", size = 0.8) +
    geom_text(data = stats_data, aes(label = label),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
              size = 1.5, check_overlap = TRUE) +
    facet_grid(as.formula(paste(treatment_var, "~", blocking_var))) +
    labs(title = paste("Fit:", indices)) +
    theme(plot.title = element_text(hjust = 0.5))
}



###################### UTILS REGRESION MODELING


##### CORRELATION


#' Calculate Correlations Between HTP Features and Agronomic Traits
#'
#' This function joins phenological features (extracted via Gaussian modeling)
#' with ground-truth data (e.g., yield, biomass) and calculates Pearson
#' correlation coefficients for every feature.
#'
#' @param df A data frame containing ground-truth agronomic data and joining keys.
#' @param htp_features A data frame containing the extracted "F" features
#'   (F1, F2, etc.) and joining keys.
#' @param trait Character. The name of the column in \code{df} to correlate
#'   against (e.g., "yield_kg_plot").
#' @param join_keys A named character vector for joining \code{df} and
#'   \code{htp_features}. Default matches "treatment" and "blocking" columns.
#'
#' @details
#' The function automatically filters out features with zero variance (constant values)
#' and requires at least 3 complete observations to perform a correlation test.
#' Statistical significance is categorized using standard p-value notation:
#' \itemize{
#'   \item \code{***}: p < 0.001
#'   \item \code{**}: p < 0.01
#'   \item \code{*}: p < 0.05
#'   \item \code{.}: p < 0.1
#'   \item \code{ns}: non-significant
#' }
#'
#' @return A data frame sorted by the absolute value of the correlation coefficient,
#'   containing the variable name, correlation value, p-value, and significance stars.
#' @export
htp_correlations <- function(traits, htp_features, trait = "yield_kg_plot",
                             join_keys = c("treatment", "blocking")) {

  combined_df <- traits %>%
    dplyr::left_join(htp_features, by = join_keys)

  feat_names <- combined_df %>%
    dplyr::select(dplyr::starts_with("F")) %>%
    names()

  feat_names <- feat_names[vapply(combined_df[feat_names], function(x) sd(x, na.rm = TRUE) > 0, logical(1))]

  if(length(feat_names) == 0) return(NULL)

  cor_results <- purrr::map_dfr(feat_names, function(f) {
    if (!f %in% names(combined_df)) return(NULL)

    test <- tryCatch({
      cor.test(combined_df[[f]], combined_df[[trait]], use = "pairwise.complete.obs")
    }, error = function(e) return(NULL))

    if(is.null(test)) return(NULL)

    data.frame(
      Variable = f,
      Correlation = round(as.numeric(test$estimate), 3),
      P_Value = round(test$p.value, 4),
      Significance = dplyr::case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01  ~ "**",
        test$p.value < 0.05  ~ "*",
        test$p.value < 0.1   ~ ".",
        TRUE                 ~ "ns"
      )
    )
  })

  return(cor_results %>% dplyr::arrange(desc(abs(Correlation))))
}



#### REGRESION

#' This function performs multiple linear regression to predict an agronomic trait
#' using extracted HTP features. It handles data cleaning, scaling, and compares
#' a Full Model against Forward and Backward stepwise selection models.
#'
#' @param df A data frame containing ground-truth agronomic data.
#' @param htp_df A data frame containing HTP-derived features (F-metrics).
#' @param trait Character. The target response variable (e.g., "PMSEM").
#'   Note: The function applies a log10(x + 1e-6) transformation to this variable.
#' @param na_threshold Numeric (0-1). The maximum proportion of missing values
#'   allowed for a feature to be included in the model. Default is 0.3.
#' @param join_keys A named character vector for joining \code{df} and \code{htp_df}.
#'
#' @details
#' The function automates several preprocessing steps:
#' \itemize{
#'   \item \bold{Feature Filtering:} Removes columns exceeding the \code{na_threshold}.
#'   \item \bold{Transformation:} Applies $log_{10}$ to the response variable to stabilize variance.
#'   \item \bold{Scaling:} Z-score scales all numeric predictors to allow for direct comparison of coefficients.
#'   \item \bold{Overfitting Guard:} Checks that the number of observations exceeds the number of predictors.
#' }
#' It uses \code{stats::step} based on AIC for model selection.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{summary}: A data frame comparing R², Adjusted R², and RSE across model types.
#'   \item \code{models}: A list containing the "Full", "Forward", and "Backward" model objects.
#'   \item \code{data}: The cleaned and scaled data frame used for modeling.
#'   \item \code{trait}: The name of the response variable.
#' }
#' @export
htp_regression <- function(traits, htp_features, trait = "yield_kg_plot", na_threshold = 0.3,
                           join_keys = c("treatment", "blocking")) {

  treatment_col <- join_keys[1]
  block_col <- join_keys[2]

  htp_clean_cols <- htp_features %>%
    dplyr::select(where(~ mean(is.na(.)) <= na_threshold))

  combined_df <- dplyr::inner_join(traits, htp_clean_cols, by = join_keys)

  feat_names <- htp_clean_cols %>%
    dplyr::select(dplyr::starts_with("F"), where(is.numeric)) %>%
    dplyr::select(-dplyr::any_of(unname(join_keys)), -dplyr::any_of(names(join_keys))) %>%
    names()

  data_clean <- combined_df %>%
    dplyr::select(dplyr::all_of(trait), dplyr::all_of(feat_names)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ is.finite(.)))

  data_clean <- data_clean %>%
    dplyr::select(dplyr::all_of(trait), where(~ sd(., na.rm = TRUE) > 0))

  feat_names <- setdiff(names(data_clean), trait)

  if (nrow(data_clean) <= (length(feat_names) + 1)) {
    warning(paste("Overfitting risk: Only", nrow(data_clean), "rows for", length(feat_names), "features."))
  }

  data_clean[[trait]] <- log10(data_clean[[trait]] + 1e-6)

  data_scaled <- data_clean %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(feat_names), ~ as.numeric(scale(.))))

  full_f <- as.formula(paste("`", trait, "` ~ .", sep = ""))
  null_f <- as.formula(paste("`", trait, "` ~ 1", sep = ""))

  m_full <- lm(full_f, data = data_scaled)

  m_fwd <- tryCatch(
    step(lm(null_f, data = data_scaled), scope = list(lower = null_f, upper = m_full), direction = "forward", trace = 0),
    error = function(e) return(NULL)
  )

  m_back <- tryCatch(
    step(m_full, direction = "backward", trace = 0),
    error = function(e) return(NULL)
  )

  extract_stats <- function(m, name) {
    if(is.null(m)) return(data.frame(Approach = name, R2 = NA, Adj_R2 = NA, RSE = NA, Num_Features = NA))
    s <- summary(m)
    data.frame(
      Approach = name,
      R2 = round(s$r.squared, 4),
      Adj_R2 = round(s$adj.r.squared, 4),
      RSE = round(s$sigma, 4),
      Num_Features = length(coef(m)) - 1,
      stringsAsFactors = FALSE
    )
  }

  summary_tab <- dplyr::bind_rows(
    extract_stats(m_full, "Full Model"),
    extract_stats(m_fwd, "Forward Selection"),
    extract_stats(m_back, "Backward Elimination")
  )

  return(list(
    summary = summary_tab,
    models = list("Full" = m_full, "Forward" = m_fwd, "Backward" = m_back),
    data = data_scaled,
    trait = trait
  ))
}


export_regression <- function(results_obj, filename = "Regression_Results.xlsx") {

  full_coef <- as.data.frame(summary(results_obj$models$Full)$coefficients)
  forward_coef <- as.data.frame(summary(results_obj$models$Forward)$coefficients)
  backward_coef <- as.data.frame(summary(results_obj$models$Backward)$coefficients)

  export_list <- list(
    "Summary_Stats" = results_obj$summary,
    "Full_Model" = full_coef,
    "Forward_Selection" = forward_coef,
    "Backward_Elimination" = backward_coef,
    "Data_Used" = results_obj$data
  )

  write.xlsx(export_list, file = filename, rowNames = TRUE, overwrite = TRUE)

  cat("File saved:", filename, "\n")
}

########## ANOVA


#' Run Mixed-Effects ANOVA and Post-Hoc Tests
#'
#' This function fits a Linear Mixed Model (LMM) to a target response variable,
#' treating "treatment" as a fixed effect and "blocking" as a random effect.
#' It performs Type III ANOVA, calculates R² for mixed models, and generates
#' Compact Letter Displays (CLD) for mean comparisons.
#'
#' @param df A data frame containing ground-truth or phenotypic data.
#' @param htp_df A data frame containing HTP-derived features (F-metrics) or
#'   other predictors to be joined.
#' @param response Character. The dependent variable to be analyzed
#'   (e.g., "yield_kg_plot" or an HTP feature like "F1").
#' @param join_keys A named character vector for joining \code{df} and \code{htp_df}.
#'   Default is \code{c("treatment" = "treatment", "bloque" = "blocking")}.
#'
#' @details
#' The function fits the following model:
#' \deqn{Response \sim Treatment + (1 | Block)}
#' It utilizes:
#' \itemize{
#'   \item \code{lme4::lmer} for model fitting.
#'   \item \code{MuMIn::r.squaredGLMM} for conditional R² calculation.
#'   \item \code{DHARMa} for residual diagnostics (dispersion testing).
#'   \item \code{emmeans} and \code{multcomp} for post-hoc means and letter grouping.
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Fixed_Effects}: Tidy summary of fixed effect estimates.
#'   \item \code{Random_Effects}: Tidy summary of random effect parameters.
#'   \item \code{Fit}: Model glance including AIC, BIC, R², and dispersion.
#'   \item \code{Anova}: Type III ANOVA table.
#'   \item \code{Means}: Data frame of estimated marginal means with CLD groupings.
#' }
#' @export
run_glmer_anova <- function(traits, htp_features, trait = "yield_kg_plot",
                            join_keys = c("treatment", "blocking")) {

  treatment_col <- join_keys[1]
  block_col <- join_keys[2]

  combined_df <- traits %>%
    dplyr::left_join(htp_features, by = join_keys)

  if (!trait %in% names(combined_df)) return(NULL)

  clean_df <- combined_df %>%
    dplyr::filter(is.finite(!!dplyr::sym(trait)),
                  !is.na(!!dplyr::sym(treatment_col)),
                  !is.na(!!dplyr::sym(block_col))) %>%
    dplyr::mutate(!!dplyr::sym(treatment_col) := as.factor(!!dplyr::sym(treatment_col)),
                  !!dplyr::sym(block_col) := as.factor(!!dplyr::sym(block_col)))

  if (nrow(clean_df) < 5 || sd(clean_df[[trait]], na.rm = TRUE) == 0) return(NULL)

  formula_str <- paste0("`", trait, "` ~ `", treatment_col, "` + (1 | `", block_col, "`)")

  m1 <- tryCatch({
    lme4::lmer(as.formula(formula_str), data = clean_df)
  }, error = function(e) return(NULL))

  if (is.null(m1)) return(NULL)

  r2_vals <- MuMIn::r.squaredGLMM(m1)
  sim_resid <- DHARMa::simulateResiduals(m1, plot = FALSE)
  dispersion <- DHARMa::testDispersion(sim_resid, plot = FALSE)
  anova_res <- car::Anova(m1, type = "III")

  em <- emmeans::emmeans(m1, as.formula(paste("~", treatment_col)))
  res_cld <- multcomp::cld(em, Letters = letters)

  list(
    "Fixed_Effects" = broom.mixed::tidy(m1, effects = "fixed"),
    "Random_Effects" = broom.mixed::tidy(m1, effects = "ran_pars"),
    "Fit" = broom.mixed::glance(m1) %>%
      dplyr::mutate(r2_marginal = r2_vals[1,1],
                    r2_conditional = r2_vals[1,2],
                    dispersion = as.numeric(dispersion$statistic)),
    "Anova" = broom::tidy(anova_res),
    "Means" = as.data.frame(res_cld)
  )
}



### PLOT REGRESION MODEL

plot_comparison_grid <- function(analysis_output, plot_path) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  target <- analysis_output$trait
  mods <- analysis_output$models
  df <- analysis_output$data

  valid_names <- names(mods)[sapply(mods, function(x) !is.null(x))]

  all_plots_df <- lapply(valid_names, function(nm) {
    m <- mods[[nm]]
    s <- summary(m)
    data.frame(
      Obs = df[[target]],
      Fit = predict(m),
      ModelName = nm,
      Stats = paste0("R² = ", round(s$r.squared, 3), "\nRSE = ", round(s$sigma, 3))
    )
  }) %>% bind_rows()

  final_plot <- ggplot(all_plots_df, aes(x = Fit, y = Obs)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 0.8) +
    geom_label(aes(label = Stats), x = Inf, y = -Inf,
               hjust = 1.05, vjust = -0.2, size = 3.5,
               fill = "white", label.size = 0.2, label.r = unit(0, "lines")) +
    facet_wrap(~ModelName, scales = "free") +
    labs(x = "predicted.values", y = "observed.values") +
    theme_grey() +
    theme(
      aspect.ratio = 1,
      strip.text = element_text(size = 10),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      axis.title = element_text(size = 12)
    )

  if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

  ggsave(file.path(plot_path, "regresion_plots.png"),
         plot = final_plot, width = 10, height = 5, dpi = 300)

  return(final_plot)
}


####
#### PLOT OBS VS PRED
####


# get_predictions <- function(model, data, source){
#
#   fitted_values <- fitted(model)
#
#   plot_data <- data.frame(
#     Observed = data$PGRA,
#     Fitted = fitted_values,
#     VAR = data$VAR  # Add the VAR column
#   ) %>% dplyr::mutate(input_data = source)
#
#   return(plot_data)
#
# }
#
# observed_vs_fitted_plot <- function(data) {
#   plot <- ggplot(data, aes(x = Observed, y = Fitted)) +
#     geom_point(color = "blue", fill = "blue", size = 3, shape = 21, alpha = 0.3) + # Blue points
#     geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", size = 1) + # 1:1 line
#     geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", size = 0.5) + # Regression line
#
#     stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
#                  formula = y ~ x, parse = TRUE,
#                  color = "black", size = 4,
#                  label.y = 0.1, label.x = 0.5) +
#
#     facet_wrap(~ input_data, scales = "fixed") + # Custom facet labels
#     theme_bw() +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_line(colour = "gray", size = 0.5),
#       panel.border = element_rect(colour = "black", fill = NA),
#       strip.background = element_blank(),
#       strip.text = element_text(size = 15),
#       plot.title = element_text(hjust = 0.5, size = 20),
#       axis.text = element_text(size = 12.5),
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       axis.title = element_text(size = 15),
#       legend.text = element_text(size = 15),
#       legend.title = element_text(size = 15)
#     ) +
#     labs(
#       title = " ",
#       x = expression(Observed ~  weight ~ (g ~ plant^{-1})),
#       y = expression(Predicted ~  weight ~ (g ~ plant^{-1}))
#     )
#   return(plot)
# }




