
calculate <- function(params) {
    output <- list()

    # working directory should be a unique, empty folder
    setwd(params$directory)
    set.seed(as.integer(params$rand_seed))

    if(params$gis){
        output <- plot_map(params)
    }
    else{

        if(params$win == "unit_circle") {
            #  actual unit circle is spatstat::disc() # defaults: radius 1, center: 0, 0
            window <- spatstat::disc(radius = 0.5, centre = c(0.5, 0.5))
        } else if(params$win == "unit_square") {
            window <- spatstat::unit.square()
        } else if(params$win == "rectangle") {
            x1 <- params$x_origin
            x2 <- x1 + params$width
            y1 <- params$y_origin
            y2 <- y1 + params$height

            window <- spatstat::owin(c(x1, x2), c(y1, y2))
        } else if(params$win == "circle") {
            x <- params$x_origin
            y <- params$y_origin

            window <- spatstat::disc(radius = params$radius, centre = c(x,y))
        }


        sp_params = list(
            win = window,
            sim_total = as.integer(params$sim_total),
            x_case = as.double(unlist(params$x_case, use.names=FALSE)),
            y_case = as.double(unlist(params$y_case, use.names=FALSE)),
            samp_case = params$samp_case,
            samp_control = params$samp_control,
            x_control = as.double(unlist(params$x_control, use.names=FALSE)),
            y_control = as.double(unlist(params$y_control, use.names=FALSE)),
            n_case = as.integer(unlist(params$n_case, use.names=FALSE)),
            n_control = as.integer(unlist(params$n_control, use.names=FALSE)),
            r_case = as.double(unlist(params$r_case, use.names=FALSE)),
            r_control = as.double(unlist(params$r_control, use.names=FALSE)),
            s_case = as.double(unlist(params$s_case, use.names=FALSE)),
            s_control = as.double(unlist(params$s_control, use.names=FALSE)),
            lower_tail = as.double(params$alpha/2),
            upper_tail = as.double(1-(params$upper_tail/2)),
            n_core = 4
        )

        if(params$sim_total == 1) {
            # additional parameters for spatial_data are ignored
            results <- do.call(sparrpowR::spatial_data, sp_params)
        }

        else if(params$sim_total > 1) {
            results <- do.call(sparrpowR::spatial_power, sp_params)
            s_stat <- t.test(results$s_obs, mu = 1, alternative = "two.sided")
            t_stat <- t.test(results$t_obs, mu = 0, alternative = "two.sided") 

            output$summary <- list(
                mean_n_con = mean(results$n_con),
                mean_n_cas = mean(results$n_cas),
                mean_bandw = mean(results$bandw),
                sd_n_con = sd(results$n_con),
                sd_n_cas = sd(results$n_cas),
                sd_bandw = sd(results$bandw),
                s_test_stat = s_stat$statistic,
                t_test_stat = t_stat$statistic,
                s_pval = s_stat$p.value,
                t_pval = t_stat$p.value
            )
        }

        # save output file
        saveRDS(results, "results.rds")

        # generate plots and return output
        output$plots <- plot_results(results, params)
        output$id <- params$id
        output
    }
}

replot <- function(params) {
    setwd(params$directory)
    results <- readRDS(params$rds)
    list(plots = plot_results(results, params))
}

plot_results <- function(results, params) {
    # cols[3] (mid_color) is only used when plot_text == TRUE, and is not actually used for legends. 
    # cols[1:3] do not match the order of colors in the documentation
    params$cols <- c(params$insuff_color, params$suff_color, params$case_color, params$control_color)
    params$chars <- as.integer(c(params$case_symbol, params$control_symbol))
    params$sizes <- as.double(c(params$case_size, params$control_size))

    # todo: specifying width and height above default makes plotting area collide with legend
    if (!'plot_format' %in% names(params)) params$plot_format <- "png"
    if (!'plot_width' %in% names(params)) params$plot_width <- 720
    if (!'plot_height' %in% names(params)) params$plot_height <- 720

    # svg files are rather large compared to other formats due to a large number of paths
    # svg width/heights are specified in inches, not pixels

    scale <- sqrt(params$plot_width ^ 2 + params$plot_height ^ 2)/sqrt(480 ^ 2 + 480 ^ 2)

    # set up graphics device
    do.call(params$plot_format, list(
        filename = paste0("plot-%d.", params$plot_format),
        width = params$plot_width, 
        height = params$plot_height
    ))

    sparrpowR::spatial_plots(results,
            p_thresh = params$p_thresh,
            chars = params$chars,
            sizes = params$sizes,
            plot_pts = params$plot_pts,
            plot_title = params$title, 
            cascon = as.logical(params$cascon),
            scale = scale,
            plot_axes = params$axes,
            plot_square = params$plot_square,
            horizontal = params$horizontal,
            cols = params$cols)
    dev.off()

    
    file.rename(paste0("plot-1.",params$plot_format),paste0("simulated-data.",params$plot_format))
    file.rename(paste0("plot-2.",params$plot_format),paste0("local-power-continuous-scale.",params$plot_format))
    file.rename(paste0("plot-3.",params$plot_format),paste0("local-power-above-threshold.",params$plot_format))

    # add generated plots
    files <- list(paste0("simulated-data.",params$plot_format),paste0("local-power-continuous-scale.",params$plot_format),paste0("local-power-above-threshold.",params$plot_format))
    files
}

plot_map <- function(params){
    output <- list()
    # Washington, D.C. boundary
    gis_path1 <- "https://opendata.arcgis.com/datasets/7241f6d500b44288ad983f0942b39663_10.geojson"
    dc <- geojsonio::geojson_read(gis_path1,  what = "sp")  

    # American Community Survey 2018 Census Tracts
    gis_path2 <- "https://opendata.arcgis.com/datasets/faea4d66e7134e57bf8566197f25b3a8_0.geojson"
    census <- geojsonio::geojson_read(gis_path2,  what = "sp")

    clipwin <- maptools::unionSpatialPolygons(census, IDs = rep(1, length(census)))
    dcc <- rgeos::gIntersection(dc, clipwin, byid = TRUE)

    dcp <- sp::spTransform(dcc, CRSobj = sp::CRS(projargs = "+init=EPSG:32618"))
    dco <- spatstat::as.owin(dcp)

    #Navy, change to user lon and lat later
    navy <- data.frame(lon = 326414.70444451, lat = 4304571.1539442)
    spf <- sp::SpatialPoints(coords = navy, proj4string = sp::CRS(projargs = "+init=EPSG:32618"))

    sim_power <- sparrpowR::spatial_power(x_case = navy[[1]], y_case = navy[[2]], # center of cluster
                           x_control = navy[[1]], y_control = navy[[2]], # center of cluster
                           n_case = 50, n_control = 950, # sample size of case/control
                           samp_case = "MVN", samp_control = "MVN", # samplers
                           s_case = 1000, s_control = 2000, # approximate size of clusters
                           lower_tail = 0.025, upper_tail = 0.975, # two-tailed alpha
                           sim_total = 2, # number of iterations
                           win = dco, # study area
                           resolution = 100, # number gridded knots on x-axis
                           edge = "diggle", # correct for edge effects
                           adapt = FALSE, # fixed-bandwidth
                           h0 = NULL, # automatically select bandwidth for each iteration
                           verbose = FALSE) # no printout

    s_stat <- t.test(sim_power$s_obs, mu = 1, alternative = "two.sided")
    t_stat <- t.test(sim_power$t_obs, mu = 0, alternative = "two.sided") 
    output$summary <- list(
                mean_n_con = mean(sim_power$n_con),
                mean_n_cas = mean(sim_power$n_cas),
                mean_bandw = mean(sim_power$bandw),
                sd_n_con = sd(sim_power$n_con),
                sd_n_cas = sd(sim_power$n_cas),
                sd_bandw = sd(sim_power$bandw),
                s_test_stat = s_stat$statistic,
                t_test_stat = t_stat$statistic,
                s_pval = s_stat$p.value,
                t_pval = t_stat$p.value)

        
    # save output file
    saveRDS(sim_power, "results.rds")

    # generate plots and return output
    output$plots <- plot_results(sim_power, params)
    output$id <- params$id
    
    expandbb <- function(bb, f) {
        x <- bb[3] - bb[1] # range of x values
        y <- bb[4] - bb[2] # range of y values
        nb <- bb # make a copy
        nb[1] <- bb[1] - (f * x) # xmin - left
        nb[3] <- bb[3] + (f * x) # xmax - right
        nb[2] <- bb[2] - (f * y) # ymin - bottom
        nb[4] <- bb[4] + (f * y) # ymax - top
        return(nb)}

        dcbb <- expandbb(bb = sp::bbox(dc), f = 0.01)
        base_map <- ggmap::get_map(location = dcbb, maptype = "roadmap", source = "osm")

        sim_pts <- sim_power$sim  # extract points from first iteration
        sim_pts <- maptools::as.SpatialPointsDataFrame.ppp(sim_pts) # convert to spatial data frame
        raster::crs(sim_pts) <- sp::proj4string(dcp) # set initial projection
        sim_pts_wgs84 <-  sp::spTransform(sim_pts, CRSobj = sp::CRS(projargs = "+init=epsg:4326")) # project to basemap
        sim_pts_df <- tibble::tibble(data.frame(sim_pts_wgs84)) # convert to tidy data frame

        dc_df <- broom::tidy(dcc) # conver to a tidy dataframe
        dcc$polyID <- sapply(slot(dcc, "polygons"), function(x) slot(x, "ID")) # preserve polygon id for merge
        dc_df <- merge(dc_df, dcc, by.x = "id", by.y="polyID") # merge data

        pvalprop <- tibble::tibble(x = sim_power$rx, y = sim_power$ry,
                           z = sim_power$pval_prop_cas) # extract proportion significant
        lrr_narm <- na.omit(pvalprop) # remove NAs
        sp::coordinates(lrr_narm) <- ~ x + y # coordinates
        sp::gridded(lrr_narm) <- TRUE # gridded
        pvalprop_raster <- raster::raster(lrr_narm) # convert to raster
        rm(pvalprop, lrr_narm) # conserve memory
        raster::crs(pvalprop_raster) <- raster::crs(dcp) # set output project (UTM 18N)
        pvalprop_raster <- raster::projectRaster(pvalprop_raster, crs = raster::crs(dc)) # unproject (WGS84)
        rtp <- raster::rasterToPolygons(pvalprop_raster) # convert to polygons
        rtp@data$id <- 1:nrow(rtp@data)   # add id column for join
        rtpFort <- broom::tidy(rtp, data = rtp@data) # convert to tibble
        rtpFortMer <- merge(rtpFort, rtp@data, by.x = 'id', by.y = 'id')  # join data
        rampcols <- grDevices::colorRampPalette(colors = c(cols[5], cols[2]), space="Lab")(length(raster::values(pvalprop_raster))) # set colorramp
        
        map1 <- ggmap::ggmap(base_map) + # basemap
            ggplot2::geom_polygon(data = dc_df, # original boundary
                        ggplot2::aes(x = long, y = lat, group = group),
                        fill = "transparent",
                        colour = "black") +
            ggplot2::geom_polygon(data = rtpFortMer, # output raster as polygons
                        ggplot2::aes(x = long, y = lat, group = group, fill = z), 
                        size = 0, 
                        alpha = 0.5) +
            ggplot2::scale_fill_gradientn(colours = rampcols) + # colors for polygons
            ggplot2::geom_point(data = sim_pts_df, # simulated point-locations
                        ggplot2::aes(x = mx, y = my, color = marks, shape = marks),
                        alpha = 0.8) + 
            ggplot2::scale_color_manual(values = cols[4:5]) + # fill of point-locations
            ggplot2::scale_shape_manual(values = chars) + # shope of point-locations
            ggplot2::labs(x = "", y = "", fill = "Power", color = "", shape = "") # legend labels

        pvalprop_reclass <- raster::reclassify(pvalprop_raster, c(-Inf, params$p_thresh-0.0000001, 1,
                                                params$p_thresh-0.0000001, Inf, 2))
        rtp <- raster::rasterToPolygons(pvalprop_reclass) # convert to polygons
        rtp@data$id <- 1:nrow(rtp@data)   # add id column for join
        rtpFort <- broom::tidy(rtp, data = rtp@data) # convert to tibble
        rtpFortMer <- merge(rtpFort, rtp@data, by.x = 'id', by.y = 'id')  # join data

        map2 <- ggmap::ggmap(base_map) + # basemap 
            ggplot2::geom_polygon(data = dc_df, # original boundary
                        ggplot2::aes(x = long, y = lat, group = group),
                        fill = "transparent",
                        colour = "black") +
            ggplot2::geom_polygon(data = rtpFortMer, # output raster as polygons
                        ggplot2::aes(x = long, y = lat, group = group, fill = as.factor(z)), 
                        size = 0, 
                        alpha = 0.5) +
            ggplot2::scale_fill_manual(values = cols[c(5,2)],
                                        labels = c("insufficient", "sufficient")) + # colors for polygons
            ggplot2::labs(x = "", y = "", fill = "Power") # legend labels

    output
}