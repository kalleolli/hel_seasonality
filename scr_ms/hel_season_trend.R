
# script for Helsinki outer archipelago phytoplankton time series seasonal and long-term decomposition

# load libraries

library(dplyr)
library(tibble)
library(tidyr) # pivot_wider, pivot_longer
library(vegan)
library(sp)
# library(rgdal)
library(parallel)
library(mgcv)
library(gratia)
library(lubridate)
library(ggside) # for marginal distribution plots
library(cowplot) # plot_grid
library(proxy) # for row.dist, col.dist



# finds where is the project folder [has to be somewhere in #homedirectory/Documents]
(hel_season_dir <- list.files(path = '~/Documents', full.names = TRUE, recursive = TRUE, pattern = 'hel_seasonality.Rproj') %>% dirname())

# Load data ####

if(T){
  
  load(paste0(hel_season_dir, "/dat_ms/hel_season.rda"))  # 1.4 mb file
  # loaded objects dat [176840 x 5] and met [5537 x 9] meta-data
  # dat = sampleID, valid_AphiaID, valid_name, UnitsL, ww
  # met = sampleID, stn, obstime [1966-2018, Jan - Dec], chla, lat, lon, Nimi, Salin, Temp
  
  # make  community matrix - species in columns, samples in rows
  # this matrix is used for ordination
} # loads data; provides objects: dat, met


# Community matrix and meta table ####

if(T){ # make community matrix datcm and meta table

  #  tmp <- group_by(dat, sampleID, valid_name) %>% summarise(ww = sum(ww, na.rm = T)) %>% pivot_wider(names_from = valid_name, values_from = ww) 
#  datcm  <- data.matrix(tmp[ , -1]); dim(datcm) # community matrix on dat 4630 samples × 619 taxa
#  rownames(datcm) <- tmp[[1]] # add sampleID as rownames
#  datcm[is.na(datcm)] <- 0 # change all NA's to zeros
  
  datcm <- dplyr::summarise(dat, ww = sum(ww, na.rm = T), .by = c('sampleID', 'valid_name')) %>% pivot_wider(names_from = valid_name, values_from = ww) %>% tibble::column_to_rownames('sampleID') %>% data.matrix() # alternative way to make datcm, same result
  datcm[is.na(datcm)] <- 0 # change all NA's to zeros
  
  
  # supplement meta table with additional time and diversity variables;  5537   13
  meta <- mutate(met,
                 obstime = as.Date(obstime), # convert obstime to Date format
                 year = as.numeric(format(obstime, '%Y')), # extract sampling year
                 mon = as.numeric(format(obstime, '%m')), # extract sampling month
                 jul = as.numeric(format(obstime, '%j')), # extract day of the year
                 week = as.numeric(format(obstime, '%W')), # extract sampling month
                 ww = rowSums(datcm[as.character(met$sampleID), ]),  # ww biomass
                 wwl = log(ww), # log biomass
                 fstn = factor(stn), # station as a factor
                 time = lubridate::decimal_date(obstime), # continuous time in year units
                 N0 = rowSums(datcm[as.character(met$sampleID), ] > 0),  # Species richness
                 H = vegan::diversity(datcm[as.character(met$sampleID), ]), # Shannon entropy (base e)
                 N1 = exp(H), # Shannon diversity (base e)
                 N2 = vegan::diversity(datcm[as.character(met$sampleID), ], "inv"), # Simpson diversity
                 rny = vegan::renyi(datcm[as.character(met$sampleID), ], scales = 0.1, hill = F),
                 J = H / log(N0), # Pielou evenness
                 E10 = N1 / N0   # Shannon evenness (Hill's ratio)
  )
  
  datcm <- datcm[as.character(meta$sampleID), ] # re-order datcm to match meta
  
} # Prepare the data, provides datcm [5537 x 581] mx, meta [5537 x 24] df


# FIG 1; map ####
# requires:
# GSHHS_h_L1.shp high res shapefile for coastlines somewhere in Documents structure
# GSHHS_f_L1.shp full res shapefile 
# to obtain: https://www.soest.hawaii.edu/pwessel/gshhg/

if(T){
  
  library(sf)
  library(ggsflabel)
  library(ggspatial) # scales and N arrows
  library(gghighlight)
  
  sf_use_s2(FALSE)
  
  
  if(file.exists('dat_ms/hel_season_maps.rda')){
    load('dat_ms/hel_season_maps.rda') # loads BS, ggm1, hel_cropped_sf, ggm2
  } # loads BS, hel_cropped_sf
  else {
    # you need GSHHS_h_L1.shp to run this
    gshhs_i <- list.files(path = '~/Documents', full.names = TRUE, recursive = TRUE, pattern = 'GSHHS_i_L1.shp')
    gshhs_f <- list.files(path = '~/Documents', full.names = TRUE, recursive = TRUE, pattern = 'GSHHS_f_L1.shp')
    
    # cut out Baltic Sea insert map
    shp_i <- st_read(gshhs_i) # load GSHHS L1 coastline shapefile, "28.5 Mb" size
    # Helsinki archipelago map
    shp_f <- st_read(gshhs_f) # load GSHHS L1 coastline shapefile, 270.3 Mb" size
 
    hel_box <- st_bbox(c(xmin = 24.5, ymin = 59.9 , xmax = 25.5, ymax = 60.3), crs = st_crs(shp_f))
    bs_box <- st_bbox(c(xmin = 10, ymin = 53 , xmax = 31, ymax = 67), crs = st_crs(shp_i))
    
    
    BS <- st_crop(shp_h, bs_box) # Baltic Sea insert coastline for ggm1
    hel_cropped_sf <- st_crop(shp_f, hel_box) # Helsinki outer archipelago shapefile for ggm1
    
    save(BS, hel_cropped_sf, file = 'dat_ms/hel_season_maps.rda') # save maps for future use
    
  } # provides BS, hel_cropped_sf
  
  # make insert map with Baltic Sea and red box for Helsinki archipelago
ggm1 <- ggplot(data = BS) + geom_sf(fill = 'whitesmoke') + theme_void() + geom_sf(data = st_as_sfc(bs_box), fill = NA, color = 'black', linewidth = 0.5 ) + geom_sf(data = st_as_sfc(hel_box), fill = NA, color = 'red', linewidth = 0.5 )
    

# station locations
metag <- group_by(meta, lat, lon, stn) %>% summarise(n = n()) %>% ungroup() # pel is mean biomass
metag <-  summarise(meta, n = n(), .by = c('lat','lon','stn'))

df <- st_as_sf(x = metag, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

ggm2 <- ggplot(data = hel_cropped_sf) + # Helsinki archipelago map
  geom_sf(fill= "whitesmoke", color = gray(.6)) + 
  coord_sf(xlim = c(24.55, 25.3), ylim = c(59.935, 60.25), expand = FALSE) +
  geom_point(data = df, aes(geometry = geometry, after_stat(X), after_stat(Y),  size = sqrt(n)), stat = StatSfCoordinates, fun.geometry = identity, show.legend = FALSE, shape = 21, fill = 'lightblue') +
  geom_sf_label_repel(data = df, aes(label = stn), label.size = NA, label.padding = unit(0.1, 'lines'), alpha = 0.7, fontface = 'bold') +
  xlab('Longitude') + ylab('Latitude') +
  annotation_scale(location = "bl", width_hint = 0.3, pad_x = unit(0.5, "in"), pad_y = unit(0.3, "in")) +
  annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.5, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.2), panel.background = element_rect(fill = 'aliceblue'), panel.border = element_rect(colour = "black", fill=NA, size=1))

rm(shp_f, shp_i) # garbage collection

  } # Helsinki outer archipelago map, provides ggm2
  
if(T){ # finalize and save Fig1
    
    Fig1 <- ggdraw() +
      draw_plot(ggm2) +
      draw_plot(ggm1, scale = 0.4, halign = 1.11, valign = 0.145)
    
    save_plot("./figs_ms/Fig1.pdf", Fig1,  ncol=1, base_height = 6.5, base_width = 7.5)
} # finalise and save Fig1
  


# FIG 2; biomass seasonality and trend global effect sizes ####
if(T){
  Fig2 <- ggplot(filter(meta, year > 1970), aes(x = time, y = wwl, col = jul)) + 
    geom_point( size = 1.2, alpha = 0.5) +
    geom_smooth(method = quantreg::rq,  se = FALSE,  method.args = list(tau = 0.05 ), col = 'darkgray', linewidth = 0.75) +
    geom_smooth(method = quantreg::rq,  se = FALSE,  method.args = list(tau = 0.95 ), col = 'darkgray', linewidth = 0.75) +
    stat_summary(mapping = aes(x = year+mon/12, y = log(ww)), fun = mean, geom = 'line', col = 'black', linewidth = 0.5) +
    geom_smooth(method = gam, formula = y ~ s(x), method.args = list(method = 'REML')) +
    xlab('Year') + ylab('Log biomass') +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) + 
    ylim(4, 10.5) + 
    theme(legend.position = "inside", legend.direction = "horizontal") + theme(legend.position.inside =  c(0.5, 0.1)) + theme(legend.key.width = unit(1.5, "cm")) + theme(legend.background = element_blank()) +
    scale_color_gradientn(colours = c("darkgray", "blue", "blue", "red","red", "darkgray"), values = c(0, 90, 140, 200, 250, 365)/365, name = "doy")
  
  cowplot::save_plot("./figs_ms/Fig2.pdf", Fig2,  ncol=1, base_height = 6, base_width = 8)
  
}  # FIGURE 2

# FIG 3; GAMM main effects ####
if(T){
  ctrl <- list(niterEM = 0, msVerbose = FALSE, optimMethod="L-BFGS-B")
  ww.gam <- gamm(wwl ~ ti(jul, bs = 'cc', k = 12) + ti(time, bs = 'tp', k = 10) + ti(jul, time, bs = c("cc","tp")), data = filter(meta, year > 1970), correlation = corCAR1(form = ~ 1|year), control = ctrl, method = 'REML')$gam
  
  if(F){
    var(predict(ww.gam, type="terms")[, "ti(jul)"])  # variance from seasonality 0.65
    var(predict(ww.gam, type="terms")[, "ti(time)"])  # variance from long term trend 0.11
    var(predict(ww.gam, type="terms")[, "ti(jul,time)"])  # variance from interaction term 0.017
  } # GAMM main effects variance partitioning
  
  Fig3a <- smooth_estimates(ww.gam, select = c("ti(jul)")) %>% draw() + labs(x='Day of year', title = 'Seasonal smooth',caption = NULL)
  Fig3b <- smooth_estimates(ww.gam, select = c("ti(time)")) %>% draw() + labs(x='Calendar year', title = 'Trend smooth', caption = NULL)
  
  Fig3 <- plot_grid(Fig3a, Fig3b, ncol = 2, labels = "AUTO")
  cowplot::save_plot("./figs_ms/Fig3.pdf", Fig3,  ncol=2, base_height = 4, base_width = 4)
  
  ww.gam %>% summary() # GAMM model statistics
  
} # FIG 3 GAM main effects; saves figs_ms/Fig3.pdf

# FIG 4; wwl GAM 3D model ####
# Spring bloom biomass peak around day 150. Is the biomass peak time-invariant?

if(T){ # 3D model of season and time interaction
 
  # bloom ww peak
  if(T){ # find the biomass peak day of the year
    
    peaks <- bloom_start <- bloom_end <- integer()
    eps <- 1e-07
    xy_gam <- gamm(wwl ~ te(jul, year, bs = c("cc","tp"), k = c(12, 10)), data = meta, correlation = corCAR1(form = ~ 1|year), control = ctrl, method = 'REML')$gam
    
    for(i in 1966:2019){
      #i <- 1970
      ww.fit <- predict(xy_gam, newdata = data.frame(jul = 1:200, year = i))
      fit2   <- predict(xy_gam, newdata = data.frame(jul = 1:200 + eps, year = i))
      der1 <- (fit2-ww.fit)/eps
      peaks <- c(peaks, which(ww.fit == max(ww.fit)))
      bloom_start <- c(bloom_start, which(der1 == max(der1)))
      bloom_end <- c(bloom_end, which(der1 == min(der1)))
    }
    
    
    peaks.df <- data.frame(peak = peaks, time = 1966.6:2019.6, ww = 0, bloom_start = bloom_start, bloom_end = bloom_end, bloom_dur = bloom_end-bloom_start)
    
    # bloom_start.df <- data.frame(jul = bloom_start, time = 1966.6:2019.6, ww = 0)
    # bloom_end.df <- data.frame(jul = bloom_end, time = 1966.6:2019.6, ww = 0)
    
    peak_lm <- lm(peak ~ time, data = peaks.df)# significant negative slope, bloom peaks get earlier 0.27 days per year
    
    # peaks.df %>% mutate(obs_juldate=as.Date(peak, start = '01/01/2020')) %>% ggplot(aes(x = time, y = obs_juldate)) + geom_point() + geom_smooth(method = 'lm')
    # bloom timing meanders quite a bit, but trend is decreasing
    
  } # find the biomass peak day of the year; provides peaks.df, bloom_start.df, bloom_end.df
  
  # make full tensor product te(jul, time) and ... 
  # marginal distribution plots for ti(jul) and time ti(time) with separate main effects and interaction ti(jul, time)
  xy_gam <- gamm(wwl ~ te(jul, time, bs = c("cc","tp"), k = c(12, 10)), data = meta, correlation = corCAR1(form = ~ 1|year), control = ctrl, method = 'REML')$gam

  chl_gam <- gamm(chla ~ te(jul, time, bs = c("cc","tp"), k = c(12, 10)), data = meta, correlation = corCAR1(form = ~ 1|year), control = ctrl, method = 'REML')$gam
  chl.gam <- gamm(wwl ~ ti(jul, bs = 'cc', k = 12) + ti(time, bs = 'tp', k = 10) + ti(jul, time, bs = c("cc","tp")), data = filter(meta, year > 1970, jul %in% 100:300), correlation = corCAR1(form = ~ 1|year), control = ctrl, method = 'REML')$gam
  
  gratia::smooth_estimates(chl.gam, newdata = newd, select = 'ti(jul,time)', dist = 0.1) %>% draw()
  sm_chl <- smooth_estimates(chl.gam, select = c('te(jul, time)'))
  
  margin_gam <- gamm(wwl ~ ti(jul, bs = 'cc', k = 12) + ti(time, bs = 'tp', k = 10) + ti(jul, time, bs = c("cc","tp")), data = meta, correlation = corCAR1(form = ~ 1|year), control = ctrl, method = 'REML')$gam
  
  # weekly grid for all julian days and years
  newd <- expand.grid(jul = seq(1, 365, by = 7), time = seq(1966, 2019)) # weekly grid 2862 x   2
  
  # use gratia::smooth_estimates to get the smooth estimates for the full tensor product (ww_xy) and marginal distributions (ww_x, ww_y)
  ww_xy <- gratia::smooth_estimates(xy_gam, newdata = newd, select = 'te(jul,time)', dist = 0.1) %>% mutate(ww = .estimate + coef(xy_gam)[1])#
  sm_chl <- gratia::smooth_estimates(chl.gam, select = 'te(jul, time)')
  
  ww_xmar <- gratia::smooth_estimates(margin_gam, newdata = newd, select = 'ti(jul)') %>% mutate(ww = .estimate + coef(margin_gam)[1])
  ww_ymar <- gratia::smooth_estimates(margin_gam, newdata = newd, select = 'ti(time)') %>% mutate(ww = .estimate + coef(margin_gam)[1])
  

  # point 2D smooth legend labels to midpoints
  mybreaksa <- seq(min(ww_xy$ww, na.rm=T), max(ww_xy$ww, na.rm=T), length.out = 20)
  mybreaksa <- round((head(mybreaksa, -1) + tail(mybreaksa, -1)) / 2 ,1)
  
Fig4 <- ggplot(data = ww_xy, aes(jul, time, z = ww)) + 
    ggside::geom_ysideline(data = ww_ymar, aes(y = time, x = ww), col = 1, orientation = 'y') + 
    ggside::geom_xsideline(data = ww_xmar, aes(y = ww, x = jul),  col = 'gray50', size = 1) + 
    geom_contour_filled(bins = 20, show.legend = TRUE, linetype = 1, color = 3, size = 0.3) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_point(data = meta, aes(x = jul, y = time), shape  = 19, size = .7, alpha = 0.3, show.legend = FALSE) +
    geom_line( data = peaks.df,      aes(x = peak, y = time), size = 1,  col = 'white', orientation = 'y') + 
    geom_line( data = peaks.df,      aes(x = bloom_start, y = time), size = 0.5,  col = 'black', orientation = 'y') + 
    geom_line( data = peaks.df,      aes(x = bloom_end, y = time), size = 0.5,  col = 'black', orientation = 'y') + 
    scale_fill_manual(values = colorRamps::matlab.like(19), labels = round(mybreaksa, 1), name = 'Biomass') +
    scale_ysidex_continuous(guide = guide_axis(angle = 90)) +
    scale_xsidey_continuous(expand = c(0, 0.02)) +
    theme(ggside.panel.scale = .15) + xlab('Season') + ylab('Calendar year') +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14), strip.text.x = element_text(size = 18))


  save_plot("./figs_ms/Fig4.pdf", Fig4, ncol = 1, base_height = 6, base_width = 6) 

  

} # 3D model of season and time interaction; saves figs_ms/Fig4.pdf


# FIG 5 wwl rbeast ####

if(T){ # beast irreg for ww seasonality breakpoints
  library(Rbeast) # for beast irreg
  # beast irreg for seasonality breakpoints
  meta_b <- filter(meta, mon %in% 4:10) # filter meta for julian days 100-300, i.e. spring and summer months)
  
  beast.ww <- beast.irreg(y = meta_b$wwl, time = meta_b$obstime, freq = 12, deltat = 1/12, period = 1, tseg.min = 12*5, sseg.min =  12*5)

quartz(width = 6, height = 8)
  plot(beast.ww, vars = c('y','s','scp','t','tcp','error'), main = "Biomass seasonality and trend detection") 
dev.copy2pdf(file='./figs_ms/Fig5.pdf', out.type='pdf')

} # FIG 5 ww rbeast; saves figs_ms/Fig5.pdf


# FIG 6; MDS block ####

if(T){  # dbRDA
  datcm_wbray <- wisconsin(sqrt(datcm)) %>% vegdist() # wisconsin standardisation and bray dist of the community matrix; "117.3 Mb" Mb
  # metric multidimensional scaling with bray distance
  if(!exists("dbrMDS")){system.time(dbrMDS <- vegan::dbrda(datcm_wbray ~ 1))} # distance based RDA, equivalent to metric MDS
  
  dat_ord <- data.frame(meta[, c('time','jul')], scores(dbrMDS))
  env_fit <- envfit(list(sites = scores(dbrMDS)), meta[, c('time','jul')], permutations = 9999) 
  # envfit with doy 0.37 and time 0.83 r2
  env_ar <- env_fit$vectors$arrows * env_fit$vectors$r * 1.5 # make the arrows a bit longer
  rownames(env_ar) <- c('Year','Doy') #
  
  Fig6a_mMDS <- ggplot(arrange(dat_ord, jul), aes(x=MDS1, y=MDS2, fill = jul)) +
    geom_point(size = 3, shape = 21) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'Doy') +
    theme(legend.position = c(0.9, 0.2), legend.background = element_rect(fill = rgb(1,1,1,0.7))) +
    scale_x_continuous(expand = c(0,0), limits = c(-1.7, 1.8)) +
    scale_y_continuous(expand = c(0,0), limits = c(-1.7, 1.8)) + 
    annotate('segment',x = 0, y = 0, xend = env_ar[, 'MDS1'], yend = env_ar[,'MDS2'], arrow = arrow(length = unit(0.5, "cm")), size = 1) +
    annotate('label', x = env_ar[, 'MDS1'], y = env_ar[,'MDS2'], label = rownames(env_ar), size = 5, fontface = 'bold',  vjust = "outward", hjust = "outward", fill = 'white', alpha =  0.7, label.size = NA, label.padding = unit(0.1, "lines")) +
    labs(x = 'dbRDA1', y = 'dbRDA2')
  
  # Fig6b season
  target_matrix <- cbind(scale(meta$jul), rep(0, length(ncol(meta))))
  dbrMDS_proc <- procrustes(target_matrix, scores(dbrMDS, display = "sites"))
  mMDS_rot <- data.frame(meta[, c('time','jul')], dbrMDS_proc$Yrot)
  
  Fig6b_mMDS <- ggplot(arrange(mMDS_rot, time), aes(x = jul, y = X1, fill = time)) + 
    geom_point(size = 3, shape = 21) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'Year') + 
    theme(legend.position = c(0.1, 0.8), legend.background = element_rect(fill = rgb(1,1,1,0.7))) +
    labs(x = 'Day of year', y = 'dbRDA')
  
  # Fig6c trend
  target_matrix <- cbind(scale(meta$jul), rep(0, length(meta$jul)))
  dbrMDS_proc <- procrustes(target_matrix, scores(dbrMDS, display = "sites"))
  mMDS_rot <- data.frame(meta[, c('time','jul')], dbrMDS_proc$Yrot)
  
  Fig6c_mMDS <- ggplot(arrange(mMDS_rot, jul), aes(x=time, y=X2, fill = jul)) + 
    geom_point(size = 3, shape = 21) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(50), name = 'doy') + 
    theme(legend.position = c(0.85, 0.7), legend.background = element_rect(fill = rgb(1,1,1,0.7))) +
    labs(x = 'Calendar year', y = 'dbRDA')
  
  
  Figure6_mMDS <- plot_grid(Fig6a_mMDS, Fig6b_mMDS, Fig6c_mMDS, ncol = 3, labels = "AUTO")
  save_plot("./figs_ms/Fig6.pdf", Figure6_mMDS, ncol = 3, base_height = 4, base_width = 4)
  
} # MMDS provides Figure6_mMDS



# FIG 7; composition beast ####

if(T){
  # if(!exists("datcm_wbray")){load('./aux/dat/dist_mx.rda')} # load the distance matrix "datcm_wbray"
  if(!exists("dbrMDS")){dbrMDS <- vegan::dbrda(datcm_wbray ~ 1)}
  
  # rotate the dbrMDS so that the first axis is aligned with time
  target_matrix <- cbind(scale(meta$time), rep(0, length(ncol(meta))))
  dbrMDS_proc <- procrustes(target_matrix, scores(dbrMDS, display = "sites"))
  # extract the rotated scores, amend with obstime for seasonality detection
  dbrMDS_rot <- data.frame(meta[, c('jul','obstime')], dbrMDS_proc$Yrot)
  # actual seasonality brakepoint detection with beast irreg
  beast.ord <- Rbeast::beast.irreg(y = dbrMDS_rot$X2, time = dbrMDS_rot$obstime, freq = 12, deltat = 1/12, period = 1, tseg.min = 12*5, sseg.min =  12*5); plot(beast.ord)
  # detects 1999 seasonality breakpoint as 1.0 probability; also 1971 at 0.966 probability
  # note that although there is trend change in the 90's, the seasonality remains, but changes in 2000 as the P rich deep anoxic water reaches GOF
  
  # save the plot
  quartz(width = 6, height = 8)
  plot(beast.ord, vars = c('y','s','scp','error'), main = "Community composition seasonality detection")
  dev.copy2pdf(file='./figs_ms/Fig7.pdf', out.type='pdf')
}

# source('./scr/hel_season_trend.R') # in one run

# end of script