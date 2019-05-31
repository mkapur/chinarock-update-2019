require(dplyr)
require(readr)
require(purrr)
require(kaputils)
require(r4ss)
devtools::source_url("https://raw.githubusercontent.com/mkapur/kaputils/master/R/SS_readforecastMK.R") ## use dev version
devtools::source_url("https://raw.githubusercontent.com/r4ss/r4ss/development/R/SS_ForeCatch.R") ## use dev version
devtools::source_url("https://raw.githubusercontent.com/mkapur/kaputils/master/R/SS_writeforecastMK.R") ## use dev version
# devtools::source_url("https://raw.githubusercontent.com/mkapur/kaputils/master/R/SS_executivesummaryMK.R")
source("./R/SS_executivesummaryMK.R")
## Execute automated forecasts

# for(r in c('North','Central','South')){
#   rootdir.temp <- paste0("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/cr",r)
#   catch_projections <- read.csv(paste0(rootdir.temp,"/cproj_",r,".csv"))
#   SS_autoForecast(rootdir = rootdir.temp,
#                   basedir = "base2015",
#                   catch_proportions = catch_projections[5,5:ncol(catch_projections)],
#                   # catch_proportions = c(0.5,0.08426184,0.4157382),
#                   forecast_start = 2021,
#                   forecast_end = 2031,
#                   fixed_catches = catch_projections[1:4,5:ncol(catch_projections)],
#                   Flimitfraction = catch_projections$PSTAR_0.45[catch_projections$YEAR >2020])
#   read.csv(paste0(rootdir,"/forecasts/decision_table_base.csv"))
#
# }

## TABLE G stitching ----
rootdir <- "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/tgfiles"
regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
regnames <- basename(regdirs)

df <- list.files(regdirs, recursive = TRUE, full.names = TRUE)[grep('decision_table_base.csv',list.files(regdirs,recursive = TRUE,full.names = TRUE))] %>%
  map_df(~ read_csv(.,
                    col_types = cols(.default = "c")),
         .id="index") %>%
  mutate( REG = gsub( "cr", "", regnames[as.numeric(index)] ))
write.csv(df,paste0(rootdir,"./",Sys.Date(),"table_g.csv"), row.names = FALSE)

## fixed catch for next decision table is average of 2019/2020 input
## I am going to hard code this -- not manual but not function
#
forecast_start <- 2021; forecast_end <- 2031; t = 11
## base M is -2.94, low is -2.99, high is -2.41
for(r in c('North','Central','South')){ ## loop regions
  for(catch in c('constant','ABC','upper')){ ## loop catch scen
    for(state in c('low','base','high')){

      if(catch == 'ABC' & state == 'base'){
        rootdir.temp <- paste0(rootdir,"/cr",r )
        lastrun <- paste0(rootdir.temp,"/forecasts/forecast2030")
        newdir.temp <- paste0(rootdir,"/cr",r,"_",catch,"_",state)
        dir.create(newdir.temp) ## make special folder and copy files
        file.copy(list.files(lastrun,
                             full.names = TRUE,
                             recursive = TRUE),
                  to = newdir.temp, overwrite = TRUE)
          next('already ran ABC with base values, copied Forecast2030 into placeholder')
      }
      # if (catch == 'upper'){
      #   next('not yet implemented upper catch stream, skipping')
      # }
      df<-data.frame()
      rootdir.temp <- paste0(rootdir,"/cr",r )

      ## specific catch values here
      catch_projections <- read.csv(paste0(rootdir.temp,"/cproj_",r,".csv"))
      Flimitfraction <- catch_projections$PSTAR_0.45[catch_projections$YEAR ==2030]
      catch_proportions <- catch_projections[5,5:ncol(catch_projections)]
      const.catch <- mean(rowSums(catch_projections[3:4,5:ncol(catch_projections)])) ## avg 2019/2020
      fixed_catches <- catch_projections[1:4,5:ncol(catch_projections)]

      replist0 <- SS_output(paste0(rootdir.temp,"/base2015")) ## get values specific to this region

      lastrun <- paste0(rootdir.temp,"/forecasts/forecast2030")
      newdir.temp <- paste0(rootdir,"/cr",r,"_",catch,"_",state)
      dir.create(newdir.temp) ## make special folder and copy files
      file.copy(list.files(lastrun,
                           full.names = TRUE,
                           recursive = TRUE),
                to = newdir.temp, overwrite = TRUE)
      setwd(newdir.temp)

      ## Update Starter to read from CTL
      if(state != 'base'){
        strt <- SS_readstarter(file = "starter.ss")
        strt$init_values_src <- 0
        strt$last_estimation_phase <- 10 ## could go as high as 20
        SS_writestarter(strt, file = "starter.ss", overwrite = TRUE)
      }
      ## update CTL file with state of nature low/base/high (all fixed, reading from par)
      mctl <- readLines(list.files(newdir.temp)[grep('_control', list.files(newdir.temp))])
      LOI <- grep("NatM_p_1_Fem_GP_1",mctl)[1] ## get line(s) containing data after natm, ignoring comment
      NewLine <- strsplit(mctl[LOI],"   ") ## split elements

      NewLine[[1]][3] <- ifelse(state == 'low', 0.05, ifelse(state == 'high', 0.08, 0.07))
      mctl[LOI][1] = paste0(NewLine[[1]], collapse = " ")

      writeLines(text=mctl, con= paste(list.files(newdir.temp)[grep('_control', list.files(newdir.temp))])) ## save it
      ## read in and replace forecast file with appropriate catches ----
      ## REPLACE WITH ORIGINAL FORECAST FILE SINCE SS_READFORECAST CAN'T DEAL WITH OPTION #2 YET
      file.copy(from = paste0(rootdir.temp,"/base2015/forecast.ss"), to = paste0(newdir.temp,"/forecast.ss"), overwrite = TRUE)

      fore <- SS_readforecast(file = './forecast.ss',
                              Nareas = replist0$nareas,
                              Nfleets = replist0$nfishfleets,
                              version = paste(replist0$SS_versionNumeric),
                              readAll = TRUE)

      fore$Nforecastyrs <- 2031-replist0$endyr
      fore$FirstYear_for_caps_and_allocations <- forecast_start+(t-1)
      fore$Ncatch <- replist0$nfishfleets*(t+forecast_start-replist0$endyr-2)
      fore$InputBasis <- 2 ## discards

      ## Now Add Catch data/projections thru the year before forecast_start.
      ## This acts similarly to SS_ForeCatch except it reads directly from your inputs.
      inityr <- max(fore$ForeCatch$Year)
      for(k in 1:(forecast_start-1-inityr)){
        term <- nrow(fore$ForeCatch) ## intital final row
        for(i in 1:replist0$nfishfleets){
          fore$ForeCatch[term+i,'Year'] <- inityr+k
          fore$ForeCatch[term+i,'Seas'] <- 1
          fore$ForeCatch[term+i,'Fleet'] <- i
          fore$ForeCatch[term+i,'Catch_or_F'] <- fixed_catches[k,i]
        } ## end nfleets
      } ## end yrs to 2020
      ## Fix forecast file to end year selectivity
      fore$Bmark_years[1:6] <- 0
      fore$Fcast_years[1:4] <- 0
      ## Fix trawl relative F to reflect proportional catch amounts by fleet in forecast.
      fore$fleet_relative_F <- 2 ## will cause original r4ss write_forecast to fail
      fore$vals_fleet_relative_f <- paste(paste0(catch_proportions, collapse = " "))
      fore$basis_for_fcast_catch_tuning <- 2 ## dead biomass

      ##  Input correct buffer fraction for this year
      fore$Flimitfraction <- Flimitfraction
      mod1 <- SS_output(paste0(rootdir.temp,"/forecasts/forecast2021"), covar = FALSE) ## just load once
      predOFLs_startForecast <-  mod1$derived_quants[grep(paste0("ForeCatch_",(forecast_start+(t-2)),collapse = "|"), mod1$derived_quants$Label),"Value"]

      if(catch == 'ABC'){
        tempForeCatch <- SS_ForeCatch(mod1,
                                      yrs = 2021:(2021+(t-2)),
                                      average = FALSE,
                                      total = predOFLs_startForecast) ## will update based on changing buffer
        fore$ForeCatch[(nrow(fore$ForeCatch)+1):(nrow(fore$ForeCatch)+nrow(tempForeCatch)),] <- tempForeCatch[,1:4]
        write.csv(tempForeCatch, file = "./tempForeCatch.csv",row.names = FALSE) ## save final year ABC catch

      } else if(catch == 'constant'){
        tempForeCatch <- SS_ForeCatch(mod1,
                                      yrs = 2021:(2021+(t-2)),
                                      average = FALSE,
                                      total = const.catch)
        fore$ForeCatch[(nrow(fore$ForeCatch)+1):(nrow(fore$ForeCatch)+nrow(tempForeCatch)),] <- tempForeCatch[,1:4]
        write.csv(tempForeCatch, file = "./tempForeCatch.csv",row.names = FALSE) ## save final year ABC catch

        } else if (catch == 'upper'){
        upperStream <- 1.5*mod1$derived_quants[grep("ForeCatch_2021", mod1$derived_quants$Label),"Value"]
        tempForeCatch <- SS_ForeCatch(mod1,
                                      yrs = 2021:(2021+(t-2)),
                                      average = FALSE,
                                      total = upperStream)
        fore$ForeCatch[(nrow(fore$ForeCatch)+1):(nrow(fore$ForeCatch)+nrow(tempForeCatch)),] <- tempForeCatch[,1:4]
        write.csv(tempForeCatch, file = "./tempForeCatch.csv",row.names = FALSE) ## save final year ABC catch

        # next('not yet implemented upper catch stream, skipping')
      }
      ## save file
      SS_writeforecastMK(fore, file = './forecast.ss', overwrite = TRUE)

      ## execute this model
      setwd(newdir.temp); system('ss3 -nohess') ## works

      ## extract values of interest ----
      modX <- SS_output(newdir.temp, covar = FALSE)

      YOI <- (replist0$endyr+1):(forecast_end-1); lYOI <- length(YOI)
      ## this will read the output of the first model and save the OFLs
      ## which will get used to compute subsequent mods
      ## https://github.com/melmonk/StockAssessment_template/blob/master/8a_Tables.Rmd
      df[1:lYOI,"Year"] <- YOI
      df[1:lYOI,"PredOFL"] <-  modX$derived_quants[grep(paste0("OFLCatch_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"]
      df[1:lYOI,"ForeCatch_ABC"] <- modX$derived_quants[grep(paste0("ForeCatch_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"]
      df[1:lYOI,"SpawnBio"] <- modX$derived_quants[grep(paste0("SSB_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"]
      df[1:lYOI,"Depletion"] <- paste0(round(modX$derived_quants[grep(paste0("Bratio_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"],3)*100,"%")

      df$PredOFL[df$Year < forecast_start] <- df$ForeCatch_ABC[df$Year < forecast_start]<- NA
      df[,2:4] <- round(df[,2:4],2)
      df %>% mutate('REG' = r, 'State of Nature' = state, 'Catch' = catch) %>%
      write.csv(.,file = paste0(newdir.temp,"/summary_table_f.csv"),row.names = FALSE)



    } ## end states of nature
  } ## end catch scenarios
} ## end regions


require(stringr)
rootdir <- "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019"
regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
#
regnames <- basename(dirname(list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep('*summary_table_f.csv',list.files(regdirs,recursive = FALSE,full.names = TRUE))]))
df <- list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep('*summary_table_f.csv',list.files(regdirs,recursive = FALSE,full.names = TRUE))] %>%
 lapply(.,read.csv)%>% bind_rows()

write.csv(df,paste0(rootdir,"./",Sys.Date(),"table_f.csv"), row.names = TRUE)




## Table I ----
## base M is -2.94, low is -2.99, high is -2.41
for(r in c('North','Central','South')){ ## loop regions
  rootdir.temp <- paste0(rootdir,"/cr",r )
  lastrun <- paste0(rootdir.temp,"/forecasts/forecast2030")
  # SS_executivesummaryMK(dir = paste0(lastrun))
  ## make sure exec w hess -- can disable if certain
  # setwd(rootdir.temp); system("ss3")
  YOI <- 2010:2020
  basemod10 <- SS_output(lastrun, covar = TRUE)

  data.frame("SPR" = 1- basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Value"],
             "sd" =  basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"StdDev"],
             "Yr" = as.numeric(as.character(substr(basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Label"],10,14)))) %>%
    filter(Yr < 2021) %>%
    mutate(lwr = SPR-1.96*sd, upper = SPR+1.96*sd ) %>%
    ggplot(., aes(x= Yr, y = 1-SPR)) +
    theme_classic()+
    geom_point()+
    geom_hline(yintercept =1, col = 'red')+
    # scale_x_continuous(limits = c(1900, 2019), breaks = seq(1900,2020,20)) +
    scale_y_continuous(limits = c(0,1)) +
    geom_errorbar(aes(ymin = 1-lwr, ymax = 1-upper), width = 0, col = 'grey22' ) +
    labs(x = 'Year', y = paste(basemod10$SPRratioLabel), title = r)
  ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_1-SPR.png"),
         height = 4, width = 6, units = 'in', dpi = 420)

  df<- data.frame()
  for(y in 1:length(YOI)){
    df["CommLandings",y] <- NA
    df["TotCatch",y] <- basemod10$timeseries[, grepl('Yr|dead[(]B', names(basemod10$timeseries))] %>% filter(Yr == YOI[y]) %>% select(-Yr) %>% rowSums(.)
    df["OFL",y] <-  ifelse(length(basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]) == 0, NA, basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"])
    df["ACL",y] <-  NA

    df["1-SPR",y] <-  1- basemod10$derived_quants[grep(paste0("SPRratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
    df["ExploitationRate",y] <-  basemod10$derived_quants[grep(paste0("F_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]

    df["A10+Biomass",y] <- subset(basemod10$timeseries[, c('Yr', 'Bio_smry')], Yr == YOI[y])$Bio_smry
    df["SpawnBiomass",y] <-  basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
    df["Spawnbio95CI",y] <-  paste0( round( as.numeric(df["SpawnBiomass",y])-1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]),"--",
                                     round(as.numeric(df["SpawnBiomass",y])+1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]))

    df["Rec",y] <- round(basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],2)
    df["Rec95CI",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
                               round(as.numeric(df["Rec",y]) + 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
    # dq - qnorm(1-(1-quant)/2)*sd chantel this is 1.96

    df["Depletion",y] <- paste0(round(basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],3)*100,"%")
    df["Depl95CI",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
                                round(as.numeric(df["Rec",y] )+ 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))


  }
  colnames(df)<- YOI
  write.csv(df,paste0(rootdir.temp,"./",Sys.Date(),r,"table_i.csv"), row.names = TRUE)


  ## Figure E TS of Harvest rate: points with SE, red line with SPR50%
  SSplotSummaryF(basemod10,  print = T, Ftgt = basemod10$derived_quants[grep("Fstd_SPRtgt", basemod10$derived_quants$Label),"Value"], plotdir = rootdir.temp)
  SSplotCatch(basemod10, subplot = 2,  print = T, plotdir = rootdir.temp)

  SSplotTimeseries(basemod10, subplot = 11, print = T, plotdir = rootdir.temp, forecastplot = FALSE) ## recdevs
  SSplotTimeseries(basemod10, subplot = 9,  print = T, plotdir = rootdir.temp, forecastplot = FALSE)
  SSplotTimeseries(basemod10, subplot = 7,  print = T, plotdir = rootdir.temp, forecastplot = FALSE)
  ## by hand rel depl
  depl <- data.frame("depl" = basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Value"],
                     "upper" = basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Value"]+1.96*basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"StdDev"],
                     "lwr" = basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Value"]-1.96* basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"StdDev"]) %>%
    mutate('Yr' = as.numeric(as.character(substr(basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Label"],8,11)))) %>% filter(Yr < 2020)
  ggplot(depl, aes(x= Yr, y = depl)) +
    theme_classic()+
    geom_line(lwd =1.1) +
    geom_hline(yintercept = 0.4, col = 'red', lwd = 1.1)+
    geom_hline(yintercept = 0.25, col = 'red',lwd = 1.1)+
    geom_text(aes(x = 1915, y = 0.42, label = "Management Target")) +
    geom_text(aes(x = 1915, y = 0.26, label = "Minimum Stock Size Threshold")) +

    geom_ribbon(aes(ymin = lwr, ymax = upper), alpha = 0.2) +
    labs(x = 'Year', y = 'Relative Depletion', title = r)
  ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_Depl.png"), height = 4, width = 6, units = 'in', dpi = 420)

  SSplotSPR(basemod10,subplot = 4, print = TRUE, plotdir = rootdir.temp, forecastplot = FALSE) ## 1-SPRraw time series
  SSplotSPR(basemod10,subplot = 2, print = T, plotdir = rootdir.temp, forecastplot = FALSE, uncertainty = TRUE) ## 1-SPR kobe type


  # SSB_sd.Yr <- basemod10$derived_quants[grep("SSB_", basemod10$derived_quants$Label),"StdDev"][1:134]
  # basemod10$timeseries %>%
  #   select(Yr, SpawnBio) %>%
  #   mutate(upper = SpawnBio + 1.96*SSB_sd.Yr, lwr = SpawnBio - 1.96*SSB_sd.Yr) %>%
  #   filter(Yr < 2020) %>%
  # ggplot(., aes(x= Yr, y = SpawnBio)) +
  #   theme_classic()+
  #   geom_line(lwd =1.1) +
  #   # scale_x_continuous(limits = c(1900, 2019), breaks = seq(1900,2020,20)) +
  #   geom_vline(xintercept = 2019) +
  #   geom_ribbon(aes(ymin = lwr, ymax = upper), alpha = 0.2) +
  #   labs(x = 'Year', y = 'Spawning Biomass (000 mt)', title = r)
  # ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_SSBtimeseries.png"), height = 4, width = 6, units = 'in', dpi = 420)




}
## now merge them all - can also be used for table b
rootdir <- "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019"
regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
# regnames <- basename(regdirs)
regnames <- basename(dirname(list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep('table_i.csv',list.files(regdirs,recursive = FALSE,full.names = TRUE))]))
df <- list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep(paste0('table_i.csv'),list.files(regdirs,recursive = FALSE,full.names = TRUE))] %>%
  lapply(.,read.csv) %>%
  bind_rows()

for(i in 1:length(basename(regnames))){
  idx <- (1:13)+c(0,13,26)[i]
  df$REG[idx] <- sub("cr","",basename(regnames)[i])
}
write.csv(df,paste0(rootdir,"/",Sys.Date(),"summary_table_i.csv"), row.names = TRUE)


## Table E summary of reference points
rootdir <- "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/tgfiles"

for(r in c('North','Central','South')){ ## loop regions
  rootdir.temp <- paste0(rootdir,"/cr",r )
  lastrun <- paste0(rootdir.temp,"/base2015") ## if running base there needs to be a space between SPR_TARGET line 370
  # basemod1  <- SS_output(lastrun, covar = TRUE)
  SS_executivesummaryMK(dir = paste0(lastrun))
  # SS_executivesummaryMK(dir = lastrun)
  # SSB_sd.Yr <- basemod1$derived_quants[grep("SSB_", basemod1$derived_quants$Label),"StdDev"][1:nrow(basemod1$timeseries)]
  # basemod1$timeseries %>%
  #   select(Yr, SpawnBio) %>%
  #   mutate(upper = SpawnBio + 1.96*SSB_sd.Yr, lwr = SpawnBio - 1.96*SSB_sd.Yr) %>%
  #   ggplot(., aes(x= Yr, y = SpawnBio)) +
  #   theme_classic()+
  #   geom_line(lwd =1.1) +
  #   # scale_x_continuous(limits = c(1900, 2019), breaks = seq(1900,2020,20)) +
  #   # scale_y_continuous(limits = c(0,), breaks = seq(0,90,15)) +
  #   geom_ribbon(aes(ymin = lwr, ymax = upper), alpha = 0.2) +
  #   labs(x = 'Year', y = 'Spawning Biomass (000 mt)', title = r)
  # ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_SSBtimeseries.png"), height = 4, width = 6, units = 'in', dpi = 420)
}
#   df<- data.frame()
#   for(i in c("Value","95CI")){
#     if(i == '95CI'){
#
#
#
#     }
#
#   df["UnfishedBio",i] <- NA
#   df["Unfished_A10",y] <- subset(basemod10$timeseries[, c('Yr', 'Bio_smry')], Yr == YOI[y])$Bio_smry
#   df["Unfished_Recr",y] <-  ifelse(length(basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]) == 0, NA, basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"])
#   df["SpBio_2015",y] <-  NA
#
#   df["ProxySB40",y] <-  1- basemod10$derived_quants[grep(paste0("SPRratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
#   df["SPR_SB40",y] <-  basemod10$derived_quants[grep(paste0("F_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
#
#   df["F_SB40",y] <-
#   df["YieldwSPR_SB40",y]
#
#   df["SpawnBio_2019",y]  <-  basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
#   df["SpawnBioCI",y] <-  paste0( round( as.numeric(df["SpawnBiomass",y])-1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]),"--",
#                                    round(as.numeric(df["SpawnBiomass",y])+1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]))
#
#   df["SB_MSY",y] <- round(basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],2)
#   df["SPR_MSY",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
#                              round(as.numeric(df["Rec",y]) + 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
#   # dq - qnorm(1-(1-quant)/2)*sd chantel this is 1.96
#
#   df["F_SPRY_MSY",y] <- paste0(round(basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],3)*100,"%")
#   df["MSY",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
#                               round(as.numeric(df["Rec",y] )+ 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
#
#
# }
sumlist<-list.dirs("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/mods", recursive = F) %>%
lapply(.,SS_output) %>%
  SSsummarize(.)

sumlist$SPRratio %>% spread(Yr,Label)
df0<- with(sumlist, cbind(melt(sumlist$SPRratio,id='Yr',value.name = 'SPRratio'),
                    melt(sumlist$SPRratioLower,id='Yr',value.name = 'SPRratioLower'),
                    melt(sumlist$SPRratioUpper,id='Yr',value.name = 'SPRratioUpper')))
df0$SPRratio <- as.numeric(as.character(df0$SPRratio))
df0$SPRratioLower <- as.numeric(as.character(df0$SPRratioLower))
df0$SPRratioUpper <- as.numeric(as.character(df0$SPRratioUpper))
names(df0)[c(2,4,5,7,8)] <- c('model','yx','vx','yy','vy')
df0 <- df0 %>%  select(-vx,-vy, -yx,-yy) %>%
  filter(model != 'Label')
levels(df0$model) <- c('Central','North','South','Label')


  # data.frame("SPR" = 1- basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Value"],
  #            "sd" =  basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"StdDev"],
  #            "Yr" = as.numeric(as.character(substr(basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Label"],10,14)))) %>%
png("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/sprratio_all.png",
    width = 8, height = 6, units = 'in', res = 520)
df0 %>%
filter(Yr < 2020) %>%
  # mutate(lwr = SPR-1.96*sd, upper = SPR+1.96*sd ) %>%
  ggplot(., aes(x= Yr, y = SPRratio, color = model)) +
  theme_classic()+
  theme(legend.position = c(0.2,0.9),
        legend.text = element_text(size = 22),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  geom_point() +
  geom_hline(yintercept =1, col = 'red')+
  scale_color_brewer(palette = 'Dark2')+
  scale_x_continuous(limits = c(1900, 2021), breaks = seq(1900,2020,20)) +
  scale_y_continuous(limits = c(0,2)) +
  geom_errorbar(aes(ymin = SPRratioLower, ymax = SPRratioUpper, col = model ), width = 0) +
  labs(x = 'Year', y = paste(basemod10$SPRratioLabel), color = '')
dev.off()# ggsave(plot = last_plot(), file = "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/sprratio_all.png")
#   SSplotComparisons(., print = TRUE, shadeForecast = TRUE, legendlabels = c('Central','North','South'),
#                     plotdir = "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/plots")
mod.cols = c("#7570B3" ,"#D95F02","#1B9E77")
png("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/yield_curve.png",
    width = 8, height = 6, units = 'in', res = 520)
SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/mods/South_forecast2030") %>%
  SSplotYield(., col= mod.cols[1], subplot=1)
grid()
SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/mods/North_forecast2030") %>%
SSplotYield(., col= mod.cols[2], subplot=1,add = TRUE)

SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/mods/Central_forecast2030") %>%
SSplotYield(., col= mod.cols[3], subplot=1, add=TRUE)
legend('topright', legend=c('South','North','Central'), col=mod.cols, lwd=3, bg='white', bty='n', cex = 3)

dev.off()
