## Script to automate catch-only forecasting updates
## Maia Sosa Kapur kapurm@uw.edu

require(r4ss)
require(dplyr)
# require(ThorsonUtilities)
nforeyrs = 3

## set your WD somewhere you'd like several models
rootdir <- "C:/Users/MKapur/Dropbox/UW/chinarock_update"
source(paste0(rootdir,"/SS_writeforecastMK.R"))
catch_projections <- read.csv("./china_north_cproj.csv")
## Step 1. First, take the base model input files from the assessment 
# This includes: .dat, .ctl, forecast, starter, ss3.exe, .par
# (I will provide these, and you will also need: Report.sso)
basedir <- "base2015"#; system("ss3") ## if you haven't executed yet
mod <- SS_output(paste0(rootdir,"/",basedir))
## Check spawning biomass and total like (there are other ways to extract)
# mod$derived_quants[mod$derived_quants$Label == "SSB_Virgin",]
# mod$derived_quants[mod$derived_quants$Label == "SSB_2013",]
# mod$likelihoods_used[[1]]


## Step 2 Change first numerical term in Starter file to 1
## I like to copy into a new folder just in case -- this will paste EVERYTHING from basedir into base_1
df <- data.frame('ForeCatch' =NA,'OFLCatch' = NA, 'Iter' = NA)
for(t in 1:nforeyrs){
  base_temp <- paste0("forecasts/forecast", (t-1)+2021)
  setwd(rootdir); dir.create(base_temp)
  file.copy(list.files(
    ifelse(t == 1, paste0(rootdir, "/", basedir), paste0("forecasts/forecast", (t-2)+2021)),
    full.names = TRUE,
    recursive = TRUE), base_temp)
  
  
  ## execute the changed model -- on a mac use "shell" instead of "system"
  setwd(base_temp)
  strt <- SS_readstarter(file = "starter.ss")
  strt$init_values_src <- 1 
  strt$last_estimation_phase <- 10 ## could go as high as 20
  SS_writestarter(strt, file = "starter.ss", overwrite = TRUE)
  # system("ss3 -nohess") ## testing OK
  
  ## Steps 3. Par File - add zeroes
  file.copy("C:/Users/mkapur/Dropbox/UW/chinarock_update/base2015_inputfiles/ss3.par",
            "ss3.par", overwrite = TRUE)
  
  mpar <- readLines("ss3.par")
  LOI <- grep("Fcast",mpar)+1 ## get line(s) containing data after fcast
  NewLine <- strsplit(mpar[LOI],"0 ") ## split elements
  length(NewLine[[1]]);length(NewLine[[2]])
  for(a in 1:length(NewLine)){
    ltemp <- length(NewLine[[a]])
    NewLine[[a]][(ltemp+1):(ltemp+5)] <- 0.000000000000 ## add 6 zeroes to end
    mpar[LOI][a] = paste0(NewLine[[a]], collapse = " ")
  }
  NewLine <- strsplit(mpar[LOI],"0 ") ## split elements
  length(NewLine[[1]]);length(NewLine[[2]])
  writeLines(text=mpar, con="ss3.par") ## save it
  system("ss3 -nohess") ## testing 
  
  
  
  ## Step 4a. Add catch/projections through 2020.
  fore <- SS_readforecast(file = './forecast.ss',
                          Nareas = mod$nareas,
                          Nfleets = mod$nfleets, 
                          version = '3.24',
                          readAll = TRUE)
  fore$Nforecastyrs <- 2031-mod$endyr
  fore$FirstYear_for_caps_and_allocations <- 2021+(t-1)
  # fore$N_forecast_loops <- 24+(3*t) ## check this
  fore$Ncatch <- 15+(3*t) ## check this
  fore$InputBasis <- 2
  
  ## Now Add Catch data/projections thru 2020 if this is base_1. 
  ## Otherwise this gets copied.
  if(t == 1){
    inityr <- max(fore$ForeCatch$Year)
    for(k in 1:(2020-inityr)){
      term <- nrow(fore$ForeCatch) ## intital final row
      for(i in 1:mod$nfleets){
        fore$ForeCatch[term+i,'Year'] <- inityr+k
        fore$ForeCatch[term+i,'Seas'] <- 1
        fore$ForeCatch[term+i,'Fleet'] <- i
        fore$ForeCatch[term+i,'Catch_or_F'] <- catch_projections[catch_projections$YEAR == inityr+k,5+(i-1)]
      } ## end nfleets
    } ## end yrs to 2020
  }
  
  # SS_ForeCatch(mod,
  #              yrs = inityr:2022, average = TRUE,
  #              avg.yrs = 2010,
  #              total = c(2,2))
  ## Steps 4.b. fix forecast file to end year selectivity
  fore$Bmark_years[1:6] <- 0
  fore$Fcast_years[1:4] <- 0
  
  ## Step 4.d. Fix trawl relative F to reflect proportional catch amounts by fleet in forecast. 
  fore$fleet_relative_F <- 2 ## THIS IS WRONG CHANGE MANUALLY
  fore$vals_fleet_relative_f <- paste(0.5,0.084262,0.4157381616)
  fore$basis_for_fcast_catch_tuning <- 2 ## deadbiomass
  
  ## Step 5a. Input correct 2021 buffer fraction and run model forecast first time.
  fore$Flimitfraction <- ifelse(is.na(catch_projections$PSTAR_0.45[catch_projections$Year == 2021+(t-1)]), 0.913,
                                catch_projections$PSTAR_0.45[catch_projections$Year == 2021+(t-1)])
  
  ## save file
  # writeLines(text=paste(fore), con="forecast.ss") ## save it
  SS_writeforecastMK(fore, file = './forecast.ss', overwrite = TRUE)
  ## execute model forecast first time
  if(t==1) system('ss3 -nohess')
  
  # Step 5b. Iterate the forecast file
  ## Find the total forecasted catch and OFL in 2021 in the Report.sso file

  ## Allocate this catch among the fleets according to the given proportions; add this to forecast file
  # catch.prop <- if(fore$N_allocation_groups == 0){rep(1/mod$nfleets,3)}
  if(t != 1){
    term <- nrow(fore$ForeCatch) ## intital final row
    for(i in 1:mod$nfleets){
      fore$ForeCatch[term+i,'Year'] <- 2021+(t-1)
      fore$ForeCatch[term+i,'Seas'] <- 1
      fore$ForeCatch[term+i,'Fleet'] <- i
      fore$ForeCatch[term+i,'Catch_or_F'] <- catch_projections[catch_projections$YEAR == 2021+(t-1),5+(i-1)]
    }
    
    cat(paste0('Added forecast catch for year',2021+(t-1),"\n"))
  }
  mod <- SS_output(base_temp)
  Fore202X <- mod$derived_quants[mod$derived_quants$Label == paste0("ForeCatch_",2021+(t-1)),"Value"]
  OFL202X <- mod$derived_quants[mod$derived_quants$Label == paste0("OFLCatch_",2021+(t-1)),"Value"]
  df['ForeCatch',t] <- Fore202X;   df['OFLCatch',t] <- OFL202X; df['Iter',t] <- t
  if(t!=1) system('ss3 -nohess -nox');    cat(paste0('Executed model with forecast thru year',2021+(t-1),"\n"))

  # Step 5c. Iterate through 2030 -- the loop will continue making a new folder each time
  
}