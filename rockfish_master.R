require(dplyr)
require(readr)
require(purrr)
require(kaputils)

## Exeute automated forecasts
fixed
for(r in c('North','Central','South')){
  rootdir.temp <- paste0("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/cr",r)
  catch_projections <- read.csv(paste0(rootdir.temp,"/cproj_",r,".csv"))
  SS_autoForecast(rootdir = rootdir.temp,
                  basedir = "base2015",
                  catch_proportions = catch_projections[5,5:ncol(catch_projections)],
                  # catch_proportions = c(0.5,0.08426184,0.4157382),
                  forecast_start = 2021,
                  forecast_end = 2031,
                  fixed_catches = catch_projections[1:4,5:ncol(catch_projections)],
                  Flimitfraction = catch_projections$PSTAR_0.45[catch_projections$YEAR >2020])
  read.csv(paste0(rootdir,"/forecasts/decision_table_base.csv"))

}

## Stich together all regions for output tabb
rootdir <- "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019"
regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
regnames <- basename(regdirs)

df <- list.files(regdirs, recursive = TRUE, full.names = TRUE)[grep('decision_table_base.csv',list.files(regdirs,recursive = TRUE,full.names = TRUE))] %>%
  map_df(~ read_csv(.,
                    col_types = cols(.default = "c")),
         .id="index") %>%
  mutate( REG = gsub( "cr", "", regnames[as.numeric(index)] ))
write.csv(df,paste0(rootdir,"./table_g.csv"), row.names = FALSE)

## fixed catch for next decision table is average of 2019/2020 input
## I am going to hard code this -- not manual but not function
const.catch <- mean(rowSums(catch_projections[3:4,5:ncol(catch_projections)]))
