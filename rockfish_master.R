require(dplyr)
require(readr)
require(purrr)
# devtools::install_deps("C:/Users/mkapur/Dropbox/kaputils")
# devtools::install_github("mkapur/kaputils", dependencies = F)
library(kaputils)
# devtools::install_github("r4ss/r4ss@2663227")
library(r4ss)
# devtools::source_url("https://raw.githubusercontent.com/mkapur/kaputils/master/R/SS_readforecastMK.R") ## use dev version
devtools::source_url("https://raw.githubusercontent.com/r4ss/r4ss/development/R/SS_ForeCatch.R") ## use dev version
# devtools::source_url("https://raw.githubusercontent.com/mkapur/kaputils/master/R/SS_writeforecastMK.R") ## use dev version
# devtools::source_url("https://raw.githubusercontent.com/mkapur/kaputils/master/R/SS_executivesummaryMK.R")
# source("./R/SS_executivesummaryMK.R")
compname <- c('mkapur',"Maia Kapur")[2]


## Run ABC x 3 States [low/base/high] ----
## Execute automated forecasts for three regions and three states of nature. Takes a bit, esp. hessian for Central
## will only invoke CTL par changes when state !=base
# for(r in c('North','Central','South')){
#   for(state in c('base','low','high')){
#     
#     rootdir.temp <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/cr",r,"_ABC_",state)
#     # newdir.temp <- paste0(rootdir,"/cr",r,"_",catch,"_",state)
#     if(state == 'base'){ ## only read in once per region
#     catch_projections <- read.csv(paste0(rootdir.temp,"/cproj_",r,".csv"))
#     }
#     kaputils:::SS_autoForecast(rootdir = rootdir.temp,
#                                basedir = "base2015",
#                                catch_proportions = catch_projections[7,5:ncol(catch_projections)],
#                                state = state,
#                                forecast_start = 2021,
#                                forecast_end = 2031,
#                                fixed_catches = catch_projections[1:4,5:ncol(catch_projections)],
#                                Flimitfraction = catch_projections$PSTAR_0.45[catch_projections$YEAR >2020])
#     # read.csv(paste0(rootdir,"/forecasts/decision_table_base.csv"))
#   }
# }

## sanity check: these should be different. the low state should be less productive.
# modb<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crNorth_ABC_base/forecasts/forecast2030")
# modh<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crNorth_ABC_high/")
# modl<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crNorth_ABC_low/")
# 
# modb$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"] ## start around 18
# modh$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"] ## start around 50s
# modl$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modl$derived_quants$Label),"Value"] ## 10s
# 
# # modh21<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crNorth_ABC_high/")
# 
# # modh2<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/dep/forecasts_0530/crNorth_ABC_high/")
# # modh221<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/dep/forecasts_0530/crNorth_ABC_high/")
# modh<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crCentral_ABC_high/")
# modl<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crCentral_ABC_low/")
# modh$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"] ## 109 
# modl$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modl$derived_quants$Label),"Value"] ## 20

# modh2$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"]

# modh21$derived_quants[grep(paste0("OFLCatch_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"]
# modh$derived_quants[grep(paste0("OFLCatch_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"]
# modh2$derived_quants[grep(paste0("OFLCatch_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"]

# modh<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crSouth_ABC_high/")
# modl<- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crSouth_ABC_low/")
# modh$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modh$derived_quants$Label),"Value"] ## 22
# modl$derived_quants[grep(paste0("SSB_",2019:2030,collapse = "|"), modl$derived_quants$Label),"Value"] ## 14


## TABLE G stitching ----
# rootdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/tgfiles")
# regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
# regnames <- basename(regdirs)
# 
# df <- list.files(regdirs, recursive = TRUE, full.names = TRUE)[grep('decision_table_base.csv',list.files(regdirs,recursive = TRUE,full.names = TRUE))] %>%
#   map_df(~ read_csv(.,
#                     col_types = cols(.default = "c")),
#          .id="index") %>%
#   mutate( REG = gsub( "cr", "", regnames[as.numeric(index)] ))
# write.csv(df,paste0(rootdir,"./",Sys.Date(),"table_g.csv"), row.names = FALSE)

## Run Const/Upper catch x 3 states ----
## Did these separately because the constant nature of the catch requires 
## non-iteration
rootdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/")
forecast_start <- 2021; forecast_end <- 2031; t = 10
## base M is -2.94, low is -2.99, high is -2.41
for(r in c('North','Central','South')){ ## loop regions
  for(catch in c('constant','upper')){ ## loop catch scen
    for(state in c('low','base','high')){

      df<-data.frame()
      # rootdir.temp <- paste0(rootdir,"cr",r,"_ABC_",state) ## copy region & state
      catch_projections <- read.csv(paste0(rootdir,"/cr",r,"_ABC_base/cproj_",r,".csv")) ## from base dir
      Flimitfraction <- catch_projections$PSTAR_0.45[catch_projections$YEAR ==2030]
      catch_proportions <- catch_projections[7,5:ncol(catch_projections)]
      const.catch <- mean(rowSums(catch_projections[catch_projections$TYPE == 'PROJECTION',5:ncol(catch_projections)])) ## avg 2019/2020
      fixed_catches <- catch_projections[catch_projections$TYPE == 'ACTUAL',5:ncol(catch_projections)]

      # replist0 <- SS_output(paste0(rootdir,"cr",r,"_ABC_base/base2015")) ## get values specific to this region
      if(state != 'base'){
      lastrun <- paste0(rootdir,"cr",r,"_ABC_",state)
      } else if(state == 'base'){
        lastrun <- paste0(rootdir,"cr",r,"_ABC_",state,"/forecasts/forecast2030")
        
      }
      mod1 <- SS_output(lastrun, covar = FALSE, hidewarn = T, verbose = F) ## just load once for structure this hasn't executed yet
      
      newdir.temp <- paste0(rootdir,"cr",r,"_",catch,"_",state)
      dir.create(newdir.temp) ## make special folder and copy files
      file.copy(list.files(lastrun,
                           full.names = TRUE,
                           recursive = TRUE),
                to = newdir.temp, overwrite = TRUE)
      setwd(newdir.temp) ## now forecast2030 appropriate to state is replicated here

      ## Update Starter to read from CTL
      # if(state != 'base'){
      #   strt <- SS_readstarter(file = "starter.ss")
      #   strt$init_values_src <- 0 !boom
      #   strt$last_estimation_phase <- 10 ## could go as high as 20
      #   SS_writestarter(strt, file = "starter.ss", overwrite = TRUE)
      # }
      # ## update CTL file with state of nature low/base/high (all fixed, reading from par)
      # mctl <- readLines(list.files(newdir.temp)[grep('_control', list.files(newdir.temp))])
      # LOI <- grep("NatM_p_1_Fem_GP_1",mctl)[1] ## get line(s) containing data after natm, ignoring comment
      # NewLine <- strsplit(mctl[LOI],"   ") ## split elements
      # 
      # NewLine[[1]][3] <- ifelse(state == 'low', 0.05, ifelse(state == 'high', 0.08, 0.07))
      # mctl[LOI][1] = paste0(NewLine[[1]], collapse = " ")
      # 
      # writeLines(text=mctl, con= paste(list.files(newdir.temp)[grep('_control', list.files(newdir.temp))])) ## save it

      ## read in and replace forecast file with appropriate catches 
      ## REPLACE WITH ORIGINAL FORECAST FILE SINCE SS_READFORECAST CAN'T DEAL WITH OPTION #2 YET
      # file.copy(from = paste0(rootdir.temp,"/base2015/forecast.ss"), to = paste0(newdir.temp,"/forecast.ss"), overwrite = TRUE)

      ## only need to change catches in forecast file
      fore <- SS_readforecastMK(file = './forecast.ss',
                                Nareas = mod1$nareas,
                                Nfleets = mod1$nfishfleets,
                                nseas = 1,
                                version = paste(mod1$SS_versionNumeric),
                                readAll = TRUE)

      # fore$Nforecastyrs <- 2031-replist0$endyr
      # fore$FirstYear_for_caps_and_allocations <- forecast_start+(t-1)
      # fore$Ncatch <- replist0$nfishfleets*(t+forecast_start-replist0$endyr-2)
      # fore$InputBasis <- 2 ## discards

      ## Now Add Catch data/projections thru the year before forecast_start.
      ## This acts similarly to SS_ForeCatch except it reads directly from your inputs.
      # inityr <- max(fore$ForeCatch$Year)
      # for(k in 1:(forecast_start-1-inityr)){
      #   term <- nrow(fore$ForeCatch) ## intital final row
      #   for(i in 1:replist0$nfishfleets){
      #     fore$ForeCatch[term+i,'Year'] <- inityr+k
      #     fore$ForeCatch[term+i,'Seas'] <- 1
      #     fore$ForeCatch[term+i,'Fleet'] <- i
      #     fore$ForeCatch[term+i,'Catch_or_F'] <- fixed_catches[k,i]
      #   } ## end nfleets
      # } ## end yrs to 2020
      ## Fix forecast file to end year selectivity
      # fore$Bmark_years[1:6] <- 0
      # fore$Fcast_years[1:4] <- 0
      # ## Fix trawl relative F to reflect proportional catch amounts by fleet in forecast.
      # fore$fleet_relative_F <- 2 ## will cause original r4ss write_forecast to fail
      fore$vals_fleet_relative_f <- paste(paste0(catch_proportions, collapse = " "))
      # fore$basis_for_fcast_catch_tuning <- 2 ## dead biomass

      ##  Input correct buffer fraction for this year
      # fore$Flimitfraction <- Flimitfraction
      # fore$ControlRuleMethod <- 3 ## 3: ramp does catch=f(SSB), buffer on catch
      
      # mod1 <- SS_output(paste0(rootdir.temp,"/forecasts/forecast2021"), covar = FALSE) ## just load once
      
      # predOFLs_startForecast <-  mod1$derived_quants[grep(paste0("ForeCatch_",(forecast_start+(t-2)),collapse = "|"), mod1$derived_quants$Label),"Value"]

      # if(catch == 'ABC'){
      #   tempForeCatch <- SS_ForeCatch(mod1,
      #                                 yrs = 2021:(2021+(t-2)),
      #                                 average = FALSE,
      #                                 total = predOFLs_startForecast) ## will update based on changing buffer
      #   fore$ForeCatch[(nrow(fore$ForeCatch)+1):(nrow(fore$ForeCatch)+nrow(tempForeCatch)),] <- tempForeCatch[,1:4]
      # 
      # 
      #   write.csv(tempForeCatch, file = "./tempForeCatch.csv",row.names = FALSE) ## save final year ABC catch
      # 
      # } else
        
      if(catch == 'constant'){
        tempForeCatch <- SS_ForeCatch(mod1,
                                      yrs = 2021:(2021+(t-2)),
                                      average = FALSE,
                                      total = const.catch)
        fore$ForeCatch[min(which((fore$ForeCatch$Year==2021))):nrow(fore$ForeCatch),]<- tempForeCatch[,1:4]
        # fore$ForeCatch[(nrow(fore$ForeCatch)+1):(nrow(fore$ForeCatch)+nrow(tempForeCatch)),] <- tempForeCatch[,1:4]
        write.csv(tempForeCatch, file = "./tempForeCatch.csv",row.names = FALSE) ## save final year ABC catch
        
      } else if (catch == 'upper'){
        upperStream <- 1.5*mod1$derived_quants[grep("ForeCatch_2021", mod1$derived_quants$Label),"Value"]
        tempForeCatch <- SS_ForeCatch(mod1,
                                      yrs = 2021:(2021+(t-2)),
                                      average = FALSE,
                                      total = upperStream)
        fore$ForeCatch[min(which((fore$ForeCatch$Year==2021))):nrow(fore$ForeCatch),]<- tempForeCatch[,1:4]
        # fore$ForeCatch[(nrow(fore$ForeCatch)+1):(nrow(fore$ForeCatch)+nrow(tempForeCatch)),] <- tempForeCatch[,1:4]
        write.csv(tempForeCatch, file = "./tempForeCatch.csv",row.names = FALSE) ## save final year ABC catch
      }
      ## save file
      SS_writeforecastMK(fore, file = './forecast.ss', overwrite = TRUE)
      ## execute this model
      setwd(newdir.temp); system('ss3 -nohess') ## works

      ## extract values of interest ----
      # modX <- SS_output(newdir.temp, covar = FALSE)

      # YOI <- (replist0$endyr+1):(forecast_end-1); lYOI <- length(YOI)
      # ## this will read the output of the first model and save the OFLs
      # ## which will get used to compute subsequent mods
      # ## https://github.com/melmonk/StockAssessment_template/blob/master/8a_Tables.Rmd
      # df[1:lYOI,"Year"] <- YOI
      # df[1:lYOI,"PredOFL"] <-  modX$derived_quants[grep(paste0("OFLCatch_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"]
      # df[1:lYOI,"ForeCatch_ABC"] <- modX$derived_quants[grep(paste0("ForeCatch_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"]
      # df[1:lYOI,"SpawnBio"] <- modX$derived_quants[grep(paste0("SSB_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"]
      # df[1:lYOI,"Depletion"] <- paste0(round(modX$derived_quants[grep(paste0("Bratio_",YOI,collapse = "|"), modX$derived_quants$Label),"Value"],3)*100,"%")
      # 
      # df$PredOFL[df$Year < forecast_start] <- df$ForeCatch_ABC[df$Year < forecast_start]<- NA
      # df[,2:4] <- round(df[,2:4],2)
      # df %>% mutate('REG' = r, 'State of Nature' = state, 'Catch' = catch) %>%
      # write.csv(.,file = paste0(newdir.temp,"/summary_table_f.csv"),row.names = FALSE)



    } ## end states of nature
  } ## end catch scenarios
} ## end regions

## Creation of objects for use in template ----
## mod RData
modN <- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crNorth_ABC_base/forecasts/forecast2030")
modC <- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crCentral_ABC_base/forecasts/forecast2030")
modS <- SS_output("C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/crSouth_ABC_base/forecasts/forecast2030")
save(modN,modC,modS, file = "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/r4ss/China_SS_output2019.RData")

## Comparison plots
load("./r4ss/China_SS_output2019.RData")
list(modN, modC, modS) %>% SSsummarize(.) %>% 
  SSplotComparisons(shadeForecast = TRUE,  
                    endyrvec = 2030,
                    legendlabels = c("North",'Central','South'),
                    png = T, print = T, 
                    plotdir = "C:/Users/Maia Kapur/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/r4ss/plots_compare")


mod.cols = c("#7570B3" ,"#D95F02","#1B9E77")

png(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/r4ss/plots_compare/yield_comparison_3_models.png"),
    width = 8, height = 6, units = 'in', res = 520)
SSplotYield(modS, col= mod.cols[1], subplot=1)
grid()
SSplotYield(modN, col= mod.cols[2], subplot=1,add = TRUE)
SSplotYield(modC, col= mod.cols[3], subplot=1, add=TRUE)
legend('topright', legend=c('South','North','Central'), col=mod.cols, lwd=3, bg='white', bty='n', cex = 1.5)
dev.off()

## Build decision table (not in SS_executive summary) ----
YOI <- 2021:2030
for(r in c('North','Central','South')){ ## loop regions
  dec_table <- matrix(NA, nrow = length(YOI)*3, ncol = 9)
  dec_table <- data.frame(dec_table)
  names(dec_table) <- c('Scenario','Year','catch',paste(c("spawnbio","depl"),rep(c('low','base','high'),each = 2)))
  dec_table$Year <- rep(YOI,3)
  idxr <- idxc <- 1
  for(catch in c('constant','ABC','upper')){ ## loop catch scen
    idxc <- 1 ## reset to initial column for new catch scenario
    for(state in c('low','base','high')){
      if(catch != 'ABC'  | state != 'base'){
        tempdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/cr",r,"_",catch,"_",state)
      } else if(catch == 'ABC'){
        tempdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/cr",r,"_",catch,"_",state,"/forecasts/forecast2030")
      }
      mod <- SS_output(tempdir, covar = F)
      if(catch == 'constant' & idxc ==2){
        catch_projections <- read.csv(paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/cr",r,"_ABC_base/cproj_",r,".csv"))
        const.catch <- mean(rowSums(catch_projections[catch_projections$TYPE == 'PROJECTION',5:ncol(catch_projections)])) ## avg 2019/2020
        dec_table$catch[idxr:(idxr+length(YOI)-1)] <- round(const.catch,2)
        
      } else if(catch == 'upper' & idxc ==2){
        upperStream <- 1.5*mod$derived_quants[grep("ForeCatch_2021", mod$derived_quants$Label),"Value"]
        
        dec_table$catch[idxr:(idxr+length(YOI)-1)] <- round(upperStream,2)
        
      } else if (catch == 'ABC' &  idxc ==2){
        catchvals <- mod$timeseries[, grepl('Yr|dead[(]B', names(mod$timeseries))] %>% 
          filter(Yr %in% YOI) %>%
          select(-Yr) %>% rowSums(.) %>% round(.,2)
        
        dec_table$catch[idxr:(idxr+length(YOI)-1)] <- round(catchvals,2)
        
      }
      # read.csv(paste0(tempdir,"/tempForeCatch.csv"))
    

      ## input what was given to forecast file
      dec_table$Scenario[idxr:(idxr+length(YOI)-1)] <- rep(catch, length(idxr:(idxr+length(YOI)-1)))

      
      dec_table[idxr:(idxr+length(YOI)-1),idxc*2+2] <-  mod$derived_quants[grep(paste0("SSB_",YOI,collapse = "|"),
                                                                                mod$derived_quants$Label),"Value"]
      dec_table[idxr:(idxr+length(YOI)-1),idxc*2+3] <-  mod$derived_quants[grep(paste0("Bratio_",YOI,collapse = "|"),
                                                                                mod$derived_quants$Label),"Value"]
      idxc <- idxc+1 ## move to next set of columns as state updates
    # idxc <- idxc+3; idxr <-
    #     df["Depletion",y] <- paste0(round(basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],3)*100,"%")
    
    
  } ## end state
  idxr <- idxr+length(YOI) ## jump down to next set of years when catch scenario updates
  } ## end catch
  ## rename to look nice
  dec_table$Scenario[dec_table$Scenario == 'constant'] <- c(rep(" ",5),'Constant (2019-2020 Average)',rep(" ",5))
  dec_table$Scenario[dec_table$Scenario == 'ABC'] <- c(rep(" ",5),'40-10 Rule',rep(" ",5))
  dec_table$Scenario[dec_table$Scenario == 'upper'] <- c(rep(" ",5),'Upper Stream',rep(" ",5))
  ## save dec_table
  write.csv(dec_table, 
            file = paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/txt_files/decision_table_",
                          r,".csv"),
                          row.names = F)
} ## end regions



# require(stringr)
# rootdir <- paste0("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019")
# regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
# #
# regnames <- basename(dirname(list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep('*summary_table_f.csv',list.files(regdirs,recursive = FALSE,full.names = TRUE))]))
# df <- list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep('*summary_table_f.csv',list.files(regdirs,recursive = FALSE,full.names = TRUE))] %>%
#  lapply(.,read.csv)%>% bind_rows()
# # 
# write.csv(df,paste0(rootdir,"./",Sys.Date(),"table_f.csv"), row.names = TRUE)
# 
# 
# 
# 
# ## Table I ----
# sumlist<-list.dirs("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/mods", recursive = F) %>%
# lapply(.,SS_output) %>%
  # SSsummarize(.)
# ## base M is -2.94, low is -2.99, high is -2.41
# for(r in c('North','Central','South')){ ## loop regions
#   rootdir.temp <- paste0(rootdir,"/cr",r )
#   lastrun <- paste0(rootdir.temp,"/forecasts/forecast2030")
#   # SS_executivesummaryMK(dir = paste0(lastrun))
#   ## make sure exec w hess -- can disable if certain
#   # setwd(rootdir.temp); system("ss3")
#   YOI <- 2010:2020
#   basemod10 <- SS_output(lastrun, covar = TRUE)
# 
#   data.frame("SPR" = 1- basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Value"],
#              "sd" =  basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"StdDev"],
#              "Yr" = as.numeric(as.character(substr(basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Label"],10,14)))) %>%
#     filter(Yr < 2021) %>%
#     mutate(lwr = SPR-1.96*sd, upper = SPR+1.96*sd ) %>%
#     ggplot(., aes(x= Yr, y = 1-SPR)) +
#     theme_classic()+
#     geom_point()+
#     geom_hline(yintercept =1, col = 'red')+
#     # scale_x_continuous(limits = c(1900, 2019), breaks = seq(1900,2020,20)) +
#     scale_y_continuous(limits = c(0,1)) +
#     geom_errorbar(aes(ymin = 1-lwr, ymax = 1-upper), width = 0, col = 'grey22' ) +
#     labs(x = 'Year', y = paste(basemod10$SPRratioLabel), title = r)
#   ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_1-SPR.png"),
#          height = 4, width = 6, units = 'in', dpi = 420)
# 
#   df<- data.frame()
#   for(y in 1:length(YOI)){
#     df["CommLandings",y] <- NA
#     df["TotCatch",y] <- basemod10$timeseries[, grepl('Yr|dead[(]B', names(basemod10$timeseries))] %>% filter(Yr == YOI[y]) %>% select(-Yr) %>% rowSums(.)
#     df["OFL",y] <-  ifelse(length(basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]) == 0, NA, basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"])
#     df["ACL",y] <-  NA
# 
#     df["1-SPR",y] <-  1- basemod10$derived_quants[grep(paste0("SPRratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
#     df["ExploitationRate",y] <-  basemod10$derived_quants[grep(paste0("F_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
# 
#     df["A10+Biomass",y] <- subset(basemod10$timeseries[, c('Yr', 'Bio_smry')], Yr == YOI[y])$Bio_smry
#     df["SpawnBiomass",y] <-  basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
#     df["Spawnbio95CI",y] <-  paste0( round( as.numeric(df["SpawnBiomass",y])-1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]),"--",
#                                      round(as.numeric(df["SpawnBiomass",y])+1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]))
# 
#     df["Rec",y] <- round(basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],2)
#     df["Rec95CI",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
#                                round(as.numeric(df["Rec",y]) + 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
#     # dq - qnorm(1-(1-quant)/2)*sd chantel this is 1.96
# 
#     df["Depletion",y] <- paste0(round(basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],3)*100,"%")
#     df["Depl95CI",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
#                                 round(as.numeric(df["Rec",y] )+ 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
# 
# 
#   }
#   colnames(df)<- YOI
#   write.csv(df,paste0(rootdir.temp,"./",Sys.Date(),r,"table_i.csv"), row.names = TRUE)
# 
# 
#   ## Figure E TS of Harvest rate: points with SE, red line with SPR50%
#   SSplotSummaryF(basemod10,  print = T, Ftgt = basemod10$derived_quants[grep("Fstd_SPRtgt", basemod10$derived_quants$Label),"Value"], plotdir = rootdir.temp)
#   SSplotCatch(basemod10, subplot = 2,  print = T, plotdir = rootdir.temp)
# 
#   SSplotTimeseries(basemod10, subplot = 11, print = T, plotdir = rootdir.temp, forecastplot = FALSE) ## recdevs
#   SSplotTimeseries(basemod10, subplot = 9,  print = T, plotdir = rootdir.temp, forecastplot = FALSE)
#   SSplotTimeseries(basemod10, subplot = 7,  print = T, plotdir = rootdir.temp, forecastplot = FALSE)
#   ## by hand rel depl
#   depl <- data.frame("depl" = basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Value"],
#                      "upper" = basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Value"]+1.96*basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"StdDev"],
#                      "lwr" = basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Value"]-1.96* basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"StdDev"]) %>%
#     mutate('Yr' = as.numeric(as.character(substr(basemod10$derived_quants[grep("Bratio_", basemod10$derived_quants$Label),"Label"],8,11)))) %>% filter(Yr < 2020)
#   ggplot(depl, aes(x= Yr, y = depl)) +
#     theme_classic()+
#     geom_line(lwd =1.1) +
#     geom_hline(yintercept = 0.4, col = 'red', lwd = 1.1)+
#     geom_hline(yintercept = 0.25, col = 'red',lwd = 1.1)+
#     geom_text(aes(x = 1915, y = 0.42, label = "Management Target")) +
#     geom_text(aes(x = 1915, y = 0.26, label = "Minimum Stock Size Threshold")) +
# 
#     geom_ribbon(aes(ymin = lwr, ymax = upper), alpha = 0.2) +
#     labs(x = 'Year', y = 'Relative Depletion', title = r)
#   ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_Depl.png"), height = 4, width = 6, units = 'in', dpi = 420)
# 
#   SSplotSPR(basemod10,subplot = 4, print = TRUE, plotdir = rootdir.temp, forecastplot = FALSE) ## 1-SPRraw time series
#   SSplotSPR(basemod10,subplot = 2, print = T, plotdir = rootdir.temp, forecastplot = FALSE, uncertainty = TRUE) ## 1-SPR kobe type
# 
# 
#   # SSB_sd.Yr <- basemod10$derived_quants[grep("SSB_", basemod10$derived_quants$Label),"StdDev"][1:134]
#   # basemod10$timeseries %>%
#   #   select(Yr, SpawnBio) %>%
#   #   mutate(upper = SpawnBio + 1.96*SSB_sd.Yr, lwr = SpawnBio - 1.96*SSB_sd.Yr) %>%
#   #   filter(Yr < 2020) %>%
#   # ggplot(., aes(x= Yr, y = SpawnBio)) +
#   #   theme_classic()+
#   #   geom_line(lwd =1.1) +
#   #   # scale_x_continuous(limits = c(1900, 2019), breaks = seq(1900,2020,20)) +
#   #   geom_vline(xintercept = 2019) +
#   #   geom_ribbon(aes(ymin = lwr, ymax = upper), alpha = 0.2) +
#   #   labs(x = 'Year', y = 'Spawning Biomass (000 mt)', title = r)
#   # ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_SSBtimeseries.png"), height = 4, width = 6, units = 'in', dpi = 420)
# 
# 
# 
# 
# }
# ## now merge them all - can also be used for table b
# rootdir <- "C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019"
# regdirs <- list.dirs(rootdir,recursive = FALSE)[grep('cr',list.dirs(rootdir,recursive = FALSE))]
# # regnames <- basename(regdirs)
# regnames <- basename(dirname(list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep('table_i.csv',list.files(regdirs,recursive = FALSE,full.names = TRUE))]))
# df <- list.files(regdirs, recursive = FALSE, full.names = TRUE)[grep(paste0('table_i.csv'),list.files(regdirs,recursive = FALSE,full.names = TRUE))] %>%
#   lapply(.,read.csv) %>%
#   bind_rows()
# 
# for(i in 1:length(basename(regnames))){
#   idx <- (1:13)+c(0,13,26)[i]
#   df$REG[idx] <- sub("cr","",basename(regnames)[i])
# }
# write.csv(df,paste0(rootdir,"/",Sys.Date(),"summary_table_i.csv"), row.names = TRUE)
# 
# 
# ## Table E summary of reference points
# rootdir <- "C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/tgfiles"
# 
# for(r in c('North','Central','South')){ ## loop regions
#   rootdir.temp <- paste0(rootdir,"/cr",r )
#   lastrun <- paste0(rootdir.temp,"/base2015") ## if running base there needs to be a space between SPR_TARGET line 370
#   # basemod1  <- SS_output(lastrun, covar = TRUE)
#   SS_executivesummaryMK(dir = paste0(lastrun))
#   # SS_executivesummaryMK(dir = lastrun)
#   # SSB_sd.Yr <- basemod1$derived_quants[grep("SSB_", basemod1$derived_quants$Label),"StdDev"][1:nrow(basemod1$timeseries)]
#   # basemod1$timeseries %>%
#   #   select(Yr, SpawnBio) %>%
#   #   mutate(upper = SpawnBio + 1.96*SSB_sd.Yr, lwr = SpawnBio - 1.96*SSB_sd.Yr) %>%
#   #   ggplot(., aes(x= Yr, y = SpawnBio)) +
#   #   theme_classic()+
#   #   geom_line(lwd =1.1) +
#   #   # scale_x_continuous(limits = c(1900, 2019), breaks = seq(1900,2020,20)) +
#   #   # scale_y_continuous(limits = c(0,), breaks = seq(0,90,15)) +
#   #   geom_ribbon(aes(ymin = lwr, ymax = upper), alpha = 0.2) +
#   #   labs(x = 'Year', y = 'Spawning Biomass (000 mt)', title = r)
#   # ggsave(last_plot(), file = paste0(rootdir.temp,"/",r,"_SSBtimeseries.png"), height = 4, width = 6, units = 'in', dpi = 420)
# }
# #   df<- data.frame()
# #   for(i in c("Value","95CI")){
# #     if(i == '95CI'){
# #
# #
# #
# #     }
# #
# #   df["UnfishedBio",i] <- NA
# #   df["Unfished_A10",y] <- subset(basemod10$timeseries[, c('Yr', 'Bio_smry')], Yr == YOI[y])$Bio_smry
# #   df["Unfished_Recr",y] <-  ifelse(length(basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]) == 0, NA, basemod10$derived_quants[grep(paste0("OFLCatch_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"])
# #   df["SpBio_2015",y] <-  NA
# #
# #   df["ProxySB40",y] <-  1- basemod10$derived_quants[grep(paste0("SPRratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
# #   df["SPR_SB40",y] <-  basemod10$derived_quants[grep(paste0("F_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
# #
# #   df["F_SB40",y] <-
# #   df["YieldwSPR_SB40",y]
# #
# #   df["SpawnBio_2019",y]  <-  basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"]
# #   df["SpawnBioCI",y] <-  paste0( round( as.numeric(df["SpawnBiomass",y])-1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]),"--",
# #                                    round(as.numeric(df["SpawnBiomass",y])+1.96*basemod10$derived_quants[grep(paste0("SSB_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"]))
# #
# #   df["SB_MSY",y] <- round(basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],2)
# #   df["SPR_MSY",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
# #                              round(as.numeric(df["Rec",y]) + 1.96*basemod10$derived_quants[grep(paste0("Recr_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
# #   # dq - qnorm(1-(1-quant)/2)*sd chantel this is 1.96
# #
# #   df["F_SPRY_MSY",y] <- paste0(round(basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"Value"],3)*100,"%")
# #   df["MSY",y] <- paste0( round(as.numeric(df["Rec",y])  - 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2),"--",
# #                               round(as.numeric(df["Rec",y] )+ 1.96*basemod10$derived_quants[grep(paste0("Bratio_",YOI[y],collapse = "|"), basemod10$derived_quants$Label),"StdDev"],2))
# #
# #
# # }
# sumlist<-list.dirs("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/mods", recursive = F) %>%
# lapply(.,SS_output) %>%
#   SSsummarize(.)
# 
# sumlist$SPRratio %>% spread(Yr,Label)
# df0<- with(sumlist, cbind(melt(sumlist$SPRratio,id='Yr',value.name = 'SPRratio'),
#                     melt(sumlist$SPRratioLower,id='Yr',value.name = 'SPRratioLower'),
#                     melt(sumlist$SPRratioUpper,id='Yr',value.name = 'SPRratioUpper')))
# df0$SPRratio <- as.numeric(as.character(df0$SPRratio))
# df0$SPRratioLower <- as.numeric(as.character(df0$SPRratioLower))
# df0$SPRratioUpper <- as.numeric(as.character(df0$SPRratioUpper))
# names(df0)[c(2,4,5,7,8)] <- c('model','yx','vx','yy','vy')
# df0 <- df0 %>%  select(-vx,-vy, -yx,-yy) %>%
#   filter(model != 'Label')
# levels(df0$model) <- c('Central','North','South','Label')
# 
# 
#   # data.frame("SPR" = 1- basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Value"],
#   #            "sd" =  basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"StdDev"],
#   #            "Yr" = as.numeric(as.character(substr(basemod10$derived_quants[grep("SPRratio_", basemod10$derived_quants$Label),"Label"],10,14)))) %>%
# png("C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/sprratio_all.png",
#     width = 8, height = 6, units = 'in', res = 520)
# df0 %>%
# filter(Yr < 2020) %>%
#   # mutate(lwr = SPR-1.96*sd, upper = SPR+1.96*sd ) %>%
#   ggplot(., aes(x= Yr, y = SPRratio, color = model)) +
#   theme_classic()+
#   theme(legend.position = c(0.2,0.9),
#         legend.text = element_text(size = 22),
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20))+
#   geom_point() +
#   geom_hline(yintercept =1, col = 'red')+
#   scale_color_brewer(palette = 'Dark2')+
#   scale_x_continuous(limits = c(1900, 2021), breaks = seq(1900,2020,20)) +
#   scale_y_continuous(limits = c(0,2)) +
#   geom_errorbar(aes(ymin = SPRratioLower, ymax = SPRratioUpper, col = model ), width = 0) +
#   labs(x = 'Year', y = paste(basemod10$SPRratioLabel), color = '')
# dev.off()# ggsave(plot = last_plot(), file = "C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/sprratio_all.png")
# #   SSplotComparisons(., print = TRUE, shadeForecast = TRUE, legendlabels = c('Central','North','South'),
# #                     plotdir = "C:/Users/",compname,"/Dropbox/UW/assessments/china_2019_update/chinarock-update-2019/sprfiles/plots")

