# Code to run bias estimation 

# Description of the script ----
## Dependencies
## model_functions.R
## modified_functions.R
## scripts/calc_functions_soil_modelling.R
## scripts/weather_data_pulling_functions.R
## Weather database should be accessible. 
## json file should have all paths, including correct path for database. 
## Data with farmer information - must have the following columns: 
## latitude, longitude 
## C inputs are taken from modelling the baseline. 
## G:\Shared drives\Climate Farmers\02_Product\01_02_Transition_Finance\Carbon_Credits\FarmerData\Francisco Alves\result\parcel_Cinput.csv


# Future amends -----
########################### AUTOMATED SOIL MODEL RUNNING SCRIPT
# Missing automated pull of DPM/RPM ratio from input (aboveground crop/pasture/root exudates) type
# Missing automated pull of bare soil
# Missing automated pull of clay content: it is not set where we will find it
# Missing automated pull of standard errors for input parameters
# Missing observed SOC values for bias analysis
# Missing change RothC parameters to semi-arid according to Farina et al. 2013
# dt should be a dataset that includes columns for: latitude, longitude, clay content
# land_use
# tree.vs has same info as dt, but only tree covered land uses

library(rjson)
library(SoilR)

# Function ----
run_soil_model <- function(soil_loc,project_loc,project_name,modelling_data_loc,weatherDB_loc, dt, tree.vs){
  
  source(file.path(soil_loc, "model_functions.R"))
  source(file.path(soil_loc, "modified_functions.R"))
  source(file.path(soil_loc, "scripts/calc_functions_soil_modelling.R"))
  source(file.path(modelling_data_loc, "scripts/weather_data_pulling_functions.R"))
  

  lat_farmer <- dt$latitude 
  lon_farmer <- dt$longitude
  clay_c     <- dt$clay_value_avg # TO BE PICKED FROM YOUR DATA
  land_use   <- dt$landcover
  # Weather inputs ---- 
  
  weather_data = data.frame(past_temperature=rep(NA,12))
  weather_data[,c("past_temperature", "future_temperature_rcp4.5")] <- get_monthly_mean_temperature(lon_farmer,lat_farmer, scenario="rcp4.5")
  weather_data[,c("past_precipitation", "future_precipitation_rcp4.5")] <- get_monthly_mean_precipitation(lon_farmer,lat_farmer,scenario="rcp4.5")
  weather_data[,c("past_pevap", "future_pevap_rcp4.5")] <- get_monthly_mean_pevap(lon_farmer,lat_farmer,scenario="rcp4.5")
  
  # Carbon inputs - from Fcos Farm -----
  # ################# Pulling calculation factors
  animal_factors <- read_csv(file.path(modelling_data_loc,"data", "carbon_share_manure.csv")) %>% filter(type=="manure") %>%
    rename(species=manure_source)
  agroforestry_factors <- read_csv(file.path(modelling_data_loc,"data", "agroforestry_factors.csv"))
  crop_data <- read_csv(file.path(modelling_data_loc,"data", "crop_factors.csv"))#, col_types =  "cdddddddd")
  pasture_data <- read_csv(file.path(modelling_data_loc,"data", "pasture_factors.csv"))
   
  # ################# Pulling calculation inputs
  animal_inputs <- read_csv(file.path(project_loc,project_name,"inputs", "animal_inputs.csv"))
  agroforestry_inputs <- read_csv(file.path(project_loc,project_name,"inputs", "agroforestry_inputs.csv"))
  crop_inputs <- read_csv(file.path(project_loc,project_name,"inputs", "crop_inputs.csv"))
  pasture_inputs <- read_csv(file.path(project_loc,project_name,"inputs", "pasture_inputs.csv"))
  parcel_inputs <- read_csv(file.path(project_loc,project_name,"inputs", "parcel_inputs.csv"))
  ################# Calculations per parcel and scenario
  # 
  parcel_Cinputs =data.frame(parcel_ID=c(),scenario=c(),agroforestry_Cinput=c(),animal_Cinput=c(),crop_Cinputs=c(),pasture_Cinputs=c())
  for(parcel in parcel_inputs$parcel_ID){
    for(scenario in c("current","future","baseline","previous_years")){
      parcel_Cinputs<-rbind(parcel_Cinputs,data.frame(parcel_ID=parcel,
                                                      scenario=scenario,
                                                      agroforestry_Cinput=get_monthly_Cinputs_agroforestry(agroforestry_inputs, agroforestry_factors,
                                                                                                           scenario, parcel, lat_farmer),
                                                      animal_Cinput=get_monthly_Cinputs_animals(animal_inputs, animal_factors, scenario, parcel),
                                                      crop_Cinputs=get_monthly_Cinputs_crop(crop_inputs, crop_data, scenario, parcel),
                                                      pasture_Cinputs=get_monthly_Cinputs_pasture(pasture_inputs, pasture_data, scenario, parcel)))
    }
  }
  parcel_Cinputs <- parcel_Cinputs %>% mutate(tot_Cinputs=agroforestry_Cinput+animal_Cinput+crop_Cinputs+pasture_Cinputs)


  # All inputs for initialisation ----  
  ################# Initialisation by making the model reach SOC of natural areas of the pedo-climatic area
  mean_input<-data.frame(run=rep(0,12), 
                         dr_ratios=c(0.67,rep(NA,11)), 
                         bare=as.factor(c(logical(12))), 
                         past_temp=weather_data$past_temperature, 
                         past_precip=weather_data$past_precipitation, 
                         past_evap=weather_data$past_pevap, 
                         soil_thick=c(25,rep(NA,11)), 
                         clay=c(clay_c,rep(NA,11)), 
                         pE=c(0.75,rep(NA,11)), 
                         tilling_factor=c(1.0,rep(NA,11)))

  # Phase 1 - model initialization ----
  # I could average all "tree cover" areas. I don't need to validate these. 
  # The baseline scenario for calibration should be ind. from the baseline sc. in validation
  # Phase 1 - a) Understanding input data ----
  tree.vs<-tree.vs[tree.vs$lower_depth>10, ]
  cstock25<-tree.vs$cstock_t.ha*25/tree.vs$lower_depth
  tree.vs<-tree.vs[ order(abs(tree.vs$clay_value_avg-clay_c))[1:10], ]
  tree.clay=median(tree.vs$clay_value_avg)
  SOC_nveg= median(cstock25)
  clay.sd = sd(tree.vs$clay_value_avg)
  cstock.sd= sd(cstock25)
  dr_ratio_forest = 0.25
  dr_ratios_savanna = 0.67
  time_horizon = 1000
  
  FallIOM=0.049*SOC_nveg^(1.139) #IOM using Falloon method
  # TO CHECK: CODE USED CLAY; i CHAGED IT TO TREE.CLAY 
  # We use pedotransfer functions for compartment equilibrium 
  # From Weihermuller et al 2013
  RPM = (0.1847*SOC_nveg+0.1555)*((tree.clay+1.2750)^-0.1158)
  HUM = (0.7148*SOC_nveg+0.5069)*((tree.clay+0.3421)^0.0184)
  BIO = (0.0140*SOC_nveg+0.0075)*((tree.clay+8.8473)^0.0567)
  # And let the model get to equilibrium 


  # Phase 1 - solve for C inputs ---- 
  # we model 10 C inputs. and plot them
  
  c_cinput_balance<-data.frame(matrix(nrow=2, ncol=2))
  names(c_cinput_balance)<-c("carboninput", "carbonstock")
  c_cinput_balance$carboninput=c(10, 20)   #Annual C inputs to soil in Mg/ha/yr
  pedotransfer_ini_soil_content=c(DPM=0, RPM=RPM, BIO=BIO, HUM=HUM, IOM=FallIOM)
  model1 <- calc_carbon_over_time(time_horizon,
                                  field_carbon_in = rep(c_cinput_balance[1,1],time_horizon),
                                  dr_ratios = rep(dr_ratio_forest,time_horizon),
                                  bare = mean_input$bare,
                                  temp = mean_input$past_temp,
                                  precip = mean_input$past_precip,
                                  evap = mean_input$past_evap,
                                  soil_thick = mean_input$soil_thick[1],
                                  clay = mean_input$clay[1],
                                  pE = mean_input$pE[1],
                                  PS = pedotransfer_ini_soil_content,
                                  tilling_factor = mean_input$tilling_factor[1])
  c_cinput_balance[1,"carbonstock"]<-tail(model1$TOT, 1)
  
  model2 <- calc_carbon_over_time(time_horizon,
                                  field_carbon_in = rep(c_cinput_balance[2,1],time_horizon),
                                  dr_ratios = rep(dr_ratio_forest,time_horizon),
                                  bare = mean_input$bare,
                                  temp = mean_input$past_temp,
                                  precip = mean_input$past_precip,
                                  evap = mean_input$past_evap,
                                  soil_thick = mean_input$soil_thick[1],
                                  clay = mean_input$clay[1],
                                  pE = mean_input$pE[1],
                                  PS = pedotransfer_ini_soil_content,
                                  tilling_factor = mean_input$tilling_factor[1])
  c_cinput_balance[2,"carbonstock"]<-tail(model2$TOT, 1)

  slope=(c_cinput_balance[2,2]-c_cinput_balance[1,2])/(c_cinput_balance[2,1]-c_cinput_balance[1,1])
  Cinput_leading_to_observed_SOC_past_land_use = SOC_nveg/slope

  # Natural area initialisation ---- 
  # Initialisation - We need to do this for grassland and arable fields. 
  # Attention: Why is the field_carbon_in is the Cinput_baseline, and not Cinput_leading_to_observed_SOC_past_land_use
  # I changed this. In the run_soil_model.R was ok. 
  mean_input$field_carbon_in <- rep(Cinput_leading_to_observed_SOC_past_land_use,12)
  # Why not start with 0? 
  starting_soil_content <- estimate_starting_soil_content(SOC=Cinput_leading_to_observed_SOC_past_land_use,clay=mean_input$clay[1])
  time_horizon = 1000
  # Attention: Cinput_leading_to_observed_SOC_past_land_use does this reflect the yearly input or monthly? 
  # Attention: Why does the model start at calculated C and not 0? 

  C0_df <- calc_carbon_over_time(time_horizon,
                                 field_carbon_in = rep(mean_input$field_carbon_in[1],time_horizon), # should this value not be Cinput_leading_to_observed_SOC_past_land_use
                                 dr_ratios = rep(dr_ratio_forest,time_horizon),
                                 bare = mean_input$bare,
                                 temp = mean_input$past_temp,
                                 precip = mean_input$past_precip,
                                 evap = mean_input$past_evap,
                                 soil_thick = mean_input$soil_thick[1],
                                 clay = mean_input$clay[1],
                                 pE = mean_input$pE[1],
                                 PS = starting_soil_content,
                                 tilling_factor = mean_input$tilling_factor[1])
  nveg_soil_content <- as.numeric(tail(C0_df,1))[c(1:5)]

  cinput_montado_bl<-parcel_Cinputs$tot_Cinputs[parcel_Cinputs$parcel_ID=="Montado" & parcel_Cinputs$scenario=="baseline"]
  
  mean_input$field_carbon_in <- 0.5*(Cinput_leading_to_observed_SOC_past_land_use+cinput_montado_bl)
  time_horizon = 350
  
  C0_df <- calc_carbon_over_time(time_horizon,
                                 field_carbon_in = rep(mean_input$field_carbon_in[1],time_horizon),
                                 dr_ratios = rep(dr_ratios_savanna,time_horizon),
                                 bare = mean_input$bare,
                                 temp = mean_input$past_temp,
                                 precip = mean_input$past_precip,
                                 evap = mean_input$past_evap,
                                 soil_thick = mean_input$soil_thick[1],
                                 clay = mean_input$clay[1],
                                 pE = mean_input$pE[1],
                                 PS = nveg_soil_content,
                                 tilling_factor = mean_input$tilling_factor[1])
  initialized_soil_content <- as.numeric(tail(C0_df,1))[c(1:5)]

  # Here is where we divide results for arable farms and grasslands    
  cinput_cropland_bl<-parcel_Cinputs$tot_Cinputs[parcel_Cinputs$parcel_ID=="Arable-crops" & parcel_Cinputs$scenario=="baseline"]
  cinput_grassland_bl<-parcel_Cinputs$tot_Cinputs[parcel_Cinputs$parcel_ID=="Montado-bale-grazing" & parcel_Cinputs$scenario=="baseline"]
  if (land_use=="Cropland"){
  time_horizon = 70
  C0_df <- calc_carbon_over_time(time_horizon,
                                 field_carbon_in = rep(cinput_cropland_bl,time_horizon),
                                 dr_ratios = rep(mean_input$dr_ratios[1],time_horizon),
                                 bare = c(rep(FALSE, 7), rep(TRUE, 4), FALSE),  # Attention: change, arable land is bare at times
                                 temp = mean_input$past_temp,
                                 precip = mean_input$past_precip,
                                 evap = mean_input$past_evap,
                                 soil_thick = mean_input$soil_thick[1],
                                 clay = mean_input$clay[1],
                                 pE = mean_input$pE[1],
                                 PS = initialized_soil_content,
                                 tilling_factor = mean_input$tilling_factor[1]) # Attention: change: tilling occurs at least 1/2 times a year. 
  print(tail(C0_df,1))  
  res <- as.numeric(tail(C0_df,1))[c(1:6)] } else if (land_use=="Grassland") {
  time_horizon = 70
  C0_df <- calc_carbon_over_time(time_horizon,
                                 field_carbon_in = rep(cinput_grassland_bl,time_horizon), # Not sure that this is the right input
                                 dr_ratios = rep(mean_input$dr_ratios[1],time_horizon),
                                 bare = mean_input$bare,
                                 temp = mean_input$past_temp,
                                 precip = mean_input$past_precip,
                                 evap = mean_input$past_evap,
                                 soil_thick = mean_input$soil_thick[1],
                                 clay = mean_input$clay[1],
                                 pE = mean_input$pE[1],
                                 PS = initialized_soil_content,
                                 tilling_factor = mean_input$tilling_factor[1]) 
  
  res <-  as.numeric(tail(C0_df,1))[c(1:6)]} else print("land_use not recognised")
  names(res)<-names(C0_df)
  ObservedC<-dt$cstock_t.ha*25/dt$lower_depth
  data.frame(PredictedC=res[6], ObservedC = ObservedC, land_use=land_use, clay=clay_c)
}

