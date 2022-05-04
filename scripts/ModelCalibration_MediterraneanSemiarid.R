# Mediterranean Semiarid validation and calibration 


# File locations
soil_loc <- "C:/Users/vazqu003/OneDrive - Wageningen University & Research/Climate farmers/GithubForks/soil-modelling"
project_loc <- "G:/Shared drives/Climate Farmers/02_Product/01_02_Transition_Finance/Carbon_Credits/FarmerData"
# Check the loc for weather DB. 
project_name<- "/Francisco Alves"
# project_name<- "/validation_temperatedryarea"
modelling_data_loc <- "C:/Users/vazqu003/OneDrive - Wageningen University & Research/Climate farmers/GithubForks/modelling-data"
weatherDB_loc <- "G:/Shared drives/Climate Farmers/07_Tech/Modelling/WeatherDB/sis-biodiversity-cmip5-regional_data"




lat<-38.9
long<-8.22
clay=15
land_use<-"Cropland" # or Grassland

LU<-read.csv(file.path(modelling_data_loc,"data", "/soil_landuse.portugal.csv"))
unique(LU$landcover)
LU<-LU[LU$latitude < 39.5,] # Roughly where dry temperate weather starts in portugal 
tree.vs     <- LU[LU$landcover=="Tree cover", ]
validation.set<-LU[LU$landcover=="Cropland" | LU$landcover=="Grassland", ]
run_soil_model(soil_loc,project_loc,project_name,modelling_data_loc,weatherDB_loc, dt=validation.set[1,], tree.vs)

val.red<-list()
for (i in 1:dim(validation.set)[1]){
  val.red[[i]]<-run_soil_model(soil_loc,project_loc,project_name,modelling_data_loc,weatherDB_loc, dt=validation.set[i,], tree.vs)
} 

val.dataset<-do.call(rbind.data.frame, val.red)
val.dataset<- val.dataset[!is.na(val.dataset$ObservedC),]
head(val.dataset)
sum(is.na(val.dataset$ObservedC))





efficiency<-function(p, o){1-(sum((p-o)^2)/sum((o-(sum(o)/length(o)))^2))}
rmse<-function(p, o){sqrt(sum(((p-o)^2)/length(o)))}
  
Obs.C_Grassland  <- val.dataset$ObservedC[val.dataset$land_use=="Grassland"]
Pred.C_Grassland <- val.dataset$PredictedC[val.dataset$land_use=="Grassland"]
RMSE.grassland <- rmse(Pred.C_Grassland, Obs.C_Grassland)
EF.grassland <- efficiency(Pred.C_Grassland, Obs.C_Grassland)

Obs.C_Cropland  <- val.dataset$ObservedC[val.dataset$land_use=="Cropland"]
Pred.C_Cropland <- val.dataset$PredictedC[val.dataset$land_use=="Cropland"]
RMSE.cropland <- rmse(Pred.C_Cropland, Obs.C_Cropland)
EF.cropland<-efficiency(Pred.C_Cropland, Obs.C_Cropland)


par(mfrow=c(1,2))
plot(Obs.C_Grassland ~ Pred.C_Grassland, xlim=c(10, 50), ylim=c(10, 50))
abline(0,1)

plot(Obs.C_Cropland ~ Pred.C_Cropland, xlim=c(10, 50), ylim=c(10, 50))
abline(0,1)




