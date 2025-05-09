####################################
## Include Me in your scripts with
## "source("include.R")
####################################

# morally always needed
library(ggplot2)
library(grDevices)
library(png)

Mode <- function(x,verbose=0) {
	ux <- unique(x)
	if(verbose){
		print(table(x))
	}
	ux[which.max(tabulate(match(x, ux)))]
}

########################
# useful aliases
########################
size = dim
end = len = lenght = length
assert = stopifnot
ln = log

# extrema = function(vec){
# 	return(c(min(vec),max(vec)))
# }
extrema <- function(...) {
	# Combine all input vectors into one
	combined <- c(...)
	
	# Return the minimum and maximum of the combined vector
	return(c(min(combined), max(combined)))
}

count_na = function(vec){
	na_vals = sum(as.numeric(is.na(vec)))
	return(na_vals)
}
na_count = quanti_na = count_na

########################
# loading procedure, with feedback
########################
library(crayon)
library(hash)

########################
# df_stat split
########################
create_df_stat = function(df_data){
	df_stat_ = hash()
	if("IDStations" %in% colnames(df_data)){
		stations_ = unique(df_data$IDStations)
		for (st in stations_){
			df_stat_[[st]] = df_data[which(df_data$IDStations == st),]
		}
		rm(st)
		rm(stations_)
		return(df_stat_)
	}
	else{
		cat(crayon::red("Error: missing IDStations field in df.\n"))
		return(NA)
	}
}

########################
# function to get colors for plotting
########################
library(RColorBrewer)
fun_colori = function(len=2, seed=33, show=1, seed_div = "Set3"){
	hcols_ = hcl.pals()
	if(seed=="rand"){
		seed = round(runif(1,0,115))
		col.ramp_ = hcl.colors(len,palette=hcols_[seed%%115+1])
	}
	if(seed=="div"){ # create a divergent palette
		col.ramp_ = brewer.pal(len, seed_div)
		# possible seed_div choices (and max len supported)
		# Set3	    12
		# Paired    12
		# Pastel1   9
		# Set1	    9
		# Accent    8
		# Dark2     8
		# Pastel2   8
		# Set2	    8
	}
	else{
		col.ramp_ = hcl.colors(len,palette=hcols_[seed%%115+1])
	}
	if(show==1){
		dati_ <- matrix(1:100, ncol = 1)
		par(mar=rep(2,4))
		image(dati_, col = col.ramp_, axes = FALSE)
		title(main=paste("palette",seed,"of",len,"colors"))
		# title(main=paste("palette",seed))
	}
	return(col.ramp_)
	
}
colori_fun = colorami = colora = fun_colori # aliases
# usage: cols = colora(7) for getting a palette with 7 colors
# usage: cols = colora(7,456) for getting a palette of seed 456 with 7 colors

# cat(crayon::cyan(h,"Created function to get color palettes. Available as"),crayon::red("colora(len, seed, show).\n"))
# cat(crayon::italic("Try for example colora(10,56,1).\n"))


########################
# function to show the hash content
########################
# what_is <- new.env(hash = TRUE, parent = emptyenv(), size = NA)
what_is = hash()
# usage: explain("AQ_pm10") show AQ_pm10 explanation
# usage: explain("pm10") yes works also with partial names
# usage: explain("pm") here explains both pm 10 and 2.5

spiegami = function(pattern) {
	matching_vars <- names(what_is)[grep(pattern, names(what_is))]
	# chat gpt help here :/	
	if (length(matching_vars) == 0) {
		cat("no matching variable name found\n")
	} else {
		for (var in matching_vars) {
			cat(crayon::red(var), ":", what_is[[var]], "\n")
		}
	}
}
spiega = tellme = tell_me = tell = explain = spiegami # aliases
# spiegami = function(colonna){
# 	cat(what_is[[colonna]])
# }


what_is[["AQ"]] = "Air quality - Air pollutants concentrations"
what_is[["WE"]] = "Wheather - Estimates of atmospheric variables"
what_is[["EM"]] = "Emission - Emission inventories"
what_is[["LI"]] = "Livestock - Livestock inventories"
what_is[["LA"]] = "Land cover -  Estimates of land cover variables"

what_is[["IDStations"]] = "Unique ID (uid) of station present in sensor network [For detail see Station_Registry.csv file]. \n[ - ]"
what_is[["Latitude"]] = "Latitude of station \n[ Degrees ]"
what_is[["Longitude"]] = "Longitude of station \n[ Degrees ]"
what_is[["Time"]] = "Date \n[ yyyy-mm-dd ]"
what_is[["Altitude"]] = "Elevations of the stations \n[ m ]"
what_is[["AQ_pm10"]] = "Concentration of particles with an aerodynamic diameter of less than 10 micrometers (µm) \n[ μg m^(-3) ]"
what_is[["AQ_pm25"]] = "Concentration of particles with an aerodynamic diameter of less than 2.5 micrometers (µm) \n[ μg m^(-3) ]"
what_is[["AQ_co"]] = "Concentration of carbon monoxide \n[ μg m^(-3) ]"
what_is[["AQ_nh3"]] = "Concentration of ammonia \n[ μg m^(-3) ]"
what_is[["AQ_nox"]] = "Concentration of several oxides of nitrogen most of which are produced in combustion and are considered to be atmospheric pollutants \n[ μg m^(-3) ]"
what_is[["AQ_no2"]] = "Concentration of nitrogen dioxide \n[ μg m^(-3) ]"
what_is[["AQ_so2"]] = "Concentration of sulfur dioxide \n[ μg m^(-3) ]"
what_is[["WE_temp_2m"]] = "Temperature of air at 2 m above the surface of land or sea or inland water \n[ Celsius ]"
what_is[["WE_wind_speed_10m_mean"]] = "Mean horizontal wind speed at a height of 10 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_wind_speed_10m_max"]] = "Maximum horizontal wind speed at a height of 10 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_mode_wind_direction_10m"]] = "Direction of horizontal wind at a height of 10 m above the surface of the Earth \n[ categorical ]"
what_is[["WE_tot_precipitation"]] = "The accumulated liquid and frozen water (comprising rain and snow) that falls to the Earth’s surface \n[ m ]"
what_is[["WE_precipication_t"]] = "The type of precipitation at the surface at the specified time. \n[ categorical ]"
what_is[["WE_surface_pressure"]] = "The pressure (force per unit area) of the atmosphere at the surface of land and sea and inland water \n[ Pa ]"
what_is[["WE_solar_radiation"]] = "Amount of solar radiation that reaches a horizontal plane at the surface of the Earth (both direct and diffuse) minus the amount reflected by the Earth’s surface \n[ J m^(-2) ]"
what_is[["WE_wind_speed_100m_mean"]] = "Mean horizontal wind speed at a height of 100 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_wind_speed_100m_max"]] = "Max horizontal wind speed at a height of 100 m above the surface of the Earth \n[ m/s ]"
what_is[["WE_mode_wind_direction_100m"]] = "Direction of horizontal wind at a height of 100 m above the surface of the Earth \n[ categorical ]"
what_is[["WE_blh_layer_max"]] = "The maximum depth of air next to the Earth’s surface which is most affected by the resistance to the transfer of momentum or heat or moisture across the surface \n[ m ]"
what_is[["WE_blh_layer_min"]] = "The minimum depth of air next to the Earth’s surface which is most affected by the resistance to the transfer of momentum or heat or moisture across the surface \n[ m ]"
what_is[["WE_rh_min"]] = "Minimum of amount of water vapour present in air expressed as a percentage of the amount needed for saturation at the same temperature \n[ % ]"
what_is[["WE_rh_mean"]] = "Mean of amount of water vapour present in air expressed as a percentage of the amount needed for saturation at the same temperature \n[ % ]"
what_is[["WE_rh_max"]] = "Maximum of amount of water vapour present in air expressed as a percentage of the amount needed for saturation at the same temperature \n[ % ]"
what_is[["EM_nh3_livestock_mm"]] = "Emissions of NH3 originating from the livestock sector within the manure management \n[ mg m^(-2) ]"
what_is[["EM_nh3_agr_soils"]] = "Emissions of NH3 originating from the agriculture soils \n[ mg m^(-2) ]"
what_is[["EM_nh3_agr_waste_burn"]] = "Emissions of NH3 originating from the burning of agriculture waste \n[ mg m^(-2) ]"
what_is[["EM_nh3_sum"]] = "Total emissions of NH3 across all sectors (anthropogenic) \n[ mg m^(-2) ]"
what_is[["EM_nox_traffic"]] = "Emissions of NOx from the on road transportation \n[ mg m^(-2) ]"
what_is[["EM_nox_sum"]] = "Emissions of NOx across all sectors (anthroprogenic) \n[ mg m^(-2) ]"
what_is[["EM_so2_sum"]] = "Total emissions of SO2 across all sectors (anthroprogenic) \n[ mg m^(-2) ]"
what_is[["LI_pigs"]] = "Density of pigs \n[ Number km^(-2) ]"
what_is[["LI_bovine"]] = "Density of bovines \n[ Number km^(-2) ]"
what_is[["LI_pigs_v2"]] = "Density of pigs \n[ Number km^(-2) ]"
what_is[["LI_bovine_v2"]] = "Density of bovines \n[ Number km^(-2) ]"
what_is[["LA_hvi"]] = "One-half of the total green leaf area per unit horizontal ground surface area for high vegetation type \n[ unitless (m^2/m^2) ]"
what_is[["LA_lvi"]] = "One-half of the total green leaf area per unit horizontal ground surface area for low vegetation type. \n[ unitless (m^2/m^2) ]"
what_is[["LA_land_use"]] = "Land use across 44 sectors \n[ categorical ]"
what_is[["LA_soil_use"]] = "Lombardy soil use across 21 sectors \n[ categorical ]"

# cat(crayon::cyan(h,"Created utility to explain covariates. Available as"),crayon::red("spiega(string).\n"))
# cat(crayon::italic("Try for example spiega(\"wind\").\n\n"))

# rm(h)

##### NAs
na_summary = function(df){
	for (c in colnames(df)){
		cat(which(colnames(df)==c),typeof(df[,c][[1]]),sprintf("%28s, #NAs = %f\n",c,sum(as.numeric(is.na(df[,c])))))
	}
}
summary_na = na_summary
interpola_NA <- function(obs) {
	if (all(is.na(obs))){
		return(obs)
	}
	non_na_indices <- which(!is.na(obs))
	na_indices <- which(is.na(obs))
	# print(non_na_indices)
	# print(na_indices)
	
	# loop on Na indices
	for (na_index in na_indices) {
		# Trova l'indice più vicino a sinistra e quello più vicino a destra
		left_index <- suppressWarnings(max(non_na_indices[non_na_indices < na_index]))
		right_index <- suppressWarnings(min(non_na_indices[non_na_indices > na_index]))
		# cat("NA idx",na_index," - LR values = ",obs[left_index],",",obs[right_index],"\n")
		
		# Esegui l'interpolazione se ci sono valori a sinistra e a destra
		if (!is.na(left_index) && !is.na(right_index)) {
			left_value <- obs[left_index]
			right_value <- obs[right_index]
			
			# Calcola il valore interpolato
			interpolated_value <- left_value + 
				((right_value - left_value) / (right_index - left_index)) * (na_index - left_index)
			
			# Sostituisci il valore NA con quello interpolato
			obs[na_index] <- interpolated_value
		}
		if(is.na(obs[na_index])){
			obs[na_index] = na.omit(c(obs[left_index],obs[right_index]))
		}
		# cat(obs[na_index],"\n")
	}
	return(obs)
}

# example
# test <- c(NA,NA,10, 15, NA, 20, 24, 20, 18, NA, NA, 10, 9, 9, NA)
# print(test)
# result <- interpola_NA(test)
# print(round(result))

