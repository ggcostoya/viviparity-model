
## Testing model predictions and sensitivity analysis ## 

## Packages ----

library(tidyverse)

## Define functions ----

# thermal performance curve function
tpc <- function(temps, ctmin, topt, ctmax){
  
  # define performance object
  perf <- rep(NA, length(temps))
  
  # define psig
  psig <- -((ctmin - topt)/4)
  
  for(i in 1:length(temps)){
    
    # if temperature is below topt
    if(temps[i] < topt){
      
      perf[i] <- exp(-((temps[i]-topt)/(2*psig))^2)
      
      # if temperature is above topt
    }else{
      
      perf[i] <- 1 - ((temps[i] - topt)/(topt - ctmax))^2
      
    }
    
  }
  
  # correct for negative values
  perf <- ifelse(perf < 0, 0, perf)
  
  return(perf)
}

# model function
model <- function(topt_a, ctmin_a, ctmax_a, topt_e, ctmin_e, ctmax_e,
                  za, ze, n, tg, tep, tap, alpha, gamma){
  
  # define sequence of possible d values
  d <- seq(0,1, by = 0.01)
  
  # estimate survival probabilities
  sag <- tpc(tg, ctmin_a, topt_a, ctmax_a)
  seg <- tpc(tg, ctmin_e, topt_e, ctmax_e)
  sap <- tpc(tap, ctmin_a, topt_a, ctmax_a)
  sep <- tpc(tep, ctmin_e, topt_e, ctmax_e)
  
  # egg carrying costs
  C <- 1 - d * alpha * ((n * ze * seg)/za)
  
  # present reproductive success
  Pr <- ((sag * C)^d) * n * (seg^d) * (sep^(1-d))
  
  # future reproductive success
  Fu <- ((sag * C)^d) * (sap^(1-d)) * gamma
  
  # get fitness
  fitness <- Pr + Fu
  
  # determine optimal d
  optimal_d <- d[which.max(fitness)]
  
  # determine fitness at optima
  fitness_optimal_d <- max(fitness)
  
  return(c(optimal_d, fitness_optimal_d))
  
}

## Load and prepare data ----

# load embryonic thermal physiology data
tphys_e <- read.csv("data/embryonic_thermal_phys.csv")

# define embryonic thermal physiology traits
topt_e <- mean(tphys_e$topt_e, na.rm = T)
ctmin_e <- mean(tphys_e$ctmin_e, na.rm = T)
ctmax_e <- mean(tphys_e$ctmax_e, na.rm = T)

# load data on adult thermal physiology and LH traits and additional info
database <- read.csv("data/adult_data.csv")

# add column to list months when embryonic development occurs
database$dev <- rep(NA, nrow(database)) 

# loop to create list of months when embryonic development occurs
for(i in 1:nrow(database)){
  
  # if eggs are seen after hatchlings list includes Dec & Jan
  if(database$eggs_seen[i] > database$hatchlings_seen[i]){
    
    database$dev[i] <- list(c(database$eggs_seen[i]:12,
                                1:database$hatchlings_seen[i]))
    
    # if eggs are seen before hatchlings do simple list
  }else{
    
    database$dev[i] <- list(
      c(database$eggs_seen[i]:database$hatchlings_seen[i]))
    
  }
}

# load ecophysiology data
eco_data <- read.csv("data/eco_data.csv")

# remove "X" columns
eco_data$X <- NULL
database$X <- NULL

# merge data for model predictions
model_data <- merge(eco_data, database, by = c("species", "lat", "lon", "elev"))

# add column to check if a month is within months when embryo dev occurs
model_data$dev_check <- rep(NA, nrow(model_data))

# loop to assign values to dev_check
for(i in 1:nrow(model_data)){
  
  model_data$dev_check[i] <- ifelse(model_data$month[i] %in% model_data$dev[[i]], 1, 0)
    
}

# add columns for sensitivity analysis
species <- unique(model_data$species)
alpha <- seq(0,1, by = 0.25)
gamma <- c(0,1,2,5,10)
combinations <- expand.grid(species, alpha, gamma) 
colnames(combinations) <- c("species", "alpha", "gamma")

# merge model data with combinations
model_data <- merge(model_data, combinations, by = c("species"), all = TRUE)

# add optimal d and fitness columns
model_data$opt_d <- NA
model_data$fitness <- NA

# filter model data to account for species for which there is information
model_data <- model_data |> 
  mutate(shade = ifelse(!is.na(nest_shade), NA, shade)) |> 
  mutate(depth = ifelse(!is.na(nest_depth), NA, depth)) |> 
  as_tibble() |>
  unique()

## Fit model ----

# loop to run model 
for(i in 1:nrow(model_data)){
  
  # determine model outcome
  outcome <- model(topt_a = model_data$topt_a[i],
                   ctmin_a = model_data$ctmin_a[i],
                   ctmax_a = model_data$ctmax_a[i],
                   topt_e = topt_e,
                   ctmin_e = ctmin_e,
                   ctmax_e = ctmax_e,
                   tg = model_data$tg[i],
                   tep = model_data$tep[i],
                   tap = model_data$tap[i],
                   n = model_data$n[i],
                   ze = model_data$ze[i],
                   za = model_data$za[i],
                   alpha = model_data$alpha[i],
                   gamma = model_data$gamma[i])
  
  # keep optimal fitness
  model_data$opt_d[i] <- outcome[1]
  model_data$fitness[i] <- outcome[2]
  
  if(i %in% seq(1000,nrow(model_data),by = 1000)){
    print(paste(round(i/nrow(model_data)*100,2), "% done"))
    }
  
}

## Save data ----

# rename
model_test_data <- model_data

# save as an .RData file to account for nested list
save(model_test_data,file = "data/model_test_data.RData")

