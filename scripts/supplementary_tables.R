
## Supplementary Tables ##

## Packages ----

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(mgcv)
library(MuMIn)

##  Table S2 ----

# read embryonic thermal physiology data
tphys_e <- read.csv("data/embryonic_thermal_phys.csv")

# summarize embryo thermal phys data
tphys_e_sum <- tphys_e |> 
  pivot_longer(cols = c("topt_e", "ctmin_e", "ctmax_e"), 
               names_to = "trait", values_to = "value") |>
  group_by(trait, parity) |> # remove parity groupping for generall averages 
  summarise(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T)) |> 
  mutate(stage = "embryo")

# read adult data 
adult_data <- read_csv("data/adult_data.csv")

# summarize adult data
tphys_a_sum <- adult_data |> 
  pivot_longer(cols = c("topt_a", "ctmin_a", "ctmax_a"), 
               names_to = "trait", values_to = "value") |>
  group_by(parity, trait) |> # remove parity grouping for general averages
  summarise(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T)) 

# species for which both adult and embryo thermal data was available
length(intersect(tphys_e$species, unique(adult_data$species)))

##  Table S3 ----

# read model_data
load("data/model_test_data.RData")

# obtain average estimates of temperatures experienced (Table S3)
model_test_data |> 
  filter(dev_check == 1) |>
  filter(if_else(!is.na(depth), depth == 5, TRUE)) |> 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) |>
  group_by(parity, species, lat, lon, elev) |> 
  summarise(tg = mean(tg), tap = mean(tap), tep = mean(tep), tenv = mean(tenv)) |> 
  ungroup()|>
  group_by(parity) |> 
  summarise(mtg = mean(tg), mtap = mean(tap), 
            mtep = mean(tep), mtenv = mean(tenv),
            sdtg = sd(tg), sdtap = sd(tap),
            sdtep = sd(tep), sdtenv = sd(tenv))

##  Table S4 ----

# read model_data
load("data/model_test_data.RData")

# prepare data for model
eco_data <- model_test_data |> 
  filter(dev_check == 1) |>
  filter(alpha == 0.5) |> # filter only species with alpha = 0.5
  filter(gamma == 2) |> 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) |> 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) |>
  group_by(parity, species, lat, lon, elev) |> 
  summarise(Dtg = mean(tg), Ctep = mean(tep), Btenv = mean(tenv)) |> 
  ungroup() |> 
  pivot_longer(cols = c("Dtg", "Ctep", "Btenv"), 
               names_to = "type", values_to = "temp") |> 
  mutate(abs_lat_scale = scale(abs(lat)),
         elev_scale = scale(elev), 
         temp_scale = scale(temp))

# run model
eco_data_model <- lmer(temp_scale ~ type * parity + abs_lat_scale + elev_scale + (1|species), 
                       data = eco_data)

# obtain model summaries
summary(eco_data_model) # get model summary
confint(eco_data_model) # get confidence intervals
MuMIn::r.squaredGLMM(eco_data_model)


##  Table S5 ----

# read model_data
load("data/model_test_data.RData")

# obtain average estimates of temperatures experienced by embryo post gestation (Table S6)
model_test_data |> 
  filter(dev_check == 1) |>
  filter(!is.na(depth), !is.na(shade)) |> 
  group_by(parity, species, lat, lon, elev, depth, shade) |> 
  summarise(tep = mean(tep)) |> 
  ungroup() |>
  group_by(depth,parity, shade) |> 
  summarise(mtep = mean(tep), sdtep = sd(tep)) |> 
  mutate(est = paste(round(mtep, 1), " ± ", round(sdtep, 2), sep = "")) |> 
  select(depth, parity, shade, est) |> 
  pivot_wider(names_from = "shade", values_from = "est") |> 
  arrange(depth)


##  Table S6 ----

# read model_data
load("data/model_test_data.RData")

# prepare data for stats
model_test_data_stats <- model_test_data |> 
  filter(dev_check == 1) |>
  filter(alpha == 0.5) |> 
  filter(gamma == 2) |> 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) |> 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) |>
  group_by(species, lat, lon, elev, parity, n, ze, za) |> 
  summarise(opt_d = mean(opt_d, na.rm = T)) |> 
  ungroup() |> 
  mutate(msp = (ze * n)/za) |> 
  mutate(z_lat = as.vector(scale(abs(lat))),
         z_elev = as.vector(scale(elev)),
         z_msp = as.vector(scale(msp)))

# run model
d_model <- glmer(opt_d ~ z_lat * z_elev * z_msp + (1|species), 
                 family = binomial(link = "logit"),
                 data = model_test_data_stats)
summary(d_model)
MuMIn::r.squaredGLMM(d_model)

##  Table S7 ----

# load model data and embryonic development data to calculate Topt_e
load("data/model_test_data.RData")
topt_e <- read.csv("data/embryonic_thermal_phys.csv") |> 
  summarise(mean(topt_e, na.rm = T)) |> as.numeric()

# prepare data for model
model_dat <- model_test_data |> 
  filter(alpha == 0.5) |> # filter only species with alpha = 0.5
  filter(gamma == 2) |> 
  filter(dev_check == 1) |> # filter only months when embryonic development occurs
  filter(if_else(!is.na(depth), depth == 5, TRUE)) |> # if no depth data use 5 cm 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) |> # if no shade data use 0.5
  mutate(abs_diff = abs(topt_e - tep)) |> # calculate absolute Topt_e to tep difference
  group_by(parity, species, lat, lon, elev) |>
  summarise(d = mean(opt_d, na.rm = T), # get average d*
            diff = mean(abs_diff, na.rm = T)) |> # get average difference
  ungroup() |>
  mutate(zabslat = as.vector(scale(abs(lat))),
         zelev = as.vector(scale(elev)),
         zdiff = as.vector(scale(diff)))

# run model
model <- glmer(d ~ diff + (1|species), family = binomial(link = "logit"),
               data = model_dat)

# obtain model summaries
summary(model)
MuMIn::r.squaredGLMM(model)

## Table S8 ----

# read model_data
load("data/model_test_data.RData")

tables8 <- model_test_data |> 
  filter(dev_check == 1) |>
  filter(alpha == 0.5) |> 
  filter(gamma == 2) |> 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) |> 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) |> 
  group_by(species, lat, lon, elev, parity, n, ze, za, nest_shade, 
           nest_depth) |> 
  summarise(opt_d = mean(opt_d, na.rm = T),
            fitness = mean(fitness, na.rm = T)) |> 
  ungroup() |> 
  mutate(msp = (ze * n)/za) |> 
  mutate(parity = ifelse(parity == "Oviparous", "O", "V")) |>
  mutate(species = gsub("_", " ", species)) |>
  select(species, parity, lat, elev, msp, nest_shade, nest_depth, opt_d, fitness)
  
  
