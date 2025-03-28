
## Summaries of data sets and statistics ##

## Packages ----

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(mgcv)

## Embryonic thermal physiology summaries ----

# read embryonic thermal physiology data

# total number of species with any data
nrow(tphys_e) # 40

# number of species with any data for each parity mode
table(tphys_e$parity) # Oviparous = 39, Viviparous = 1

# average values of thermal physiology traits 
mean(tphys_e$topt_e, na.rm = T) # Topt_e = 27.57
mean(tphys_e$ctmin_e, na.rm = T) # CTmin_e = 19.41
mean(tphys_e$ctmax_e, na.rm = T) # CTmax_e = 33.83

# number of species with thermal physiology data at embryo and adult stage
length(intersect(tphys_e$species, unique(model_data$species))) # 6

## Adult thermal physiology summaries ----

# read adult thermal physiology, life history and additional info

# load adult data
adult_data <- read_csv("data/adult_data.csv")

# determine number of populations
length(unique(paste(adult_data$species, adult_data$lat, adult_data$elev))) # 89

# determine number of species
length(unique(adult_data$species)) # 79

# determine number of oviparous and viviparous species and populatuions
adult_data %>% select(parity, species, lat, elev) %>% unique() %>% 
  group_by(parity) %>% summarise(n = n(), nspecies = length(unique(species)))

# number of populations where gamma could be calculated 
adult_data %>% filter(!is.na(gamma_ex)) %>% 
  select(parity, species, lat, elev) %>% unique() %>% 
  summarise(n = n(), nspecies = length(unique(species))) # 35 pop, 27 spp

# number of populations with data on the nesting conditions
adult_data %>% filter(!is.na(nest_depth) |!is.na(nest_shade)) %>% 
  select(parity, species, lat, elev) %>% unique() %>% 
  summarise(n = n(), nspecies = length(unique(species))) # 35 pop, 27 spp

# load complete adult thermal physiology data
tphys_a <- read_csv("data/adult_thermal_phys.csv")

# determine number of populations
length(unique(paste(tphys_a$species, tphys_a$lat, tphys_a$elev))) # 89

# determine number of species
length(unique(tphys_a$species)) # 113

# determine number of species with life history trait info
length(unique(lh_traits$species))


## Eco physiology summaries ----

# read model_data
load("data/model_test_data.RData")

# obtain average estimates of temperatures experienced (Table S3)
model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(parity, species, lat, lon, elev) %>% 
  summarise(tg = mean(tg), tap = mean(tap), tep = mean(tep), tenv = mean(tenv)) %>% 
  ungroup()%>%
  group_by(parity) %>% 
  summarise(mtg = mean(tg), mtap = mean(tap), 
            mtep = mean(tep), mtenv = mean(tenv),
            sdtg = sd(tg), sdtap = sd(tap),
            sdtep = sd(tep), sdtenv = sd(tenv))

# obtain average estimates of temperatures experienced by embryo post gestation (Table S6)
model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(!is.na(depth), !is.na(shade)) %>% 
  group_by(parity, species, lat, lon, elev, depth, shade) %>% 
  summarise(tep = mean(tep)) %>% 
  ungroup() %>%
  group_by(parity, depth, shade) %>% 
  summarise(mtep = mean(tep), sdtep = sd(tep)) %>% 
  mutate(est = paste(round(mtep, 1), " ± ", round(sdtep, 2), sep = "")) %>% 
  select(depth, parity, shade, est) %>% 
  pivot_wider(names_from = "shade", values_from = "est") %>% 
  arrange(depth)

## Eco physiology model ----

# read model_data
load("data/model_test_data.RData")

# prepare data for model
eco_data <- model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(alpha == 0.5) |> # filter only species with alpha = 0.5
  filter(if_else(!is.na(gamma), gamma == 2, TRUE)) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(parity, species, lat, lon, elev) %>% 
  summarise(Dtg = mean(tg), Ctep = mean(tep), Btenv = mean(tenv)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c("Dtg", "Ctep", "Btenv"), 
               names_to = "type", values_to = "temp") %>% 
  mutate(abs_lat_scale = scale(abs(lat)),
         elev_scale = scale(elev), 
         temp_scale = scale(temp))

# run model
eco_data_model <- lmer(temp_scale ~ type * parity + abs_lat_scale + elev_scale + (1|species), 
                       data = eco_data)

# obtain model summaries
summary(eco_data_model)
confint(eco_data_model)

## d ~ |Topt - Tep| model ----

# load model data and embryonic development data to calculate Topt_e
load("data/model_test_data.RData")
topt_e <- read.csv("data/embryonic_thermal_phys.csv") |> 
  summarise(mean(topt_e, na.rm = T)) |> as.numeric()

# prepare data for model
model_dat <- model_test_data |> 
  filter(alpha == 0.5) |> # filter only species with alpha = 0.5
  filter(gamma == 2) %>% 
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
confint(model)

# prepare data for model
model_dat <- model_test_data |> 
  filter(dev_check == 1) |> # filter only months when embryonic development occurs
  filter(if_else(!is.na(gamma), gamma == 2, TRUE)) %>% 
  filter(is.na(nest_depth)) |>
  filter(is.na(nest_shade)) |>
  mutate(abs_diff = abs(topt_e - tep)) |> # calculate absolute Topt_e to tep difference
  group_by(parity, species, lat, lon, elev, shade, depth) |>
  summarise(d = mean(opt_d, na.rm = T), # get average d*
            diff = mean(abs_diff, na.rm = T)) |> # get average difference
  ungroup() |>
  mutate(zabslat = as.vector(scale(abs(lat))),
         zelev = as.vector(scale(elev)),
         zdiff = as.vector(scale(diff)))

# run model
model <- glmer(d ~ diff * shade * depth + (1|species), family = binomial(link = "logit"),
               data = model_dat)

# obtain model summaries
summary(model)


## Model results summaries ----

# read model_data
load("data/model_test_data.RData")

# average d* across all species by parity mode
model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(alpha == 0.5) %>% 
  filter(gamma  == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(opt_d = mean(opt_d, na.rm = T)) %>% 
  group_by(parity) %>%
  summarise(mopt_d = mean(opt_d, na.rm = T),
            sdopt_d = sd(opt_d, na.rm = T))

# % of oviparous and viviparous species with d == 1, 0.5 < d < 1, d <= 0.5
model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(alpha == 0.5) %>% 
  filter(gamma  == 2) %>%  
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(opt_d = mean(opt_d, na.rm = T)) %>% 
  #filter(parity == "Oviparous") |> 
  #filter(opt_d > 0.5) |>
  group_by(parity) %>%
  summarise(d1 = sum(opt_d == 1, na.rm = T),
            d05 = sum(opt_d > 0.5 & opt_d < 1, na.rm = T),
            d0 = sum(opt_d <= 0.5, na.rm = T))

## Model results statistics ----

# read model_data
load("data/model_test_data.RData")

# prepare data for stats
model_test_data_stats <- model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity, n, ze, za) %>% 
  summarise(opt_d = mean(opt_d, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(msp = (ze * n)/za) %>% 
  mutate(z_lat = as.vector(scale(abs(lat))),
         z_elev = as.vector(scale(elev)),
         z_msp = as.vector(scale(msp)))

# run model
d_model <- glmer(opt_d ~ z_lat * z_elev * z_msp + (1|species), 
                 family = binomial(link = "logit"),
                 data = model_test_data_stats)
summary(d_model)
confint(d_model)
  
