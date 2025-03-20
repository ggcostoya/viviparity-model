
set.seed(123)

# plot conceptual figure

library(tidyverse)
library(ggpubr)

# define function to simulate soil temperature
soil_temp_fluct <- function(tmin, tmax, noise){
  time <- 0:23 # define time
  amp <- (tmax-tmin)/2 # define amplitude
  mean <- (tmax+tmin)/2 # define mean temperature
  temp <- mean + amp * sin(2 * pi * (time - 10) / 24) # calculate temp
  temp <- temp + rnorm(length(seq(temp)), 0, noise) # add noise
  return(temp)
}

# define function to simulate gestation temperature
gest_temp <- function(topt, ctmin, correction, noise){
  time <- 0:23
  temp <- soil_temp_fluct(tmin = ctmin + 1, tmax = ctmin + correction, noise = 0)
  temp <- ifelse(time > 5 & time < 18, topt, temp)
  temp <- temp + rnorm(length(seq(temp)), 0, noise)
  return(temp)
}

# generate mock data for optimal climate
optimal_climate_temp <- data.frame(
  time = 0:23,
  gest = gest_temp(topt = topt_e, ctmin = ctmin_e + 2, noise = 0.1, correction = 5),
  post = soil_temp_fluct(tmin = 19, tmax = 32, noise = 0.25)
) |> 
  pivot_longer(cols = c("gest", "post"), values_to = "temp", names_to = "stage") |> 
  mutate(stage = ifelse(stage == "gest", "Gestation", "Post-gestation")) |>
  mutate(climate = "Optimal climate")

# generate mock data for sub-optimal climate
suboptimal_climate_temp <- data.frame(
  time = 0:23,
  gest = gest_temp(topt = topt_e, ctmin = ctmin_e, noise = 0.1, correction = 5),
  post = soil_temp_fluct(tmin = 16.5, tmax = 26.5, noise = 0.25)
) |> 
  pivot_longer(cols = c("gest", "post"), values_to = "temp", names_to = "stage") |> 
  mutate(stage = ifelse(stage == "gest", "Gestation", "Post-gestation")) |>
  mutate(climate = "Suboptimal climate") 

# combine data
climate_temp <- rbind(optimal_climate_temp, suboptimal_climate_temp)

# get averages per climate and stage
climate_temp_avg <- climate_temp |> group_by(climate, stage) |> 
  summarize(mean_temp = mean(temp)) |> ungroup()

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

# calculate performance
climate_temp_avg$surv <- tpc(climate_temp_avg$mean_temp, ctmin_e, topt_e, ctmax_e)

# generate mock_tpc data
tpc_data <- data.frame(
  temp = seq(15, 35, 0.1),
  surv = tpc(seq(15, 35, 0.1), ctmin_e, topt_e, ctmax_e)
)

tgopt <- climate_temp_avg$mean_temp[1] - 3
tgsub <- climate_temp_avg$mean_temp[3] 
tepopt <- climate_temp_avg$mean_temp[2] - 3
tepsub <- climate_temp_avg$mean_temp[4] 

# prep test data
figure_pred_data <- data.frame(
  climate = c("Optimal climate", "Suboptimal climate"),
  topt_a = 33.6, ctmin_a = 13.1, ctmax_a = 42.5,
  za = 17.2, ze = 0.85, n = 9.5,
  tg = c(tgopt, tgsub), tep = c(tepopt, tepsub), tap = c(29, 23), 
  alpha = 0.5, gamma = 2
)

# modified model function
model <- function(species, topt_a, ctmin_a, ctmax_a, topt_e, ctmin_e, ctmax_e,
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
  
  # return data
  pred_data <- data.frame(d = d, fitness = fitness, species = species)
  
  return(pred_data)
  
}

# holder dataset
pred_holder <- data.frame(d = c(), fitness = c(), species = c())

# loop to predict fitness
for(i in 1:nrow(figure_pred_data)){
  
  # run model
  pred_data <- model(species = figure_pred_data$climate[i],
                     topt_a = figure_pred_data$topt_a[i],
                     ctmin_a = figure_pred_data$ctmin_a[i],
                     ctmax_a = figure_pred_data$ctmax_a[i],
                     topt_e = topt_e,
                     ctmin_e = ctmin_e,
                     ctmax_e = ctmax_e,
                     za = figure_pred_data$za[i],
                     ze = figure_pred_data$ze[i],
                     n = figure_pred_data$n[i],
                     tg = figure_pred_data$tg[i],
                     tep = figure_pred_data$tep[i],
                     tap = figure_pred_data$tap[i],
                     alpha = figure_pred_data$alpha[i],
                     gamma = figure_pred_data$gamma[i])
  
  # bind to holder data set
  pred_holder <- rbind(pred_holder, pred_data)
  
}

# plot panels A and D
labels_top <- data.frame(letter = c("A)", "D)"), 
                         climate = c("Optimal climate", "Suboptimal climate"))

label_traits <- data.frame(time = 1, climate = "Optimal climate")

top <- climate_temp |> 
  ggplot(aes(x = time, y = temp)) +
  geom_hline(yintercept = topt_e, alpha = 0.5, col = "lightgray") +
  geom_hline(yintercept = ctmin_e, alpha = 0.5, col = "lightgray") +
  geom_hline(yintercept = ctmax_e, alpha = 0.5, col = "lightgray") +
  geom_line(aes(color = stage), linewidth = 0.75) +
  geom_hline(data = climate_temp_avg, aes(yintercept = mean_temp, col = stage),
             linetype = 5, linewidth = 0.75) +
  geom_text(data = labels_top, aes(x = 21.5, y = 35, label = letter),
            size = 5) +
  annotate(geom = "text", x = 0.5, y = topt_e - 1, hjust = 0, 
           label = expression(T[opt[e]]), size = 2.75) +
  annotate(geom = "text", x = 0.5, y = ctmin_e - 1, hjust = 0, 
           label = expression(CT[min[e]]), size = 2.75) +
  annotate(geom = "text", x = 0.5, y = ctmax_e - 1, hjust = 0, 
           label = expression(CT[max[e]]), size = 2.75) +
  facet_wrap(~climate, axes = "all") +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("tomato", "darkslategray4")) +
  scale_alpha_manual(values = c(1,0)) +
  xlab("Hour of the day") +
  ylab("Embryo body temperature (°C)") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 11)),
    panel.spacing.x = unit(0.75, "lines")
  )

# plot panels B and E
climate_temp_avg_mod <- climate_temp_avg 
climate_temp_avg_mod$mean_temp[2] <- climate_temp_avg$mean_temp[2]
climate_temp_avg_mod$surv <- tpc(climate_temp_avg_mod$mean_temp, ctmin_e, topt_e, ctmax_e)


labels_middle <- data.frame(letter = c("B)", "E)"), 
                            climate = c("Optimal climate", "Suboptimal climate"))

middle <- ggplot() +
  geom_line(data = tpc_data, aes(x = temp, y = surv), linewidth = 0.75) +
  geom_segment(data = climate_temp_avg_mod, 
               aes(x = mean_temp, xend = mean_temp, y = 0, yend = surv, col = stage),
               linetype = 5, linewidth = 0.75) +
  geom_segment(data = climate_temp_avg_mod, 
               aes(x = 15, xend = mean_temp, y = surv, yend = surv, col = stage),
               linetype = 5, linewidth = 0.75) +
  geom_text(data = labels_middle, aes(x = 34, y = 0.97, label = letter),
            size = 5) +
  theme_classic() +
  scale_y_continuous(limits = c(0,1.01), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("tomato", "darkslategray4")) +
  facet_wrap(~climate, axes = "all") +
  xlab("Embryo body temperature (°C)") +
  ylab("Embryo survival probability") +
  theme(
    legend.position = "none",
    strip.text = element_blank()
  )

# calculate place holder optimas
pred_holder_optimas <- pred_holder |> group_by(species) |> 
  filter(fitness == max(fitness)) |> ungroup()
pred_holder_optimas$fitness[2] <- pred_holder_optimas$fitness[2] - 1

labels_bottom <- data.frame(letter = c("C)", "F)"), 
                            species = c("Optimal climate", "Suboptimal climate"))

# plot panels C and F
bottom <- pred_holder |> 
  mutate(fitness = ifelse(species == "Suboptimal climate", fitness - 1, fitness)) |> 
  ggplot(aes(x = d, y = fitness)) +
  geom_hline(yintercept = 2, linetype = 2, col = "lightgray") +
  geom_line(linewidth = 0.75) +
  geom_point(data = pred_holder_optimas, size = 2) +
  geom_segment(data = pred_holder_optimas, 
               aes(x = d, xend = d, y = 0, yend = fitness),
               linetype = 2) +
  geom_text(data = pred_holder_optimas,
            aes(y = fitness + 0.5, x = d), label = expression(paste(italic(d),"*"))) +
  geom_text(data = labels_bottom, aes(x = 0.98, y = 4.77, label = letter),
            size = 5) +
  facet_wrap(~species, axes = "all") +
  scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  theme_classic()+
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    axis.title.y = element_text(margin = margin(r = 14.5)),
    panel.spacing.x = unit(1.25, "lines")
  ) +
  xlab(expression(paste("Gestation length (", italic(d),")"))) +
  ylab("Fitness") 

plot <- ggarrange(top, middle, bottom, ncol = 1, nrow = 3, heights = c(1.5, 1, 1))

ggsave(plot, filename = "test.png", dpi = "retina")



