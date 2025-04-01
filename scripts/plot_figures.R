
## Plotting figures ## 

## Packages ----

library(ggrepel)
library(ggpubr)
library(mgcv)
library(ggh4x)
library(tidyverse)

## Figure 2 ----

# load adult database
adult_data <- read.csv("data/adult_data.csv")

# plot world data
world_plot <- ggplot() +
  geom_polygon(data =  map_data("world") %>% filter(region != "Antarctica"), 
               aes(x = long, y = lat, group = group),
               alpha = 0.4) +
  geom_point(data = adult_data, 
             aes(x = lon, y = lat, shape = parity, fill = parity), 
             size = 2) +
  scale_fill_manual(values = c("darkorange", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0),) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.6, 0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.direction = "horizontal",
        panel.grid = element_blank(),
        plot.margin = margin(r = 0)) +
  guides(fill = guide_legend(override.aes = list(size=4))) 

# latitude plot
lat_plot <-  ggplot(data = adult_data) +
  geom_histogram(aes(y = lat, fill = parity, x = ..count..), binwidth = 5,
                 position = position_dodge(preserve = "single"), col = "black",
                 size = 0.1) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Latitude (°)") +
  xlab("Number of populations") +
  scale_x_reverse(expand = c(0, 0), position = "top") + 
  scale_fill_manual(values = c("darkorange", "yellowgreen")) +
  scale_y_continuous( position = "right") +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        plot.margin = margin(l = -25, t = 50, b = 20, r = 25)) 

# elevation plot
elev_plot <-  ggplot(data = adult_data) +
  geom_histogram(aes(y = elev, fill = parity), binwidth = 250,
                 position = position_dodge(preserve = "single"), col = "black",
                 size = 0.1) +
  ylab("Elevation (m a.s.l.)") +
  xlab("Number of populations") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("darkorange", "yellowgreen")) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6),
        plot.margin = margin(r = -100, t = 50, b = 10, l = 100))

# combine plots
figure2 <- ggarrange(elev_plot, world_plot, lat_plot,
                     widths = c(0.15, 0.6, 0.1),
                     ncol = 3)

# save figure
ggsave(figure2, file = "figures/figure2.png", dpi = "retina")

## Figure 3 ----

# load model data
load("data/model_test_data.RData")

# load embryonic thermal physiology data
tphys_e <- read.csv("data/embryonic_thermal_phys.csv")

# define tpc function
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

embryo_tpc <- tibble(temp = seq(5,35,by = 0.1),
                     perf = tpc(temp, ctmin = mean(tphys_e$ctmin_e, na.rm = T),
                                topt = mean(tphys_e$topt_e, na.rm = T),
                                ctmax = mean(tphys_e$ctmax_e, na.rm = T)))
# plot 
figure3 <- model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>% 
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(tg = mean(tg, na.rm = T), 
            tep = mean(tep, na.rm = T),
            tenv = mean(tenv, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(tg, tep, tenv), names_to = "st", values_to = "temp") %>%
  ggplot(aes(x = temp)) +
  geom_point(aes(y = parity, col = st, fill = st, shape = st), 
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05),
             alpha = 0.5, size = 2) +
  stat_summary(aes(y = parity, shape = st, fill = st),
               position = position_dodge(width = 0.5),
               size = 0.75, linewidth = 0.75, col = "black") +
  geom_line(data = embryo_tpc, aes(y = 2.5*perf), linetype = 2) +
  scale_shape_manual(values = c(23,22,21),
                     labels = c("Environment", "Post-gestation", "Gestation")) +
  scale_color_manual(values = c("darkgray" ,"darkslategray", "tomato"), 
                     labels = c("Environment","Post-gestation", "Gestation")) +
  scale_fill_manual(values = c("darkgray", "darkslategray", "tomato"), 
                    labels = c("Environment", "Post-gestation", "Gestation")) +
  theme_minimal() +
  theme(
    legend.position= c(0.15,0.13),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    panel.grid = element_blank()
  ) +
  xlab("Temperature (°C)") +
  ylab("Parity mode")

ggsave(figure3, file = "figures/figure3.png", dpi = "retina")

## Figure 4 ----

figure4 <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T)) %>% 
  ggplot(aes(x = parity, y = mean_opt_d, fill = parity, shape = parity)) +
  geom_violin(alpha = 0.25, color = NA) +
  geom_jitter(col = "black", width = 0.1, size = 2, alpha = 0.6) +
  stat_summary(size = 1, linewidth = 0.75) +
  scale_fill_manual(values = c("orange2", "yellowgreen")) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, colour = "black"),
    legend.position = "none"
  ) +
  ylab("Predicted optimal gestation length (d*)")

ggsave(figure4, file = "figures/figure4.png", dpi = "retina")

## Figure 5 ----

# Panel A 
f5a <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T)) %>% 
  ggplot(aes(x = abs(lat), y = mean_opt_d)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              col = "black", alpha = 0.25) +
  geom_point(
    aes(col = parity, fill = parity, shape = parity),
    size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("orange2", "yellowgreen")) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab("Absolute latitude (°)")+
  theme(plot.title = element_text(size = 12, hjust = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))  +
  guides(colour = guide_legend(override.aes = list(size=4)))

# Panel B
f5b <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T)) %>% 
  ggplot(aes(x = elev, y = mean_opt_d)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              col = "black", alpha = 0.25) +
  geom_point(
    aes(col = parity, fill = parity, shape = parity),
    size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("orange2", "yellowgreen")) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.title.y = element_text(colour = "white")
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab("Elevation (m a.s.l.)")+
  theme(plot.title = element_text(size = 12, hjust = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))  +
  guides(colour = guide_legend(override.aes = list(size=4)))

# Panel C
f5c <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity, n, za, ze) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T)) %>% 
  ggplot(aes(x = (n*ze)/za, y = mean_opt_d)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              col = "black", alpha = 0.25) +
  geom_point(
    aes(col = parity, fill = parity, shape = parity),
    size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("orange2", "yellowgreen")) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab("Mass-specific production") +
  theme(plot.title = element_text(size = 12, hjust = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))  +
  guides(colour = guide_legend(override.aes = list(size=4)))

# panel D
f5d <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>%
  group_by(species, lat, lon, elev, parity) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T),
            tep = mean(tep, na.rm = T)) %>% 
  ggplot(aes(x = abs(tep - 27.57), y = mean_opt_d)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              col = "black", alpha = 0.25) +
  geom_point(
    aes(col = parity, fill = parity, shape = parity),
    size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("orange2", "yellowgreen")) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    plot.title = element_text(size = 12, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(colour = "white")
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab(expression(paste("| ",T[opt[e]] - T[ep], " |"))) +
  coord_cartesian(ylim = c(0.01, 0.97)) +
  guides(colour = guide_legend(override.aes = list(size=4)))

# panel E
f5e <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(!is.na(depth)) %>% 
  filter(!is.na(shade)) %>%
  group_by(species, lat, lon, elev, parity, tenv, depth, shade) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T),
            tep = mean(tep, na.rm = T),
            tenv = mean(tenv, na.rm = T)) %>% 
  ungroup() %>% 
  filter(tenv < mean(tenv)) %>% 
  mutate(shade_factor = as.factor(paste(shade*100, "%"))) %>% 
  mutate(shade_factor = fct_reorder(shade_factor, shade)) %>%
  ggplot(aes(x = depth, y = mean_opt_d, col = shade_factor)) +
  geom_jitter(width = 0.1, height = 0, size = 1.5, alpha= 0.25) +
  geom_smooth(aes(col = shade_factor, fill = shade_factor),
                  method = "glm", method.args = list(family = "binomial")) +
  scale_color_manual(values = c("#FFCC99", "#806680", "#000066"), name = "Nest shade") +
  scale_fill_manual(values = c("#FFCC99", "#806680", "#000066"), name = "Nest shade") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "top"
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab("Nest depth (cm)")

# panel F
f5f <- model_test_data %>% 
  as_tibble() %>%
  filter(dev_check == 1) %>% 
  filter(alpha == 0.5) %>% 
  filter(gamma == 2) %>% 
  filter(!is.na(depth)) %>% 
  filter(!is.na(shade)) %>%
  group_by(species, lat, lon, elev, parity, tenv, depth, shade) %>% 
  summarise(mean_opt_d = mean(opt_d, na.rm = T),
            tep = mean(tep, na.rm = T),
            tenv = mean(tenv, na.rm = T)) %>% 
  ungroup() %>% 
  filter(tenv > mean(tenv)) %>% 
  mutate(shade_factor = as.factor(paste(shade*100, "%"))) %>% 
  mutate(shade_factor = fct_reorder(shade_factor, shade)) %>%
  ggplot(aes(x = depth, y = mean_opt_d, col = shade_factor)) +
  geom_jitter(width = 0.1, height = 0, size = 1.5, alpha= 0.25) +
  geom_smooth(aes(col = shade_factor, fill = shade_factor),
              method = "glm", method.args = list(family = "binomial")) +
  scale_color_manual(values = c("#FFCC99", "#806680", "#000066"), name = "Nest shade") +
  scale_fill_manual(values = c("#FFCC99", "#806680", "#000066"), name = "Nest shade") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "top",
    axis.title.y = element_text(colour = "white")
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab("Nest depth (cm)")

# combine top panels
fig5top <- ggarrange(f5a, f5b, f5c, f5d, ncol = 2, nrow = 2, 
          common.legend = TRUE, legend = "top",
          labels = c("A)", "B)", "C)", "D)"),
          font.label = list(face = "plain"))

# combine bottom panels
fig5bottom <- ggarrange(f5e, f5f, ncol = 2, nrow = 1, 
                        common.legend = TRUE, legend = "top",
                        labels = c("E)", "F)"),
                        font.label = list(face = "plain"))

# combine all panels
figure5 <- ggarrange(fig5top, fig5bottom, ncol = 1, nrow = 2, heights = c(0.65,0.35))

# save figure
ggsave(figure5, file = "figures/figure5.png", dpi = "retina")

## Figure S1 ----

# load embryonic thermal physiology data
tphys_e <- read.csv("data/embryonic_thermal_phys.csv")
tphys_e$X <- NULL

# re-structure thphys_e data
tphys_e <- tphys_e |>
  pivot_longer(cols = c("ctmin_e", "topt_e", "ctmax_e"), 
              names_to = "trait", values_to = "temp") |> 
  mutate(trait = ifelse(trait == "ctmin_e", "ctmin",
                        ifelse(trait == "topt_e", "topt", "ctmax"))) |> 
  mutate(stage = "Embryo")
  
# load adult thermal physiology data
adult_data <- read.csv("data/adult_data.csv")

# filter columns of interest and re-structure tphys_a data
tphys_a <- adult_data |> select(species, parity, ctmin_a, topt_a, ctmax_a) |> 
  pivot_longer(cols = c("ctmin_a", "topt_a", "ctmax_a"), 
              names_to = "trait", values_to = "temp") |> 
  mutate(trait = ifelse(trait == "ctmin_a", "ctmin",
                        ifelse(trait == "topt_a", "topt", "ctmax"))) |> 
  mutate(stage = "Adult")

# combine adult and embryo data
tphys <- rbind(tphys_a, tphys_e)

# get tphys_summary data 
tphys_summary <- tphys %>% 
  group_by(parity, trait, stage) %>% 
  summarise(mean = mean(temp, na.rm = T), 
            sd = sd(temp, na.rm = T)) %>% 
  ungroup()

# plot Figure S1
figures1 <- tphys %>% 
  mutate(trait = factor(trait, levels = c("ctmin", "topt", "ctmax"))) %>% 
  ggplot() +
  geom_point(aes(x = trait, y = temp, fill = parity, col = parity, shape = stage),
             position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1), 
             size = 2, alpha = 0.5) +
  geom_pointrange(data = tphys_summary,
                  aes(x = trait, y = mean, ymin = mean - sd, 
                      ymax = mean + sd, fill = parity, col = parity, shape = stage),
                  position = position_dodge(width = 1), size = 1, linewidth = 0.75) +
  geom_pointrange(data = tphys_summary,
                  aes(x = trait, y = mean, ymin = mean - sd, 
                      ymax = mean + sd, fill = parity, shape = stage),
                  position = position_dodge(width = 1), 
                  col = "black", size = 1, linewidth = 0.75, show.legend = FALSE) +
  scale_x_discrete(labels = c(bquote(CT[min]), bquote(T[opt]), bquote(CT[max]))) +
  scale_fill_manual(values = c("orange2", "yellowgreen")) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, colour = "black"),
  ) +
  ylab("Temperature (°C)") +
  guides(fill = guide_legend(title = "Parity mode"),
         color = guide_legend(title = "Parity mode"),
         shape = guide_legend(title = "Life stage"))

# save figure s1
ggsave(figures1, file = "figures/figureS1.png", dpi = "retina")

## Figure S2 ----

# load embryonic thermal physiology data
tphys_e <- read.csv("data/embryonic_thermal_phys.csv")
tphys_e$X <- NULL

# load adult thermal physiology data
adult_data <- read.csv("data/adult_data.csv")

# filter columns of interest and re-structure tphys_a data
tphys_a <- adult_data |> select(species, parity, ctmin_a, topt_a, ctmax_a) 

# merge data
tphys <- merge(tphys_e, tphys_a, by = c("species","parity"), all = TRUE)

# plot figure S2 A 
s2_A <- tphys %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  select(parity, species, ctmin_a, ctmin_e) %>%
  filter(!is.na(ctmin_e)) %>%
  filter(!is.na(ctmin_a)) %>%
  unique() %>%
  ggplot(aes(x = ctmin_a, y = ctmin_e)) +
  geom_point(aes(col = species), size = 3) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  xlab(bquote(CT[min[a]]~"("~degree~C~")")) +
  ylab(bquote(CT[min[e]]~"("~degree~C~")")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    legend.text = element_text(face = "italic"),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line()
  ) +
  ggtitle("A")

# plot figure S2 B
s2_B <- tphys %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  select(parity, species, topt_a, topt_e) %>%
  filter(!is.na(topt_e)) %>%
  filter(!is.na(topt_a)) %>%
  unique() %>%
  ggplot(aes(x = topt_a, y = topt_e)) +
  geom_point(aes(col = species), size = 3) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  xlab(bquote(T[opt[a]]~"("~degree~C~")")) +
  ylab(bquote(T[opt[e]]~"("~degree~C~")")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    legend.text = element_text(face = "italic"),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line()
  )+
  ggtitle("B")

# plot figure S2 C
s2_C <- tphys %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  select(parity, species, ctmax_a, ctmax_e) %>%
  filter(!is.na(ctmax_e)) %>%
  filter(!is.na(ctmax_a)) %>%
  unique() %>%
  ggplot(aes(x = ctmax_a, y = ctmax_e)) +
  geom_point(aes(col = species), size = 3) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  xlab(bquote(CT[max[a]]~"("~degree~C~")")) +
  ylab(bquote(CT[max[e]]~"("~degree~C~")")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    legend.text = element_text(face = "italic"),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line()
  ) +
  ggtitle("C")

# save 
figS2 <- ggarrange(s2_A, s2_B, s2_C, ncol = 1, common.legend = FALSE)
ggsave(figS2, file = "figures/figureS2.png", dpi = "retina")
  
## Figure S3 ----

# load model data
load("data/model_test_data.RData")

# plot figure
figS3 <- model_test_data %>% 
  select(species, parity, lat, elev, dev) %>% 
  unique() %>% 
  mutate(hem = ifelse(lat >0,"Northern hemisphere", "Southern hemisphere")) %>% 
  group_by(parity, hem) %>%
  mutate(n_pop = n()) %>% 
  as_tibble() %>% 
  unnest(dev) %>%
  group_by(parity, hem, dev, n_pop) %>% 
  summarise(howmany = n()) %>%
  mutate(freq = howmany/n_pop) %>%
  ungroup() %>%
  ggplot(aes(x = dev, y = freq*100, fill = parity)) +
  geom_col(col = "white", position = position_dodge(width = 1, preserve = "single")) +
  scale_fill_manual(values = c("darkorange", "yellowgreen")) +
  scale_x_continuous(breaks = 1:12, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,25,50,75,100)) +
  facet_wrap(vars(hem), axes = "all") +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        strip.text = element_text(size = 10),
        axis.line = element_line(),
        axis.ticks = element_line(),
        legend.text = element_text(size = 10),
        panel.grid = element_blank()) +
  xlab("Month of the year") +
  ylab("% of populations") 

ggsave(figS3, file = "figures/figureS3.png", dpi = "retina")

## Figure S4 ----

# load data
load("data/model_test_data.RData")

# plot
figs4 <- model_test_data %>% 
  select(species, parity, lat, elev, eggs_seen, hatchlings_seen, dev) %>% 
  unique()  %>% 
  as_tibble() %>%
  unnest(dev) %>% 
  mutate(hemisphere = ifelse(lat < 0, "S", "N")) %>%
  mutate(species = gsub("_", " ", species)) %>%
  mutate(population = paste(species, " (", round(lat, digits = 1), 
                            "°", " ,", elev, " m)", sep = "")) %>%
  group_by(hemisphere) %>% 
  arrange(desc(population)) %>% 
  ggplot(aes(y = fct_inorder(population), x = dev, fill = parity, group = hemisphere)) +
  geom_tile(alpha = 0.75, col = NA) +
  geom_tile(aes(x = eggs_seen))+
  facet_grid(rows = vars(hemisphere), scales = "free_y", space = "free") +
  xlab("Month of the year") +
  scale_x_continuous(breaks = 1:12, expand = c(0,0)) +
  scale_fill_manual(values = c("darkorange", "yellowgreen")) +
  scale_color_manual(values = c("darkorange", "yellowgreen")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(), 
        legend.position = "top",
        panel.border = element_rect(fill = NA, linewidth = 0.1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(),
        axis.text.y = element_text(face = "italic", size = 6),
        panel.grid = element_line(linewidth = 0.1),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.title = element_text(hjust = -0.38, vjust=2.12))

# save figure
ggsave(figs4, file = "figures/figureS4.png", dpi = "retina")

## Figure S5 ----

# load data
load("data/model_test_data.RData")

# plot
figs5 <- model_test_data %>% 
  filter(dev_check == 1) %>%
  filter(!is.na(depth), !is.na(shade)) %>% 
  group_by(parity, species, lat, lon, elev, depth, shade) %>% 
  summarise(tep = mean(tep)) %>% 
  ungroup() %>%
  ggplot(aes(x = depth, y = tep, col = as.factor(shade*100))) +
  geom_jitter(width = 0.25, alpha = 0.5) +
  stat_summary(geom = "line", linewidth = 1) +
  scale_color_manual(values = c("#FFCC99", "#806680", "#000066")) +
  scale_y_continuous(limits = c(10,40)) +
  facet_wrap(~parity, scales = "free") +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.ticks = element_line(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
  ) +
  xlab("Nest depth (cm belowground)") +
  ylab(bquote(~T[ep]~"(°C)")) +
  guides(color = guide_legend(title = "Nest shade (%)"))

# save figure
ggsave(figs5, file = "figures/figureS5.png", dpi = "retina")
  
## Figure S6 ----

# load data
load("data/model_test_data.RData")

# plot
figures6 <- model_test_data %>% 
  as_tibble() %>%
  filter(if_else(!is.na(depth), depth == 5, TRUE)) %>% 
  filter(if_else(!is.na(shade), shade == 0.5, TRUE)) %>% 
  group_by(parity, species, lat, lon, elev, alpha, gamma) %>% 
  summarise(opt_d = mean(opt_d, na.rm = T)) %>% 
  ggplot(aes(x = gamma, y = opt_d, col = parity, linetype = as.factor(alpha))) +
  stat_summary(geom = "line", linewidth = 1) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,1,2,5,10)) +
  scale_color_manual(values = c("orange2", "yellowgreen")) +
  theme_minimal() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.title.x = element_text(size = 16),
  ) +
  ylab(expression(paste(Model~prediction~(d~"*")))) +
  xlab(expression(gamma)) +
  guides(color = guide_legend(title = "Parity mode"),
         linetype = guide_legend(title = expression(alpha),
                                 theme = theme(legend.title = element_text(size = 16))))

# save
ggsave(figures6, file = "figures/figures6.png", dpi = "retina")


