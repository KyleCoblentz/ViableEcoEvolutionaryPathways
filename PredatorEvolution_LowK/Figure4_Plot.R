################################################################################
### Figure 4 code -- evolution on the landscape
################################################################################

### load libraries

library(dplyr); library(ggplot2); library(cowplot);library(ggnewscale); library(viridis)

### load data for the lanscape

ExtData<- read.csv('ProportionExt_HighIntComp.csv')

################################################################################
### make heatmap
################################################################################

ggplot(data = ExtData, aes(x = HandlingTime, y = AttackRate, fill = NumberPersist)) + geom_tile() + 
  theme_cowplot()

### want to add lines with analytical stability and feasibility boundaries

q <- 0.005; r <- 3.5; e <- 0.3; m <- 0.6;

feasibility <- function(h, q, e, m) {
  
  hout <- h
  
  aout <- q*m/((e-h*m))
  
  out <- data.frame(h = hout,
                    a = aout)
  
  return(out)
  
}

stability <- function(h, q, e, m) {
  
  hout <- h
  
  aout <- q*(e + h*m)/(h*(e - h*m))
  
  out <- data.frame(h = hout,
                    a = aout)
  
  return(out)
  
}



f_high <- feasibility(h = seq(0.01,0.65, length.out = 1000), q = q, e = e, m =m)

s_high <- stability(h = seq(0.01, 0.65, length.out = 1000), q = q, e = e, m = m)

### put them all together

HeatMap_PopPersist <- ggplot(data = ExtData, aes(x = HandlingTime, y = AttackRate)) + geom_tile(aes(fill = NumberPersist)) + 
  theme_cowplot() + geom_line(data = f_high, aes(x = h, y = a), color = 'white')  + 
  geom_line(data = s_high, aes(x = h, y = a), color = 'white') + xlab('Handling Time') + 
  ylab('Space Clearance Rate') + labs(fill = 'Percent Persist') + 
  ylim(c(0,0.15)) + xlim(c(0, 0.5)) + theme_cowplot() +
  annotate(geom = 'text', x = 0.225, y = 0.1, color = 'white', label = 'Unstable') + 
  annotate(geom = 'text', x = 0.4, y = 0.015, color = 'white', label = 'Infeasible')

HeatMap_PopPersist2 <- ggplot(data = ExtData, aes(x = HandlingTime, y = AttackRate)) + geom_tile(aes(fill = NumberPersist)) + 
  theme_cowplot() + geom_line(data = f_high, aes(x = h, y = a), color = 'white')  + 
  geom_line(data = s_high, aes(x = h, y = a), color = 'white') + xlab('Handling Time') + 
  ylab('Space Clearance Rate') + labs(fill = 'Percent Persist') + 
  ylim(c(0,0.15)) + xlim(c(0, 0.5)) + theme_cowplot() 

#########################################################################################
### Individual evolution runs for plot
#########################################################################################

### load individual run data

inddata <- read.csv('indrundata.csv', na.strings = 'NaN')

### add maximum attack rate points for each run to the heatmap plot

### create dataset that has row containing the mean end of the runs

SCRMeanEnd <- inddata %>% filter(!is.na(SCRate)) %>% group_by(Pop) %>%  slice_tail(., n = 1) %>% group_by(Pop) %>% summarise(SCRate = mean(SCRate, na.rm = TRUE), HandlingTime = mean(HandlingTime, na.rm = TRUE), Extinct = mean(Extinct, na.rm = TRUE))

SCRMeanEnd <- SCRMeanEnd %>% mutate(ExtinctFactor = ifelse(Extinct == 1, 'Extinct', 'Extant'))

SCRStart <- inddata %>% filter(!is.na(SCRate)) %>% group_by(Pop) %>% slice_head(., n = 1) 

### add points to the heat map that are different colors whether the population went extinct or not

ExtinctionHeatMap <- HeatMap_PopPersist2 + new_scale_color() + 
  geom_point(data = SCRMeanEnd, aes(x = HandlingTime, y = SCRate, shape = ExtinctFactor, color = ExtinctFactor ), size = 1, stroke = 1.25) + 
  geom_point(aes(x = 0.25, y = 0.04), size = 3, color = 'white')   + 
  scale_color_manual(values = c('magenta', 'green'), breaks = c('Extinct', 'Extant'), name = 'Extinct?') + 
  scale_shape_manual(values = c(4, 16), breaks = c('Extinct', 'Extant'), name = 'Extinct?')

### average trajectories of extinct and extant populations

ExtinctTrajectories <- inddata %>% filter(!is.na(SCRate)) %>% mutate(ExtinctFactor = ifelse(Extinct == 1, "Extinct", "Extant")) %>% 
  filter(Extinct == 1) %>% group_by(Pop) %>% slice_tail(., n = 1) %>% summarise(SCRate = mean(SCRate), HandlingTime = mean(HandlingTime))

ExtinctStart <- SCRStart %>% filter(Extinct == 1) %>% ungroup() %>% summarise(SCRate = mean(SCRate),  HandlingTime = mean(HandlingTime))

ExtinctTrajectory <- data.frame(SCRate = c(ExtinctStart$SCRate, mean(ExtinctTrajectories$SCRate)),
                                HandlingTime = c(ExtinctStart$HandlingTime, mean(ExtinctTrajectories$HandlingTime)))

### Extant Trajectories

ExtantTrajectories <- inddata %>% filter(!is.na(SCRate)) %>% mutate(ExtinctFactor = ifelse(Extinct == 1, "Extinct", "Extant")) %>% 
  filter(Extinct == 0) %>% group_by(Pop) %>% slice_tail(., n = 1) %>% summarise(SCRate = mean(SCRate), HandlingTime = mean(HandlingTime))

ExtantStart <- SCRStart %>% filter(Extinct == 0) %>%  ungroup() %>% summarise(SCRate = mean(SCRate), HandlingTime = mean(HandlingTime))

ExtantTrajectory <- data.frame(SCRate = c(ExtantStart$SCRate, mean(ExtantTrajectories$SCRate)),
                               HandlingTime = c(ExtantStart$HandlingTime, mean(ExtantTrajectories$HandlingTime)))

### Quantitative Genetics trajectory

### load data

QG_data <- read.csv('TraitDynamics_QG.csv')

##################################################################################################
### add trajectories to the plots
##################################################################################################

AverageTrajectoryPlot <- ExtinctionHeatMap +   
  geom_line(data = QG_data, aes(x =HandlingTime, y = SpaceClearanceRate), color = 'white', 
               arrow = arrow(length = unit(0.03, "npc"), ends = 'first'))

##################################################################################################
### add mean trait values to the plot
##################################################################################################

### load trait data

inddata <- inddata %>% mutate(ExtinctFactor = ifelse(Extinct == 1, 'Extinct', 'Extant'))

AvgTrajectory <- inddata  %>% group_by(Time, ExtinctFactor) %>% summarize(SCRate = mean(SCRate, na.rm = TRUE), HandlingTime = mean(HandlingTime, na.rm = TRUE))

SCRatePopulations <- ggplot(data = inddata, aes(x = Time, y = SCRate, group = Pop, color = ExtinctFactor)) + geom_line(alpha = 0.5) + 
  theme_cowplot() + scale_color_manual(values = c('Extinct' = 'magenta', 'Extant' = 'green'), name = "Extinct?") + ylab('Space Clearance\n Rate') + 
  geom_line(data = QG_data, aes(x = Time, y = SpaceClearanceRate), inherit.aes = FALSE)
  
SCRateAvgPopulations <- ggplot(data = AvgTrajectory, aes(x = Time, y = SCRate, color = ExtinctFactor)) + geom_line(size = 1) +
  theme_cowplot() + scale_color_manual(values = c('Extinct' = 'magenta', 'Extant' = 'green'), name = 'Extinct?') + ylab('Average Space\n Clearance Rate') + 
  geom_line(data = QG_data, aes(x = Time, y = SpaceClearanceRate), inherit.aes = FALSE)

HandlingTimePopulations <- ggplot(data = inddata, aes(x = Time, y = HandlingTime, group = Pop, color = ExtinctFactor)) + geom_line(alpha = 0.5) + 
  theme_cowplot() + scale_color_manual(values = c('Extinct' = 'magenta', 'Extant' = 'green'), name = "Extinct?") + ylab('Handling Time') + 
  geom_line(data = QG_data, aes(x = Time, y = HandlingTime), inherit.aes = FALSE)

HandlingTimeAvgPopulations <- ggplot(data = AvgTrajectory, aes(x = Time, y = HandlingTime, color = ExtinctFactor)) + geom_line(size = 1) +
  theme_cowplot() + scale_color_manual(values = c('Extinct' = 'magenta', 'Extant' = 'green'), name = 'Extinct?') + ylab('Average \n Handling Time') +
  geom_line(data = QG_data, aes(x = Time, y = HandlingTime), inherit.aes = FALSE)

# which populations went extinct?

which_extinct <- inddata %>% group_by(Pop) %>% filter(Extinct == 1)

which_extinct <- unique(which_extinct$Pop)

rand_extinct <- sample(x = which_extinct, size = 10)

extinct_ex_data <- inddata %>% filter(Pop %in% rand_extinct)

### same for extant populations

which_extant <- inddata %>% group_by(Pop) %>% filter(Extinct == 0)

which_extant <- unique(which_extant$Pop)

rand_extant <- sample(x = which_extant, size = 10)

extant_ex_data <- inddata %>% filter(Pop %in% rand_extant)

average_extinct_traj <- extinct_ex_data %>% group_by(Time) %>% summarize(SCRate = mean(SCRate, na.rm = TRUE), HandlingTime = mean(HandlingTime, na.rm = TRUE)) 

average_extant_traj <- extant_ex_data %>% group_by(Time) %>% summarize(SCRate = mean(SCRate, na.rm = TRUE), HandlingTime = mean(HandlingTime, na.rm = TRUE))

Trajectory_Zoomed <- ExtinctionHeatMap + new_scale_color() + geom_line(data = extant_ex_data, aes(x = HandlingTime, y = SCRate, group = Pop), color = 'green', alpha = 0.4) +
  geom_line(data = extinct_ex_data, aes(x = HandlingTime, y = SCRate, group = Pop), color = 'magenta', alpha = 0.4) + 
  geom_line(data = average_extinct_traj, aes(x = HandlingTime, y = SCRate), color = 'magenta', size = 1) + 
  geom_line(data = average_extant_traj, aes(x = HandlingTime, y = SCRate), color = 'green', size = 1) +
  xlim(c(0, 0.3)) + ylim(c(0,0.1)) +  geom_line(data = QG_data, aes(x =HandlingTime, y = SpaceClearanceRate), color = 'white', 
                                                arrow = arrow(length = unit(0.03, "npc"), ends = 'first'))

### Put plots together

Plots_Together <- plot_grid(AverageTrajectoryPlot, Trajectory_Zoomed, SCRatePopulations, SCRateAvgPopulations,
          HandlingTimePopulations, HandlingTimeAvgPopulations, ncol = 2, nrow = 3, rel_heights = c(1,0.8,0.8),
          align = 'hv', labels = 'AUTO')

Title <- ggplot() + labs(title = 'K = 200') + theme(plot.title = element_text(hjust = 0.5, size = 18))

Plots_Together <- plot_grid(Title, Plots_Together, nrow = 2, ncol = 1, rel_heights = c(0.05,1))


save_plot(filename = 'Figure_4_LowK.png', plot = Plots_Together, nrow = 2.5, ncol = 2)

### plot of space clearance rate and handling time variance

traitdata <- read.csv('mediantraitdata.csv')

SCRVarPlot <- ggplot(data = traitdata, aes(x = time, y = medavar)) + geom_line(size = 1) + geom_ribbon(aes(ymin = loweravar, ymax = upperavar), alpha = 0.3) + 
  theme_cowplot() + xlab('Time') + ylab('Space Clearance Rate Variance')

save_plot(filename = 'SCRVarPlot.png', plot = SCRVarPlot)

HandlingVarPlot <- ggplot(data = traitdata, aes(x = time, y = medhvar)) + geom_line(size = 1) + geom_ribbon(aes(ymin = lowerhvar, ymax = upperhvar), alpha = 0.3) + 
  theme_cowplot() + xlab("Time") + ylab("Handling Time Variance")

save_plot(filename = 'HandlingVariationPlot.png', HandlingVarPlot)

### put plots together

VarPlots <- plot_grid(SCRVarPlot, HandlingVarPlot, nrow = 1, ncol = 2)

save_plot(filename = 'VarPlots_LowK.png', VarPlots, ncol = 2, nrow = 1, bg = 'white')




