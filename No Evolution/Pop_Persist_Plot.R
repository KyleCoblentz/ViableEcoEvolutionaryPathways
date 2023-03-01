################################################################################
### Population persistence graph
################################################################################

### load packages

library(tidyverse); library(cowplot); library(ggnewscale); library(viridis);

### load data

HighIntData <- read.csv('ProportionExt_HighIntComp.csv')

ggplot(data = HighIntData, aes(x = HandlingTime, y = AttackRate, fill = NumberPersist)) + geom_tile() + 
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



f_high <- feasibility(h = seq(0.01,0.5, length.out = 1000), q = q, e = e, m =m)

s_high <- stability(h = seq(0.01, 0.5, length.out = 1000), q = q, e = e, m = m)

### put them all together

HeatMap_PopPersist_HighInt <- ggplot(data = HighIntData, aes(x = HandlingTime, y = AttackRate)) + geom_tile(aes(fill = NumberPersist)) +
  theme_cowplot() + geom_line(data = f_high, aes(x = h, y = a), color = 'white')  +
  geom_line(data = s_high, aes(x = h, y = a), color = 'white') + xlab('Handling Time') +
  ylab('Space Clearance Rate') + labs(fill = 'Percent Persist') +
  ylim(c(0,0.15)) + xlim(c(0, 0.5)) + theme_cowplot() +
  annotate(geom = 'text', x = 0.225, y = 0.125, color = 'white', label = 'Unstable') +
  annotate(geom = 'text', x = 0.4, y = 0.015, color = 'white', label = 'Infeasible') 

#############################################################################################
### Qualitative stability results
#############################################################################################

### load data for stability on the a,h plane

stab_data <- read.csv('Stability_ahPlane.csv')

### change QualResult column

stab_data$QualResult <- ifelse(stab_data$QualResult == 'NoPosSolution', 'Infeasible',
                               ifelse(stab_data$QualResult == "False", "Stable",
                                      ifelse(stab_data$QualResult == "True", "Damped Oscillations", 
                                             "Unstable")))

### reorder the QualResult variable to change the order of the variables in the plot

stab_data$QualResult <- factor(stab_data$QualResult, levels = c('Stable', 'Damped Oscillations', 'Infeasible', 'Unstable'))

### make a plot showing the qualitative results

QualStabPlot <- ggplot(data = stab_data, aes(x = HandlingTime, y = SpaceClearanceRate, fill = QualResult)) + geom_tile() + 
  scale_fill_viridis(discrete = TRUE) + geom_line(data = f_high, aes(x = h, y = a), color = 'white', inherit.aes = FALSE)  + 
  geom_line(data = s_high, aes(x = h, y = a), color = 'white', inherit.aes = FALSE) + xlim(c(0, 0.5)) + ylim(c(0, 0.15)) +
  labs(fill = "Qualitative \nDynamics") + theme_cowplot() + ylab('Space Clearance Rate') + xlab('Handling Time')


###########################################################################################################
### Plot of maximum eigenvalues
###########################################################################################################


stab_data <- stab_data %>% rowwise() %>% mutate(MaxEigenvalue = max(Eigenvalue1, Eigenvalue2))

EigenvaluePlot <- ggplot(data = stab_data, aes(x = HandlingTime, y = SpaceClearanceRate, fill = MaxEigenvalue)) + geom_tile() + 
  scale_fill_viridis() + geom_line(data = f_high, aes(x = h, y = a), color = 'white', inherit.aes = FALSE)  + 
  geom_line(data = s_high, aes(x = h, y = a), color = 'white', inherit.aes = FALSE) + xlim(c(0, 0.5)) + ylim(c(0, 0.15)) + 
  xlab('Handling Time') + ylab('Space Clearance Rate') + labs(fill = "Maximum \nEigenvalue") + theme_cowplot() 

############################################################################################
### plot of minimum population sizes
############################################################################################

ExtData <- HighIntData %>% rowwise() %>% mutate(MinPop = min(MinPrey, MinPred))

MinPopPlot <- ggplot(data = ExtData, aes(x = HandlingTime, y = AttackRate)) + geom_tile(aes(fill = MinPop)) + scale_fill_viridis() + theme_cowplot() + 
  xlab('Handling Time') + ylab('Space Clearance Rate') + labs(fill = 'Minimum \nPopulation \nSize') 


MinPopPlot <-  MinPopPlot + geom_line(data = f_high, aes(x = h, y = a), color = 'white')  +
  geom_line(data = s_high, aes(x = h, y = a), color = 'white') + xlim(c(0, 0.5)) + ylim(c(0, 0.15))



############################################################################################
### put all of the plots together
############################################################################################

together_plot <- plot_grid(HeatMap_PopPersist_HighInt, QualStabPlot, EigenvaluePlot, MinPopPlot,
                           ncol = 2, nrow = 2, labels = 'AUTO', align = 'hv')

title <- ggplot() + labs(title = 'K = 200') + theme(plot.title = element_text(hjust = 0.5, size = 18))

together_plot <- plot_grid(title, together_plot, rel_heights = c(0.05, 1), nrow = 2)

save_plot('Figure2.png', together_plot, nrow = 2, ncol = 2)

