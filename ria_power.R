################################################################################
#                   Author: Joshua Thompson
#   O__  ----       Email:  joshuajamesdavidthompson@gmail.com
#  c/ /'_ ---
# (*) \(*) --
# ======================== Script  Information =================================
# PURPOSE: Power analysis for no-mow study 
#
# PROJECT INFORMATION:
#   Name:  Power analysis for no-mow study 
#
# HISTORY:----
#   Date		        Remarks
#	-----------	   ---------------------------------------------------------------
#	 04/26/2023   Created script                                   JThompson (JT)
#===============================  Environment Setup  ===========================
#==========================================================================================

library(dplyr)
#library(truncnorm)
library(tidyverse)

# RIA Funtcion
RIA = function(df, control, impact, before, after, simName, iterations) {
  
  # df should be a data frame of observations 
  # control and impact are columns of the dependent variable of interest, one control (control) and one impact (impact)
  # before and after are the names which identify the periods - these should be stored in a column with the field name 'period'
  # par name is the field name of the column that contains the name of the dependent variable
  # iterations are how many random samples you want to conduct
  df <- df %>% mutate(Diff = log10(1+ pull(df[,impact])) - log10(1 + pull(df[,control])))
  
  # Real difference
  before.df = df[df$period == before,]
  after.df = df[df$period == after,]
  
  beforeMean = mean(before.df$Diff, na.rm = TRUE)
  afterMean = mean(after.df$Diff, na.rm = TRUE)
  MeanDiffReal = abs(afterMean - beforeMean)
  
  
  # File for shuffling
  
  DiffShufMeans = matrix(nrow = iterations,ncol = 1) #matrix for shuffled before and after means
  
  for(x in 1: iterations){
    df$RanDiffs = sample(df$Diff, size = length(df$Diff), replace = FALSE)
    beforeMeanRanDiffs = mean(df[df$period == before,]$RanDiffs, na.rm = TRUE)
    afterMeanRanDiffs = mean(df[df$period == after,]$RanDiffs, na.rm = TRUE)
    DiffMeanRanDiffs = abs(afterMeanRanDiffs - beforeMeanRanDiffs)
    DiffShufMeans[x,] = DiffMeanRanDiffs
  }
  DiffShufMeans
  
  
  pvalue = length(DiffShufMeans[DiffShufMeans[,1] > MeanDiffReal,])/iterations
  
  PValList = tibble(simName,pvalue)
  return(PValList)
}


simulate_data <- function(n_sites, n_samples, effect_size,site_variance,mean_conc,conc_variance) {
  # Generate site IDs
  site_ids <- rep(1:n_sites, each = n_samples)
  
  # Generate time period (before/after treatment)
  time_period <- rep(rep(c("Before", "After"), each = (n_samples/2)),n_sites)
  
  # Generate treatment (control/treatment)
  treatment <- c(rep("Control", (n_sites/2)*n_samples), rep("Treatment", (n_sites/2)*n_samples))
  
  n=1
  sims<-NULL
  for (i in 1:(n_sites/2)) {
    siteID <- n
    control_mean <- rnorm(1, mean = mean_conc, sd = sqrt(site_variance))
    treatment_mean <- control_mean * (1 - effect_size)
    control_load <- tibble("site_ids" =siteID,pollutant_load=rnorm(n_samples, mean = control_mean, sd = sqrt(conc_variance)),"week_no"=1:n_samples)
    beforetreatment_load <- rnorm(n_samples/2, mean = control_mean, sd = sqrt(conc_variance))
    aftertreatment_load <- rnorm(n_samples/2, mean = treatment_mean, sd = sqrt(conc_variance))
    treatment_load <- tibble("site_ids" =siteID+5,"pollutant_load"=c(beforetreatment_load,aftertreatment_load),"week_no"=1:n_samples)
    sims <- bind_rows(sims,bind_rows(control_load,treatment_load)) %>% arrange(site_ids)
    n<-n+1
  }
  
  df <- tibble(site_ids,time_period,treatment,sims %>% select(-c(site_ids)))
  minload<- min(df$pollutant_load)
  if(minload<0){
    df$pollutant_load <- ifelse(df$pollutant_load<0,df$pollutant_load+abs(minload),df$pollutant_load)
  }
  return(df)
}
  

l.n_sites= c(2,4,6,8,10,12,14,16,18,20)
l.n_samples=c(1*52*2,2*52*2,3*52*2)
l.effect_size=c(0.1,0.2,0.3,0.4,0.5)
l.site_variance=3
l.mean_conc=5
l.conc_variance=2
ria_res = NULL
it=0
tot_it=length(l.n_sites)*length(l.n_samples)*length(l.effect_size)
simulated_data_store <- NULL
for(nsi in l.n_sites){
  for(nsa in l.n_samples){
    for(efsa in l.effect_size){
      ria_compare <- NULL
      for(i in 1:500){
        simulated_data <- simulate_data(n_sites=nsi, n_samples=nsa, effect_size=efsa,site_variance=l.site_variance,mean_conc=l.mean_conc,conc_variance=l.conc_variance)
        simulated_data_store <- bind_rows(simulated_data_store,simulated_data)
        av.Sims <- simulated_data %>% group_by(treatment,week_no) %>% 
          summarize(time_period= unique(time_period),
                    mean_pollutant_load=mean(pollutant_load)) %>%
          tidyr::pivot_wider(names_from=c(treatment), values_from=c(mean_pollutant_load)) %>% 
          select(-c(week_no)) %>% 
          rename(period=time_period)
        ria_compare <- bind_rows(ria_compare,RIA(av.Sims, control="Control", impact="Treatment", before="Before", after="After", simName=i, iterations=100)) 
      } 
      it = it+1
      ria_res <- bind_rows(ria_res,bind_cols("simPower"=pull(ria_compare %>% summarize(sigSim=sum(pvalue<=0.05)))/pull(ria_compare %>% summarize(n=n())),
                           "numSites" = nsi,"numsampledYears"=(nsa/52)/2,"reductionPerc"=efsa,"SiteVariance_mgL"=l.site_variance,"meanConc_mgL"=l.mean_conc,"concVariance_mg"=l.conc_variance))
    message(paste0("Finished iteration ",it+1," of a total ",tot_it," simulations."))
      }
  }
}


#ria_res

ggplot(ria_res %>% mutate(Title=paste0("Expected Load Reduction: ",reductionPerc*100,"%")),
       aes(x= numsampledYears, y=simPower*100, fill = numSites, colour = numSites), size = 5) + 
  facet_wrap(~Title,ncol = 2) + 
  labs(fill='No. of Sites',col='No. of Sites') +
  scale_fill_gradient(low = "red", high = "blue", na.value = NA)+
  scale_color_gradient(low = "red", high = "blue", na.value = NA)+
  geom_hline(yintercept=80, linetype="dashed", color = "black",size=1.2)+
  geom_point()+ xlab("Number of Before/After Years") + ylab("Statistical Power %")+
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5,colour="black"),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white',face="bold",size = 8),
        panel.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                          colour = "lightgrey"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "lightgrey"),
        panel.border=element_rect(colour="black",size=1,fill=NA),
        axis.text.y=element_text(colour="black"),
        legend.background = element_rect(colour = NA, fill = 'white', linetype='solid'),
        legend.spacing.x = unit(0, "pt"),
        legend.margin=margin(t=-0.5,l=0.05,b=0.05,r=0.05, unit='cm'))
ggsave("riaPower.png",height = 7*0.8,width = 7*0.8)



