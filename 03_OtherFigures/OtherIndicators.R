###############

library(ggplot2)
library(readxl)
library(stringr)  

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

indicators <- read_excel("./OtherIndicators.xlsx")

my_colors <- c("#CB4C15","#004550")

indicators$type = factor(indicators$type, levels = c("Posterior","Prior") )
indicators$meas = factor(indicators$meas, levels = c("subtoclin","mintoinfec","mortality","duration") )

plothoriz <-    ggplot(indicators) +
  geom_point(aes(x=median,y=meas,color=type),size=6,position=position_dodge(.5)) +
  geom_errorbar(aes(y=meas,xmin=min,xmax=max,color=type),width=.4,size=1.5,position=position_dodge(.5)) +
  scale_y_discrete("",labels=c("Subclinical-to-Clinical Prevalence Ratio",
                               "Minimal-to-Infectious Prevalence Ratio",
                               "Duration of Infectious Disease"                               )) +
                               #str_wrap("Subclinical-to-Clinical Prevalence Ratio",width=16))) + 
  scale_x_continuous("Median (95% CI or CrI)",expand=c(0, 0),lim=c(0,4.5)) +
  theme_classic(base_size=35) +
  scale_color_manual(values = c("#CB4C15","#004550"),name=NULL) +
  guides(color = guide_legend(reverse=TRUE)) +
  theme(legend.position="right")
png(paste("./Output/OtherIndicators.png"),width=2000,height=400)
print(plothoriz)
dev.off()

