# Clear environment
rm(list = ls())

# Load required libraries
library(tidyr)
library(tibble)
library(rbi.helpers)
library(ggplot2)
library(patchwork)
library(here)
library(rriskDistributions)
library(ggplot2)
library(PropCIs)
library(bayesplot)
library(binom)
library(beepr)
library(dplyr)
library(ggforce)
library(ggimage)
library(ggpubr)

# Set working directory (to this file's location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
dir.create(paste(getwd(),"Output",sep = "/"))

source("./InfDisNTI_DataPrep.R")
libbi_file <- "./InfDisNTI_Model50Clin.bi"

model <- rbi::bi_model(file = libbi_file)

model <- fix(model)

initial_fit <- rbi::sample(model, target = "posterior", nsamples = 10000, nparticles = 1,
                           input = input, obs = obs, verbose = FALSE,
                           start_time = 0, end_time = 10, noutputs = 100)

adapted <- initial_fit %>%
    adapt_proposal(min = 0.2, max = 0.3, adapt = "both", max_iter = 25, truncate = TRUE, verbose = FALSE)

posterior <- adapted %>%
    sample(nsamples = 150000, verbose = FALSE)

#x <- bi_read(posterior)

post_out = summary(posterior,type = "state", quantiles = c(0.025, 0.975))
post_out[,2:8] = lapply(post_out[,2:8],as.numeric)

fitting_traj <- post_out[which(post_out$var %in% c("sub_over_clin_p",
                                                   "min_over_infect_p",
                                                   "duration") &
                                   post_out$time == 2),]

params = summary(posterior, type = "param", quantiles = c(0.025, 0.975))

params <- rbind(params, fitting_traj[,-2])

dic <- c("DIC",DIC(posterior), rep(0,5))
params <- rbind(params, dic)

duration = -2*log(2)/log(post_out[which(post_out$var == "I"),"Median"])[21]
-2*log(2)/log(post_out[which(post_out$var == "I"),"2.5%"])[21]
-2*log(2)/log(post_out[which(post_out$var == "I"),"97.5%"])[21]
source("InfDisNTI_ResultsPlots.R")
traces = get_traces(posterior)

trace_plot = mcmc_trace(traces)
corr_plot = mcmc_pairs(traces, diag_fun = "dens",
                       off_diag_fun = "hex",
                       off_diag_args = list(size = 0.5, alpha = 0.2, shape = 23))
hist_plot = mcmc_hist(traces)

fit_plots = total_fit_plots(data, post_out)

#write.csv(x, paste(filename,"biread.csv",sep = "/"))
write.csv(post_out,file="./Output/posterior-state.csv")
write.csv(params,file="./Output/posterior-params.csv",)
write.csv(traces,file="./Output/param_traces.csv")
png(filename = "./Output/fit-plots.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
fit_plots
dev.off()
png(filename = "./Output/trace-plots.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
trace_plot
dev.off()
png(filename = "./Output/corr-plots.png",
    width = 2900/3,
    height = 1800/3,
    units = "mm",
    res = 600)
corr_plot
dev.off()
png(filename = "./Output/hist-plots.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
hist_plot
dev.off()

param_round <- params
param_round[,2] <- round(as.numeric(params[,2]),2)
param_round[,3] <- round(as.numeric(params[,3]),2)
param_round[,4] <- round(as.numeric(params[,4]),2)
param_round[,5] <- round(as.numeric(params[,5]),2)
param_round[,6] <- round(as.numeric(params[,6]),2)
param_round[,7] <- round(as.numeric(params[,7]),2)

write.csv(param_round, file="./Output/params_round.csv")

## INCIDENCE FOLLOWING INFECTION ##

colnames(post_out) <- c("var","time","min","lower","median","mean","upper","max")

fitmin  <- post_out[ which(post_out$var=='Z_min' & post_out$time %in% c(1,2,3,4,5,6,7,8,9,10)), ]
fitsub  <- post_out[ which(post_out$var=='Z_sub' & post_out$time %in% c(1,2,3,4,5,6,7,8,9,10)), ]
fitclin <- post_out[ which(post_out$var=='Z_clin'& post_out$time %in% c(1,2,3,4,5,6,7,8,9,10)), ]
write.csv(fitmin ,file="./Output/fitmin.csv")
write.csv(fitsub ,file="./Output/fitsub.csv")
write.csv(fitclin,file="./Output/fitclin.csv")

cummin  <- cuminc(fitmin)
cumsub  <- cuminc(fitsub)
cumclin <- cuminc(fitclin)
write.csv(cummin ,file="./Output/cummin.csv")
write.csv(cumsub ,file="./Output/cumsub.csv")
write.csv(cumclin,file="./Output/cumclin.csv")

# Incidence of minimal disease

mintargets <- merge(Minobs,Mininput,by=c("time","series"))
mintargets$value <- mintargets$value.x/mintargets$value.y
mintargets$lower <- as.double(unlist(binom.confint(mintargets$value.x,mintargets$value.y,conf.level = 0.95, method = "exact")["lower"]))
mintargets$upper <- as.double(unlist(binom.confint(mintargets$value.x,mintargets$value.y,conf.level = 0.95, method = "exact")["upper"]))
for (i in 1:nrow(mintargets)){
  if (mintargets$series[i] == "daniels") mintargets$time[i] <- mintargets$time[i] - 0.05
  if (mintargets$series[i] == "madsen")  mintargets$time[i] <- mintargets$time[i] + 0.05
}
write.csv(mintargets,file="./Output/mintargets.csv")

incminplot <- incplotmin(fitmin)
png(paste("./Output/IncMinPlot.png",sep=''),width=1200, height=800)
print(incminplot)
dev.off()
inccumminplot <- inccumplot(cummin)
inccumminplot <- inccumminplot + 
  scale_y_continuous(str_wrap("Cumulative incidence of minimal TB",width=24),expand=c(0, 0),lim=c(0,0.12))
png(paste("./Output/IncMinPlotCum.png",sep=''),width=1200, height=800)
print(inccumminplot)
dev.off()

# Incidence of subclinical disease

subtargets <- merge(Subobs,Subinput,by=c("time","series"))
subtargets$value <- subtargets$value.x/subtargets$value.y
subtargets$lower <- as.double(unlist(binom.confint(subtargets$value.x,subtargets$value.y,conf.level = 0.95, method = "exact")["lower"]))
subtargets$upper <- as.double(unlist(binom.confint(subtargets$value.x,subtargets$value.y,conf.level = 0.95, method = "exact")["upper"]))
for (i in 1:nrow(subtargets)){
  if (subtargets$series[i] == "NTI1") subtargets$time[i] <- subtargets$time[i] - 0.05
  if (subtargets$series[i] == "NTI3")  subtargets$time[i] <- subtargets$time[i] + 0.05
}
write.csv(subtargets,file="./Output/subtargets.csv")

incsubplot <- incplotsub(fitsub)
png(paste("./Output/IncSubPlot.png",sep=''),width=1200, height=800)
print(incsubplot)
dev.off()
inccumsubplot <- inccumplot(cumsub)
inccumsubplot <- inccumsubplot + 
  scale_y_continuous(str_wrap("Cumulative incidence of subclinical TB",width=24),expand=c(0, 0),lim=c(0,0.12))
png(paste("./Output/IncSubPlotCum.png",sep=''),width=1200, height=800)
print(inccumsubplot)
dev.off()

# Incidence of clinical disease

clintargets <- merge(Clinobs,Clininput,by=c("time","series"))
clintargets$value <- clintargets$value.x/clintargets$value.y
clintargets$lower <- as.double(unlist(binom.confint(clintargets$value.x,clintargets$value.y,conf.level = 0.95, method = "exact")["lower"]))
clintargets$upper <- as.double(unlist(binom.confint(clintargets$value.x,clintargets$value.y,conf.level = 0.95, method = "exact")["upper"]))
for (i in 1:nrow(clintargets)){
  if (clintargets$series[i] == "NTI1") clintargets$time[i] <- clintargets$time[i] - 0.05
  if (clintargets$series[i] == "NTI3") clintargets$time[i] <- clintargets$time[i] + 0.05
}
write.csv(clintargets,file="./Output/clintargets.csv")

incclinplot <- incplotclin(fitclin)
png(paste("./Output/IncClinPlot.png",sep=''),width=1200, height=800)
print(incclinplot)
dev.off()
inccumclinplot <- inccumplot(cumclin)
inccumclinplot <- inccumclinplot + 
  scale_y_continuous(str_wrap("Cumulative incidence of clinical TB",width=24),expand=c(0, 0),lim=c(0,0.12))
png(paste("./Output/IncClinPlotCum.png",sep=''),width=1200, height=800)
print(inccumclinplot)
dev.off()
beep()
