# Clear environment
rm(list = ls())

# Load required libraries
library(reshape2)
library(ggplot2)
library(patchwork)
library(Rmisc)
library(ggbreak)
library(dplyr)
library(stringr)
library(binom)
library(networkD3)
library(htmlwidgets)
library(gridExtra)

# Set working directory (to this file's location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Source model and results functions
source(file="./InfDisCohort_IBM.R")
source(file="./InfDisCohort_ResultsPlots.R")

# Read parameter values from calibrated model
traces <- read.csv(file = "./param_traces.csv")
traces <- traces[,c(2:10)]
colnames(traces) <- c("mort","inf_clr","inf_min","inf_sub","min_rec","min_sub","sub_min","sub_clin","clin_sub")

# Set cohort parameter values 
nreps   <- 1000            # Number of repetitions
nsteps  <- 121             # Number of time steps (months)
nyears  <- (nsteps-1)/12   # Number of time steps (years)
npeople <- 10000           # Number of people in each cohort
treatment <- 0             # Treatment on/off
threshold_months  <- 9     # Threshold duration to define transitional state
threshold_changes <- 3     # Threshold transitions to define transitional state

## Run cohort simulation 
#for (i in 1:nreps){
#  IBM_temp <- full_simulation_inf(nsteps, npeople, 1, traces, "random", treatment)
#  save(IBM_temp,file=paste("./IBM_i/IBM_",i,".Rdata",sep=""))
#}
## Subset to only those who develop any form of disease
#for (i in 1:nreps){
#  disease <- as.data.frame(array(0,dim=c(nsteps,0)))
#  load(paste("./IBM_i/IBM_",i,".Rdata",sep=""))
#  for (j in 1:npeople){
#    if (length(which(IBM_temp[,j] %in% c('m','s'))) > 0)
#      disease <- cbind(disease,IBM_temp[,j])
#  }
#  save(disease,file=paste("./disease_i/disease_",i,".Rdata",sep=""))
#  rm("IBM_temp","disease")
#}

### CLEARANCE FOLLOWING INFECTION ###

# Calculate proportion who cleared within 10 years
clr <- array(0,dim=c(nreps,2))
for (i in 1:nreps){
  load(paste("./disease_i/disease_",i,".Rdata",sep=""))
  clr[i,1] <- npeople-ncol(disease)
  clr[i,2] <- (npeople-ncol(disease))/npeople
}
apply(clr,2,median)
apply(clr,2,quantile,probs=c(0.025,0.975))

# Calculate proportion who cleared within two years
xtwo <- array(0,dim=c(nreps,2))
for (i in 1:nreps){
  load(paste("./IBM_i/IBM_",i,".Rdata",sep=""))
  IBM_temp <- IBM_temp[1:25,]
  x <- 0
  for (j in 1:ncol(IBM_temp)){
    if (any(IBM_temp[,j] == 'x') == TRUE)  x <- x+1
  }
  xtwo[i,1] <- x
  xtwo[i,2] <- x/ncol(IBM_temp)
}
apply(xtwo,2,median)
apply(xtwo,2,quantile,probs=c(0.025,0.975))

# Calculate proportion who developed disease in ten years
disten <- array(0,dim=c(nreps,2))
for (i in 1:nreps){
  load(paste("./disease_i/disease_",i,".Rdata",sep=""))
  x <- 0
  for (j in 1:ncol(disease)){
    if (any(disease[,j] %in% c('m','s')) == TRUE)  x <- x+1
  }
  disten[i,1] <- x
  disten[i,2] <- x/npeople
}
apply(disten,2,median)
apply(disten,2,quantile,probs=c(0.025,0.975))

# Calculate proportion who developed disease in two years
distwo <- array(0,dim=c(nreps,2))
for (i in 1:nreps){
  load(paste("./disease_i/disease_",i,".Rdata",sep=""))
  disease <- disease[1:25,]
  x <- 0
  for (j in 1:ncol(disease)){
    if (any(disease[,j] %in% c('m','s')) == TRUE)  x <- x+1
  }
  distwo[i,1] <- x
  distwo[i,2] <- x/nreps
}
apply(distwo,2,median)
apply(distwo,2,quantile,probs=c(0.025,0.975))

### INCIDENCE FOLLOWING INFECTION ###

# Incidence of minimal disease
incmin       <- incidence(nreps, nyears, nsteps, npeople, 1, traces, "random", "min")
save(incmin,file="./Output/incmin.Rdata")
incminsum    <- incsum(incmin)
write.csv(incminsum,file="./Output/incmin.csv")
inccummin    <- inccum(incmin, nyears, nreps)
inccumminsum <- incsum(inccummin)
write.csv(inccumminsum,file="./Output/cummin.csv")
incminplot   <- incplotmin(incminsum[-1,])
png(paste("./Output/IncMinPlot.png",sep=''),width=1200, height=800)
print(incminplot)
dev.off()
inccumminplot <- inccumplot(inccumminsum) + 
  scale_y_continuous(str_wrap("Cumulative incidence of minimal TB",width=24),expand=c(0, 0),lim=c(0,0.11))
png(paste("./Output/IncMinPlotCum.png",sep=''),width=1200, height=800)
print(inccumminplot)
dev.off()

# Incidence of subclinical disease
incsub     <- incidence(nreps, nyears, nsteps, npeople, 1, traces, "random", "sub")
save(incsub,file="./Output/incsub.Rdata")
incsubsum  <- incsum(incsub)
write.csv(incsubsum,file="./Output/incsub.csv")
inccumsub    <- inccum(incsub, nyears, nreps)
inccumsubsum <- incsum(inccumsub)
write.csv(inccumsubsum,file="./Output/cumsub.csv")
incsubplot <- incplotsub(incsubsum[-1,])
png(paste("./Output/IncSubPlot.png",sep=''),width=1200, height=800)
print(incsubplot)
dev.off()
inccumsubplot <- inccumplot(inccumsubsum)
inccumsubplot <- inccumsubplot + 
  scale_y_continuous(str_wrap("Cumulative incidence of subclinical TB",width=24),expand=c(0, 0),lim=c(0,0.11))
png(paste("./Output/IncSubPlotCum.png",sep=''),width=1200, height=800)
print(inccumsubplot)
dev.off()

# Incidence of clinical disease
incclin     <- incidence(nreps, nyears, nsteps, npeople, 1, traces,"random","clin")
save(incclin,file="./Output/incclin.Rdata")
incclinsum  <- incsum(incclin)
write.csv(incclinsum,file="./Output/incclin.csv")
inccumclin    <- inccum(incclin, nyears, nreps)
inccumclinsum <- incsum(inccumclin)
write.csv(inccumclinsum,file="./Output/cumclin.csv")
incclinplot <- incplotclin(incclinsum[-1,])
png(paste("./Output/IncClinPlot.png",sep=''),width=1200, height=800)
print(incclinplot)
dev.off()
inccumclinplot <- inccumplot(inccumclinsum)
inccumclinplot <- inccumclinplot + 
  scale_y_continuous(str_wrap("Cumulative incidence of clinical TB",width=24),expand=c(0, 0),lim=c(0,0.11))
png(paste("./Output/IncClinPlotCum.png",sep=''),width=1200, height=800)
print(inccumclinplot)
dev.off()

#### PATHWAYS FOLLOWING INFECTION ###
#
#states <- array(0,dim=c(nreps,8))
#for (i in 1:nreps){
#  load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#  tm <- 0
#  ts <- 0
#  ns <- 0
#  s <- 0
#  nc <- 0
#  c <- 0
#  d <- 0
#  r <- 0
#  for (j in 1:ncol(disease)){
#    str <- paste(disease[,j],collapse="")
#    str <- str_remove_all(str,"-")
#    im <- str_count(str,"im")
#    is <- str_count(str,"is")
#    mr <- str_count(str,"mr")
#    ms <- str_count(str,"ms")
#    sm <- str_count(str,"sm")
#    sc <- str_count(str,"sc")
#    cs <- str_count(str,"cs")
#    cd <- str_count(str,"d")
#    if (im > 0)                   tm <- tm + 1
#    if (is > 0)                   ts <- ts + 1
#    if (im > 0 && ms == 0)        ns <- ns + 1
#    if ((is+ms) > 0)               s <- s + 1
#    if ((is+ms) > 0 && sc == 0)   nc <- nc + 1
#    if (sc > 0)                    c <- c + 1
#    if (cd > 0)                    d <- d + 1
#    if (mr > 0)                    r <- r + 1
#  }
#  disnum <- ncol(disease)
#  states[i,1] <- tm/disnum
#  states[i,2] <- ts/disnum
#  states[i,3] <- ns/disnum
#  states[i,4] <-  s/disnum
#  states[i,5] <- nc/disnum
#  states[i,6] <-  c/disnum
#  states[i,7] <-  d/disnum
#  states[i,8] <-  r/disnum
#}
#apply(states,2,median)
#apply(states,2,quantile,probs=c(0.025,0.975))
#
#### PATHWAY TRAJECTORIES FOLLOWING INFECTION ###
#
#category <- array(0,dim=c(nreps,5))
#class <- array(0,dim=c(nreps,npeople))
#for (i in 1:nreps){
#  load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#  infect <- 0
#  for (j in 1:ncol(disease)){
#    str <- paste(disease[,j],collapse="")
#    str <- str_remove_all(str,"-")
#    is <- str_count(str,"is")
#    mr <- str_count(str,"mr")
#    ms <- str_count(str,"ms")
#    sm <- str_count(str,"sm")
#    sc <- str_count(str,"sc")
#    cs <- str_count(str,"cs")
#    s <- str_count(str,"s")
#    if (mr == 0 && sm == 0 && cs == 0 && sc == 0)      class[i,j] <- "PRO"
#    else if (mr == 0 && sm == 0 && cs == 0 && sc > 0)  class[i,j] <- "PRC"
#    else if ((is+ms) >= 2 || sc >= 2)                  class[i,j] <- "UND"
#    else if (ms <= 1 && sc <= 1 && sm <= 1 && cs <= 1) class[i,j] <- "REC"
#    if (s > 0) infect <- infect + 1
#  }
#  category[i,1] <-  sum(class[i,] == "REC")/sum(class[i,]>0)
#  category[i,2] <-  sum(class[i,] == "UND")/sum(class[i,]>0)
#  category[i,3] <-  sum(class[i,] == "UND")/infect
#  category[i,4] <- (sum(class[i,] == "PRO")+sum(class[i,] == "PRC"))/sum(class[i,]>0)
#  category[i,5] <-  sum(class[i,] == "PRC")/sum(class[i,]>0)
#}
#apply(category,2,median)
#apply(category,2,quantile,probs=c(0.025,0.975))
#
#transnum <- c()
#for (i in 1:nreps){
#  load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#  for (j in 1:ncol(disease)){
#    str <- paste(disease[,j],collapse="")
#    str <- str_remove_all(str,"-")
#    is <- str_count(str,"is")
#    ms <- str_count(str,"ms")
#    sc <- str_count(str,"sc")
#    if ((is+ms) >= 2 || sc >= 2){
#      sm <- str_count(str,"sm")
#      cs <- str_count(str,"cs")
#      transcount <- ms+sm+sc+cs
#      transnum <- c(transnum,transcount)
#    }
#  } 
#}
#median(transnum)
#quantile(transnum,probs=c(0.25,0.75,0.025,0.975,0,1))
#
### DOMINANT DISEASE STATES ##
#
## Summarise annual distributions across states
#traj_sum_inf <- array(0,dim=c(nreps,nyears+1,10))
#for (i in 1:nreps){
#  load(paste("./IBM_i/IBM_",i,".Rdata",sep=""))
#  traj_inf_temp <- report_trajectories_inf(IBM_temp, npeople, 10, threshold_months, threshold_changes)
#  traj_sum_inf[i,,]  <- summarise_trajectories_inf(traj_inf_temp, 10) 
#  rm("IBM_temp")
#}
#
## Calculate median annual distribution across states
#traj_sum_med <- array(0,dim=c(nyears+1,10))
#for (i in 1:nyears+1){
#  for (j in 1:10){
#    traj_sum_med[i,j] <- median(traj_sum_inf[,i,j])
#  }
#}
#
## Convert counts to percentages
#traj_pct_med <- array(0,dim=c(nyears+1,10))
#traj_pct_med[,1] <- c(0:10)
#traj_pct_med[,2:10]<- round((traj_sum_med[,-1]/rowSums(traj_sum_med[,-1]))*100,2)
#colnames(traj_pct_med) <- c("year","death","treat","clinical","transitional","subclinical","minimal","infection","recover","cleared")
#
## Plot median annual distribution across states
#traj_plot_inf <- plot_trajectories_inf(traj_pct_med, 10, "")
#png(paste("./Output/Plot_AnnualDistribution.png",sep=''),width=1200, height=800)
#print(traj_plot_inf)
#dev.off()
#
### SIMPLIFIED PATHWAYS ##
#
## Generate data for Sankey plot
#pathways <- array(data=0,dim=c(14,3))
#colnames(pathways) <- c("source","target","value")
#
#  # Progression from infection
#  infec <- array(0,dim=c(nreps,2))
#  for (i in 1:nreps){
#    load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#    tm <- 0
#    ts <- 0
#    for (j in 1:ncol(disease)){
#      str <- paste(disease[,j],collapse="")
#      str <- str_remove_all(str,"-")
#      im <- str_count(str,"im")
#      is <- str_count(str,"is")
#      if (im > 0)   tm <- tm + 1
#      if (is > 0)   ts <- ts + 1
#    }
#    infec[i,1] <- tm
#    infec[i,2] <- ts
#  }
#  pathways[ 1,] <- c("infection","initmin",round(median(infec[,1]),digits=0))
#  pathways[ 2,] <- c("infection","initsub",round(median(infec[,2]),digits=0))
#
#  # End points for those who never progressed beyond minimal
#  minonly <- array(0,dim=c(nreps,2))
#  for (i in 1:nreps){
#    load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#    endmin <- 0
#    endrec <- 0
#    for (j in 1:length(disease)){
#      str <- paste(disease[,j],collapse="")
#      str <- str_remove_all(str,"-")
#      sub <- str_count(str,"s")
#      rec <- str_count(str,"r") 
#      if (sub == 0 && rec > 0)      endrec = endrec + 1
#      if (sub == 0 && rec == 0)     endmin = endmin + 1
#    }
#    minonly[i,1] <- endrec
#    minonly[i,2] <- endmin
#  }
#  pathways[ 3,] <- c("initmin","recovery",round(median(minonly[,1]),digits=0))
#  pathways[ 4,] <- c("initmin","endmin"  ,round(median(minonly[,2]),digits=0))
#
#  # Progression from infection to minimal to subclinical
#  infminsub <- array(0,dim=c(nreps,1))
#  for (i in 1:nreps){
#    load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#    minsub <- 0
#    for (j in 1:length(disease)){
#      str <- paste(disease[,j],collapse="")
#      str <- str_remove_all(str,"-")
#      im <- str_count(str,"im")
#      ms <- str_count(str,"ms") 
#      if (im > 0 && ms > 0)      minsub = minsub + 1
#    }
#    infminsub[i,1] <- minsub
#  }
#  pathways[ 5,] <- c("initmin","initsub",round(median(infminsub[,1]),digits=0))
#  
#  # End points for those who develop subclinical
#  subends <- array(0,dim=c(nreps,4))
#  for (i in 1:nreps){
#    load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#    subrec <- 0
#    submin <- 0
#    subsub <- 0
#    subcli <- 0
#    for (j in 1:length(disease)){
#      str <- paste(disease[,j],collapse="")
#      str <- str_remove_all(str,"-")
#      s <- str_count(str,"s")
#      c <- str_count(str,"c")
#      sr <- str_count(str_sub(str,-1,-1),"r")
#      sm <- str_count(str_sub(str,-1,-1),"m")
#      ss <- str_count(str_sub(str,-1,-1),"s")
#      if (s > 0 && c == 0 && sr > 0)      subrec = subrec + 1
#      if (s > 0 && c == 0 && sm > 0)      submin = submin + 1
#      if (s > 0 && c == 0 && ss > 0)      subsub = subsub + 1
#      if (s > 0 && c > 0)                 subcli = subcli + 1
#    }
#    subends[i,1] <- subrec
#    subends[i,2] <- submin
#    subends[i,3] <- subsub
#    subends[i,4] <- subcli
#  }
#  pathways[ 6,] <- c("initsub"  ,"recovery",round(median(subends[,1]),digits=0))
#  pathways[ 7,] <- c("initsub"  ,"endmin"  ,round(median(subends[,2]),digits=0))
#  pathways[ 8,] <- c("initsub"  ,"endsub"  ,round(median(subends[,3]),digits=0))
#  pathways[ 9,] <- c("initsub"  ,"clinical",round(median(subends[,4]),digits=0))
#
#  # End points for those who develop clinical
#  clinends <- array(0,dim=c(nreps,5))
#  for (i in 1:nreps){
#    load(paste("./disease_i/disease_",i,".Rdata",sep=""))
#    clirec   <- 0
#    climin   <- 0
#    clisub   <- 0
#    clicli   <- 0
#    clideath <- 0
#    for (j in 1:length(disease)){
#      str <- paste(disease[,j],collapse="")
#      str <- str_remove_all(str,"-")
#      c <- str_count(str,"c")
#      cr <- str_count(str_sub(str,-1,-1),"r")
#      cm <- str_count(str_sub(str,-1,-1),"m")
#      cs <- str_count(str_sub(str,-1,-1),"s")
#      cc <- str_count(str_sub(str,-1,-1),"c")
#      cd <- str_count(str_sub(str,-1,-1),"d")
#      if (c > 0 && cr > 0)      clirec   = clirec + 1
#      if (c > 0 && cm > 0)      climin   = climin + 1
#      if (c > 0 && cs > 0)      clisub   = clisub + 1
#      if (c > 0 && cc > 0)      clicli   = clicli + 1
#      if (c > 0 && cd > 0)      clideath = clideath + 1
#    }
#    clinends[i,1] <- clirec  
#    clinends[i,2] <- climin  
#    clinends[i,3] <- clisub  
#    clinends[i,4] <- clicli  
#    clinends[i,5] <- clideath
#  }
#  pathways[10,] <- c("clinical" ,"recovery",round(median(clinends[,1]),digits=0))
#  pathways[11,] <- c("clinical" ,"endmin"  ,round(median(clinends[,2]),digits=0))
#  pathways[12,] <- c("clinical" ,"endsub"  ,round(median(clinends[,3]),digits=0))
#  pathways[13,] <- c("clinical" ,"endclin" ,round(median(clinends[,4]),digits=0))
#  pathways[14,] <- c("clinical" ,"death"   ,round(median(clinends[,5]),digits=0))
#
#write.csv(as.data.frame(pathways)[],"./Output/sankeypathways.csv")
#
## Generate Sankey plot
#paths <- read.csv(file = "./Output/sankeypathways.csv")
#paths <- paths[,-1]
#nodes <- data.frame(
#  name=c(as.character(paths$source), 
#         as.character(paths$target)) %>% unique()
#)
#paths$IDsource <- match(paths$source,nodes$name)-1 
#paths$IDtarget <- match(paths$target,nodes$name)-1
#my_color <- 'd3.scaleOrdinal() .domain(["infection","initmin","initsub","clinical","recovery","endmin","endsub","endclin","death"]) .range(["#537D8D","#CBA715","#CB4C15","#AE0D0A","#8A8A8A","#CBA715","#CB4C15","#AE0D0A","#545454"])'
#
#sn <- sankeyNetwork(Links = paths, Nodes = nodes,
#                    Source = "IDsource", Target = "IDtarget",
#                    Value = "value", NodeID = "name", 
#                    nodeWidth = 40,colourScale = my_color,
#                    fontSize = 0, sinksRight = TRUE)
#
#sn # Plot can be manipulated and saved by hand#