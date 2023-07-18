rm(list = ls())

library(rriskDistributions)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

data = read.csv("included_studies_shrunk1.csv", stringsAsFactors = FALSE)

data_use <- data[which(data$d.of.y != 0),]
studies <- unique(data_use$record.id)
n_studies <- length(studies)
n_datapoints <- nrow(data_use)
n_people <- sum(data_use$d.of.y/data_use$repeats)
n_transitions <- 0
for (i in seq(1,n_studies)) {
    temp <- data_use[which(data_use$record.id == studies[i]),]
    trans <- unique(temp$Transition)
    for (j in seq(1,length(trans))) {
        temp2 <- temp[which(temp$Transition == trans[j]),]
        n_transitions <- n_transitions + temp2$n.of.y[nrow(temp2)]
    }
}
n_min_sub <- length(which(data_use$model.transition == "Min-Sub"))
n_min_clin <- length(which(data_use$model.transition == "Min-Clin"))
n_clin_min <- length(which(data_use$model.transition == "Clin-Min"))
n_min_inf <- length(which(data_use$model.transition == "Min-Inf"))
n_inf_min <- length(which(data_use$model.transition == "Inf-Min"))

data$times <- round(data$x.months/12,1)
data$d.of.y[which(data$Transition %in% c("Submin","Clinmin","Clinmin_clin","Submin_sub","InfminI","Infmin_infI"))] <- round(data$d.of.y[which(data$Transition %in% c("Submin","Clinmin","Clinmin_clin","Submin_sub","InfminI","Infmin_infI"))] * 3/4 )
data$n.of.y_weight <- round(data$n.of.y / data$repeats)
data$d.of.y_weight <- round(data$d.of.y / data$repeats)

################################################################################
# DEFINE DISTRIBUTION BASED ON POULSEN 1957

# Based on cumulative proportions from Poulsen 1957 (Table 52)
#get.cauchy.par(p=c(0.17,0.52,0.70,0.87,0.93,0.97,0.98,1.00),q=c(1,2,3,4,5,6,8,10))
#dcauchy(c(1,2,3,4,5,6,7,8,9,10),2.0600957,0.8282028)

# Define function for left- and right-truncated gamma distribution
rcauchyt <- function(n, range, location, scale) {
    F.a <- pcauchy(min(range), location=location, scale=scale)
    F.b <- pcauchy(max(range), location=location, scale=scale)
    u <- runif(n, min=F.a, max=F.b)
    qcauchy(u, location=location, scale=scale)
}

################################################################################
# PROCESS DATA FOR CALIBRATION

sim_num <- 1000 

## Daniels (infection to minimal) ##

# Reported case numbers:           # Adjusted case numbers (reduced by 25%)
# 19 disease in  1-12 months       # 14 disease in  1-12 months
#  4 disease in 13-24 months       #  3 disease in 13-24 months
#  3 disease in 25-36 months       #  2 disease in 25-36 months
#  1 disease in 37-48 months       #  1 disease in 37-48 months
#  0 disease in 49-60 months       #  0 disease in 49-60 months

# Reported at risk:                # Adjusted at risk (add back in cases)
# 248 at risk in  1-12 months      # 248 at risk in  1-12 months
# 198 at risk in 13-24 months      # 217 at risk in 13-24 months
# 134 at risk in 25-36 months      # 138 at risk in 25-36 months
#  58 at risk in 37-48 months      #  61 at risk in 37-48 months
#  16 at risk in 49-60 months      #  17 at risk in 49-60 months

daniels <- array(0,dim=c(sim_num,5)) 

for (i in 1:sim_num){

    temp <- array(0,dim=c(20,5))
    
    # Select months from conversion to conversion detection
    temp[,1] <- runif(20,1,12)

    # Select months from conversion detection to disease onset
    for (j in 1:14){
        truncmin <- (temp[j,1]+0)/12
        truncmax <- (temp[j,1]+12)/12
        temp[j,2] <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }

    for (j in 15:17){
        truncmin <- (temp[j,1]+12)/12
        truncmax <- (temp[j,1]+24)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 18:19){
        truncmin <- (temp[j,1]+24)/12
        truncmax <- (temp[j,1]+36)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 20){
        truncmin <- (temp[j,1]+36)/12
        truncmax <- (temp[j,1]+48)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }

    # Convert months to years
    temp <- as.data.frame(temp)
    daniels[i,1] <- sum(temp$V2 <= 12)
    daniels[i,2] <- sum(temp$V2 > 12 & temp$V2 <= 24)
    daniels[i,3] <- sum(temp$V2 > 24 & temp$V2 <= 36)
    daniels[i,4] <- sum(temp$V2 > 36 & temp$V2 <= 48)
    daniels[i,5] <- sum(temp$V2 > 48 & temp$V2 <= 60)

}

# Select median cases and at risk population
daniels_summary <- array(0,dim=c(5,3))
daniels_summary[,1] <- c(1:5)
for (i in 1:5){
    daniels_summary[i,2] <- median(daniels[,i])
}
 daniels_summary[, 3] <- c(248,
                           217-daniels_summary[1,2],
                           138-daniels_summary[2,2],
                           61-daniels_summary[3,2],
                           17-daniels_summary[4,2])

# Convert to data frame
data_daniels <- as.data.frame(daniels_summary)
names(data_daniels)[names(data_daniels)=="V1"] <- "year"
names(data_daniels)[names(data_daniels)=="V2"] <- "cases"
names(data_daniels)[names(data_daniels)=="V3"] <- "atrisk"

write.csv(data_daniels,file="./Output/data_daniels_preweight.csv")

################################################################################

## Madsen (infection to minimal) ##

# Reported case numbers:           # Adjusted case numbers (reduced by 25%)
# 41 disease at     0 months       # 31 disease at     0 months
#  9 disease in  1-12 months       #  7 disease in  1-12 months
#  2 disease in 13-24 months       #  2 disease in 13-24 months
#  0 disease in 25-36 months       #  0 disease in 25-36 months
#  0 disease in 37-48 months       #  0 disease in 37-48 months
#  0 disease in 49-60 months       #  0 disease in 49-60 months

# Reported at risk:                # Adjusted at risk (add back in cases)
# 208 at risk at     0 months      # 208 at risk at     0 months
# 162 at risk in  1-12 months      # 203 at risk in  1-12 months
# 109 at risk in 13-24 months      # 118 at risk in 13-24 months
#  63 at risk in 25-36 months      #  65 at risk in 25-36 months
#  26 at risk in 37-48 months      #  26 at risk in 37-48 months
#  16 at risk in 49-60 months      #  16 at risk in 49-60 months

madsen <- array(0,dim=c(sim_num,6))

for (i in 1:sim_num){

    temp <- array(0,dim=c(40,6))

    # Select months from conversion to conversion detection
    temp[,1] <- runif(40,1,15)

    for (j in 1:31){
        truncmin <- 0
        truncmax <- (temp[j,1]+0)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 32:38){
        truncmin <- (temp[j,1]+0)/12
        truncmax <- (temp[j,1]+12)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 39:40){
        truncmin <- (temp[j,1]+12)/12
        truncmax <- (temp[j,1]+24)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }

    # Convert months to years
    temp <- as.data.frame(temp)
    madsen[i,1] <- sum(temp$V2 <= 12)
    madsen[i,2] <- sum(temp$V2 > 12 & temp$V2 <= 24)
    madsen[i,3] <- sum(temp$V2 > 24 & temp$V2 <= 36)
    madsen[i,4] <- sum(temp$V2 > 36 & temp$V2 <= 48)
    madsen[i,5] <- sum(temp$V2 > 48 & temp$V2 <= 60)
    madsen[i,6] <- sum(temp$V2 > 60 & temp$V2 <= 72)
}

# Select median cases and at risk pop
madsen_summary <- array(0,dim=c(6,3))
madsen_summary[,1] <- c(1:6)
for (i in 1:6){
    madsen_summary[i,2] <- median(madsen[,i])
}
 madsen_summary[, 3] <- c(208,
                          203-madsen_summary[1,2],
                          118-madsen_summary[2,2],
                          65-madsen_summary[3,2],
                          26-madsen_summary[4,2],
                          16-madsen_summary[5,2])

# Convert to data frame
data_madsen <- as.data.frame(madsen_summary)
names(data_madsen)[names(data_madsen)=="V1"] <- "year"
names(data_madsen)[names(data_madsen)=="V2"] <- "cases"
names(data_madsen)[names(data_madsen)=="V3"] <- "atrisk"

write.csv(data_madsen,file="./Output/data_madsen_preweight.csv")

################################################################################

## NTI-1 (infection to infectious disease) ##

# Infected between surveys 1 & 2. No report of loss to follow-up.

# Reported case numbers:
# 19 disease in  1-18 months
#  2 disease in 19-36 months
#  4 disease in 37-60 months

# Reported at risk:
# 674 at risk at 0 months

nti1 <- array(0,dim=c(sim_num,5))

for (i in 1:sim_num){

    temp <- array(0,dim=c(25,5))

    # Select months from conversion to conversion detection
    temp[,1] <- runif(25,1,18)

    # Select months from conversion detection to disease onset
    for (j in 1:19){
        truncmin <- 0
        truncmax <- (temp[j,1]+0)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 20:21){
        truncmin <- (temp[j,1]+0)/12
        truncmax <- (temp[j,1]+18)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 22:25){
        truncmin <- (temp[j,1]+18)/12
        truncmax <- (temp[j,1]+42)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }

    # Convert months to years
    temp <- as.data.frame(temp)
    nti1[i,1] <- sum(temp$V2 <= 12)
    nti1[i,2] <- sum(temp$V2 > 12 & temp$V2 <= 24)
    nti1[i,3] <- sum(temp$V2 > 24 & temp$V2 <= 36)
    nti1[i,4] <- sum(temp$V2 > 36 & temp$V2 <= 48)
    nti1[i,5] <- sum(temp$V2 > 48 & temp$V2 <= 60)
    
}

# Select median cases and at risk pop
nti1_summary <- array(0,dim=c(5,3))
nti1_summary[,1] <- c(1:5)
for (i in 1:5){
    nti1_summary[i,2] <- median(nti1[,i])
}
nti1_summary[,3] <- c(674,
                      674-nti1_summary[1,2],
                      674-(nti1_summary[1,2]+nti1_summary[2,2]),
                      674-(nti1_summary[1,2]+nti1_summary[2,2]+nti1_summary[3,2]),
                      674-(nti1_summary[1,2]+nti1_summary[2,2]+nti1_summary[3,2]+nti1_summary[4,2]))

# Convert to data frame
data_nti1 <- as.data.frame(nti1_summary)
names(data_nti1)[names(data_nti1)=="V1"] <- "year"
names(data_nti1)[names(data_nti1)=="V2"] <- "cases"
names(data_nti1)[names(data_nti1)=="V3"] <- "atrisk"

write.csv(data_nti1,file="./Output/data_nti1_preweight.csv")

################################################################################

## NTI-2 (infection to infectious disease) ##

# Infected between surveys 2 & 3. No report of loss to follow-up.

# Reported case numbers:
#  7 disease in  1-18 months
#  2 disease in 19-42 months

# Reported at risk:
# 451 at risk at 0 months

nti2 <- array(0,dim=c(sim_num,4))

for (i in 1:sim_num){

    temp <- array(0,dim=c(9,4))

    # Select months from conversion to conversion detection
    temp[,1] <- runif(9,1,18)

    for (j in 1:7){
        truncmin <- 0
        truncmax <- (temp[j,1]+0)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }
    for (j in 8:9){
        truncmin <- (temp[j,1]+0)/12
        truncmax <- (temp[j,1]+24)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }

    # Convert months to years
    temp <- as.data.frame(temp)

    nti2[i,1] <- sum(temp$V2 <= 12)
    nti2[i,2] <- sum(temp$V2 > 12 & temp$V2 <= 24)
    nti2[i,3] <- sum(temp$V2 > 24 & temp$V2 <= 36)
    nti2[i,4] <- sum(temp$V2 > 36 & temp$V2 <= 48)

}

# Select median cases and at risk pop
nti2_summary <- array(0,dim=c(4,3))
nti2_summary[,1] <- c(1:4)
for (i in 1:4){
    nti2_summary[i,2] <- median(nti2[,i])
}
nti2_summary[,3] <- c(451,
                     451-nti2_summary[1,2],
                     451-(nti2_summary[1,2]+nti2_summary[2,2]),
                     451-(nti2_summary[1,2]+nti2_summary[2,2]+nti2_summary[3,2]))

# Convert to data frame
data_nti2 <- as.data.frame(nti2_summary)
names(data_nti2)[names(data_nti2)=="V1"] <- "year"
names(data_nti2)[names(data_nti2)=="V2"] <- "cases"
names(data_nti2)[names(data_nti2)=="V3"] <- "atrisk"

write.csv(data_nti2,file="./Output/data_nti2_preweight.csv")

################################################################################

## NTI-3 (infection to infectious disease) ##

# Infected between surveys 2 & 3. No report of loss to follow-up.

# Reported case numbers:
#  8 disease in  1-24 months

# Reported at risk:
# 413 at risk at 0 months

nti3 <- array(0,dim=c(sim_num,6))

for (i in 1:sim_num){

    temp <- array(0,dim=c(8,6))

    # Select months from conversion to conversion detection
    temp[,1] <- runif(8,1,24)

    for (j in 1:8){
        truncmin <- 0
        truncmax <- (temp[j,1]+0)/12
        temp[j,2]  <- (rcauchyt(1,c(truncmin,truncmax),2.0600957,0.8282028))*12
    }

    # Convert months to years
    temp <- as.data.frame(temp)

    nti3[i,1] <- sum(temp$V2 <= 12)
    nti3[i,2] <- sum(temp$V2 > 12 & temp$V2 <= 24)

}

# Select median cases and at risk pop
nti3_summary <- array(0,dim=c(2,3))
nti3_summary[,1] <- c(1:2)
for (i in 1:2){
    nti3_summary[i,2] <- median(nti3[,i])
}
nti3_summary[,3] <- c(413,
                      413-nti3_summary[1,2])

# Convert to data frame
data_nti3 <- as.data.frame(nti3_summary)
names(data_nti3)[names(data_nti3)=="V1"] <- "year"
names(data_nti3)[names(data_nti3)=="V2"] <- "cases"
names(data_nti3)[names(data_nti3)=="V3"] <- "atrisk"

write.csv(data_nti3,file="./Output/data_nti3_preweight.csv")

data_daniels$series = "daniels"
data_daniels$cases = round(data_daniels$cases / 6)
data_daniels$atrisk = round(data_daniels$atrisk / 6)
data_madsen$series = "madsen"
data_madsen$cases = round(data_madsen$cases / 6)
data_madsen$atrisk = round(data_madsen$atrisk / 6)

data_nti1$series = "NTI1"
data_nti2$series = "NTI2"
data_nti3$series = "NTI3"

data_nti1$cases = round(data_nti1$cases / 5)
data_nti2$cases = round(data_nti2$cases / 4)
data_nti3$cases = round(data_nti3$cases / 3)

data_nti1$atrisk = round(data_nti1$atrisk / 5)
data_nti2$atrisk = round(data_nti2$atrisk / 4)
data_nti3$atrisk = round(data_nti3$atrisk / 3)

    data_daniels$inc_med =  data_daniels$cases/data_daniels$atrisk
    data_daniels$inc_low = data_daniels$inc_med
    data_daniels$inc_upp = data_daniels$inc_med
    
    for(i in seq(1,dim(data_daniels)[1])){
        data_daniels$inc_low[i] = unlist(exactci(data_daniels$cases[i], data_daniels$atrisk[i], conf.level=0.95))[[1]]
        data_daniels$inc_upp[i] = unlist(exactci(data_daniels$cases[i], data_daniels$atrisk[i], conf.level=0.95))[[2]]
    }
    
    data_madsen$inc_med =  data_madsen$cases/data_madsen$atrisk
    data_madsen$inc_low = data_madsen$inc_med
    data_madsen$inc_upp = data_madsen$inc_med
    
    for(i in seq(1,dim(data_madsen)[1])){
        data_madsen$inc_low[i] = unlist(exactci(data_madsen$cases[i], data_madsen$atrisk[i], conf.level=0.95))[[1]]
        data_madsen$inc_upp[i] = unlist(exactci(data_madsen$cases[i], data_madsen$atrisk[i], conf.level=0.95))[[2]]
    }
    
    data_nti1$inc_med =  data_nti1$cases/data_nti1$atrisk
    data_nti1$inc_low = data_nti1$inc_med
    data_nti1$inc_upp = data_nti1$inc_med
    
    for(i in seq(1,dim(data_nti1)[1])){
        data_nti1$inc_low[i] = unlist(exactci(data_nti1$cases[i], data_nti1$atrisk[i], conf.level=0.95))[[1]]
        data_nti1$inc_upp[i] = unlist(exactci(data_nti1$cases[i], data_nti1$atrisk[i], conf.level=0.95))[[2]]
    }
    
    data_nti2$inc_med =  data_nti2$cases/data_nti2$atrisk
    data_nti2$inc_low = data_nti2$inc_med
    data_nti2$inc_upp = data_nti2$inc_med
    
    for(i in seq(1,dim(data_nti2)[1])){
        data_nti2$inc_low[i] = unlist(exactci(data_nti2$cases[i], data_nti2$atrisk[i], conf.level=0.95))[[1]]
        data_nti2$inc_upp[i] = unlist(exactci(data_nti2$cases[i], data_nti2$atrisk[i], conf.level=0.95))[[2]]
    }
    
    data_nti3$inc_med =  data_nti3$cases/data_nti3$atrisk
    data_nti3$inc_low = data_nti3$inc_med
    data_nti3$inc_upp = data_nti3$inc_med
    
    for(i in seq(1,dim(data_nti3)[1])){
        data_nti3$inc_low[i] = unlist(exactci(data_nti3$cases[i], data_nti3$atrisk[i], conf.level=0.95))[[1]]
        data_nti3$inc_upp[i] = unlist(exactci(data_nti3$cases[i], data_nti3$atrisk[i], conf.level=0.95))[[2]]
    }

Minobs = rbind(data_daniels,data_madsen)[,c("year","cases","series")]
names(Minobs) = c("time","value","series")
save(Minobs,file="./Output/Minobs.Rdata")

Subobs = rbind(data_nti1,data_nti2,data_nti3)[,c("year","cases","series")]
names(Subobs) = c("time","value","series")
save(Subobs,file="./Output/Subobs.Rdata")

Mininput = rbind(data_daniels,data_madsen)[,c("year","atrisk","series")]
names(Mininput) = c("time","value","series")
save(Mininput,file="./Output/Mininput.Rdata")

Subinput = rbind(data_nti1,data_nti2,data_nti3)[,c("year","atrisk","series")]
names(Subinput) = c("time","value","series")
save(Subinput,file="./Output/Subinput.Rdata")

#load(Minobs.R)
#load(Subobs.R)
#load(Mininput.R)
#load(Subinput.R)

obs <- list(
    Minobs = Minobs,
    Subobs = Subobs,
    Subminobs = data[which(data$Transition == "Submin"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),
    Minclinobs = data[which(data$Transition == "Minclin"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),
    Clinminobs = data[which(data$Transition == "Clinmin"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),

    Clinmin_clinobs = data[which(data$Transition == "Clinmin_clin"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),
    Submin_subobs = data[which(data$Transition == "Submin_sub"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),
    Minclin_minobs = data[which(data$Transition == "Minclin_min"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),

    InfminIobs = data[which(data$Transition == "InfminI"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),
    MininfIobs = data[which(data$Transition == "MininfI"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),
    Infmin_infIobs = data[which(data$Transition == "Infmin_infI"),] %>%
        dplyr::select(time = times, value = n.of.y_weight, series = record.id),


    t_med_obs = data.frame(time = 2, value = 2, series = "#666"),

    sub_over_clin_obs = data.frame(time = 2, value = 1, series = "sub:clin ratio"),
    min_over_infect_obs = data.frame(time = 2, value = 2.5, series = "min:infect ratio"),
    duration_obs = data.frame(time = 2, value = 2, series = "#666")
)

input <- list(
    Mininput = Mininput,
    Subinput = Subinput,
    Submininput = data[which(data$Transition == "Submin"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),
    Minclininput = data[which(data$Transition == "Minclin"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),
    Clinmininput = data[which(data$Transition == "Clinmin"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),

    Clinmin_clininput = data[which(data$Transition == "Clinmin_clin"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),
    Submin_subinput = data[which(data$Transition == "Submin_sub"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),
    Minclin_mininput = data[which(data$Transition == "Minclin_min"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),

    InfminIinput = data[which(data$Transition == "InfminI"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),
    MininfIinput = data[which(data$Transition == "MininfI"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id),
    Infmin_infIinput = data[which(data$Transition == "Infmin_infI"),] %>%
        dplyr::select(time = times, value = d.of.y_weight, series = record.id)
)


