## ---Load Libraries------------------------------------------------------------------ ##
rm(list = ls())

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

load("Minobs.Rdata")
load("Subobs.Rdata")
load("Mininput.Rdata")
load("Subinput.Rdata")

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


