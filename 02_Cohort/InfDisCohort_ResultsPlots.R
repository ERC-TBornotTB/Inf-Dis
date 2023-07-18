report_trajectories_inf <- function(simulation, n_people, years, threshold_months = 9, threshold_changes = 2) {
  if (n_people > ncol(simulation)) {
    n_people = ncol(simulation)
    warning("fewer people in simulation than requested, analysing on full simulation size")
  }
  if (years > (12 * nrow(simulation))) {
    years = floor(nrow(simulation))
    warning("fewer years in simulation than requested, analysing on full years in simulation")
  }
  if (threshold_months < 0 || threshold_months > 12) stop("provide threshold_months in range 0 - 12, or use default of 9")
  if (threshold_changes < 0 || threshold_changes > 12) stop("provide threshold_changes in range 0 - 12, or use default of 2")
  Tr <- matrix(data = NA ,nrow = years + 1, ncol = n_people)
  for (p in seq(1,n_people)) {
    if (simulation[1,p] == "i") {
      Tr[1,p] = "infection"
    } else if (simulation[1,p] == "m") {
      Tr[1,p] = "minimal"
    } else if (simulation[1,p] == "s") {
      Tr[1,p] = "subclinical"
    } else if (simulation[1,p] == "c") {
      Tr[1,p] = "clinical"
    }
    for (y in seq(1,years)) {
      start <- (y * 12) - 10
      end <- (y * 12) + 1
      if (Tr[y,p] %in% c("treat","death","recover","cleared")) {
        Tr[y + 1,p] = Tr[y,p]
      } else {
        Tr[y + 1,p] <- transitional_inf(simulation[(start:end),p], threshold_months, threshold_changes)
      }
    }
  }
  return(Tr)
}

transitional_inf <- function(data, threshold_months, threshold_changes) {
  treat   <- length(which(data == "t"))
  death   <- length(which(data == "d"))
  clear   <- length(which(data == "x"))
  recover <- length(which(data == "r"))
  clin    <- length(which(data == "c"))
  sub     <- length(which(data == "s"))
  min     <- length(which(data == "m"))
  infect  <- length(which(data == "i"))
  changes <- count_changes(data)
  if (treat == 1) {
    return("treat")
  } else if (death == 1) {
    return("death")
  } else if (clear == 1) {
    return("cleared")
  } else if (recover == 1) {
    return("recover")
  } else if (clin >= threshold_months && changes <= threshold_changes) {
    return("clinical")
  } else if (sub >= threshold_months && changes <= threshold_changes) {
    return("subclinical")
  } else if (min >= threshold_months && changes <= threshold_changes) {
    return("minimal")
  } else if (infect >= threshold_months && changes <= threshold_changes) {
    return("infection")
  } else {
    return("transitional")
  }
}

count_changes <- function(data){
  changes <- 0
  old <- data[1]
  for (i in seq(2,12)) {
    new <- data[i]
    if (new != old) {
      changes <-  changes + 1
    }
    old <- new
  }
  return(changes)
}

summarise_trajectories_inf <- function(trajectories, years) {
  Tr <- matrix(data = NA, nrow = years + 1, ncol = 10)
  colnames(Tr) <- c("year","death","treat","clinical","transitional","subclinical","minimal","infection","recover","cleared")
  for (y in seq(1,years + 1)) {
    Tr[y,"year"] = y - 1
    Tr[y,"death"] = length(which(trajectories[y,] == "death"))
    Tr[y,"treat"] = length(which(trajectories[y,] == "treat"))
    Tr[y,"clinical"] = length(which(trajectories[y,] == "clinical"))
    Tr[y,"transitional"] = length(which(trajectories[y,] == "transitional"))
    Tr[y,"subclinical"] = length(which(trajectories[y,] == "subclinical"))
    Tr[y,"minimal"] = length(which(trajectories[y,] == "minimal"))
    Tr[y,"infection"] = length(which(trajectories[y,] == "infection"))
    Tr[y,"recover"] = length(which(trajectories[y,] == "recover"))
    Tr[y,"cleared"] = length(which(trajectories[y,] == "cleared"))
  }
  return(Tr)
}

plot_trajectories_inf <- function(summary, years, title) {
  pal <- c("death"        = '#545454',
           "clinical"     = '#AE0D0A', 
           "transitional" = '#754668',
           "subclinical"  = '#CB4C15',
           "minimal"      = '#CBA715',
           "infection"    = '#537D8D',
           "recover"      = '#8A8A8A',
           "cleared"      = '#C0C0C0')
  summary_melt <- melt(as.data.frame(summary), id.vars = c("year"), variable.name = "state", value.name = "freq")
  p <- ggplot(data = summary_melt, aes(x = year, y = freq, fill = state)) +
       geom_bar(stat = "identity") +
       scale_fill_manual(labels = c("Death","Clinical",
                                    "Transitional","Subclinical",
                                    "Minimal",
                                    "Infection","Recovered","Cleared"),
                         values = pal) +
       scale_y_break(c(0,5)) +
       scale_y_break(c(10,80)) +
       scale_x_continuous(labels = as.character(c(0,1,2,3,4,5,6,7,8,9,10)), breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
       scale_y_continuous(labels = as.character(c(0,10,20,30,40,50,60,65,70,75,80,85,90,95,100)),
                          breaks = c(0,10,20,30,40,50,60,65,70,75,80,85,90,95,100),expand=c(0, 0)) +
       xlab("Years since infection") +
       ylab("Percent of cohort") +
       theme_classic(base_size=35)
  return(p)
}

# Incidence of each state

incidence <- function(n_reps, n_years, n_steps, n_people, prop_inf, params, method, endpoint, treatment=0, prop_min=0, prop_sub=0){
  res <- as.data.frame(c(0:n_years))
  colnames(res) <- c("year")
  for (i in (1:n_reps)){
    IBM_sim_inc  <- full_simulation_inc(n_steps, n_people, prop_inf, params, method, endpoint, treatment=0, prop_min=0, prop_sub=0)
    if (endpoint=="min"){
      nstate <- 3
    } else if (endpoint=="sub"){
      nstate <- 4
    } else if (endpoint=="clin"){
      nstate <- 5
    }
    counts <- matrix(NA,n_people,7)
    for (j in 1:n_people){
      counts[j,1] <- length(which(IBM_sim_inc[,j]=='i'))-1 # i
      counts[j,2] <- length(which(IBM_sim_inc[,j]=='r'))   # r
      counts[j,3] <- length(which(IBM_sim_inc[,j]=='m'))   # m
      counts[j,4] <- length(which(IBM_sim_inc[,j]=='s'))   # s
      counts[j,5] <- length(which(IBM_sim_inc[,j]=='c'))   # c
    }
    inccounts <- subset(counts, counts[,nstate]>0)
    if (dim(inccounts)[1] != 0) {
      for (k in 1:dim(inccounts)[1]){
        inccounts[k,6] <- sum(inccounts[k,1:(nstate-1)])
        for (m in 1:n_years){
          if (inccounts[k,6] <= (n_years-m+1)*12){
            inccounts[k,7] <- (n_years-m+1)
          }
        }
      }
      temp <- as.data.frame(table(inccounts[,7])/n_people)
    } else {
      temp <- as.data.frame(matrix(0,1,2))
    }
    colnames(temp) <- c("year","")
    res <- merge(res,temp,by="year",all.x = TRUE)
    res[is.na(res)] <- 0
  }
  return(res)
}



inccum <- function(inc,n_years,n_reps){
  cummul <- as.data.frame(c(0:n_years))
  colnames(cummul) <- c("year")
  for (i in 1:n_years){
    for (j in 1:n_reps){
      cummul[1  ,j+1] <- 0
      cummul[i+1,j+1] <- cummul[i,j+1] + inc[i+1,j+1]
    }
  }
  return(cummul)
}

incsum <- function(incres){
  results_summary <- apply(incres[,-1], 1, function(x) CI_min_max(x))
  results_summary <- as.data.frame(t(results_summary))
  colnames(results_summary) <- c("med","lowci","uppci","min","max")
  results_summary$year <- c(0:10)
  results_summary <- results_summary[,c(6,1,2,3,4,5)]
  return(results_summary)
}

CI_min_max <- function(data){
  ci <- unname(CI(data, ci = 0.95))
  ci <- unname(quantile(data,c(0.025,0.5,0.975)))
  return(c(ci[2], ci[1], ci[3], min(data), max(data)))
}

incplotmin <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_ribbon(aes(x=year,ymin=lowci,ymax=uppci),size=2,fill='#004550',alpha=.15) +
    geom_line(aes(x=year,y=med),size=2,color='#004550') +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous(str_wrap("Incidence of minimal TB",width=13),expand=c(0, 0),lim=c(0,0.11)) +
    theme_classic(base_size=50) +
    theme(legend.position=c(0.835, 0.92))
  return(plot)
}

incplotsub <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_ribbon(aes(x=year,ymin=lowci,ymax=uppci),size=2,fill='#004550',alpha=.15) +
    geom_line(aes(x=year,y=med),size=2,color='#004550') +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous(str_wrap("Incidence of subclinical TB",width=15),expand=c(0, 0),lim=c(0,0.11)) +
    theme_classic(base_size=50) +
    theme(legend.position=c(0.68, 0.9))
  return(plot)
}

incplotclin <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_ribbon(aes(x=year,ymin=lowci,ymax=uppci),size=2,fill='#004550',alpha=.15) +
    geom_line(aes(x=year,y=med),size=2,color='#004550') +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous(str_wrap("Incidence of clinical TB",width=15),expand=c(0, 0),lim=c(0,0.11)) +
    theme_classic(base_size=50)+
    theme(legend.position=c(0.57, 0.9))
  return(plot)
}

inccumplot <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_ribbon(aes(x=year,ymin=lowci,ymax=uppci),size=1,fill='#004550',alpha=.15) +
    geom_line(aes(x=year,y=med),size=2,color='#004550') +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous("Incidence rate",expand=c(0, 0),lim=c(0,0.11)) +
    theme_classic(base_size=50)
  return(plot)
}
