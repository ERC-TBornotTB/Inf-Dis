
fit_plot <- function(transition, posterior, data, title){

    data <-  data[which(data$Transition == transition),]

    posterior <-  posterior[which(posterior$var == transition),]
    colour <- ifelse(transition %in% (c("Submin","Minsub","Clinmin","Minclin","Clinsub","Subclin","InfminI","MininfI")),
                     "dodgerblue4",
                     "tomato4")

    plot <- ggplot() +
        geom_ribbon(data = posterior, aes(x = time, ymin = `2.5%`, ymax = `97.5%`), fill = colour, alpha = 0.5) +
        geom_line(data = posterior, aes(x = time, y = Median), colour = colour) +
        xlim(c(0,15)) +
        ylim(c(0,1)) +
        ggtitle(title) +
        xlab("") +
        ylab("") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))

    if (nrow(data) > 0) {

        plot <- plot +
            geom_pointrange(data = data, aes(x = times, y = value,ymin = error_low, ymax = error_upp), colour = colour, shape = 20, size = 0.25)#  +
        #theme(legend.position="bottom") +
        #geom_text(data = data, aes(label = record.id,x = times, y = value),hjust = 0, vjust = 0, size = 2.5)


    }

    return(plot)
}


total_fit_plots <- function(data, posterior){

data$value = data$n.of.y_weight/data$d.of.y_weight
data$error_low = as.double(unlist(binom.confint(data$n.of.y_weight,data$d.of.y_weight,conf.level = 0.95, method = "wilson")["lower"]))
data$error_upp = as.double(unlist(binom.confint(data$n.of.y_weight,data$d.of.y_weight,conf.level = 0.95, method = "wilson")["upper"]))


min_sub_cs <-  fit_plot("Submin", posterior, data, "Minimal to Subclinical")
sub_min_cs <-  fit_plot("Minsub", posterior, data, "Subclinical to Minimal")
sub_clin_cs <-  fit_plot("Clinsub", posterior, data, "Subclinical to Clinical")
clin_sub_cs <-  fit_plot("Subclin", posterior, data, "Clinical to Subclinical")
min_clin_cs <-  fit_plot("Clinmin", posterior, data, "Minimal to Clinical")
clin_min_cs <-  fit_plot("Minclin", posterior, data, "Clinical to Minimal")
min_inf_cs <- fit_plot("InfminI", posterior, data, "Minimal to Infectious")
inf_min_cs <- fit_plot("MininfI",posterior, data, "Infectious to Minimal")

min_sub_tte <-  fit_plot("Submin_sub", posterior, data, "Minimal to Subclinical")
sub_min_tte <-  fit_plot("Minsub_min", posterior, data, "Subclinical to Minimal")
sub_clin_tte <-  fit_plot("Clinsub_clin", posterior, data, "Subclinical to Clinical")
clin_sub_tte <-  fit_plot("Subclin_sub", posterior, data, "Clinical to Subclinical")
min_clin_tte <-  fit_plot("Clinmin_clin", posterior, data, "Minimal to Clinical")
clin_min_tte <-  fit_plot("Minclin_min", posterior, data, "Clinical to Minimal")
min_inf_tte <- fit_plot("Infmin_infI", posterior, data, "Minimal to Infectious")
inf_min_tte <- fit_plot("Mininf_minI",posterior, data, "Infectious to Minimal")


layout <- "
            AAAABBBBCCCCDDDD
            EEEEFFFFGGGGHHHH
            IIIIJJJJKKKKLLLL
            MMMMNNNNOOOOPPPP"

fit_plot <- (
    min_sub_tte + sub_min_tte +
    min_sub_cs + sub_min_cs +
    sub_clin_tte + clin_sub_tte +
    sub_clin_cs + clin_sub_cs +
    min_clin_tte + clin_min_tte +
    min_clin_cs + clin_min_cs +
    min_inf_tte + inf_min_tte +
    min_inf_cs + inf_min_cs +
    plot_layout(design = layout)) +
    theme_minimal()

x_lab <- "Time (years)"
y_lab <- "Proportion moved from initial to final state"

x_lab <-
    ggplot() +
    annotate(geom = "text", x = 1, y = 1, label = x_lab) +
    coord_cartesian(clip = "off") +
    theme_void()

y_lab <-
    ggplot() +
    annotate(geom = "text", x = 1, y = 1, label = y_lab, angle = 90) +
    coord_cartesian(clip = "off") +
    theme_void()

fit_plot <- (fit_plot / x_lab + plot_layout(heights = c(1, 0.01)))
fit_plot <- (y_lab + fit_plot + plot_layout(widths = c(1, 100)))

return(fit_plot)
}

alt_fit_plot <- function(transition, posterior, data, title){

    data <-  data[which(data$Transition == transition),]

    posterior_new <-  posterior[which(posterior$var == transition),] %>%
        filter(time %% 1 == 0)



    plot <- ggplot() +
        geom_ribbon(data = posterior_new, aes(x = time, ymin = `2.5%`, ymax = `97.5%`), fill = "dodgerblue4", alpha = 0.5) +
        geom_line(data = posterior_new, aes(x = time, y = Median), colour = "dodgerblue4") +
        ggtitle(title) +
        xlab("") +
        ylab("") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))

    if (nrow(data) > 0) {

        plot <- plot +
            geom_pointrange(data = data, aes(x = times, y = value,ymin = error_low, ymax = error_upp), colour = "tomato", shape = 20, size = 0.25)#  +
        #theme(legend.position="bottom") +
        #geom_text(data = data, aes(label = record.id,x = times, y = value),hjust = 0, vjust = 0, size = 2.5)


    }

    return(plot)
}

incplotmin <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_errorbar(data=mintargets,aes(x=time,ymin=lower,ymax=upper,color=series),size=2,width=0.5) +
    geom_point(data=mintargets,aes(x=time,y=value,color=series,shape=series),size=8) +
    #geom_segment(x=0.95
    #geom_segment(x=1.95
    #geom_segment(x=2.95
    geom_segment(x=3.95,xend=3.95,y=0           ,yend=0.2,color='#AE0D0A',size=2,arrow=arrow(length=unit(1.5, "cm"))) +
    geom_segment(x=4.95,xend=4.95,y=0           ,yend=0.2,color='#AE0D0A',size=2,arrow=arrow(length=unit(1.5, "cm"))) +
    geom_segment(x=1.05,xend=1.05,y=0.0480607784,yend=0.2,color='#754668',size=2,arrow=arrow(length=unit(1.5, "cm"))) +
    #geom_segment(x=2.05
    #geom_segment(x=3.05
    geom_segment(x=4.05,xend=4.05,y=0,yend=0.2,color='#754668',size=2,arrow=arrow(length=unit(1.5, "cm"))) +
    geom_segment(x=5.05,xend=5.05,y=0,yend=0.2,color='#754668',size=2,arrow=arrow(length=unit(1.5, "cm"))) +
    geom_segment(x=6.05,xend=6.05,y=0,yend=0.2,color='#754668',size=2,arrow=arrow(length=unit(1.5, "cm"))) +
    geom_ribbon(aes(x=time,ymin=lower,ymax=upper),size=2,fill='#004550',alpha=.15) +
    geom_line(aes(x=time,y=median),size=2,color='#004550') +
    scale_colour_manual(name=NULL,values =c('daniels'='#AE0D0A','madsen'='#754668'), labels = c('Daniels 1944','Madsen 1942')) +
    scale_shape_manual(name=NULL,values =c('daniels'=17,'madsen'=19), labels = c('Daniels 1944','Madsen 1942')) +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous(str_wrap("Incidence of minimal TB",width=13),expand=c(0, 0),lim=c(0,0.2)) +
    theme_classic(base_size=50) +
    theme(legend.position=c(0.835, 0.92))
  return(plot)
}

incplotsub <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_errorbar(data=subtargets,aes(x=time,ymin=lower,ymax=upper,color=series),size=2,width=0.5) +
    geom_point(data=subtargets,aes(x=time,y=value,color=series,shape=series),size=8) +
    geom_ribbon(aes(x=time,ymin=lower,ymax=upper),size=2,fill='#004550',alpha=.15) +
    geom_line(aes(x=time,y=median),size=2,color='#004550') +
    scale_colour_manual(name=NULL,values =c('NTI1'='#AE0D0A','NTI2'='#754668','NTI3'='#CB4C15'), labels = c('NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976')) +
    scale_shape_manual(name=NULL,values =c('NTI1'=17,'NTI2'=19,'NTI3'=15), labels = c('NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976')) +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous(str_wrap("Incidence of subclinical TB",width=15),expand=c(0, 0),lim=c(0,0.07)) +
    theme_classic(base_size=50) +
    theme(legend.position=c(0.68, 0.9))
  return(plot)
}

incplotclin <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_errorbar(data=clintargets,aes(x=time,ymin=lower,ymax=upper,color=series),size=2,width=0.5) +
    geom_point(data=clintargets,aes(x=time,y=value,color=series,shape=series),size=8) +
    geom_ribbon(aes(x=time,ymin=lower,ymax=upper),size=2,fill='#004550',alpha=.15) +
    geom_line(aes(x=time,y=median),size=2,color='#004550') +
    scale_colour_manual(name=NULL,values =c('NTI1'='#AE0D0A','NTI2'='#754668','NTI3'='#CB4C15'), labels = c('NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976')) +
    scale_shape_manual(name=NULL,values =c('NTI1'=17,'NTI2'=19,'NTI3'=15), labels = c('NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976','NTI 1974, Krishnamurthy 1976')) +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous(str_wrap("Incidence of clinical TB",width=15),expand=c(0, 0),lim=c(0,0.05)) +
    theme_classic(base_size=50)+
    theme(legend.position=c(0.68, 0.9))
  return(plot)
}

cuminc <- function(anninc){
  cuminc <- anninc[1,]
  for (i in 2:10){
    cuminc[i,1:3] <- anninc[i,1:3]
    for (j in 4:length(anninc)){
      cuminc[i,j] <- cuminc[i-1,j] + anninc[i,j]
    }
  }
  return(cuminc)
}

inccumplot <- function(summaryres){
  plot <- ggplot(summaryres) +
    geom_ribbon(aes(x=time,ymin=lower,ymax=upper),size=1,fill='#004550',alpha=.15) +
    geom_line(aes(x=time,y=median),size=2,color='#004550') +
    scale_x_continuous("Years since infection",expand=c(0, 0),breaks=c(0,1,2,3,4,5,6,7,8,9,10)) + 
    scale_y_continuous("Incidence rate",expand=c(0, 0),lim=c(0,0.2)) +
    theme_classic(base_size=50)
  return(plot)
}
