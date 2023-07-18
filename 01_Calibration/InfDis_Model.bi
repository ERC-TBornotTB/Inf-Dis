model TB_model {
    param mort

    param infclr
    param infmin
    param infsub
    param minrec

    param min_sub
    param sub_min
    param sub_clin
    param clin_sub

    param phi

    state I
    state t_med                         (has_output = 0)

    state Inf
    state Clr
    state Rec
    state Min1
    state Min2                          (has_output = 0)
    state Min3                          (has_output = 0)
    state Sub1                           
    state Sub2
    state Clin                          (has_output = 0)

    state Z_min
    state Z_sub
    state Z_clin

    state Minmin                        (has_output = 0)
    state Submin
    state Clinmin

    state Minsub
    state Subsub                        (has_output = 0)
    state Clinsub

    state Minclin
    state Subclin
    state Clinclin                      (has_output = 0)

    state Minmin_sub                    (has_output = 0)
    state Submin_sub

    state Minsub_min
    state Subsub_min                    (has_output = 0)
    state Clinsub_min                   (has_output = 0)

    state Minsub_clin                   (has_output = 0)
    state Subsub_clin                   (has_output = 0)
    state Clinsub_clin

    state Subclin_sub
    state Clinclin_sub                  (has_output = 0)

    state Minmin_clin                   (has_output = 0)
    state Submin_clin                   (has_output = 0)
    state Clinmin_clin

    state Minclin_min
    state Subclin_min                   (has_output = 0)
    state Clinclin_min                  (has_output = 0)

    state MinminI                       (has_output = 0)
    state SubminI                       (has_output = 0)
    state ClinminI                      (has_output = 0)
    state InfminI

    state MininfI
    state SubinfI                       (has_output = 0)
    state ClininfI                      (has_output = 0)
    state InfinfI                       (has_output = 0)

    state Minmin_infI                   (has_output = 0)
    state Submin_infI                   (has_output = 0)
    state Clinmin_infI                  (has_output = 0)
    state Infmin_infI

    state Mininf_minI
    state Subinf_minI                   (has_output = 0)
    state Clininf_minI                  (has_output = 0)
    state Infinf_minI                   (has_output = 0)

    state sub_over_clin_p
    state min_over_infect_p
    state duration                           

    dim series(36)

    obs t_med_obs[series]
    obs sub_over_clin_obs[series]
    obs min_over_infect_obs[series]
    obs duration_obs[series]

    obs Minobs[series]
    obs Subobs[series]

    obs Subminobs[series]
    obs Minsubobs[series]
    obs Clinminobs[series]
    obs Minclinobs[series]

    obs Submin_subobs[series]
    obs Clinmin_clinobs[series]
    obs Minclin_minobs[series]

    obs InfminIobs[series]
    obs MininfIobs[series]

    obs Infmin_infIobs[series]
    obs Mininf_minIobs[series]

    input Mininput[series]
    input Subinput[series]

    input Submininput[series]
    input Minsubinput[series]
    input Clinmininput[series]
    input Minclininput[series]

    input Submin_subinput[series]
    input Clinmin_clininput[series]
    input Minclin_mininput[series]

    input InfminIinput[series]
    input MininfIinput[series]

    input Infmin_infIinput[series]
    input Mininf_minIinput[series]

    sub initial {
        Inf <- 1
	Clr <- 0
        Rec <- 0
        Min1 <- 0
        Min2 <- 0
        Min3 <- 0
        Sub1 <- 0
        Sub2 <- 0
        Clin <- 0
        Z_min <- 0
        Z_sub <- 0
        Z_clin <- 0

        Minmin <- 1
        Submin <- 0
        Clinmin <- 0

        Minsub <- 0
        Subsub <- 1
        Clinsub <- 0

        Minclin <- 0
        Subclin <- 0
        Clinclin <- 1

        Minmin_sub <- 1
        Submin_sub <- 0

        Minsub_min <- 0
        Subsub_min <- 1
        Clinsub_min <- 0

        Minsub_clin <- 0
        Subsub_clin <- 1
        Clinsub_clin <- 0

        Subclin_sub <- 0
        Clinclin_sub <- 1

        Minmin_clin <- 1
        Submin_clin <- 0
        Clinmin_clin <- 0

        Minclin_min <- 0
        Subclin_min <- 0
        Clinclin_min <- 1

        MinminI <- 1
        SubminI <- 0
        ClinminI <- 0
        InfminI <- 0

        MininfI <- 0
        SubinfI <- phi
        ClininfI <- (1 - phi)
        InfinfI <- 1

        Minmin_infI <- 1
        Submin_infI <- 0
        Clinmin_infI <- 0
        Infmin_infI <- 0

        Mininf_minI <- 0
        Subinf_minI <- phi
        Clininf_minI <- (1 - phi)
        Infinf_minI <- 1

        I <- 1
    }

    sub parameter {

        infclr ~ uniform(0,6)
        infmin ~ uniform(0,3)
        infsub ~ uniform(0,3)
        minrec ~ uniform(0,3)

        min_sub ~ uniform(0,3)
        sub_min ~ uniform(0,3)
        sub_clin ~ uniform(0,3)
        clin_sub ~ uniform(0,3)
        phi ~ uniform(0,1)

        mort ~ gaussian(0.389, 0.028)
    }

    sub transition (delta = 0.1) {

        Z_min <- (t_now % 1 == 0 ? 0 : Z_min)
        Z_sub <- (t_now % 1 == 0 ? 0 : Z_sub)
        Z_clin <- (t_now % 1 == 0 ? 0 : Z_clin)

        ode {

            dInf/dt  = - infmin*Inf  - infclr*Inf - infsub*Inf;
            dClr/dt  =   infclr*Inf;
	    dRec/dt  =   minrec*Min1 + minrec*Min2 + minrec*Min3;
            dMin1/dt =   infmin*Inf - min_sub*Min1 - minrec*Min1;
            dMin2/dt =   sub_min*Sub1 - min_sub*Min2 - minrec*Min2;
            dMin3/dt =   sub_min*Sub2 - min_sub*Min3 - minrec*Min3;
            dSub1/dt =   min_sub*Min1 + min_sub*Min2 - sub_min*Sub1 - sub_clin*Sub1 + infsub*Inf;
            dSub2/dt =   min_sub*Min3 - sub_min*Sub2 + clin_sub*Clin - sub_clin*Sub2;
            dClin/dt =   sub_clin*Sub1 + sub_clin*Sub2 - clin_sub*Clin - mort*Clin;
            dZ_min/dt  =  infmin*Inf;
            dZ_sub/dt  =  infsub*Inf + min_sub*Min1;
            dZ_clin/dt =  sub_clin*Sub1;

            dMinmin/dt  = - min_sub * Minmin + sub_min * Submin - minrec * Minmin
            dSubmin/dt  =   min_sub * Minmin - sub_min * Submin - sub_clin * Submin + clin_sub * Clinmin
            dClinmin/dt =   sub_clin * Submin - clin_sub * Clinmin - mort * Clinmin

            dMinsub/dt  = - min_sub * Minsub + sub_min * Subsub - minrec * Minsub
            dSubsub/dt  =   min_sub * Minsub - sub_min * Subsub - sub_clin * Subsub + clin_sub * Clinsub
            dClinsub/dt =   sub_clin * Subsub - clin_sub * Clinsub - mort * Clinsub

            dMinclin/dt  = - min_sub * Minclin + sub_min * Subclin - minrec * Minclin
            dSubclin/dt  =   min_sub * Minclin - sub_min * Subclin - sub_clin * Subclin + clin_sub * Clinclin
            dClinclin/dt =   sub_clin * Subclin - clin_sub * Clinclin - mort * Clinclin

            dMinmin_sub/dt = - min_sub * Minmin_sub - minrec * Minmin_sub
            dSubmin_sub/dt =   min_sub * Minmin_sub

            dMinsub_min/dt  =   sub_min * Subsub_min
            dSubsub_min/dt  = - sub_min * Subsub_min - sub_clin * Subsub_min + clin_sub * Clinsub_min
            dClinsub_min/dt =   sub_clin * Subsub_min - clin_sub * Clinsub_min - (mort * Clinsub_min)

            dMinsub_clin/dt  = - min_sub * Minsub_clin + sub_min * Subsub_clin - minrec * Minsub_clin
            dSubsub_clin/dt  =   min_sub * Minsub_clin - sub_min * Subsub_clin - sub_clin * Subsub_clin
            dClinsub_clin/dt =   sub_clin * Subsub_clin

            dSubclin_sub/dt  =   clin_sub * Clinclin_sub
            dClinclin_sub/dt = - clin_sub * Clinclin_sub - (mort * Clinclin_sub)

            dMinmin_clin/dt  = - min_sub * Minmin_clin + sub_min * Submin_clin - minrec * Minmin_clin
            dSubmin_clin/dt  =   min_sub * Minmin_clin - sub_min * Submin_clin - sub_clin * Submin_clin
            dClinmin_clin/dt =   sub_clin * Submin_clin

            dMinclin_min/dt  =   sub_min * Subclin_min
            dSubclin_min/dt  = - sub_min * Subclin_min - sub_clin * Subclin_min + clin_sub * Clinclin_min
            dClinclin_min/dt =   sub_clin * Subclin_min - clin_sub * Clinclin_min - (mort * Clinclin_min)


            dMinminI/dt  = - min_sub * MinminI + sub_min * SubminI - minrec * MinminI
            dSubminI/dt  =   min_sub * MinminI - sub_min * SubminI - sub_clin * SubminI + clin_sub * ClinminI
            dClinminI/dt =   sub_clin * SubminI - clin_sub * ClinminI - mort * ClinminI

            dMininfI/dt  = - min_sub * MininfI + sub_min * SubinfI - minrec * MininfI
            dSubinfI/dt  =   min_sub * MininfI - sub_min * SubinfI - sub_clin * SubinfI + clin_sub * ClininfI
            dClininfI/dt =   sub_clin * SubinfI - clin_sub * ClininfI - mort * ClininfI

            dMinmin_infI/dt  = - min_sub * Minmin_infI - minrec * Minmin_infI
            dSubmin_infI/dt  =   min_sub * Minmin_infI - sub_clin * Submin_infI + clin_sub * Clinmin_infI
            dClinmin_infI/dt =   sub_clin * Submin_infI - clin_sub * Clinmin_infI

            dMininf_minI/dt  =   sub_min * Subinf_minI
            dSubinf_minI/dt  = - sub_min * Subinf_minI - sub_clin * Subinf_minI + clin_sub * Clininf_minI
            dClininf_minI/dt =   sub_clin * Subinf_minI - clin_sub * Clininf_minI - (mort * Clininf_minI)
        }

        I <- SubinfI + ClininfI
        t_med <- -2*log(2)/log(I)

        InfminI <- SubminI + ClinminI
        InfinfI <- SubinfI + ClininfI
        Infmin_infI <- Submin_infI + Clinmin_infI
        Infinf_minI <- Subinf_minI + Clininf_minI
        
        sub_over_clin_p <- (clin_sub + mort) / sub_clin
        min_over_infect_p <- (
        ((sub_min * infsub + (infmin) * (sub_clin + sub_min)) * (clin_sub + mort)) - (infmin * clin_sub * sub_clin)
        ) / (
        (infsub * (min_sub + minrec) + (infmin) * min_sub) * (clin_sub + mort + sub_clin)
        )
       duration <- -2 * log(2) / log(SubinfI + ClininfI)
    }

    sub observation {

        Minobs[series] ~ binomial(size = Mininput[series], prob = Z_min)
        Subobs[series] ~ binomial(size = Subinput[series], prob = Z_sub)

        Subminobs[series] ~ binomial(size = Submininput[series], prob = Submin)
        Minsubobs[series] ~ binomial(size = Minsubinput[series], prob = Minsub)
        Clinminobs[series] ~ binomial(size = Clinmininput[series], prob = Clinmin)
        Minclinobs[series] ~ binomial(size = Minclininput[series], prob = Minclin)

        Submin_subobs[series] ~ binomial(size = Submin_subinput[series], prob = Submin_sub)
        Clinmin_clinobs[series] ~ binomial(size = Clinmin_clininput[series], prob = Clinmin_clin)
        Minclin_minobs[series] ~ binomial(size = Minclin_mininput[series], prob = Minclin_min)

        InfminIobs[series] ~ binomial(size = InfminIinput[series], prob = InfminI)
        MininfIobs[series] ~ binomial(size = MininfIinput[series], prob = MininfI)

        Infmin_infIobs[series] ~ binomial(size = Infmin_infIinput[series], prob = Infmin_infI)
        Mininf_minIobs[series] ~ binomial(size = Mininf_minIinput[series], prob = Mininf_minI)

        
        sub_over_clin_obs[series] ~ gaussian(mean = sub_over_clin_p, std = 0.25)
        min_over_infect_obs[series] ~ gaussian(mean = min_over_infect_p, std = 0.5)
        duration_obs[series] ~ gaussian(mean = duration, std = 0.5)

    }

    sub proposal_parameter {
        infclr ~ truncated_gaussian(mean = infclr, std = 0.01, lower = 0, upper = 6)
        infmin ~ truncated_gaussian(mean = infmin, std = 0.01, lower = 0, upper = 3)
        infsub ~ truncated_gaussian(mean = infsub, std = 0.01, lower = 0, upper = 3)
        minrec ~ truncated_gaussian(mean = minrec, std = 0.01, lower = 0, upper = 3)

        min_sub  ~ truncated_gaussian(mean = min_sub,  std = 0.01, lower = 0, upper = 3)
        sub_min  ~ truncated_gaussian(mean = sub_min,  std = 0.01, lower = 0, upper = 3)
        sub_clin ~ truncated_gaussian(mean = sub_clin, std = 0.01, lower = 0, upper = 3)
        clin_sub ~ truncated_gaussian(mean = clin_sub, std = 0.01, lower = 0, upper = 3)
        phi ~ truncated_gaussian(mean = phi, std = 0.01, lower = 0, upper = 1)

        mort ~ truncated_gaussian(mean = mort, std = 0.01, lower = 0)
    }
}
