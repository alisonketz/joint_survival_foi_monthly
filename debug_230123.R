
dNegCapPosMort <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        x = integer(),
        e = integer(0), #e, age of entry
        r = integer(0), #r, age of last known alive
        s = integer(0), #s, age of known mortality
        dn1 = integer(0), #last interval of test negative
        dn = integer(0), #interval of test positive (could be interval of mortality)
        sex = integer(0),
        age2date = integer(0),
        beta_sex = double(0),
        beta0_sus = double(0),
        beta0_inf = double(0),
        age_effect_surv = double(1),
        period_effect_surv = double(1),
        f_age_foi = double(1),
        m_age_foi = double(1),
        age_lookup_f = double(1),
        age_lookup_m = double(1),
        period_lookup=double(1),
        f_period_foi=double(1),
        m_period_foi=double(1),
        space = double(0),
        log = double()
        ) {

    llik <- 0 #intialize log-likelihood
    lam_inf <- nimNumeric(s)
    lam_foi <- nimNumeric(s)
    lam_sus <- nimNumeric(s)
    lik_temp <- nimNumeric(s)

    #############################################
    # preliminary hazards for the likelihood
    #############################################
    
    #survival hazard for susceptible deer
    n_indx_sus <- length(e:(dn - 1))
    lam_sus[e:(dn - 1)] <- exp(rep(beta0_sus, n_indx_sus)  + 
                    age_effect_surv[e:(dn - 1)] +
                    period_effect_surv[(e + age2date):(dn - 1 + age2date)] +
                    rep(beta_sex * sex, n_indx_sus))
    
    #survival hazard while infected
    n_indx_inf <- length(dn1:(s - 1))
    lam_inf[dn1:(s - 1)] <- exp(rep(beta0_inf, n_indx_inf) +
        age_effect_surv[dn1:(s - 1)] +
        period_effect_surv[(dn1 + age2date):(s - 1 + age2date)] +
        rep(beta_sex * sex, n_indx_inf)
        )

    #force of infection infection hazard
    lam_foi[1:(dn - 1)] <- exp(rep(space, dn - 1) +
        sex * (f_age_foi[age_lookup_f[1:(dn - 1)]] +
            f_period_foi[period_lookup[(1 + age2date):(dn - 1 + age2date)]]) +
        (1 - sex) * (m_age_foi[age_lookup_m[1:(dn - 1)]] +
            m_period_foi[period_lookup[(1 + age2date):(dn - 1 + age2date)]])
        )
    #######################################
    ### calculating the joint likelihood
    #######################################

    lik_temp[dn1 + 1] <- lam_foi[dn1 + 1] *
                      exp(-sum(lam_inf[(dn1 + 1):(r - 1)])) *
                      (1 - exp(-sum(lam_inf[r:(s - 1)])))

    for (k in (dn1 + 2):(r - 2)) {
        lik_temp[k] <- lam_foi[k] *
                      exp(-sum(lam_foi[(dn1 + 1):k])) *
                      exp(-sum(lam_sus[(dn1 + 1):(k - 1)])) *
                      exp(-sum(lam_inf[k:(r - 1)])) *
                      (1 - exp(-sum(lam_inf[r:(s - 1)])))
    }
    for (k in (r - 1):(s - 2)) {
        lik_temp[k] <- lam_foi[k] *
                      exp(-sum(lam_foi[(dn1 + 1):k])) *
                      exp(-sum(lam_sus[(dn1 + 1):(k - 1)])) *
                      (1 - exp(-sum(lam_inf[k:(s - 1)])))
    }
    lik_temp[(s - 1)] <- lam_foi[(s - 1)] *
                        exp(-sum(lam_foi[(dn1 + 1):(s - 1)])) *
                        exp(-sum(lam_sus[(dn1 + 1):((s - 1) - 1)])) *
                        lam_inf[(s - 1)]
    lik <- exp(-sum(lam_sus[e:dn1])) *
           exp(-sum(lam_foi[1:dn1])) *
           sum(lik_temp[(dn1 + 1):(s - 1)])

    llik <- log(lik)

    returnType(double(0))
    if(log) return(llik) else return(exp(llik))    ## return log-likelihood
  })



#83,100,102,107 are returning -inf log lik

i=107
dNegCapPosMort(
        x = 1,
        e = d_fit_idead$left_age_month[i],
        r = d_fit_idead$right_age_rmonth[i],
        s = d_fit_idead$right_age_smonth[i],
        dn1 = d_fit_idead$left_age_month[i],
        dn = d_fit_idead$right_age_smonth[i],
        sex = d_fit_idead$sex[i],
        age2date = idead_age2date[i],
        beta_sex = beta_sex,
        beta0_sus = -12,#beta0_sus
        beta0_inf = -10,#beta0_inf
        age_effect_surv = age_effect_survival_test,
        period_effect_surv = period_effect_survival_test,
        f_age_foi = f_age_foi,
        m_age_foi = m_age_foi,
        age_lookup_f = age_lookup_col_f,
        age_lookup_m = age_lookup_col_m,
        period_lookup = period_lookup,
        f_period_foi = f_period_foi,
        m_period_foi = m_period_foi,
        space = 0,
        log = TRUE
        )

# test <- c()
# for(i in 1:nrow(d_fit_idead)){
# test[i] <- dNegCapPosMort(
#         x = 1,
#         e = d_fit_idead$left_age_month[i],
#         r = d_fit_idead$right_age_rmonth[i],
#         s = d_fit_idead$right_age_smonth[i],
#         dn1 = d_fit_idead$left_age_month[i],
#         dn = d_fit_idead$right_age_smonth[i],
#         sex = d_fit_idead$sex[i],
#         age2date = idead_age2date[i],
#         beta_sex = beta_sex,
#         beta0_sus = beta0_sus,
#         beta0_inf = beta0_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_col_f,
#         age_lookup_m = age_lookup_col_m,
#         period_lookup = period_lookup,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,
#         log = TRUE
#         )
#  }
# test



