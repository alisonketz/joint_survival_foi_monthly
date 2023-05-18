


##################################################################################################################################
###
### Preliminary constants for running in the model
###
##################################################################################################################################


##########################################
###
### age class indexing and period indexing
###
##########################################

age_lookup_col <- c(age_lookup_f,rep(n_age,max(d_surv$right_age_smonth,na.rm=TRUE) - n_age_lookup))
age_lookup_col_f <- age_lookup_col
age_lookup_col_m <- age_lookup_col
age_lookup_col_m[age_lookup_col == 7] <- 6

n_age_lookup_col <- length(age_lookup_col)
period_lookup_col <- c(rep(16:20, each=12), rep(21,floor(as.duration( ymd("2022-01-01") %--% death_end)/dmonths(1))+1))
#January 2017 - May 2022
n_period_lookup_col <- length(period_lookup_col)
##########################################
###
### age class indexing and period indexing
### for infected at capture
###
##########################################

age_lookup_col_inf_f  <- age_lookup_col_inf_m <- age_lookup_col_inf <- age_lookup_col
age_lookup_col_inf_m[age_lookup_col_inf_m == 7] <- 6
n_age_lookup_col_inf <- length(age_lookup_col_inf)

# pre_study_weeks <- max(d_fit_icap_mort$left_age_e) - d_fit_icap_mort$left_period_e[which.max(d_fit_icap_mort$left_age_e)] + 1
# difftime("2002-01-01","2019-02-15",units="weeks")
# 187 months before "2017-01-01"  
#= 15 years and 6 months prior to start of collar study
# 2002-01-01 is when CWD was first detected in wisconsin surveillance data
# as.POSIXct("2019-02-15") - years(17)
# period_lookup_col_inf <- c(rep(1, 6),rep(1:20, each = 12), rep(21,5))
# n_period_lookup_col_inf <- length(period_lookup_col_inf)

#repeat for longer due to the really old ass individual
# ceiling(difftime("2002-01-01","2001-05-25",units="weeks"))
period_lookup_col_inf <- c(rep(1,floor(as.duration( ymd("1994-05-15") %--%  ymd("2002-01-01"))/dmonths(1))+1),
                        #    rep(1,tot_pre_study-round(difftime("2017-01-09","2002-01-01",units="weeks"))),
                           rep(1:20, each = 12), rep(21,5))# through the end of may 2022
n_period_lookup_col_inf <- length(period_lookup_col_inf)


##############################################
###
### To run survival period effects prior to
###  the start of the study
###
##############################################


period_effect_survival <- c(rep(-6,4),#may1992-sep1992
                            rep(c(rep(-4.5,4),rep(-6,8)),27),#sep1992 - sep2021
                            rep(-4.5,5),#sep2021-jan2022
                            rep(-6,5))#jan2022-May15,2022
length(period_effect_survival)
nT_overall

period_effect_survival[(nT_period_presurv+1):(nT_overall)] <- NA
length(is.na(period_effect_survival))