#### Prepare data ----

# 1. YSA breeding biology data
# data/Life History 2006-2014.xlsx - this is YSA breeding biology data from Bonaire
# R/YSA_life_history_data.R - script which read and summarises the YSA breeding biology data

source("R/YSA_life_history_data.R")

# 2. Combine with estimates of survival and stage duration

# Manually creating data frame for now
stage <- c("1", "2", "3")

# Stage classes
class <- c("fledgling", "juvenile", "adult")

# Stage durations
# To age 12 months, Age 13-36 months,
# Age 37 months+ (as 10 years -> needs updating and a reference)
duration <- c(1, (24/12), 10)

# 0.71 from Salinas-Melgoza & Renton 2007 (CHECK), 0.875 from Rodriguez et al 2004
survival <- c(0.71, 0.875, 0.875)

# SD/standard error estimates from YSA_demog_data_master.csv (CHECK), 0.2 from
# Salinas-Melgoza & Renton 2007, 0.075 from Rodriguez et al 2004
survival_SD <- c(0.2, 0.075, 0.075)

# Reproductive output/fecundity
# Clutch size taken as half of 3.2 as sex ratio assumed 1:1 (Sam's thesis), data
# from YSA_demog_data_master_csv
fecundity <- c(0, 0, 1.6*total_summary$mean_hatch[1]*total_summary$mean_nestling_surv[1])

# Reproductive output/fecundity SEs
# SE Sam's thesis (halved) TODO modify to incorporate SEs for hatch and nestling surv:
# total_summary$se_hatch[1], total_summary$se_nestling_surv[1]
fecundity_SD <- c(0, 0, 0.1)

#### Build matrix model ----

source("R/make_projection_matrix.R")

ysa <- make_projection_matrix(survival, fecundity, duration)

lambda <- eigen(ysa)$values[1]

#### Life stage simulation analysis ----

set.seed(2021)

# Sample from vital rates

