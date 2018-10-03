################################################################################
#                                 PRINT RESULTS                                #
################################################################################

rm(list=ls())

# Load in libraries
library(dplyr)
library(ggplot2)
library(lmerTest)
library(lme4)
library(Hmisc)

# Load functions
source("R/mixture_model.R")
source("R/TZ_p_multiple.bites.R")

Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

################
# 1. Load data #
################

CT_totals <- read.csv("output/CT_totals.csv", stringsAsFactors=FALSE)
CT_data <- read.csv("output/CT_data.csv", stringsAsFactors=FALSE)
CT_data_rabid <- read.csv("output/CT_data_rabid.csv", stringsAsFactors=FALSE)
HC_bite_patient_summary <- read.csv("output/HC_bite_summary_patients.csv", stringsAsFactors=FALSE)
HC_rab_exposure_summary <- read.csv("output/HC_bite_summary_exposures.csv", stringsAsFactors=FALSE)
HC_data <- read.csv("output/HC_data.csv", stringsAsFactors=FALSE)
HC_outside_study_area <- read.csv("output/HC_outside_study.csv", stringsAsFactors=FALSE)
HC_distance_data <- read.csv("output/HC_dist_data.csv", stringsAsFactors=FALSE)
infect_results <- read.csv("output/p_infect.csv", stringsAsFactors=FALSE)
model_results <- read.csv("output/p_infect_model.csv", stringsAsFactors=FALSE)

# Create a subset of CT data, containing only Serengeti and Ngorongoro (NTz)
ser_ngor <- CT_data[which(CT_data$District %in% c("Serengeti", "Ngorongoro")),]
ser_ngor_rabid <- ser_ngor[which(ser_ngor$Rabid=="Yes"),]

###############
# 2. Abstract #
###############

paste0("Number of probably rabies exposures identified: ", format(CT_totals$n_rabid, big.mark=","))
paste0("Number of MS records: ", format(HC_outside_study_area$n_total_HC, big.mark=","))

# Proportion of CT rabid dog exposures that do not seek healthcare
n_rabid = nrow(ser_ngor_rabid)
n_seek = length(which(ser_ngor_rabid$Sought==1))
paste0("Percentage of rabid bite exposures (CT) that do not seek healthcare: ",
       round((1-(n_seek/n_rabid))*100, digits=0), "% (", n_rabid-n_seek, "/", n_rabid, ")")

# Percentage of CT individuals that sought care and did not receive PEP
n_receive = length(which(ser_ngor_rabid$n.PEP.recieved>0))
paste0("Percentage of rabid bite exposures (CT) that do not receieve PEP having sought it: ",
       round((1-(n_receive/n_seek))*100, digits=0), "% (", n_receive, "/", n_seek, ")")

n_complete = length(which(ser_ngor_rabid$n.PEP.recieved>2))
paste0("Percentage of rabid bite exposures (CT) that did not complete having initiated: ",
       round((1-(n_complete/n_receive))*100, digits=0), "% (", n_complete, "/", n_receive, ")")

##############
# 3. Methods #
##############

paste0("Timeframe for CT in NTZ: ", min(CT_data$Year.bitten[which(CT_data$Study.Area=="NTz")]), "-", max(CT_data$Year.bitten[which(CT_data$Study.Area=="NTz")]))

unique(CT_data$District[which(CT_data$Study.Area %in% c("STz", "Pemba"))])
paste0("Number of south TZ districts in CT: ", length(unique(CT_data$District[which(CT_data$Study.Area %in% c("STz", "Pemba"))]))+3) # Pemba districts grouped as 1
paste0("Timeframe for CT in STZ: ", min(CT_data$Year.bitten[which(CT_data$Study.Area %in% c("STz", "Pemba"))]), "-", max(CT_data$Year.bitten[which(CT_data$Study.Area %in% c("STz", "Pemba"))]))

# DATA COLLECTION
paste0("HC: Number of districts: ", length(unique(HC_data$Correct_clinic_district))+3) # Pemba districts grouped
paste0("HC: Number of regions: ", length(unique(HC_data$Correct_clinic_region))+1) # Pemba regions combined
paste0("HC: Timeframe: ", min(HC_data$year), "-", max(HC_data$year))

paste0("HC dates: ", min(HC_data$month_year), "-", max(HC_data$month_year))
paste0("HC number of records: ", format(HC_outside_study_area$n_total_HC, big.mark=","), " (all records), ", format(nrow(HC_data), big.mark=","), " (study area only).")

paste0("CT Number of people presenting to health facility: ", format(CT_totals$n_hosp_presentations, big.mark=","))
paste0("CT Number of people bitten by a suspect rabid animal: ", format(CT_totals$n_rabid, big.mark=","))
paste0("CT Number of people bitten by a suspect rabid animal that did not seek PEP: ", format(CT_totals$n_non_report_rabid, big.mark=","))

paste0("CT years taken for Bite incidence: ", min(HC_bite_patient_summary$year[which(HC_bite_patient_summary$District %in% c("Serengeti", "Ngorongoro"))]),
       "-", max(HC_bite_patient_summary$year[which(HC_bite_patient_summary$District %in% c("Serengeti", "Ngorongoro"))]))

paste0("CT Number of people bitten by a suspect rabid animal: ", format(CT_totals$n_rabid, big.mark=","))
paste("Number of people that received RIG: ", length(which(CT_data$Immunoglobulin==TRUE)))
paste0("Number of people that died from tetanus (n=", length(which(CT_data$Rabid=="Yes" & CT_data$Cause.of.death=="Tetanus")), ") or injury (n=", length(which(CT_data$Rabid=="Yes" & CT_data$Cause.of.death=="Injuries")), ")")

# Run function to estimate number of individuals with multiple bite wounds
multiple.bite.df <- TZ_p_multiple.bites(CT_data_rabid)
paste0("Number of people with multiple bite wounds: ", multiple.bite.df$n_multiple.bites, "/", multiple.bite.df$n_bite.victims)

paste0("Sample size for delays: ", length(which(!is.na(CT_data_rabid$late.PEP))))
paste0("Sample size completion in SD/ND: ", nrow(ser_ngor_rabid))
paste0("Sample size completion in STz: ", nrow(CT_data_rabid[which(CT_data_rabid$Study.Area %in% c("STz", "Pemba")),]))

##############
# 4. Results #
##############

#########################################################
# Incidence of bite-injury patients and probable rabies exposures
incidence_aver <- HC_bite_patient_summary %>%
  group_by(District) %>%
  summarise(mean_inc = mean(bites_per_pop_per_100000))

paste0("Average range of Incidence of bite patients presenting to health facilities across districts in Tanzania: ",
       round(min(incidence_aver$mean_inc), digits=2), "/100,000 to ",
       round(max(incidence_aver$mean_inc), digits=2), "/100,000.")

# Need to rescale (standardise) the continuous explanatory variables
HC_model_df <- cbind(HC_bite_patient_summary,
                     "R_estimatedDogs"=(HC_bite_patient_summary$est_dogs-mean(HC_bite_patient_summary$est_dogs))/sd(HC_bite_patient_summary$est_dogs),
                     "R_dog_density"=(HC_bite_patient_summary$dog_density-mean(HC_bite_patient_summary$dog_density))/sd(HC_bite_patient_summary$dog_density),
                     "R_human_dog_ratio"=(HC_bite_patient_summary$HDR-mean(HC_bite_patient_summary$HDR))/sd(HC_bite_patient_summary$HDR))

model_b <- lmer(bites_per_pop_per_100000 ~ R_human_dog_ratio + Setting +
                  (1|year) + (1|District), data=HC_model_df)
AIC(model_b)
model_coef <- coef(summary(model_b))
paste0("Summary coefficients for HDR: Estimate=", round(model_coef[2,1], digits=3),
       "; Pr=", round(model_coef[2,5], digits=5), "; Std. error=", round(model_coef[2,2], digits=3))
paste0("Summary coefficients for Setting: Estimate=", round(model_coef[3,1], digits=3),
       "; Pr=", round(model_coef[3,5], digits=5), "; Std. error=", round(model_coef[3,2], digits=3))

# Print probable rabies exposures for NTz
paste0("Probable rabies exposures for Serengeti: ", round(mean(HC_rab_exposure_summary$bites_per_pop_per_100000[which(HC_rab_exposure_summary$District=="Serengeti")]), digits=1), "/100,000 (min=",
       round(min(HC_rab_exposure_summary$bites_per_pop_per_100000[which(HC_rab_exposure_summary$District=="Serengeti")]), digits=1), "/100,000; max=",
       round(max(HC_rab_exposure_summary$bites_per_pop_per_100000[which(HC_rab_exposure_summary$District=="Serengeti")]), digits=1), "/100,000)")
paste0("Probable rabies exposures for Ngorongoro: ", round(mean(HC_rab_exposure_summary$bites_per_pop_per_100000[which(HC_rab_exposure_summary$District=="Ngorongoro")]), digits=1), "/100,000 (min=",
       round(min(HC_rab_exposure_summary$bites_per_pop_per_100000[which(HC_rab_exposure_summary$District=="Ngorongoro")]), digits=1), "/100,000; max=",
       round(max(HC_rab_exposure_summary$bites_per_pop_per_100000[which(HC_rab_exposure_summary$District=="Ngorongoro")]), digits=1), "/100,000)")

# Print HDRs for NTZ
paste0("HDR for Serengeti: ", round(mean(HC_bite_patient_summary$HDR[which(HC_bite_patient_summary$District=="Serengeti")]), digits=1))
paste0("HDR for Ngorongoro: ", round(mean(HC_bite_patient_summary$HDR[which(HC_bite_patient_summary$District=="Ngorongoro")]), digits=1))

# Percentage of CT bite victims that presented to health facility following probable exposure
n_sought = length(which(CT_data$Sought==1))
n_rabid_sought = length(which(CT_data$Sought==1 & CT_data$Rabid=="Yes"))
paste0("ALL DISTRICTS: Percentage of health seeking bite victims bitten by a probable rabid animal: ",
       round((n_rabid_sought/n_sought)*100, digits=0), "% (", n_rabid_sought, "/", n_sought, ")")

# Percentage of CT rabid dog exposures that do not seek healthcare (NTz)
n_rabid = nrow(ser_ngor_rabid)
n_seek = length(which(ser_ngor_rabid$Sought==1))
paste0("Percentage of rabid bite exposures (CT) that do not seek healthcare: ",
       round(100-((n_seek/n_rabid)*100), digits=0), "% (", n_rabid-n_seek, "/", n_rabid, ")")

#########################################################
# Rabies risk and effectiveness of PEP

# Calculate overall risk of infection following a probable rabid bite
quantiles_list <- get_quantiles(model_results)
paste0("Overall probability of infection: ", quantiles_list$p_infect_overall)

# Sites with greatest risk of rabies infection
infect_results <- infect_results[order(infect_results$p_deaths, decreasing=TRUE),]
paste0("Highest risk of infection: ", infect_results$Bite.loc[1], ", p=", round(infect_results$p_deaths[1], digits=3),
       " (95% CI: ", round(infect_results$P_death_LCI[1], digits=3), "-", round(infect_results$P_death_UCI[1], digits=3), ")")
paste0("Second highest risk of infection: ", infect_results$Bite.loc[2], ", p=", round(infect_results$p_deaths[2], digits=3),
       " (95% CI: ", round(infect_results$P_death_LCI[2], digits=3), "-", round(infect_results$P_death_UCI[2], digits=3), ")")

# Limit to individuals that have had traceback completed
traced_CT_data <- CT_data_rabid[which(CT_data_rabid$Traceback=="Yes"),]

# Find individuals that completed timely and complete PEP
time_compl_pep <- traced_CT_data[which(traced_CT_data$late.PEP=="timely" & traced_CT_data$complete=="complete"),]
paste0("N bite victims with effective PEP administration (timely & complete): ", nrow(time_compl_pep))
p_prevent_compl <- binconf(nrow(time_compl_pep), nrow(time_compl_pep))
paste0("Probability of preventing Rabies given complete and timely PEP: p=", round(p_prevent_compl[1], digits=3), " (", round(p_prevent_compl[2], digits=3), "-", round(p_prevent_compl[3], digits=3), ")")

# Find individuals that receieved late or incomplete PEP
late_incomp_pep <- traced_CT_data[which(traced_CT_data$late.PEP=="late" | traced_CT_data$complete=="incomplete"),]
paste0("N bite victims with incomplete and/or late PEP administration: ", format(nrow(late_incomp_pep), big.mark=","))
paste0("Number that died with incomplete and/or late PEP: ", length(which(late_incomp_pep$Patient.outcome=="Died")))
p_prevent_incompl <- binconf(length(which(late_incomp_pep$Patient.outcome=="Died")), nrow(late_incomp_pep))
paste0("Probability of developing Rabies given incomplete and/or late PEP: p=", round(p_prevent_incompl[1], digits=3), " (", round(p_prevent_incompl[2], digits=3), "-", round(p_prevent_incompl[3], digits=3), ")")
paste0("Probability of preventing Rabies given incomplete and/or late PEP: p=", round(1-p_prevent_incompl[1], digits=3), " (", round(1-p_prevent_incompl[3], digits=3), "-", round(1-p_prevent_incompl[2], digits=3), ")")

paste0("Number of deaths attributable to delays in PEP: ", length(which(late_incomp_pep$Patient.outcome=="Died" & late_incomp_pep$late.PEP=="late")))
paste0("Number of deaths attributable to incomplete PEP with no delay: ", length(which(late_incomp_pep$Patient.outcome=="Died" & late_incomp_pep$late.PEP=="timely" & late_incomp_pep$complete=="incomplete")))

#########################################################
# Health seeking, PEP access and provision

# Proportion of CT rabid dog exposures that do not seek healthcare
n_rabid = nrow(ser_ngor_rabid)
n_seek = length(which(ser_ngor_rabid$Sought==1))
paste0("Percentage of rabid bite exposures (CT) that do not seek healthcare: ",
       round((1-(n_seek/n_rabid))*100, digits=1), "% (", n_rabid-n_seek, "/", n_rabid, ")")

# Percentage of CT individuals that sought care and did not receive PEP
n_receive = length(which(ser_ngor_rabid$PEP.1==1))
paste0("Percentage of rabid bite exposures (CT) that do not receieve PEP having sought it: ",
       round((1-(n_receive/n_seek))*100, digits=1), "% (", n_seek-n_receive, "/", n_seek, ")")

# Average delay in PEP administration
STz_rabid <- CT_data_rabid[which(CT_data_rabid$Study.Area %in% c("STz", "Pemba")),]
paste0("Most common delay length: ", Mode(STz_rabid$delay.pep.1))
quantile(ser_ngor_rabid$delay.pep.1[which(ser_ngor_rabid$n.PEP.recieved>0)], na.rm=TRUE)[3]
quantile(STz_rabid$delay.pep.1[which(STz_rabid$n.PEP.recieved>0)], na.rm=TRUE)[3]

# Probability of patients initiating PEP
n_sought = length(which(HC_data$visit_status==1))
n_receive = length(which(HC_data$is_antirabies_available=="yes" & HC_data$visit_status==1))
paste0("HC: Proportion of health seeking bite victims that initiated PEP: ",
       round(n_receive/n_sought, digits=3), " (", n_receive, "/", n_sought, ")")
n_sought = length(which(ser_ngor_rabid$Sought==1 & ser_ngor_rabid$Patient.outcome!="Died"))
n_receive = length(which(ser_ngor_rabid$Sought==1 & ser_ngor_rabid$PEP.1==1))
paste0("CT: Proportion of probable rabid bite victims that initiated PEP: ",
       round(n_receive/n_sought, digits=3), " (", n_receive, "/", n_sought, ")")

# PEP not obtained due to shortages - for HC, ignore 0 PEP doses (positive clinical signs)
n_sought = length(which(HC_data$visit_status>0))
n_shortage = length(which(HC_data$is_antirabies_available=="no" & HC_data$visit_status>0))
paste0("HC: Percentage of health seeking bite victims that did not obtain PEP due to shortage: ",
       round((n_shortage/n_sought)*100, digits=0), "% (", n_shortage, "/", n_sought, ")")

# Proportion of individuals that completed PEP having initiated
n_receive = length(which(HC_data$is_antirabies_available=="yes" & HC_data$visit_status==1))
n_complete = length(which(HC_data$is_antirabies_available=="yes" & HC_data$visit_status==3))
paste0("HC: Proportion of health seeking bite victims that completed PEP: ",
       round(n_complete/n_receive, digits=3), " (", n_complete, "/", n_receive, ")")
n_receive = length(which(ser_ngor_rabid$Sought==1 & ser_ngor_rabid$PEP.1==1))
n_complete = length(which(ser_ngor_rabid$Sought==1 & ser_ngor_rabid$n.PEP.recieved>2))
paste0("CT: Proportion of probable rabid bite victims that completed PEP: ",
       round(n_complete/n_receive, digits=3), " (", n_complete, "/", n_receive, ")")

# Remove "sub system admin" - a Clinic name that does has not been matched
HC_sub <- HC_data[which(HC_data$clinic_name!="sub system admin"),]

# Remove patients that did not receive PEP
HC_sub <- HC_sub[which(HC_sub$visit_status!=0 & HC_sub$is_antirabies_available!="no"),]

# Collect unique patient trips i.e. 1 row for a patient that visited 1 facility multiple times,
# vs. 3 rows for a patient that visited 3 facilities
unique_trips <- unique(dplyr::select(HC_sub, Patient_id,
                                    Correct_patient_region, Correct_patient_district,
                                    Correct_clinic_region, Correct_clinic_district, clinic_name))

# Summarise number of facilities visited by each individuals
n_facilities_visited <- unique_trips %>%
  group_by(Patient_id) %>%
  summarise(n_facilities=length(Patient_id))

# Collect total number of visits where home and clinic location info is available
length(which(is.na(HC_sub$Correct_patient_district) | is.na(HC_sub$Correct_patient_region) | is.na(HC_sub$Correct_clinic_district) | is.na(HC_sub$Correct_clinic_region)))
n_total_visits = nrow(HC_sub)

# For individuals that travelled to one clinic
one_clinic_ind <- n_facilities_visited$Patient_id[which(n_facilities_visited$n_facilities==1)]
one_clinic <- HC_sub[which(HC_sub$Patient_id %in% one_clinic_ind),]

# For individuals that travelled to multiple clinics
multiple_clin_ind <- n_facilities_visited$Patient_id[which(n_facilities_visited$n_facilities>1)]
multiple_clin <- HC_sub[which(HC_sub$Patient_id %in% multiple_clin_ind),]

# Patients that visited one facility
paste0("Percentage of patient visits to a single facility: ", round((nrow(one_clinic)/n_total_visits)*100, digits=1), "% (", nrow(one_clinic), "/", n_total_visits, ")")

n_travelled_outside_dis = sum(length(which(one_clinic$Correct_patient_district!=one_clinic$Correct_clinic_district & one_clinic$Correct_patient_region==one_clinic$Correct_clinic_region)),
                              length(which(multiple_clin$Correct_patient_district!=multiple_clin$Correct_clinic_district & multiple_clin$Correct_patient_region==multiple_clin$Correct_clinic_region)))
n_travelled_outside_reg = sum(length(which(one_clinic$Correct_patient_region!=one_clinic$Correct_clinic_region)),
                              length(which(multiple_clin$Correct_patient_region!=multiple_clin$Correct_clinic_region)))
paste0("Percentage of patient visits outside of home district: ", round((n_travelled_outside_dis/n_total_visits)*100, digits=1), "% (", format(n_travelled_outside_dis, big.mark=","), "/", format(n_total_visits, big.mark=","), ")")
paste0("Percentage of patient visits outside of home region: ", round((n_travelled_outside_reg/n_total_visits)*100, digits=1), "% (", format(n_travelled_outside_reg, big.mark=","), "/", format(n_total_visits, big.mark=","), ")")

# Patients from outside study area
paste0("Percentage of MS patient records that came from outside of study districts: ",
       round((HC_outside_study_area$n_outside_study_area/HC_outside_study_area$n_total_HC)*100, digits=2), " (",
       format(HC_outside_study_area$n_outside_study_area, big.mark=","), "/", format(HC_outside_study_area$n_total_HC, big.mark=","), ")")

# Calculate actual distance travelled
actual_dist_quant <- round(quantile(HC_distance_data$dist_travelled_km, c(0.025, 0.5, 0.975), na.rm=TRUE), digits=0)
paste0("Patients travelled an average of ", round(mean(HC_distance_data$dist_travelled_km, na.rm=TRUE), digits=0), "km (", actual_dist_quant[1], "-", actual_dist_quant[3], ")")
paste0("Range of distance Patients travelled ", round(min(HC_distance_data$dist_travelled_km, na.rm=TRUE), digits=0), " to ",
       round(max(HC_distance_data$dist_travelled_km, na.rm=TRUE), digits=0), "km")

# Calculate ideal distance travelled
ideal_dist_quant <- round(quantile(HC_distance_data$ideal_dist_travelled_km, c(0.025, 0.5, 0.975), na.rm=TRUE), digits=0)
paste0("Lowest ideal average travelling distance ", round(mean(HC_distance_data$ideal_dist_travelled_km, na.rm=TRUE), digits=0), "km (", ideal_dist_quant[1], "-", ideal_dist_quant[3], ")")

########################
# Caption for figure 2 #
########################

# Remove any remaining delays outside of range
summary(traced_CT_data$delay.pep.1)
hist_cap <- traced_CT_data[which(traced_CT_data$delay.pep.1>=0 & traced_CT_data$delay.pep.1<91),]

paste0("Figure 2: Of ", length(which(hist_cap$late.PEP=="late")),
       " patients with delayed PEP, ", length(which(hist_cap$late.PEP=="late" & hist_cap$Patient.outcome=="Died")),
       " deaths occured.")

########################
# Caption for figure 3 #
########################

HC_barplot_cap <- HC_data[which(HC_data$visit_status>0 & HC_data$is_antirabies_available=="yes"),]
CT_barplot_cap <- ser_ngor_rabid[which(ser_ngor_rabid$Sought==1 & ser_ngor_rabid$Patient.outcome!="Died"),]

paste0("Dark grey indicates patients recorded in mobile phone-based surveillance data from Southern Tanzania (",
       format(nrow(HC_barplot_cap), big.mark=","), ") and light grey indicates rabies exposed patients from Serengeti and Ngorongoro districts (",
       format(nrow(CT_barplot_cap), big.mark=","), ")")
