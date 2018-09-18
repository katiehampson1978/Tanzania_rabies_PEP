################################################################################
#                            PRODUCE FIGURES/TABLES                            #
################################################################################

rm(list=ls())

# Load in libraries
library(dplyr)
library(ggplot2)
library(lmerTest)
library(lme4)

# Load functions
source("R/mixture_model.R")
source("R/TZ_p_infect.from.bite.loc.R")
source("R/TZ_p_multiple.bites.R")

# Set palette for all figures
fig_palette <- c("#0072B2", "#E69F00")

################
# 1. Load data #
################

CT_data <- read.csv("output/CT_data.csv", stringsAsFactors=FALSE)
CT_data_rabid <- read.csv("output/CT_data_rabid.csv", stringsAsFactors=FALSE)
HC_bite_summary <- read.csv("output/HC_bite_summary.csv", stringsAsFactors=FALSE)
infect_results <- read.csv("output/p_infect.csv", stringsAsFactors=FALSE)
model_results <- read.csv("output/p_infect_model.csv", stringsAsFactors=FALSE)
HC_data <- read.csv("output/HC_data.csv", stringsAsFactors=FALSE)

# Limit CT to Serengeti & Ngorongoro
ser_ngor <- CT_data[which(CT_data$District %in% c("Serengeti", "Ngorongoro")),]
ser_ngor_rabid <- ser_ngor[which(ser_ngor$Rabid=="Yes"),]

###################################################################
# Figure 1: Variation in annual incidence of patient presentation #
###################################################################

# Set the df order and District factor by dog population
boxplot_data <- HC_bite_summary[order(HC_bite_summary$Dogs),]
row.names(boxplot_data) <- NULL
levsord <- unique(boxplot_data$District)
boxplot_data$District <- factor(boxplot_data$District, levels = levsord, ordered = TRUE)
boxplot_data <- boxplot_data[with(boxplot_data, order(HDR)), ]
boxplot_data$District <- factor(boxplot_data$District, levels = unique(boxplot_data$District))

# Produce Figure 1 - Option 1 (District on x-axis)
pdf("figs/Bite_incidence_Boxplot_1.pdf", width=10, height=8)
ggplot(boxplot_data, aes(x=District, y=bites_per_pop_per_100000, fill=Setting)) +
  theme_classic()+
  geom_boxplot(position=position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(shape=21, size=2, position=position_dodge(0.9)) +
  stat_summary(fun.y=mean, geom="point", size = 2, colour="black",position=position_dodge(0.9)) +
  labs(x = "District", y="Incidence of bite patients per 100,000 & \n probable rabid bites (Serengeti & Ngorongoro Districts) per 100,000") +
  scale_fill_manual(name="Setting", values=fig_palette) + ylim(0,165) +
  geom_segment(aes(x=1, y=160, xend=length(unique(boxplot_data$District)), yend=160), arrow = arrow(ends="both", length=unit(0.5, "cm"))) +
  annotate("text", x=3.5, y=155, label=paste0("Low HDR (", round(min(boxplot_data$HDR), digits=1), ")")) +
  annotate("text", x=length(unique(boxplot_data$District))-2.5, y=155, label=paste0("High HDR (", round(max(boxplot_data$HDR), digits=1), ")")) #+
dev.off()

###############################################################
# Table 1: Probable rabies exposures and deaths for bite site #
###############################################################

# Use mixture model to calculate overall probability of infection
quantiles_list <- get_quantiles(model_results)

# Produce table
final_infect_df <- data.frame(Bite_loc = infect_results$Bite.loc,
                              Order = infect_results$p_deaths,
                              P_death = paste0(round(infect_results$p_deaths, digits=3), " (", round(infect_results$P_death_LCI, digits=3), "-", round(infect_results$P_death_UCI, digits=3), ")"),
                              N_death = infect_results$n_deaths,
                              N_no_pep = infect_results$n_no.pep,
                              P_bite = quantiles_list$p_bite_site,
                              N_rabid_bites = infect_results$n_bite.from.rabid.animal,
                              stringsAsFactors = FALSE)
quantiles_list <- get_quantiles(model_results)
new_row <- data.frame(Bite_loc = "", Order = "", P_death = quantiles_list$p_infect_overall,
                      N_death = sum(final_infect_df$N_death),
                      N_no_pep = sum(final_infect_df$N_no_pep), P_bite = "",
                      N_rabid_bites = sum(final_infect_df$N_rabid_bites),
                      stringsAsFactors = FALSE)
final_infect_df <- rbind(final_infect_df, new_row)
final_infect_df <- final_infect_df[order(final_infect_df$Order, decreasing=TRUE),]
final_infect_df$Order <- NULL
colnames(final_infect_df) <- c("Bite location",
                               "Probability of death (95% CI)", "Number of deaths",
                               "Number of bite victims that did not receive PEP",
                               "Probability of bite site (95% CI)", "Number of suspect rabies exposures")

# Save output
write.csv(final_infect_df, "output/paper/Probablity_of_death.csv", row.names=FALSE)

###########################################################################
# Table 2: Details of human deaths where late/incomplete PEP was received #
###########################################################################

# Create df of details on which bite victims died following late/incomplete pep
deaths_details <- CT_data_rabid[which(CT_data_rabid$Patient.outcome=="Died"),]
deaths_details <- deaths_details[which(deaths_details$complete=="incomplete" |  deaths_details$late.PEP=="late"),]

deaths_details <- data.frame(PEP_Failure=NA, Region=deaths_details$Region, Age=deaths_details$Age..in.years.,
                             Sex=deaths_details$Sex, PEP=deaths_details$PEP.administration,
                             n_pep_received=deaths_details$n.PEP.recieved, delay_to_PEP.1=deaths_details$delay.pep.1,
                             n_bites=deaths_details$n.bite.locs, Head=deaths_details$Head, Arm=deaths_details$Arm,
                             Hand=deaths_details$Hand, Trunk=deaths_details$Trunk, Leg=deaths_details$Leg,
                             Foot=deaths_details$Foot, Bite_loc=NA)
for(i in 1:nrow(deaths_details)){
  bite_loc_string <- c()
  if(deaths_details$Head[i]==1){
    bite_loc_string <- c(bite_loc_string, "Head")
  }
  if(deaths_details$Arm[i]==1){
    bite_loc_string <- c(bite_loc_string, "Arm")
  }
  if(deaths_details$Hand[i]==1){
    bite_loc_string <- c(bite_loc_string, "Hand")
  }
  if(deaths_details$Trunk[i]==1){
    bite_loc_string <- c(bite_loc_string, "Trunk")
  }
  if(deaths_details$Leg[i]==1){
    bite_loc_string <- c(bite_loc_string, "Leg")
  }
  if(deaths_details$Foot[i]==1){
    bite_loc_string <- c(bite_loc_string, "Foot")
  }
  deaths_details$Bite_loc[i] <- toString(bite_loc_string)
}

# Set PEP Failure & PEP route
deaths_details$PEP <- as.character(deaths_details$PEP)
for(i in 1:nrow(deaths_details)){
  if(deaths_details$n_pep_received[i]>2 & deaths_details$delay_to_PEP.1[i]>0){
    deaths_details$PEP_Failure[i] <- "Delayed"
  } else if(deaths_details$n_pep_received[i]<3 & deaths_details$delay_to_PEP.1[i]==0){
    deaths_details$PEP_Failure[i] <- "Incomplete"
  } else if(deaths_details$n_pep_received[i]<3 & deaths_details$delay_to_PEP.1[i]>0){
    deaths_details$PEP_Failure[i] <- "Delayed & Incomplete"
  }

  if(deaths_details$Region[i] %in% c("Mara", "Arusha")){
    deaths_details$PEP[i] <- "IM"
  } else if(is.na(deaths_details$PEP[i]) | deaths_details$PEP[i]=="IV"){
    deaths_details$PEP[i] <- "ID"
  }
}

deaths_details <- dplyr::select(deaths_details, PEP_Failure, Age, Sex, n_pep_received,
                                PEP, delay_to_PEP.1, n_bites, Bite_loc)

# Sort by Number of PEP doses, then number of days delay
deaths_details <- deaths_details[with(deaths_details, order(-n_pep_received, -delay_to_PEP.1)), ]

# Save output
write.csv(deaths_details, "output/paper/Death_details.csv", row.names=FALSE)

##################################################
# Figure 2A: Delay between bite date & PEP 1 Hist #
##################################################

# Remove any remaining delays outside of range
summary(CT_data_rabid$delay.pep.1)
barplot_df <- CT_data_rabid[which(CT_data_rabid$delay.pep.1>=0 & CT_data_rabid$delay.pep.1<91),]

# Read any "Other" patient outcomes as "Fine"
barplot_df$Patient.outcome[which(barplot_df$Patient.outcome=="Other")] <- "Fine"

# Seperate off NTz and STz
NTz_barplot <- barplot_df[which(barplot_df$Study.Area=="NTz"),]; nrow(NTz_barplot)
STz_barplot <- barplot_df[which(barplot_df$Study.Area %in% c("STz", "Pemba")),]; nrow(STz_barplot)

# NTZ: Create dataframe and transform to matrix
NTz_delays <- hist(NTz_barplot$delay.pep.1, breaks=c(seq(-1, 14, by=1), 21, 35, 64, 84))
NTz_barplot.df <- data.frame(Delay=NTz_delays$breaks[2:20],
                             N_Fine=NTz_delays$counts,
                             Location=rep("NTz", 19))
NTz_barplot.df$Prop = NTz_barplot.df$N_Fine/sum(NTz_barplot.df$N_Fine)
NTz_barplot.df$Delay <- factor(NTz_barplot.df$Delay, levels=sort(unique(NTz_barplot.df$Delay)))

STz_delays <- hist(STz_barplot$delay.pep.1, breaks=c(seq(-1, 14, by=1), 21, 35, 64, 84))
STz_barplot.df <- data.frame(Delay=STz_delays$breaks[2:20],
                             N_Fine=STz_delays$counts,
                             Location=rep("STz", 19))
STz_barplot.df$Prop = STz_barplot.df$N_Fine/sum(STz_barplot.df$N_Fine)
STz_barplot.df$Delay <- factor(STz_barplot.df$Delay, levels=sort(unique(STz_barplot.df$Delay)))

# Bind together
barplot.df <- rbind(NTz_barplot.df, STz_barplot.df)

pdf("figs/delay_in_pep_V3.pdf", width=10, height=6)
ggplot(data=barplot.df, aes(x=Delay, y=Prop, fill=Location)) +
  geom_col(position=position_dodge(), color="black") + theme_classic() + ylim(0,1) +
  labs(x="Delay between exposure and initiation of PEP",
       y="Proportion that presented at a health facility") +
  scale_fill_manual(name="", breaks=c("NTz", "STz"), labels=c("Serengeti & Ngorongoro", "Southern Tanzania"),
                    values=fig_palette) +
  theme(legend.position="top",
        text = element_text(size=12),
        axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_x_discrete(breaks=c(0,7,14,21,35,64,84),
                   labels=c("0 days", "7 days", "14 days", "2-3 weeks", "3-4 weeks", "1-2 months", ">2 months")) +
  geom_text(label="A)", aes(x=-Inf, y=Inf, hjust=-1, vjust=1), size=7)
dev.off()

#####################################################
# Figure 2B - Barplot of proportion doses receieved #
#####################################################

# Barplots of patients that received PEP dose 1 (as a proportion, 95%; include
# clinical signs as no PEP), then dose 2, then 3, then 4, then 5 - for Mobile
# surveillance (panel 1) and for contact tracing rabies exposures (SD/ND where
# people pay), panel 2.
HC_data_sub <- HC_data[which(HC_data$visit_status>0),]
hc_visits <- HC_data_sub %>%
  group_by(visit_status) %>%
  summarise(n_bite_victims = length(which(is_antirabies_available=="yes")))
hc_visits$p_receive <- hc_visits$n_bite_victims/length(which(HC_data_sub$visit_status==1))
hc_visits$PEP_available <- "yes"
hc_visits$data <- "MS"

CT_barplot <- ser_ngor_rabid[which(ser_ngor_rabid$Sought==1),]
ct_doses <- data.frame(visit_status=1:5,
                       n_bite_victims=c(length(which(CT_barplot$n.PEP.recieved>0)),
                                        length(which(CT_barplot$n.PEP.recieved>1 & CT_barplot$PEP.1==1)),
                                        length(which(CT_barplot$n.PEP.recieved>2 & CT_barplot$PEP.1==1)),
                                        length(which(CT_barplot$n.PEP.recieved>3 & CT_barplot$PEP.1==1)),
                                        length(which(CT_barplot$n.PEP.recieved>4 & CT_barplot$PEP.1==1))))
ct_doses$p_receive <- ct_doses$n_bite_victims/nrow(CT_barplot)
ct_doses$PEP_available <- "yes"
ct_doses$data <- "CT"

barplot_data <- rbind(hc_visits, ct_doses)
barplot_data$data <- factor(barplot_data$data, levels = c("CT", "MS"))

pdf("figs/Proportion_PEP_completion.pdf", width=10, height=6)
ggplot(data=barplot_data, aes(x=visit_status, y=p_receive, fill=data)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_classic() + ylim(0,1) +
  scale_fill_manual(name="", values=fig_palette,
                    breaks=c("CT", "MS"), labels=c("Serengeti & Ngorongoro", "Southern Tanzania")) +
  labs(x="PEP doses received", y="Proportion of patients that received PEP") +
  theme(legend.position="top",
        text = element_text(size=12)) +
  geom_text(label="B)", aes(x=-Inf, y=Inf, hjust=-1, vjust=1), size=7)
dev.off()

###############################
# Create output for schematic #
###############################

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

# Create dataframe
one_n_same_district = length(which(one_clinic$Correct_patient_district==one_clinic$Correct_clinic_district & one_clinic$Correct_patient_region==one_clinic$Correct_clinic_region))
one_n_same_region = length(which(one_clinic$Correct_patient_district!=one_clinic$Correct_clinic_district & one_clinic$Correct_patient_region==one_clinic$Correct_clinic_region))
one_n_diff_region = length(which(one_clinic$Correct_patient_region!=one_clinic$Correct_clinic_region))
one_clinic_summary <- data.frame(category="1 clinic", n_total=nrow(one_clinic),
                                 n_inside_home_dis=one_n_same_district, per_inside_home_dis=round((one_n_same_district/n_total_visits)*100, digits=1),
                                 n_outside_home_dis=one_n_same_region, per_outside_home_dis=round((one_n_same_region/n_total_visits)*100, digits=1),
                                 n_outisde_home_reg=one_n_diff_region, per_outside_home_reg=round((one_n_diff_region/n_total_visits)*100, digits=1))
mult_n_same_district = length(which(multiple_clin$Correct_patient_district==multiple_clin$Correct_clinic_district & multiple_clin$Correct_patient_region==multiple_clin$Correct_clinic_region))
mult_n_same_region = length(which(multiple_clin$Correct_patient_district!=multiple_clin$Correct_clinic_district & multiple_clin$Correct_patient_region==multiple_clin$Correct_clinic_region))
mult_n_diff_region =  length(which(multiple_clin$Correct_patient_region!=multiple_clin$Correct_clinic_region))
mult_clinic_summary <- data.frame(category=">1 clinic", n_total=nrow(multiple_clin),
                                  n_inside_home_dis=mult_n_same_district, per_inside_home_dis=round((mult_n_same_district/n_total_visits)*100, digits=1),
                                  n_outside_home_dis=mult_n_same_region, per_outside_home_dis=round((mult_n_same_region/n_total_visits)*100, digits=1),
                                  n_outisde_home_reg=mult_n_diff_region, per_outside_home_reg=round((mult_n_diff_region/n_total_visits)*100, digits=1))
overall_clinic_summary <- rbind(one_clinic_summary, mult_clinic_summary)

# Save output
write.csv(overall_clinic_summary, "output/paper/clinic_visits_summary.csv", row.names=FALSE)
