TZ_p_infect.from.bite.loc <- function(dataframe, by_dis=FALSE, dis_label="District",
                                      rabid_label="Rabid", rabid_true="Yes",
                                      PEP.1_label="PEP.1",
                                      COD_label="Cause.of.death",
                                      head_col="Head", arm_col="Arm",
                                      hand_col="Hand", trunk_col="Trunk",
                                      leg_col="Leg", foot_col="Foot"){

  # load in required libraries
  library(dplyr)

  # Read columns as characters
  dataframe[[COD_label]] <- as.character(dataframe[[COD_label]])

  # Take only selected columns from dataframe
  sub_dataframe <- dataframe[,c(dis_label, rabid_label, PEP.1_label, COD_label, head_col,
                                arm_col, hand_col, trunk_col, leg_col, foot_col)]

  ################################################################################
  #                         GROUP DATA COLUMNS TOGETHER                          #
  ################################################################################

  # Group "arm" and "hand" columns together under "arm"
  for(i in 1:nrow(sub_dataframe)){
    if(sub_dataframe[[hand_col]][i]==1)
      sub_dataframe[[arm_col]][i] <- 1
  }
  sub_dataframe[[hand_col]] <- NULL

  # Group "leg" and "foot" columns together under "leg"
  for(i in 1:nrow(sub_dataframe)){
    if(sub_dataframe[[foot_col]][i]==1)
      sub_dataframe[[leg_col]][i] <- 1
  }
  sub_dataframe[[foot_col]] <- NULL

  #############################################################################
  #     Ignore secondary bites, only taking into account highest risk bite    #
  #############################################################################

  sub_dataframe$n.bites <- NA

  # Calculate which rows have more than one bite location
  for(j in 1:nrow(sub_dataframe)){
    sub_dataframe$n.bites[j] <- sum(sub_dataframe[[head_col]][j], sub_dataframe[[arm_col]][j],
                                    sub_dataframe[[trunk_col]][j], sub_dataframe[[leg_col]][j], na.rm=T)
  }

  # Seperate off rows with 0 or 1 bites and mutliple bites
  normal.bite_rows <- subset(sub_dataframe, n.bites<2)
  multiple.bite_rows <- subset(sub_dataframe, n.bites>1)

  # Heriachy of risk according to Shim et al (2009): head [1], arms[2], legs[5], trunk[4]
  # print(paste0("Method: Shim et al 2009"))
  # print(paste0("Heirarchy: Head > Arms/Hands > Legs/Feet > Trunk"))
  # if(nrow(multiple.bite_rows)>0){
  # 
  #   sub.head <- multiple.bite_rows[multiple.bite_rows[[head_col]]==1,]
  #   multiple.bite_rows <- multiple.bite_rows[multiple.bite_rows[[head_col]]!=1,]
  #   sub.head[[arm_col]] <- 0
  #   sub.head[[trunk_col]] <- 0
  #   sub.head[[leg_col]] <- 0
  #   multiple.bite_rows <- rbind(multiple.bite_rows, sub.head)
  # 
  #   sub.arms <- multiple.bite_rows[multiple.bite_rows[[arm_col]]==1,]
  #   multiple.bite_rows <- multiple.bite_rows[multiple.bite_rows[[arm_col]]!=1,]
  #   sub.arms[[trunk_col]] <- 0
  #   sub.arms[[leg_col]] <- 0
  #   multiple.bite_rows <- rbind(multiple.bite_rows, sub.arms)
  # 
  #   sub.legs <- multiple.bite_rows[multiple.bite_rows[[leg_col]]==1,]
  #   multiple.bite_rows <- multiple.bite_rows[multiple.bite_rows[[leg_col]]!=1,]
  #   sub.legs[[leg_col]] <- 0
  #   multiple.bite_rows <- rbind(multiple.bite_rows, sub.legs)
  # 
  # }

  # Heirachy of risk according to first simulation: head[1], trunk[4], arms[2], legs[5]
  print(paste0("Method: Re-calculated Heirarchy"))
  print(paste0("Heirarchy: Head > Trunk > Arms/Hands > Legs/Feet"))
  if(nrow(multiple.bite_rows)>0){

    sub.head <- multiple.bite_rows[multiple.bite_rows[[head_col]]==1,]
    multiple.bite_rows <- multiple.bite_rows[multiple.bite_rows[[head_col]]!=1,]
    sub.head[[arm_col]] <- 0 # Make all remaining "Arms" false
    sub.head[[trunk_col]] <- 0 # Make all remaining "Trunk" false
    sub.head[[leg_col]] <- 0 # Make all remaining "Legs" false
    multiple.bite_rows <- rbind(multiple.bite_rows, sub.head)

    sub.trunk <- multiple.bite_rows[multiple.bite_rows[[trunk_col]]==1,]
    multiple.bite_rows <- multiple.bite_rows[multiple.bite_rows[[trunk_col]]!=1,]
    sub.trunk[[arm_col]] <- 0 # Make all remaining "Arms" false
    sub.trunk[[leg_col]] <- 0 # Make all remaining "Legs" false
    multiple.bite_rows <- rbind(multiple.bite_rows, sub.trunk)

    sub.arms <- multiple.bite_rows[multiple.bite_rows[[arm_col]]==1,]
    multiple.bite_rows <- multiple.bite_rows[multiple.bite_rows[[arm_col]]!=1,]
    sub.arms[[leg_col]] <- 0 # Make all remaining "Legs" false
    multiple.bite_rows <- rbind(multiple.bite_rows, sub.arms)

  }

  # Bind both dataframes together
  sub_df <- rbind(normal.bite_rows, multiple.bite_rows)

  # Check only one bite is present per row
  for(j in 1:nrow(sub_df)){
    sub_df$n.bites[j] <- sum(sub_df[[head_col]][j], sub_df[[arm_col]][j],
                             sub_df[[trunk_col]][j], sub_df[[leg_col]][j], na.rm=T)
  }

  # If any individual still has more than one bite, stop the function
  if(any(sub_df$n.bites)>1){
    stop("Error combining bites by severity")
  }

  #############################################################################
  #                     Create column stating location of bite                #
  #############################################################################

  # Create new column populated with NAs
  sub_df$Bite.loc <- NA

  # Add information on bite location to new column
  for(m in 1:nrow(sub_df)){
    if(sub_df[[head_col]][m]==1){
      sub_df$Bite.loc[m] <- "Head"
    } else if(sub_df[[arm_col]][m]==1){
      sub_df$Bite.loc[m] <- "Arm"
    } else if(sub_df[[trunk_col]][m]==1){
      sub_df$Bite.loc[m] <- "Trunk"
    } else if(sub_df[[leg_col]][m]==1){
      sub_df$Bite.loc[m] <- "Leg"
    } else {
      sub_df$Bite.loc[m] <- NA
    }
  }

  print(paste0(length(which(is.na(sub_df$Bite.loc)&sub_df[[rabid_label]]==rabid_true)), " individuals did not provide a bite location when bitten by a rabid animal"))

  #############################################################################
  #                               Calculate estimates                         #
  #############################################################################

  if(by_dis==TRUE){
    by_bite.loc <- sub_df %>%
      group_by_(dis_label, "Bite.loc") %>%
      summarise_(n_bite.victims= ~n())

    # Remove row of NA bite loc
    by_bite.loc <- by_bite.loc[!is.na(by_bite.loc$Bite.loc),]

    # Add probability for each bite loc - ADD CONFIDENCE INTERVALS
    by_bite.loc$p_bite.victims <- by_bite.loc$n_bite.victims/sum(by_bite.loc$n_bite.victims, na.rm=T)

    # Add in information regarding # bites from rabid animals
    indx_r <- which(sub_df[[rabid_label]]==rabid_true)
    sub_dataframe <- sub_df[indx_r,]
    by_bite.loc_rabid <- sub_dataframe %>%
      group_by_(dis_label, "Bite.loc") %>%
      summarise_(n_bite.from.rabid.animal= ~n())
    prop_infect <- left_join(by_bite.loc, by_bite.loc_rabid)
    prop_infect$p_bite.from.rabid.animal <- prop_infect$n_bite.from.rabid.animal/sum(prop_infect$n_bite.from.rabid.animal, na.rm=T)

    # Add in information regarding # people that didn't recieved PEP
    sub_dataframe$PEP.0 <- NA
    for(n in 1:nrow(sub_dataframe)){
      if(is.na(sub_dataframe[[PEP.1_label]][n])){
        sub_dataframe$PEP.0[n] <- NA
      } else if(sub_dataframe[[PEP.1_label]][n]==0){
        sub_dataframe$PEP.0[n] <- 1
      } else if(sub_dataframe[[PEP.1_label]][n]==1){
        sub_dataframe$PEP.0[n] <- 0
      }
    }
    indx_nopep <- which(sub_dataframe[["PEP.0"]]==1)
    subbed_df <- sub_dataframe[indx_nopep,]
    by_bite.loc_nopep <- subbed_df %>%
      group_by_(dis_label, "Bite.loc") %>%
      summarise_(n_no.pep= ~n())
    prop_infect <- left_join(prop_infect, by_bite.loc_nopep)
    # prop_infect$p_no.pep <- prop_infect$n_no.pep/prop_infect$n_bite.from.rabid.animal

    # Add in information regarding # people died after not recieving PEP
    indx_d <- which(subbed_df[[COD_label]]=="Rabies")
    subbed_dataframe <- subbed_df[indx_d,]
    by_bite.loc_death <- subbed_dataframe %>%
      group_by_(dis_label, "Bite.loc") %>%
      summarise_(n_deaths= ~n())
    final_df <- left_join(prop_infect, by_bite.loc_death)
    final_df$p_deaths <- final_df$n_deaths/final_df$n_no.pep

    # Add row for totals
    new.row <- data.frame(District="TOTAL", Bite.loc=NA,
                          n_bite.victims=sum(final_df$n_bite.victims, na.rm=T),
                          p_bite.victims=NA,
                          n_bite.from.rabid.animal=sum(final_df$n_bite.from.rabid.animal, na.rm=T),
                          p_bite.from.rabid.animal=NA,
                          n_no.pep=sum(final_df$n_no.pep, na.rm=T),
                          n_deaths=sum(final_df$n_deaths, na.rm=T), p_deaths=NA)
    final_df <- rbind(as.data.frame(final_df), new.row)

  } else if(by_dis==FALSE){

    by_bite.loc <- sub_df %>%
      group_by_("Bite.loc") %>%
      summarise_(n_bite.victims= ~n())

    # Remove row of NA bite loc
    by_bite.loc <- by_bite.loc[!is.na(by_bite.loc$Bite.loc),]

    # Add probability for each bite loc - ADD CONFIDENCE INTERVALS
    by_bite.loc$p_bite.victims <- by_bite.loc$n_bite.victims/sum(by_bite.loc$n_bite.victims, na.rm=T)

    # Add in information regarding # bites from rabid animals
    indx_r <- which(sub_df[[rabid_label]]==rabid_true)
    sub_dataframe <- sub_df[indx_r,]
    by_bite.loc_rabid <- sub_dataframe %>%
      group_by_("Bite.loc") %>%
      summarise_(n_bite.from.rabid.animal= ~n())
    prop_infect <- left_join(by_bite.loc, by_bite.loc_rabid)
    prop_infect$p_bite.from.rabid.animal <- prop_infect$n_bite.from.rabid.animal/sum(prop_infect$n_bite.from.rabid.animal, na.rm=T)

    # Calculate 95% CI
    prop_infect$P_rabid.bite_LCI <- NA
    prop_infect$P_rabid.bite_UCI <- NA

    for(i in 1:nrow(prop_infect)){
      p_rabid.bite_CIs <- binom.test(prop_infect$n_bite.from.rabid.animal[i], sum(prop_infect$n_bite.from.rabid.animal))
      prop_infect$P_rabid.bite_LCI[i] <- p_rabid.bite_CIs$conf.int[1]
      prop_infect$P_rabid.bite_UCI[i] <- p_rabid.bite_CIs$conf.int[2]
    }

    # Add in information regarding # people that didn't recieved PEP
    sub_dataframe$PEP.0 <- NA
    for(n in 1:nrow(sub_dataframe)){
      if(is.na(sub_dataframe[[PEP.1_label]][n])){
        sub_dataframe$PEP.0[n] <- NA
      } else if(sub_dataframe[[PEP.1_label]][n]==0){
        sub_dataframe$PEP.0[n] <- 1
      } else if(sub_dataframe[[PEP.1_label]][n]==1){
        sub_dataframe$PEP.0[n] <- 0
      }
    }
    indx_nopep <- which(sub_dataframe[["PEP.0"]]==1)
    subbed_df <- sub_dataframe[indx_nopep,]
    by_bite.loc_nopep <- subbed_df %>%
      group_by_("Bite.loc") %>%
      summarise_(n_no.pep= ~n())
    prop_infect <- left_join(prop_infect, by_bite.loc_nopep)
    # prop_infect$p_no.pep <- prop_infect$n_no.pep/prop_infect$n_bite.from.rabid.animal

    # Add in information regarding # people died after not recieving PEP
    indx_d <- which(subbed_df[[COD_label]]=="Rabies")
    subbed_dataframe <- subbed_df[indx_d,]
    by_bite.loc_death <- subbed_dataframe %>%
      group_by_("Bite.loc") %>%
      summarise_(n_deaths= ~n())
    final_df <- left_join(prop_infect, by_bite.loc_death)
    final_df$p_deaths <- final_df$n_deaths/final_df$n_no.pep

    # Add in CI columns
    final_df$P_death_LCI <- NA
    final_df$P_death_UCI <- NA

    # Calculate 95% CI
    for(i in 1:nrow(final_df)){
      p_death_CIs <- binom.test(final_df$n_deaths[i], final_df$n_no.pep[i])
      final_df$P_death_LCI[i] <- p_death_CIs$conf.int[1]
      final_df$P_death_UCI[i] <- p_death_CIs$conf.int[2]
    }

    # Add row for totals
    # new.row <- data.frame(Bite.loc="Total",
    #                       n_bite.victims=sum(final_df$n_bite.victims, na.rm=T),
    #                       p_bite.victims=NA,
    #                       n_bite.from.rabid.animal=sum(final_df$n_bite.from.rabid.animal, na.rm=T),
    #                       p_bite.from.rabid.animal=NA,
    #                       n_no.pep=sum(final_df$n_no.pep, na.rm=T),
    #                       n_deaths=sum(final_df$n_deaths, na.rm=T), p_deaths=NA,
    #                       P_death_LCI=NA, P_death_UCI=NA)
    # final_df <- rbind(final_df, new.row)
    #
    # new.df <- sub_df[sub_df$Rabid=="Yes",]
    # new.df <- new.df[!is.na(new.df$Bite.loc),]
    # print((sum(new.df$PEP.1==1, na.rm=T))/(nrow(new.df)))
  }

  # Create final df to go to decision tree
  final.dataframe <- data.frame(prop_rabid.bites=(sum(final_df$n_bite.from.rabid.animal))/(sum(final_df$n_bite.victims)),
                                n_rabid.bites=sum(final_df$n_bite.from.rabid.animal, na.rm=TRUE),
                                n_bite.vics=sum(final_df$n_bite.victims, na.rm=TRUE),
                                prob_head=as.character(round(final_df[2,5], digits=8)),
                                n_head.bites=final_df$n_bite.from.rabid.animal[2],
                                prob_arm=as.character(round(final_df[1,5], digits=8)),
                                n_arm.bites=final_df$n_bite.from.rabid.animal[1],
                                prob_trunk=as.character(round(final_df[4,5], digits=8)),
                                n_trunk.bites=final_df$n_bite.from.rabid.animal[4],
                                prob_leg=as.character(round(final_df[3,5], digits=8)),
                                n_leg.bites=final_df$n_bite.from.rabid.animal[3],
                                prob_death_given_head=as.character(round(final_df[2,10], digits=8)),
                                n_no.pep.head=final_df$n_no.pep[2],
                                n_deaths.head=final_df$n_deaths[2],
                                prob_death_given_arm=as.character(round(final_df[1,10], digits=8)),
                                n_no.pep.arm=final_df$n_no.pep[1],
                                n_deaths.arm=final_df$n_deaths[1],
                                prob_death_given_trunk=as.character(round(final_df[4,10], digits=8)),
                                n_no.pep.trunk=final_df$n_no.pep[4],
                                n_deaths.trunk=final_df$n_deaths[4],
                                prob_death_given_leg=as.character(round(final_df[3,10], digits=8)),
                                n_no.pep.leg=final_df$n_no.pep[3],
                                n_deaths.leg=final_df$n_deaths[3]
                                )

  return(list(props=final_df, dec_tree=final.dataframe))
}