#' @title Calculate proportion of bite victims that recieved multiple bites
#' @description This function produces a dataframe, giving the proportion of
#' individuals that recieved mutliple bites during the bite incident, the numerator
#' and denominator for the proportion, as well as 95% confidence intervals.
#' @param dataframe the dataframe from which the proportion is calculated
#' @param head_col (default="Head") the column stating a bite wound to the head
#' @param arm_col (default="Arm") the column stating a bite wound to the arm
#' @param hand_col (default="Hand") the column stating a bite wound to the hand
#' @param trunk_col (default="Trunk") the column stating a bite wound to the trunk
#' @param leg_col (default="Leg") the column stating a bite wound to the leg
#' @param foot_col (default="Foot") the column stating a bite wound to the foot
#' @return a dataframe providing the results of the calculations
#' @export
TZ_p_multiple.bites <- function(dataframe, head_col="Head", arm_col="Arm", hand_col="Hand",
                                trunk_col="Trunk", leg_col="Leg", foot_col="Foot"){

  # Load in required libraries
  require(dplyr)

  # Take only selected columns from dataframe
  sub_dataframe <- dataframe[,c(head_col, arm_col, hand_col, trunk_col, leg_col, foot_col)]

  ################################################################################
  #                CALCULATE PROPORTION OF MULTIPLE BITE VICTIMS                 #
  ################################################################################

  sub_dataframe$n.bites <- NA

  # Calculate which rows have more than one bite location
  for(j in 1:nrow(sub_dataframe)){
    sub_dataframe$n.bites[j] <- sum(sub_dataframe[[head_col]][j], sub_dataframe[[arm_col]][j],
                                    sub_dataframe[[hand_col]][j], sub_dataframe[[trunk_col]][j],
                                    sub_dataframe[[leg_col]][j], sub_dataframe[[foot_col]][j],na.rm=T)
  }

  # Seperate off rows with mutliple bites
  multiple.bite_rows <- subset(sub_dataframe, n.bites>1)

  # Calculate 95% CI
  CIs <- binom.test(nrow(multiple.bite_rows), nrow(dataframe))

  return(data.frame(p_multiple.bites=nrow(multiple.bite_rows)/nrow(dataframe),
                    n_multiple.bites=nrow(multiple.bite_rows),
                    n_bite.victims=nrow(dataframe),
                    LCI=CIs$conf.int[1], UCI=CIs$conf.int[2],
                    n_2.bites=length(which(multiple.bite_rows$n.bites==2)),
                    n_3.bites=length(which(multiple.bite_rows$n.bites==3)),
                    n_4.bites=length(which(multiple.bite_rows$n.bites==4))))
}