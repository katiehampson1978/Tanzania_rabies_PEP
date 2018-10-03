# Load in data
transmission <- read.csv("output/transmission.csv", stringsAsFactors=FALSE)

# This function uses information on number of bite incidents, expoosures, no-PEP,
# deaths and location on the body of bites.
mixture_model<-function(n){ # Cleaveland 2002 mixture model
  
  N = transmission$n_rabid.bites # N BITES
  
  # BITES by site
  n_head = transmission$n_head.bites
  n_arm = transmission$n_arm.bites
  n_trunk = transmission$n_trunk.bites
  n_leg = transmission$n_leg.bites
  
  # DEATHS by site
  death_head = transmission$n_deaths.head
  death_arm = transmission$n_deaths.arm
  death_trunk = transmission$n_deaths.trunk
  death_leg = transmission$n_deaths.leg
  
  # No PEP by site
  n_no.pep.head = transmission$n_no.pep.head
  n_no.pep.arm = transmission$n_no.pep.arm
  n_no.pep.trunk = transmission$n_no.pep.trunk
  n_no.pep.leg = transmission$n_no.pep.leg
  
  # random draw of probability of bite site
  prob_head <- rbinom(n=n, size=N, prob=n_head/N)/N
  prob_arm <- rbinom(n=n, size=N, prob=n_arm/N)/N
  prob_trunk <- rbinom(n=n, size=N, prob=n_trunk/N)/N
  prob_leg <- rbinom(n=n, size=N, prob=n_leg/N)/N
  
  # random draw of probability of death given bite site
  prob_death_given_head <- rbinom(n=n, size=n_no.pep.head, prob=death_head/n_no.pep.head)/n_no.pep.head
  prob_death_given_arm <- rbinom(n=n, size=n_no.pep.arm, prob=death_arm/n_no.pep.arm)/n_no.pep.arm
  prob_death_given_trunk <- rbinom(n=n, size=n_no.pep.trunk, prob=death_trunk/n_no.pep.trunk)/n_no.pep.trunk
  prob_death_given_leg <-  rbinom(n=n, size=n_no.pep.leg, prob=death_leg/n_no.pep.leg)/n_no.pep.leg
  
  prob_death_given_rabid_bite <- prob_head*prob_death_given_head +
    prob_arm*prob_death_given_arm+
    prob_trunk*prob_death_given_trunk+
    prob_leg*prob_death_given_leg
  
  #return(prob_death_given_rabid_bite)
  return(data.frame(p_death=prob_death_given_rabid_bite,
                    p_head=prob_head,
                    p_death_head=prob_death_given_head*prob_head,
                    p_arm=prob_arm,
                    p_death_arm=prob_death_given_arm*prob_arm,
                    p_trunk=prob_trunk,
                    p_death_trunk=prob_death_given_trunk*prob_trunk,
                    p_leg=prob_leg,
                    p_death_leg=prob_death_given_leg*prob_leg))
}

# This function calculates quantiles for confidence of bite location and deaths
get_quantiles <- function(model_results){
  
  # Create a vector of probability of bite site for each bite location (+ CIs)
  p_head <- round(quantile(model_results$p_head, c(0.025, 0.5, 0.975)), digits=3)
  p_arm <- round(quantile(model_results$p_arm, c(0.025, 0.5, 0.975)), digits=3)
  p_trunk <- round(quantile(model_results$p_trunk, c(0.025, 0.5, 0.975)), digits=3)
  p_leg <- round(quantile(model_results$p_leg, c(0.025, 0.5, 0.975)), digits=3)
  p_bite_site <- c(paste0(p_arm[2], " (", p_arm[1], "-", p_arm[3], ")"),
                    paste0(p_head[2], " (", p_head[1], "-", p_head[3], ")"),
                    paste0(p_leg[2], " (", p_leg[1], "-", p_leg[3], ")"),
                    paste0(p_trunk[2], " (", p_trunk[1], "-", p_trunk[3], ")"))
  
  # Create a vector of probability of infection for each bite location (+ CIs)
  p_infect_head <- round(quantile(model_results$p_death_head, c(0.025, 0.5, 0.975)), digits=3)
  p_infect_arm <- round(quantile(model_results$p_death_arm, c(0.025, 0.5, 0.975)), digits=3)
  p_infect_trunk <- round(quantile(model_results$p_death_trunk, c(0.025, 0.5, 0.975)), digits=3)
  p_infect_leg <- round(quantile(model_results$p_death_leg, c(0.025, 0.5, 0.975)), digits=3)
  p_infect <- c(paste0(p_infect_arm[2], " (", p_infect_arm[1], "-", p_infect_arm[3], ")"),
                    paste0(p_infect_head[2], " (", p_infect_head[1], "-", p_infect_head[3], ")"),
                    paste0(p_infect_leg[2], " (", p_infect_leg[1], "-", p_infect_leg[3], ")"),
                    paste0(p_infect_trunk[2], " (", p_infect_trunk[1], "-", p_infect_trunk[3], ")"))
  
  # Create string of overall probability of infection (+ CIs)
  p_infect_overall <- round(quantile(model_results$p_death, c(0.025, 0.5, 0.975)), digits=3)
  p_infect_overall <- paste0(p_infect_overall[2], " (", p_infect_overall[1], "-", p_infect_overall[3], ")")
  
  return(list("p_bite_site"=p_bite_site, "p_infect"=p_infect, "p_infect_overall"=p_infect_overall))
}