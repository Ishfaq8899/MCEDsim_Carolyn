library(msm)
library(dplyr)
###########################################################################################
#get_init
#Author: JL
#Get the initial distribution of all states, assuming no one is clinically diagnosed at the start of the trial. 
#Sample from this distribution to get an intial state

#INPUTS:
#        rate.matrix=rate matrix for natural history model
#        a1= start age at the beginning of trial
#OUTPUTS: A vector of size k that provides the initial distribution
#################################################################################################
get_init<-function(rate.matrix,a1){
  
  k=dim(rate.matrix)[1]
  init=vector()
  init[k]=0
  init[k-1]=0
  
  prob_mat_a1=MatrixExp(t=a1,mat=rate.matrix)
  denom=1-prob_mat_a1[1,(k-1)]-prob_mat_a1[1,(k)]
  for(j in 1:(k-2)){
    
    init[j]=prob_mat_a1[1,j]/denom
  }
  
  init_state=sample(1:k,size=1,prob=init)
  return(init_state)
}

#######################################################################
#gettime
#Author: JL
#  The function generates a uniform random number and uses it to interpolate
#  the time at which the event occurs based on the cumulative distribution function
#INPUTS: (time, surv)
#
#OUTPUTS: A data from with two columns:
#        -time: The simulated event time based on the given survival probabilities.
#        -status: The event status (1 indicating the event occured).
########################################################################
gettime<-function(time,surv){
  
  unif1=runif(n=1)
  cdf=1-surv
  
  eventtime=approx(x=cdf, y=time, xout=unif1, rule=2)$y
  
  eventstatus=1
 
  return(data.frame(time = eventtime, status = eventstatus))
}
###########################################################################################
# sim_cancer_death 
# Author: 
# Function to simulate the time of cancer-specific death
# INPUTS: 
#    age_clin_dx: Age at clinical diagnosis
#    stage_dx: Stage at diagnosis
#    cancer_type: cancer type
#    cancer_survival_dist: Data frame containing the cancer survival distributions.
#
# OUTPUTS: A numeric value representing the simulated time of death due to a specific cancer
#################################################################################################
sim_cancer_death <- function(age_clin_dx, stage_dx, cancer_type, cancer_survival_dist, model){
  cancer_type=the_cancer_type
  # Filter the survival distribution based on the type and stage
  # browser()
  
  survival_dist_indiv = filter(ovarian_survival_dist, HGSC==HGSCtype,type == "weibull",
                               stage == stage_dx,stage4cat==FALSE)
  
  # Get the survival time based on the distribution
  death_info = gettime(time = survival_dist_indiv$time, surv = survival_dist_indiv$surv)
  death_time = death_info$time
  
  return(death_time)
}

###########################################################################################
#generate_stageshift_death
#Author: 
# Function to generate the simulated time of ovarian cancer death based on the diagnosis stage and type
#
# INPUTS: 
#     screen_diagnosis_time: Time of diagnosis via screening 
#     screen_diagnosis_stage: Stage at diagnosis via screening
#     clinical_diagnosis_time: Time of diagnosis via clinical symptoms
#     clinical_diagnosis_stage: Stage at diagnosis via clinical symptoms
#     cancer_site: Cancer site (HGSC or nonHGSC)
#     ovarian_survival_dist: Data frame containing ovarian survival distribution
#
#OUTPUTS: A numeric value representing the simulated time of death due to ovarian cancer
#################################################################################################
generate_stageshift_death<-function(screen_diagnosis_time,screen_diagnosis_stage,clinical_diagnosis_time,
                                clinical_diagnosis_stage,cancer_site,ovarian_survival_dist){
  #  browser()
  
  if(!is.na(clinical_diagnosis_time)){
    if(!is.na(screen_diagnosis_stage)&screen_diagnosis_stage=="early"){
      the_stage="early"
    }else{
      the_stage=clinical_diagnosis_stage
    }
    ovarian_death=sim_ovarian_death(age_clin_dx=clinical_diagnosis_time, stage_dx=the_stage, 
                                    type=cancer_site, ovarian_survival_dist=ovarian_survival_dist)
  }else{
    ovarian_death=NA
  }
  
  return(ovarian_death=ovarian_death)
}

###########################################################################################

###########################################################################################
#sim_individual_trial
#Author: 
# Function to simulate an individual age at diagnosis of HGSC and non HGSC based on the number
#  screen in UKCTOCS trials
#
# INPUTS: 
#     ID: individual
#     age_entry: Age of the individual at entry into the trial
#     arm: (M, U, or control)
#     number_of_screens: Number of screening events the individual undergoes
#     censoring_dist: Data frame containing the censoring distribution
#     rate_HGSC/nonHGSC: Rate matrix for HGSC/nonHGSC
#     sens_early/late_HGSC_M/nonHGSC_M: Sensitivity for early/late HGSC/nonHGSC in arm M
#     sens_early/late_HGSC_U/nonHGSC_U: Sensitivity for early/late HGSC/nonHGSC in arm U
#     age_origin: origin age for the trial
#     ovarian_survival_dist: Data frame containing ovarian survival distribution
#
#OUTPUTS: A data frame with the simulated data for the individual.
#################################################################################################
sim_individual_trial<-function( ID, age_entry, arm, number_of_screens, 
                                  censoring_dist, rate_HGSC, rate_nonHGSC,
                                  sens_early_HGSC_M, sens_late_HGSC_M, 
                                  sens_early_nonHGSC_M, sens_late_nonHGSC_M,
                                  sens_early_HGSC_U, sens_late_HGSC_U, 
                                  sens_early_nonHGSC_U, sens_late_nonHGSC_U,
                                  age_origin, ovarian_survival_dist){
 
  # simulate censoring time  
  age_cat = cut(age_entry,breaks=seq(50,75,by=5),right=T)
  if(age_entry==50){
   #  browser()
    age_cat="(50,55]"
  }
  
  censoringtime = sim_censoring_time(group = arm, age_at_entry = age_cat,
                                     number_of_screens = number_of_screens,
                                     censoring_dist = censoring_dist)
  
 
   # Get the starting state
  the_start_HGSC = get_init(rate_HGSC, age_entry-age_origin)
  the_start_nonHGSC = get_init(rate_nonHGSC, age_entry-age_origin)
  

  start_states = c(the_start_HGSC, the_start_nonHGSC) 

  
  if(arm=="M"|arm=="U"){ 
    # Set sensitivities based on the arm
    if (arm == "M"){
      sens_early_HGSC <-  sens_early_HGSC_M
      sens_late_HGSC <- sens_late_HGSC_M
      sens_early_nonHGSC <- sens_early_nonHGSC_M
      sens_late_nonHGSC <- sens_late_nonHGSC_M
    } else if (arm == "U"){
      sens_early_HGSC <-  sens_early_HGSC_U
      sens_late_HGSC <- sens_late_HGSC_U
      sens_early_nonHGSC <- sens_early_nonHGSC_U
      sens_late_nonHGSC <- sens_late_nonHGSC_U
    }
    
    ### get the screening times
   # browser()
    screen_times = seq((age_entry-age_origin),(number_of_screens + age_entry-age_origin-1), by=1)
    
    result <- sim_multiple_cancer_indiv(ID = ID,
                                        cancer_sites = c("HGSC","nonHGSC"),
                                        rate_matrices = list(rate_HGSC, rate_nonHGSC),
                                        early_sensitivities = c(sens_early_HGSC, sens_early_nonHGSC), 
                                        late_sensitivities =  c(sens_late_HGSC, sens_late_nonHGSC),
                                        specificities = c(1, 1),
                                        obs.times = screen_times,
                                        start.time = age_entry-age_origin,
                                        end.time = 100, 
                                        start.states =  start_states) 
    
  }else{
    #HGSC
    out1 = get.obs.data.individual.control(ID, rate.matrix=rate_HGSC, end.time = 10000, 
                                           start.time = age_entry - age_origin, start.state = the_start_HGSC)
    out1$cancer_site = "HGSC"
    #nonHGSC
    out2 = get.obs.data.individual.control(ID, rate.matrix = rate_nonHGSC, end.time = 10000, 
                                           start.time = age_entry - age_origin, start.state = the_start_nonHGSC)
    out2$cancer_site = "nonHGSC"
    result = rbind(out1,out2)
  }
  
  
    result$death_time_ovarian_cancer=mapply(FUN="generate_ovarian_death",result$screen_diagnosis_time,
                                          result$screen_diagnosis_stage,result$clinical_diagnosis_time,
                                          result$clinical_diagnosis_stage,
                                          result$cancer_site,MoreArgs=list(ovarian_survival_dist))
  
  result$arm = arm
  result$censoring_time=censoringtime
  return(result)
  
}  

##############################################################################################
#sim_multi_individuals_UKCTOCS
#Author: 
# Function to simulate multiple individuals age at diagnosis of HGSC and non HGSC based on the number
#  screen in UKCTOCS trials
#
# INPUTS: 
#     num_individuals, age_entry: Age of the individual at entry into the trial
#     arm: Trail arm (M, U, or control)
#     number_of_screens: Number of screening events the individual undergoes
#     censoring_dist: Data frame containing the censoring distribution
#     rate_HGSC/nonHGSC: Rate matrix for HGSC/nonHGSC
#     sens_early/late_HGSC_M/nonHGSC_M: Sensitivity for early/late HGSC/nonHGSC in arm M
#     sens_early/late_HGSC_U/nonHGSC_U: Sensitivity for early/late HGSC/nonHGSC in arm U
#     age_origin: origin age for the trial
#     ovarian_survival_dist: Data frame containing ovarian survival distribution
#
#OUTPUTS: A data frame with the combined simulated data for multiple individuals.
############################################################################################
sim_multi_individuals_UKCTOCS <- function(num_individuals,age_entry, arm, number_of_screens, 
                                          censoring_dist, rate_HGSC, rate_nonHGSC,
                                          sens_early_HGSC_M, sens_late_HGSC_M, 
                                          sens_early_nonHGSC_M, sens_late_nonHGSC_M,
                                          sens_early_HGSC_U, sens_late_HGSC_U, 
                                          sens_early_nonHGSC_U, sens_late_nonHGSC_U,
                                          age_origin, ovarian_survival_dist){
  
  # Creat a vector of IDs
  IDs <- 1:num_individuals
  
  # Use mapply to apply the sim_individual_UKCTOCS function to each ID
  results_list <- mapply(sim_individual_UKCTOCS,
                         ID = IDs,
                         MoreArgs = list(age_entry = age_entry,
                                         arm = arm,
                                         number_of_screens = number_of_screens,
                                         censoring_dist = censoring_dist,
                                         rate_HGSC = rate_HGSC,
                                         rate_nonHGSC = rate_nonHGSC,
                                         sens_early_HGSC_M = sens_early_HGSC_M,
                                         sens_late_HGSC_M = sens_late_HGSC_M, 
                                         sens_early_nonHGSC_M = sens_early_nonHGSC_M, 
                                         sens_late_nonHGSC_M = sens_late_nonHGSC_M,
                                         sens_early_HGSC_U = sens_early_HGSC_U,
                                         sens_late_HGSC_U = sens_late_HGSC_U, 
                                         sens_early_nonHGSC_U = sens_early_nonHGSC_U, 
                                         sens_late_nonHGSC_U = sens_late_nonHGSC_U,
                                         age_origin = age_origin, 
                                         ovarian_survival_dist = ovarian_survival_dist),
                         SIMPLIFY = FALSE)
  
  # Combine all individual results into one data frame
  combined_results <- do.call(rbind, results_list)
  combined_results =combined_results %>% mutate(age_at_entry=age_entry,num_screens=number_of_screens)
  
  return(combined_results)
}

###########################################################################################
#simulate_from_input
#Author: 
# Function to simulate individuals from input data
#
# INPUTS: 
#     input_data: data containing input data (num_individuals, age_entry, arm, number_of_screens)
#     censoring_dist: Data frame containing the censoring distribution
#     rate_HGSC/nonHGSC: Rate matrix for HGSC/nonHGSC
#     sens_early/late_HGSC_M/nonHGSC_M: Sensitivity for early/late HGSC/nonHGSC in arm M
#     sens_early/late_HGSC_U/nonHGSC_U: Sensitivity for early/late HGSC/nonHGSC in arm U
#     age_origin: origin age for the trial
#     ovarian_survival_dist: Data frame containing ovarian survival distribution
#
#OUTPUTS: A data frame with the combined simulated data for all individuals.
############################################################################################
simulate_from_input <- function(input_data, 
                                censoring_dist, rate_HGSC, rate_nonHGSC,
                                sens_early_HGSC_M, sens_late_HGSC_M, 
                                sens_early_nonHGSC_M, sens_late_nonHGSC_M,
                                sens_early_HGSC_U, sens_late_HGSC_U, 
                                sens_early_nonHGSC_U, sens_late_nonHGSC_U,
                                age_origin, ovarian_survival_dist){
  
 
  
  # Use mapply to simulate individuals based on input data 
  results_list <- mapply(function(num_individuals, age_entry, arm, number_of_screens){
    
    sim_multi_individuals_UKCTOCS(num_individuals = num_individuals,
                                  age_entry = age_entry,
                                  arm = arm,
                                  number_of_screens = number_of_screens,
                                  censoring_dist = censoring_dist,
                                  rate_HGSC = rate_HGSC,
                                  rate_nonHGSC = rate_nonHGSC,
                                  sens_early_HGSC_M = sens_early_HGSC_M,
                                  sens_late_HGSC_M = sens_late_HGSC_M, 
                                  sens_early_nonHGSC_M = sens_early_nonHGSC_M, 
                                  sens_late_nonHGSC_M = sens_late_nonHGSC_M,
                                  sens_early_HGSC_U = sens_early_HGSC_U,
                                  sens_late_HGSC_U = sens_late_HGSC_U, 
                                  sens_early_nonHGSC_U = sens_early_nonHGSC_U, 
                                  sens_late_nonHGSC_U = sens_late_nonHGSC_U,
                                  age_origin = age_origin, 
                                  ovarian_survival_dist = ovarian_survival_dist)
  }, num_individuals = input_data$count,
     age_entry = input_data$Age_at_entry,
     arm = input_data$group,
     number_of_screens = input_data$num_screens,
     SIMPLIFY = FALSE)
    
  # Combine all results into one data frame 
  combined_results <- do.call(rbind, results_list)
  
  return(combined_results)
  }
  
  
  
  
  
  
  
  











