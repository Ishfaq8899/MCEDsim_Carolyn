
###########################################################################################
#sim_individual_MCED
#Author: 
#
# INPUTS: 
#
#OUTPUTS: A data frame with the simulated data for the individual.
#################################################################################################
sim_individual_MCED<-function( ID, 
                               cancer_sites,
                               rates_list,
                               test_performance,
                               other_cause_death_dist,
                               starting_age,
                               num_screens,
                               screen_interval
                               #cause_specific_mortality
                               ){
  
  # simulate time of other cause death 

  other_cause_death_time = sim_othercause_death(other_cause_death_dist)
  
  
  # Get the starting states for each cancer
  start_states= sapply(rates_list, FUN="get_init",a1=starting_age)

  ### get the screening times
  screen_times = seq((starting_age),(num_screens + starting_age-1), by=screen_interval)
  
 
  result <- sim_multiple_cancer_indiv(ID = ID,
                                        cancer_sites = cancer_sites,
                                        rate_matrices = rates_list,
                                        early_sensitivities = test_performance$early_sens,
                                        late_sensitivities =  test_performance$late_sens,
                                        specificities = test_performance$specificities,
                                        obs.times = screen_times,
                                        start.time = starting_age,
                                        end.time = 500, 
                                        start.states =  start_states) 
  
  # result$death_time_ovarian_cancer=mapply(FUN="generate_ovarian_death",result$screen_diagnosis_time,
  #                                         result$screen_diagnosis_stage,result$clinical_diagnosis_time,
  #                                         result$clinical_diagnosis_stage,
  #                                         result$cancer_site,MoreArgs=list(ovarian_survival_dist))
  # 
  # result$arm = arm

#-----------------------  
# Ishfaq modification   
#-----------------------
# Identify first cancer by onset time
if ("onset_time" %in% names(result)) {
  first_cancer_row <- result %>%
    filter(!is.na(onset_time)) %>%
# selects the single row with the earliest (smallest) onset_time.    
  slice_min(order_by = onset_time, n = 1, with_ties = FALSE)
  
  result <- result %>%
    mutate(is_first_cancer = (cancer_site == first_cancer_row$cancer_site[1] &
                                onset_time == first_cancer_row$onset_time[1]))
} else {
  result$is_first_cancer <- FALSE
}
#-----------------------  

# Add other-cause death info
result$other_cause_death_status <- other_cause_death_time$status
result$other_cause_death_time <- other_cause_death_time$time

return(result)
}

###########################################################################################
#sim_multi_individuals_MCED
#Author: 
# Function to simulate multiple individuals 
# INPUTS: 
#
#OUTPUTS: A data frame with the combined simulated data for multiple individuals.
############################################################################################
sim_multiple_individuals_MCED_parallel_universe <- function(  cancer_sites,
                                         LMST_vec, 
                                         OMST_vec, 
                                         test_performance_dataframe, 
                                         starting_age,
                                         num_screens,
                                         screen_interval,
                                         num_males, 
                                         num_females, 
                                         all_rates_male, 
                                         all_rates_female, 
                                         all_meta_data_female, 
                                         all_meta_data_male, 
                                         cdc_data, 
                                         hmd_data,
                                         MCED_cdc
                                         #survival tables
                                         ){
  
  
  # Creat a vector of IDs
  IDs_male <- 1:num_males
  IDs_female <-(num_males+1):(num_females+num_males)
  
  
  # Extract rate matrices matrices based on OMST and LMST specs (Male)
  rates_list_male = get_filtered_rates(the_omsts = OMST_vec, the_lmsts = LMST_vec,
                                   all_meta_data = all_meta_data_male,
                                   all_rates = all_rates_male, the_cancer_site = cancer_sites)

  sites_male = rates_list_male$cancer_sites 
  rates_list_male = rates_list_male$rates_list

  # Extract rate matrices (Female)
  rates_list_female = get_filtered_rates(the_omsts = OMST_vec, the_lmsts = LMST_vec, 
                                     all_meta_data = all_meta_data_female, 
                                     all_rates = all_rates_female, the_cancer_site = cancer_sites)
  sites_female = rates_list_female$cancer_sites 
  rates_list_female = rates_list_female$rates_list                                    
  
  # Extract sensitivities and specificity based on selected cancer sites
  test_performance_male = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_male))
  test_performance_female = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_female))
  
 # browser()
  #Get the other-cause death tables for men and women
  other_cause_death_male=make_othercause_death_table(cdc_data=cdc_data,
                                                     MCED_cdc=MCED_cdc, 
                                                     hmd_data=hmd_data,
                                                     the_starting_age = starting_age,
                                                     the_sex="Male",
                                                     selected_cancers=sites_male,
                                                     the_year=2022)
  
  other_cause_death_female=make_othercause_death_table(cdc_data=cdc_data,
                                                       MCED_cdc=MCED_cdc, 
                                                       hmd_data=hmd_data,
                                                       the_starting_age = starting_age,
                                                       the_sex="Female",
                                                   selected_cancers=sites_female,
                                                   the_year=2022)

  # Use mapply to apply the sim_individual_UKCTOCS function to each ID
  results_list_male <- mapply(sim_individual_MCED,
                         ID = IDs_male,
                         MoreArgs = list(rates_list=rates_list_male,
                                         cancer_sites=sites_male,
                                         test_performance=test_performance_male,
                                         other_cause_death_dist=other_cause_death_male,
                                         starting_age=starting_age,
                                         num_screens=num_screens,
                                         screen_interval=screen_interval),
                         SIMPLIFY = FALSE)

  
  results_list_female <- mapply(sim_individual_MCED,
                                ID = IDs_female,
                                MoreArgs = list(rates_list=rates_list_female,
                                                cancer_sites=sites_female,
                                                test_performance=test_performance_female,
                                                other_cause_death_dist=other_cause_death_female,
                                                starting_age=starting_age,
                                                num_screens=num_screens,
                                                screen_interval=screen_interval),
                                SIMPLIFY = FALSE)
  
  # Combine all individual results into one data frame
  combined_results_males <- do.call(rbind, results_list_male)
  combined_results_females <- do.call(rbind, results_list_female)

  combined_results =bind_rows(combined_results_males,combined_results_females)%>%
                   mutate(start_age=starting_age)
  
  return(combined_results)
}

###########################################################################################
#generate_ovarian_death
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
generate_ovarian_death<-function(screen_diagnosis_time,screen_diagnosis_stage,clinical_diagnosis_time,
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





