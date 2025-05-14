
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
                               screen_interval,
                               end_time,
                               sex,
                               surv_param_table
                               ){
  
  # simulate time of other cause death 

  other_cause_death = sim_othercause_death(other_cause_death_dist)
  
  
  # Get the starting states for each cancer
  start_states= sapply(rates_list, FUN="get_init",a1=starting_age)

  ### get the screening times
  screen_times = seq((starting_age),(num_screens + starting_age-1), by=screen_interval)
  
  #number of cancer sites
  num_sites=length(cancer_sites)
 
  result <- sim_multiple_cancer_indiv(ID = ID,
                                        cancer_sites = cancer_sites,
                                        rate_matrices = rates_list,
                                        early_sensitivities = test_performance$early_sens,
                                        late_sensitivities =  test_performance$late_sens,
                                        specificities = test_performance$specificities/num_sites,
                                        obs.times = screen_times,
                                        start.time = starting_age,
                                        end.time = end_time, 
                                        start.states =  start_states) 
  
  # Add other-cause death info
  result$other_cause_death_status <- other_cause_death$status
  result$other_cause_death_time <- other_cause_death$time
  
  #-----------------------  
  # Identify first cancer by onset time
  #-----------------------
  result$is_first_cancer=FALSE
  result$cancer_death_time_no_screen=NA

  if (!is.na(min(result$onset_time,na.rm=T))) {
    
    first_cancer_row <- result %>%
      filter(!is.na(onset_time)) %>%
      # selects the single row with the earliest (smallest) onset_time.    
      slice_min(order_by = onset_time, n = 1, with_ties = FALSE)%>%mutate(is_first_cancer=TRUE)
  
    #Simulate time of death for first cancer without screening and with screening
    #if stage at clinical diagnosis is same as stage at screen diagnosis, then cancer_death_time_screen=cancer_death_time_no_screen
    #if it is different, then cancer_death_time_screen=clinical_diagnosis_time+sim_cancer_death_param(the_stage=screen_diagnosis_stage,
                                                       # the_cancer_site=cancer_site,
                                                       #the_sex=sex,
                                                      #the_model_type="Loglogistic",
                                                       #param_table=surv_param_table)
    
    first_cancer_row <-first_cancer_row  %>%mutate(cancer_death_time_no_screen=ifelse(is_first_cancer,clinical_diagnosis_time+sim_cancer_death_param(the_stage=clinical_diagnosis_stage,
                                                                                                       the_cancer_site=cancer_site,
                                                                                                       the_sex=sex,
                                                                                                       the_model_type="Loglogistic",
                                                                                                       param_table=surv_param_table),NA))%>%
                    mutate(cancer_death_time_screen=ifelse(screen_diagnosis_stage==clinical_diagnosis_stage|is.na(screen_diagnosis_time),cancer_death_time_no_screen,
                                                           clinical_diagnosis_time+sim_cancer_death_param(the_stage=screen_diagnosis_stage,
                                                                                                          the_cancer_site=cancer_site,
                                                                                                          the_sex=sex,
                                                                                                          the_model_type="Loglogistic",
                                                                                                          param_table=surv_param_table)))
    
   
      result <- result %>%
      mutate(is_first_cancer = (cancer_site == first_cancer_row$cancer_site[1] &
                                  onset_time == first_cancer_row$onset_time[1]))%>%
      mutate(cancer_death_time_no_screen=ifelse(is_first_cancer,first_cancer_row$cancer_death_time_no_screen,NA),
             cancer_death_time_screen=ifelse(is_first_cancer,first_cancer_row$cancer_death_time_screen,NA))
    
  
  } 
  

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
                                         ending_age,
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
                                         MCED_cdc,
                                         surv_param_table
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
                                         screen_interval=screen_interval,
                                         end_time=ending_age,
                                         surv_param_table=surv_param_table,
                                         sex="Male"),
                         SIMPLIFY = FALSE)

  results_list_female <- mapply(sim_individual_MCED,
                                ID = IDs_female,
                                MoreArgs = list(rates_list=rates_list_female,
                                                cancer_sites=sites_female,
                                                test_performance=test_performance_female,
                                                other_cause_death_dist=other_cause_death_female,
                                                starting_age=starting_age,
                                                num_screens=num_screens,
                                                screen_interval=screen_interval,
                                                end_time=ending_age,
                                                surv_param_table=surv_param_table,
                                                sex="Female"),
                                SIMPLIFY = FALSE)
  
  # Combine all individual results into one data frame
  combined_results_males <- do.call(rbind, results_list_male)%>%mutate(sex="Male")
  combined_results_females <- do.call(rbind, results_list_female)%>%mutate(sex="Female")

  combined_results =bind_rows(combined_results_males,combined_results_females)%>%
                   mutate(start_age=starting_age)
  
  return(combined_results)
}






