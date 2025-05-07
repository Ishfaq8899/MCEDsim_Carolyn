
# This function extracts rate matrices and cancer sites based on OMST and LMST specifications.
get_rates_list<-function(the_omst, the_lmst, all_meta_data_fits,all_fits,the_cancer_site){
  the_indices=all_meta_data_fits %>% filter(OMST%in% the_omst&LMST%in%the_lmst&cancer_site%in%the_cancer_site)%>% select("index")
  rates_list=lapply(all_fits[unlist(the_indices)],"[[","rate.matrix")
  cancer_sites=  the_indices=all_meta_data_fits %>% filter(OMST%in% the_omst&LMST%in%the_lmst&cancer_site%in%the_cancer_site)%>% select("cancer_site")

  return(list(rates_list=rates_list,cancer_sites=cancer_sites))
}

# Function to initialize the state based on the rate matrix and a given time a1
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

# Function to run the simulation
run_simulation_parallel_universe <- function(cancer_sites,
                                             LMST_vec, 
                                             OMST_vec, 
                                             test_performance_dataframe, 
                                             starting_age, 
                                             num_males, 
                                             num_females, 
                                             screening_times, 
                                             all_fits_male, 
                                             all_fits_female, 
                                             all_meta_data_fits_female, 
                                             all_meta_data_fits_male, 
                                             cdc_data, 
                                             hmd_data,
                                             MCED_cdc
                                             ){
  
  # Extract the rate matrices based on OMST and LMST specs
  rates_list_male = get_rates_list(the_omst = OMST_vec, the_lmst = LMST_vec,
                                   all_meta_data_fits = all_meta_data_fits_male,
                                   all_fits = all_fits_male, the_cancer_site = cancer_sites)
  sites_male = rates_list_male$cancer_sites %>% rename(the_cancer_site = cancer_site)
  sites_male = sites_male$the_cancer_site
  rates_list_male = rates_list_male$rates_list
  
  rates_list_female = get_rates_list(the_omst = OMST_vec, the_lmst = LMST_vec, all_meta_data_fits = all_meta_data_fits_female, all_fits = all_fits_female, the_cancer_site = cancer_sites)
  sites_female = rates_list_female$cancer_sites %>% rename(the_cancer_site = cancer_site)
  sites_female = sites_female$the_cancer_site
  rates_list_female = rates_list_female$rates_list                                    
                                             
                                             
  # Extract sensitivities and specificity based on selected cancer sites
  test_performance_male = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_male))
  test_performance_female = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_female))
  
                                             
  # Simulate males
  sim_males = sim_multiple_cancer_multiple_individuals(num_individuals = num_males,
                                                       cancer_sites = sites_male, 
                                                       rate_matrices = rates_list_male, 
                                                       early_sensitivities = test_performance_male$early_sens,
                                                       late_sensitivities = test_performance_male$late_sens, 
                                                       specificities = test_performance_male$specificities,
                                                       obs.times = screening_times, 
                                                       end.time = 500, 
                                                       start.time = starting_age, 
                                                       start.state = 1)
# start.state = sapply(rates_list_female, get_init, a1 = starting_age))
# The sapply function applies the get_init function to each rate matrix in rates_list_female, using starting_age to determine 
#   the initial state for each individual based on the model's probabilities.
  
  
  # Simulate females
  sim_females = sim_multiple_cancer_multiple_individuals(num_individuals = num_females, 
                                                         cancer_sites = sites_female, 
                                                         rate_matrices = rates_list_female, 
                                                         early_sensitivities = test_performance_female$early_sens, 
                                                         late_sensitivities = test_performance_female$late_sens, 
                                                         specificities = test_performance_female$specificities, 
                                                         obs.times = screening_times, 
                                                         end.time = 500, 
                                                         start.time = starting_age, 
                                                         start.state = 1)
                                                        
########                                           
  #Extract the rate matrices corresponding to selected cancer sites and sojourn time inputs 
  #simulate multiple individuals 
########
 
  # Combine results
  results = list(males = sim_males, females = sim_females)
  return(results)
}












