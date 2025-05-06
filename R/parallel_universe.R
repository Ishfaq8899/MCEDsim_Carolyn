
#' Get selected rate matrices for all cancer sites based on common OMST, and LMST specifications (given list of all rate matrices) 
#'
#' This function simulates the progression of multiple cancer sites for a single individual over a specified time period.
#'
#' @param the_omst          Overall mean sojourn time 
#' @param the_lmst          Late mean sojourn time
#' @param all_meta_data_fits Dataframe with site, OMST, DMST, and index corresponding to all_fits
#' @param all_fits      list of all rate matrices whose indices correspond to all_meta_data_fitsf
#' @param the_cancer_site specified cancer site
#' @return List of selected rate matrices. 
get_rates_list<-function(the_omst, the_lmst, all_meta_data_fits,all_fits,the_cancer_site){
  the_indices=all_meta_data_fits %>% filter(OMST%in% the_omst&LMST%in%the_lmst&cancer_site%in%the_cancer_site)%>% select("index")
  rates_list=lapply(all_fits[unlist(the_indices)],"[[","rate.matrix")
  cancer_sites=  the_indices=all_meta_data_fits %>% filter(OMST%in% the_omst&LMST%in%the_lmst&cancer_site%in%the_cancer_site)%>% select("cancer_site")

  return(list(rates_list=rates_list,cancer_sites=cancer_sites))
}


run_simulation_parallel_universe(cancer_sites,  
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

  #Extract the rate matrices based on OMST and LMST specs
  rates_list_male=get_rates_list(the_omst=OMST_vec,the_lmst=LMST_vec,all_meta_data_fits=all_meta_data_fits_male,
                                 all_fits=all_fits_male,the_cancer_site=cancer_sites_vec)
  sites_male=rates_list_male$cancer_sites%>%rename(the_cancer_site=cancer_site)
  sites_male=sites_male$the_cancer_site
  rates_list_male=rates_list_male$rates_list
  
  rates_list_female=get_rates_list(the_omst=OMST_vec,the_lmst=LMST_vec,all_meta_data_fits=all_meta_data_fits_female,
                                   all_fits=all_fits_female,the_cancer_site=cancer_sites_vec)
  sites_female=rates_list_female$cancer_sites%>%rename(the_cancer_site=cancer_site)
  sites_female=sites_female$the_cancer_site
  rates_list_female=rates_list_female$rates_list
  
  #Extract sensitivities and specificity based on selected cancer sites

  test_performance_male=test_performance_dataframe%>%filter(cancer_site%in%as.vector(sites_male))
  test_performance_female=test_performance_dataframe%>%filter(cancer_site%in%as.vector(sites_female))

  
  #simulate males
  sim_multiple_cancer_multiple_individuals(num_individuals=num_males, 
                                           cancer_sites=sites_male, 
                                           rate_matrices=rates_list_male,
                                           early_sensitivities=test_performance_male$early_sens,
                                           late_sensitivities = test_performance_male$late_sens,
                                           specificities = test_performance_male$specificities,
                                           obs.times, 
                                           end.time,
                                           start.time, start.state)
    
  #simulate females
  
  
  #Extract the rate matrices corresponding to selected cancer sites and sojourn time inputs 
  #simulate multiple individuals 
  #Compute different in life years with and without screening 
  #Difference in late stage diagnosis with and without screening
} 











