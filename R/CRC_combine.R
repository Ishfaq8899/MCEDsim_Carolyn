###########################################################
#' Simulate a cohort of individuals with and without Multicancer Early Detection (MCED) screening
#'
#' This function simulates cancer outcomes in a population with a designated starting age using a "parallel universe" approach.
#' That is, cancer outcomes in the population are simulated with and without MCED screening.  The natural history (i.e., the times of cancer onset and clinical diagnosis) for 
#' each individual is the same in both screening and no-screening scenarios.  
#' The user specifies cancer sites in the MCED screening test.  The user also provides sensitivity of the tests for early and late-stage disease, where early refers to AJCC 7 stages I-II and late, III-IV, except for
#' pancreas cancer, where early is stage I and late, II-IV. 
#' Natural history models are based on built-in fitted models that are calibrated to SEER 2015-2021 data by age and sex and can be specified based on 
#' user-provided inputs about the overall mean sojourn time (OMST) and the late mean sojourn time (LMST) for each cancer site. 
#' The function tracks only the first cancer diagnosis based on pre-clinical onset. 
#' 
#' Other-cause mortality is based on all cause mortality tables from the Human Mortality Database that have been adjusted to remove the mortality due to
#' cancers included in the MCED tests.  Cancer-specific mortality in the screen arm assumes a stage-shift benefit of screening. That is, individuals
#' who are diagnosed in early stage under screening but who would have been diagnosed in late stage clinically, are assume to remain in early stage post lead-time. 
#' For both screened and unscreened individuals, cancer mortality is projected from the point of clinical diagnosis (i.e., post lead-time) to prevent
#' lead-time bias.
#'  
#' @param cancer_sites Vector of cancer sites (allowable values include "Anus", ” Bladder" ,  ”Breast",  ”Colorectal" , "Esophagus"  , ”Gastric", "Headandneck", 
#' "Liver" ,  "Lung",   "Pancreas", "Prostate", "Renal", "Uterine","Ovary")
#' @param LMST_vec Vector of late mean sojourn times for each cancer site.
#' @param OMST_vec Vector of overall mean sojourn times for each cancer site.
#' @param test_performance_dataframe Data frame with test sensitivity/specificity info.
#' @param MCED_specificity Overall specificity for the MCED test (not specific to cancer site)
#' @param starting_age Starting age for simulation.
#' @param ending_age Ending age for simulation.
#' @param num_screens Number of screening rounds.
#' @param screen_interval Interval between screening rounds.
#' @param num_males Number of male individuals to simulate.
#' @param num_females Number of female individuals to simulate.
#' @param all_rates_male List of transition matrices for males.
#' @param all_rates_female List of transition matrices for females.
#' @param all_meta_data_female Metadata for female cancer sites.
#' @param all_meta_data_male Metadata for male cancer sites.
#' @param cdc_data CDC mortality data.
#' @param hmd_data Human Mortality Database data.
#' @param MCED_cdc CDC data for MCED.
#' @param surv_param_table Survival parameters table.
#' @param CRC_data CRC data.
#' @export
#'
#' @return A data frame with combined simulated results for all individuals.
#' The function returns each individual's first cancer site,  age and stage of clinical diagnosis in 
#' absence of screening, age and stage at screen diagnosis, the time of other-cause mortality, and the time of cancer-specific mortality in absence
#' and presence of screening.   Cancer diagnosis and death times are presented both with and without competing other-cause mortality.

#' @examples
#' library(MCEDsim)
#' #Load the other-cause mortality tables
#' data("cdc_hmd_data")
#' #Load the prefitted natural history models
#' data("combined_fits")
#' #Load the prefitted cause-specific survival models
#' data("parametric_surv_fits")
#'
#' #Specify inputs
#' #Cancer sites in the MCED test
#' cancer_sites_vec=c("Anus",  "Bladder",  "Breast", "Esophagus",
#'                   "Gastric", "Headandneck","Liver" , 
#'                   "Lung",   "Pancreas", "Prostate", "Renal",  "Uterine","Ovary")
#'
#' #overall mean sojourn time=2 years for all sites
#' OMST_vec=rep(2,times=13)
#' #late mean sojourn time=.5 years for all sites
#' LMST_vec=rep(.5,times=13)
#' #early-stage sensitivity
#' early_sens=rep(.65,times=13)
#' #late-stage sensitivity
#' late_sens=rep(.95,times=13)
#' #dataframe with all of the test performance inputs
#' test_performance_dataframe = data.frame(early_sens, late_sens,cancer_site = cancer_sites_vec)
#'
#' results <- sim_multiple_individuals_MCED_parallel_universe(cancer_sites = cancer_sites_vec,
#'                                                           LMST_vec = LMST_vec,
#'                                                           OMST_vec = OMST_vec,
#'                                                           test_performance_dataframe = test_performance_dataframe,
#'                                                           starting_age = 45,
#'                                                           ending_age=500,
#'                                                           num_screens = 30,
#'                                                           screen_interval = 1,
#'                                                           num_males = 20,
#'                                                           num_females = 20,
#'                                                           all_rates_male = all_rates_male,
#'                                                           all_rates_female = all_rates_female,
#'                                                           all_meta_data_female = all_meta_data_female,
#'                                                           all_meta_data_male = all_meta_data_male,
#'                                                           cdc_data = all_cause_cdc,
#'                                                           hmd_data = hmd_data,
#'                                                           MCED_cdc = MCED_cdc,
#'                                                           surv_param_table=param_table,
#'                                                           MCED_specificity = .995)
sim_MCED_parallel_universe_before_CRC <- function(cancer_sites,
                                                            LMST_vec, 
                                                            OMST_vec, 
                                                            test_performance_dataframe, 
                                                            MCED_specificity,
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
                                                            surv_param_table)
{
  
  #browser()
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
  
  #Get the other-cause death tables for men and women
  other_cause_death_male=make_othercause_death_table(cdc_data=cdc_data,
                                                     MCED_cdc=MCED_cdc, 
                                                     hmd_data=hmd_data,
                                                     the_starting_age = starting_age,
                                                     the_sex="Male",
                                                     selected_cancers=sites_male,
                                                     the_year=2018)
  
  other_cause_death_female=make_othercause_death_table(cdc_data=cdc_data,
                                                       MCED_cdc=MCED_cdc, 
                                                       hmd_data=hmd_data,
                                                       the_starting_age = starting_age,
                                                       the_sex="Female",
                                                       selected_cancers=sites_female,
                                                       the_year=2018)
  
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
                                              sex="Male",MCED_specificity=MCED_specificity),
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
                                                sex="Female",
                                                MCED_specificity=MCED_specificity),
                                SIMPLIFY = FALSE)
  
  
  
 
  return(list(results_list_female=results_list_female,results_list_male=results_list_male))
}
  
  


#' @export
combine_MCED_CRC<- function(results_list_female,
                            results_list_male,
                            starting_age,
                            ending_age,
                            CRC_data,
                            surv_param_table){
  
  #get the first cancer and additional cancers for all individuals
  first_site_female=lapply(results_list_female,"[[","first_result")
  additional_sites_female=lapply(results_list_female,"[[","stored_result")
  
  first_site_male=lapply(results_list_male,"[[","first_result")
  additional_sites_male=lapply(results_list_male,"[[","stored_result")
  
  # Combine all individual results (first cancers) 
  combined_first_results_males <- do.call(rbind, first_site_male)%>%mutate(sex="Male")
  combined_first_results_females <- do.call(rbind, first_site_female)%>%mutate(sex="Female")
  combined_first_results=bind_rows(combined_first_results_males,combined_first_results_females)%>%
    mutate(start_age=starting_age,end_time=ending_age)
  
  # Combine all individual results (additional cancers) 
  combined_additional_results_males <- do.call(rbind, additional_sites_male)%>%mutate(sex="Male")
  combined_additional_results_females <- do.call(rbind, additional_sites_female)%>%mutate(sex="Female")
  combined_additional_results=bind_rows(combined_additional_results_males,combined_additional_results_females)%>%
    mutate(start_age=starting_age,end_time=ending_age)
  
  
  #Note: consider if we want to use cancer onset or clinical diagnosis for additional cancers for purpose of calculating over diagnosis. 
  #If we decide this is important, change clinical_diagnosis_time to cancer_onset_time in the subsequent code
  
  #Filter additional cancers for those whose clinical diagnosis is prior to other cause death
  combined_additional_results <-combined_additional_results %>% filter(clinical_diagnosis_time <=other_cause_death_time)
  
  #simulate cancer-specific deaths for additional cancers
  addtl_cancer_deaths=mapply(FUN="sim_cancer_deaths_screen_no_screen",clinical_diagnosis_time=combined_additional_results$clinical_diagnosis_time,
                             clinical_diagnosis_stage=combined_additional_results$clinical_diagnosis_stage,
                             cancer_site=combined_additional_results$cancer_site,
                             sex=combined_additional_results$sex,
                             ID=combined_additional_results$ID,
                             screen_diagnosis_stage=combined_additional_results$screen_diagnosis_stage,
                             MoreArgs=list(surv_param_table=surv_param_table),SIMPLIFY=F)  
  
  #Join cancer-specific deaths with cancer diagnoses for additional cancers
  combined_additional_results=data.frame(do.call(rbind,addtl_cancer_deaths))%>%inner_join(combined_additional_results, by=c("ID","cancer_site"))
  
  #combine the CRC data with the additional cancers for reassignment.  CRC diagnoses that occur after other cause death do not 
  #need to be reassigned so these people are removed from combined_additional_results.
  combined_additional_results <- bind_rows(combined_additional_results, CRC_data)%>%
    mutate(age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0,150,by=5)))%>% 
    filter(clinical_diagnosis_time <=other_cause_death_time)
  
  #Identify people who have clinical diagnosis of first cancer prior to other cause death  
  primary_cancer <- combined_first_results %>% filter(clinical_diagnosis_time<=other_cause_death_time)
  
  #Identify people who do not have clinical diagnosis of first cancer prior to other cause death.  
  #These people are eligible for reassignment of additional cancers based on matching age at OC death. 
  #Define other cause death strata based on five year age groups. 
  no_primary_cancer <- combined_first_results %>% filter(clinical_diagnosis_time>other_cause_death_time)%>%
    mutate(age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0,150,by=5)))
  
  
  #Counts of number of multiple cancers, primary cancers in lifeteime, and no primary cancers in lifetime
  N_multiple_cancers=nrow(combined_additional_results)
  N_primary_cancers=nrow(primary_cancer)
  N_no_primary_cancer=nrow(no_primary_cancer)
  
  #set index as ID to identify specific individual without primary cancer (used in bookkeeping in next step)
  no_primary_cancer <- no_primary_cancer %>% mutate(index=seq(1,nrow(no_primary_cancer)))
  
  #This loop goes through the combined_additional_results rows and attempts to match in an individual in
  #no_primary_cancer within the same OC death strata.  Then it adds the combined_additional_results row to primary_cancer
  #and removes the matching row from no_primary cancer, since they are not eligible to matched again. 
  #Finally, it removes the row from combined_additional_results and updates the no_match_counter.  
  no_match_counter=0
  while(nrow(combined_additional_results)>0){
    test=match_individual(i=1,combined_additional_results=combined_additional_results, no_primary_cancer=no_primary_cancer)
    
    if(length(test)>1){
      primary_cancer=bind_rows(combined_additional_results[1,],primary_cancer)
      no_primary_cancer=no_primary_cancer%>%filter(index!=test$index)
    }else{
      no_match_counter=no_match_counter+1
    }
    combined_additional_results=combined_additional_results[-1,]
    
  }
  
  
  #check to see if the added and subtracted rows match expectations based on previous loop. 
  N_multiple_cancers_2=nrow(combined_additional_results)
  N_primary_cancers_2=nrow(primary_cancer)
  N_no_primary_cancer_2=nrow(no_primary_cancer)
  
  #Final data with all cancers combined (first cancers and reassigned cancers)
  combined_results=bind_rows(primary_cancer,no_primary_cancer) 
  
  #Process data with other cause death as a censoring event
  
  #Ascertain age at at screen and clinical diagnosis in presence of other cause death
  #Ascertain age at death under screening scenarios and no screening scenarios in presence of other cause death
  #Ascertain age at diagnosis under screening scenario (can be either screen or clinical)
  #Ascertain mode of diagnosis under screening scenario
  #Ascertain stage at diagnsosis under screening scenario
  #Ascertain if individual was overdiagnosed (screen detected but died due to other causes prior to clinical diagnosis)
  combined_results=combined_results %>% mutate(clin_dx_age = pmin(other_cause_death_time,clinical_diagnosis_time,end_time,na.rm = T),
                                               clin_dx_event = case_when(
                                                 clin_dx_age == other_cause_death_time ~ "other_cause_death",
                                                 clin_dx_age == end_time ~ "censor",
                                                 clin_dx_age == clinical_diagnosis_time ~ "clin_cancer_diagnosis",
                                                 .default = NA
                                               ),
                                               clin_dx_event_stage = case_when(clin_dx_event == "clin_cancer_diagnosis" & clinical_diagnosis_stage == "Early"~1,
                                                                               clin_dx_event == "clin_cancer_diagnosis" & clinical_diagnosis_stage == "Late"~2,
                                                                               .default = 3),
                                               screen_dx_age = pmin(other_cause_death_time,screen_diagnosis_time,end_time,na.rm = T),
                                               screen_dx_event = case_when(
                                                 screen_dx_age == other_cause_death_time ~ "other_cause_death",
                                                 screen_dx_age == end_time ~ "censor",
                                                 screen_dx_age == screen_diagnosis_time ~ "screen_cancer_diagnosis",
                                                 .default = NA
                                               ),
                                               screen_dx_event_stage = case_when(screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Early"~1,
                                                                                 screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Late"~2,
                                                                                 .default = 3),
                                               death_age_no_screen=pmin(other_cause_death_time,cancer_death_time_no_screen,end_time,na.rm = T),
                                               death_age_screen=pmin(other_cause_death_time,cancer_death_time_screen,end_time,na.rm = T),
                                               death_event_no_screen=case_when(
                                                 death_age_no_screen == other_cause_death_time ~ "other_cause_death",
                                                 death_age_no_screen == end_time ~ "censor",
                                                 death_age_no_screen == cancer_death_time_no_screen ~ "cancer_death",
                                                 .default = NA
                                               ),
                                               death_event_screen=case_when(
                                                 death_age_screen == other_cause_death_time ~ "other_cause_death",
                                                 death_age_screen == end_time ~ "censor",
                                                 death_age_screen == cancer_death_time_screen ~ "cancer_death",
                                                 .default = NA
                                               ),
                                               diagnosis_age_screen_scenario=pmin(clin_dx_age,screen_dx_age,na.rm=T),
                                               diagnosis_event_screen_scenario=ifelse(screen_dx_age<=clin_dx_age, screen_dx_event,
                                                                                      clin_dx_event),
                                               diagnosis_event_stage_screen_scenario=case_when(screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Early"~1,
                                                                                               screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Late" ~2,
                                                                                               (screen_dx_event!="screen_cancer_diagnosis" & clin_dx_event=="clin_cancer_diagnosis") & clinical_diagnosis_stage=="Early"~1,
                                                                                               (screen_dx_event !="screen_cancer_diagnosis" & clin_dx_event=="clin_cancer_diagnosis") & clinical_diagnosis_stage=="Late"~2,
                                                                                               .default = 3),
                                               life_years_diff=death_age_screen-death_age_no_screen,
                                               overdiagnosis=ifelse(screen_dx_event=="screen_cancer_diagnosis"&clin_dx_event=="other_cause_death",1,0)
                                               
                                               
  )
  
  return(list(results = combined_results,no_match_counter = no_match_counter))
}
