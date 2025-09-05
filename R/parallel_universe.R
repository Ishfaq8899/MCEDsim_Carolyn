#####################################################
#' Simulate a single individual with and without Multicancer Early Detection (MCED) screening
#'
#' This function simulates cancer outcomes in a single individual with a given starting age using a "parallel universe" approach.
#' That is, cancer outcomes in the population are simulated with and without MCED screening.  The natural history (i.e., the times of cancer onset and clinical diagnosis) for 
#' is the same in both screening and no-screening scenarios.  
#' The user specifies cancer sites in the MCED screening test.  The user also provides sensitivity of the tests for early and late-stage disease, where early refers to AJCC 7 stages I-II and late, III-IV, except for
#' pancreas cancer, where early is stage I and late, II-IV. 
#' Natural history models for each cancer site are specified based on a list of transition rate matrices for each site.  
#' Cancer-specific mortality is specified based on a data frame providing the parametric survival parameters for each cancer site.  Built-in tables are included in the "param_table"
#' dataframe that loads with the built-in "parametric_surv_fits" data sets. 
#' Other-cause mortality is based on other-cause life tables that can be extracted with the function "make_othercause_death_table."
#'
#' The function tracks only the first cancer diagnosis based on pre-clinical onset. 
#' Cancer-specific mortality in the screen arm assumes a stage-shift benefit of screening. That is, individuals
#' who are diagnosed in early stage under screening but who would have been diagnosed in late stage clinically, are assume to remain in early stage post lead-time. 
#' For both screened and unscreened individuals, cancer mortality is projected from the point of clinical diagnosis (i.e., post lead-time) to prevent
#' lead-time bias.
#'  
#' @param ID Numeric ID for the individual.
#' @param cancer_sites Vector of selected cancer sites.
#' @param rates_list List of transition rate matrices for each cancer site.
#' @param test_performance A list with early_sens, late_sens, and specificities for the test.
#' @param MCED_specificity Overall specificity for the MCED test (not specific to cancer site)
#' @param other_cause_death_dist A table representing other-cause mortality.
#' @param starting_age Starting age for simulation.
#' @param num_screens Number of screening rounds.
#' @param screen_interval Time interval between screenings.
#' @param end_time Ending time (age) of simulation.
#' @param sex "Male" or "Female".
#' @param surv_param_table Survival parameters table (for simulating cancer death).
#' @export
#'
#' @return A data frame with simulation results for the individual, including cancer events, screening detection, and death info.
#' 
#' @examples
#' result <- sim_individual_MCED(
#'   ID = 1,
#'   cancer_sites = c("Lung", "Colorectal"),
#'   rates_list = example_rates_list,
#'   test_performance = example_test_perf,
#'   other_cause_death_dist = example_death_dist,
#'   starting_age = 50,
#'   num_screens = 5,
#'   screen_interval = 1,
#'   end_time = 90,
#'   sex = "Male",
#'   surv_param_table = example_surv_table
#' )
sim_individual_MCED<-function( ID, 
                               cancer_sites,
                               rates_list,
                               test_performance,
                               MCED_specificity,
                               other_cause_death_dist,
                               starting_age,
                               num_screens,
                               screen_interval,
                               end_time,
                               sex,
                               surv_param_table){
  
  # simulate time of other cause death 

  other_cause_death = sim_othercause_death(other_cause_death_dist,ID=ID)
  
  
  # Get the starting states for each cancer
  start_states= unlist(lapply(rates_list, FUN="get_init",a1=starting_age))

  ### get the screening times
  screen_times = seq((starting_age),(num_screens + starting_age-1), by=screen_interval)
  
  #number of cancer sites
  num_sites=length(cancer_sites)
  
  #set cancer site specificity=1 because we are not tracking FP on a per cancer basis
  the_specificities=rep(1,times=num_sites)
  
  result <- sim_multiple_cancer_indiv(ID = ID,
                                        cancer_sites = cancer_sites,
                                        rate_matrices = rates_list,
                                        early_sensitivities = test_performance$early_sens,
                                        late_sensitivities =  test_performance$late_sens,
                                        specificities = the_specificities,
                                        obs.times = screen_times,
                                        start.time = starting_age,
                                        end.time = end_time, 
                                        start.states =  start_states) 
  
  # Add other-cause death info
  result$other_cause_death_status <- other_cause_death$status
  result$other_cause_death_time <- other_cause_death$time
  result$sex=sex
  #-----------------------  
  # Identify first cancer by onset time
  #-----------------------
  result$is_first_cancer=FALSE
  result$cancer_death_time_no_screen=NA
  result$cancer_death_time_screen=NA
  
   if (sum(!is.na(result$onset_time)>=1)) {
    
    first_cancer_row <- result %>%
      filter(!is.na(onset_time)) %>%
      # selects the single row with the earliest (smallest) onset_time.    
      slice_min(order_by = onset_time, n = 1, with_ties = FALSE)%>%mutate(is_first_cancer=TRUE)
  
    #Simulate time of death for first cancer without screening and with screening
    #if stage at clinical diagnosis is same as stage at screen diagnosis OR there is is no screen diagnosis, then cancer_death_time_screen=cancer_death_time_no_screen
    #otherwise, then cancer_death_time_screen=clinical_diagnosis_time+sim_cancer_death_param(the_stage=screen_diagnosis_stage,
                                                       # the_cancer_site=cancer_site,
                                                       #the_sex=sex,
                                                      #the_model_type="Loglogistic",
                                                       #param_table=surv_param_table)

    if(!is.na(first_cancer_row$clinical_diagnosis_stage)){
     
    first_cancer_row <-first_cancer_row  %>%mutate(cancer_death_time_no_screen=clinical_diagnosis_time+sim_cancer_death_param(the_stage=clinical_diagnosis_stage,
                                                                                                       the_cancer_site=cancer_site,
                                                                                                       the_sex=sex,
                                                                                                       the_model_type="Loglogistic",
                                                                                                       param_table=surv_param_table,ID=ID))%>%
                                    mutate(cancer_death_time_screen=ifelse((screen_diagnosis_stage!=clinical_diagnosis_stage&!is.na(screen_diagnosis_stage))&!
                                                                             is.na(clinical_diagnosis_stage),
                                                                           clinical_diagnosis_time+
                                                                             sim_cancer_death_param(the_stage="Early",
                                                                                                    the_cancer_site=cancer_site,
                                                                                                    the_sex=sex,
                                                                                                    the_model_type="Loglogistic",
                                                                                                    param_table=surv_param_table,
                                                                                                    ID=ID),
                                                                           cancer_death_time_no_screen))
     
      
                 
      
    }
    result<-first_cancer_row
 
  
  }else{
    result<-slice_head(result,n=1)%>%mutate(cancer_site=NA)
  } 

  set.seed(ID)
  result<-result %>% mutate(FP_tot=rbinom(n(),size=total_no_canc_screens,prob=1-MCED_specificity))

  return(result)
}


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
#' @param cancer_sites Vector of cancer sites (allowable values include "Anus",  ” Bladder" ,  ”Breast",  ”Colorectal" , "Esophagus"  , ”Gastric", "Headandneck", 
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
sim_multiple_individuals_MCED_parallel_universe <- function(cancer_sites,
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
  
 # browser()
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
  
  
  
  
  
  
  # Combine all individual results into one data frame
  combined_results_males <- do.call(rbind, results_list_male)%>%mutate(sex="Male")
  combined_results_females <- do.call(rbind, results_list_female)%>%mutate(sex="Female")

    combined_results=bind_rows(combined_results_males,combined_results_females)%>%
                   mutate(start_age=starting_age,end_time=ending_age)
  
  
  #Simulate times of cancer_death_time screen for individuals diagnosed early via screening but late via clinical diagnosis
   # temp=filter(combined_results,(screen_diagnosis_stage!=clinical_diagnosis_stage&!is.na(screen_diagnosis_stage))&!is.na(clinical_diagnosis_stage))
   # temp$cancer_death_time_screen=temp$clinical_diagnosis_time+mapply(temp$screen_diagnosis_stage,temp$cancer_site,temp$sex,temp$ID*8,FUN="sim_cancer_death_param",
   #             MoreArgs=list(the_model_type="Loglogistic",
   #                           param_table=surv_param_table))
   # 
  # combined_results=bind_rows(temp,subset(combined_results,!ID%in%temp$ID))

    
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
    
  return(combined_results)
}






