#' Simulate a time based on a provided survival function.
#' 
#' @param time vector of times corresponding to the survival function
#' @param surv vector of survival probabilities at each time. 
#' @return A data frame with "time" and "status" elements, where status=1 indicates the event occurred. 
#' @export
gettime<-function(time,surv){
  
  unif1=runif(n=1)
  cdf=1-surv
  
  eventtime=approx(x=cdf, y=time, xout=unif1, rule=2)$y
  
  eventstatus=1
 
  return(data.frame(time = eventtime, status = eventstatus))
}

###########################################################################################
#' Simulate a time of cancer death for a specific cancer
#' Simulate a time of cancer death for a specific cancer in a specific stage based on provided survival distributions.
#' Uses the inverse CDF method.
#' @param the_stage stage at clinical diagnosis
#' @param the_cancer_site cancer site (one of "Anus","Breast","Bladder","Colorectal","Esophagus","Headandneck","Gastric","Liver" ,"Lung","Pancreas", "Prostate", "Renal", "Ovary", "Uterine")
#' @param the_sex Male or Female
#' @param the_model_type Weibull or  Loglogistic
#' @param cancer_survival_dist: a data frame with columns: surv, time, site, stage, sex, model_type 
#
# OUTPUTS: A numeric value representing the simulated time of death due to a specific cancer
#################################################################################################
sim_cancer_death <- function(the_stage, the_cancer_site, the_sex,the_model_type, cancer_survival_dist){
  
  # Filter the survival distribution based on the type and stage
  survival_dist_indiv = filter(cancer_survival_dist,cancer_site == paste(the_cancer_site), stage==paste(the_stage),
                               sex==paste(the_sex),model_type==paste(the_model_type))

 
  # Get the survival time based on the distribution
  death_info = gettime(time = survival_dist_indiv$time, surv = survival_dist_indiv$surv)
  death_time = death_info$time
  
  return(death_time)
}

#' Simulate a time of cancer death for a specific cancer in a specific stage based on provided survival distributions.
#' 
#' @param the_stage stage at clinical diagnosis
#' @param the_cancer_site cancer site (one of "Anus","Breast","Bladder","Colorectal","Esophagus","Headandneck","Gastric","Liver" ,"Lung","Pancreas", "Prostate", "Renal", "Ovary", "Uterine")
#' @param the_sex Male or Female
#' @param the_model_type Weibull or  Loglogistic
#' @param param_table_dist a data frame with columns: intercept (based on survreg fit), scale (based on survreg fit) site, stage, sex, model_type 
#' @return A data frame with "time" and "status" elements, where status=1 indicates the event occurred. 
#' @export
sim_cancer_death_param <- function(the_stage, the_cancer_site, the_sex,the_model_type, param_table){
  
 
  # Filter the survival distribution based on the type and stage
  survival_dist_indiv = filter(param_table,cancer_site == paste(the_cancer_site), stage==paste(the_stage),
                               sex==paste(the_sex),model_type==paste(the_model_type))
  
  if(survival_dist_indiv$model_type=="Loglogistic"){
    the_survobj=s_loglogistic(intercept = survival_dist_indiv$intercept, scale =survival_dist_indiv$scale)
  }
  
  if(survival_dist_indiv$model_type=="Weibull"){
    the_survobj=s_weibull(intercept = survival_dist_indiv$intercept, scale =survival_dist_indiv$scale)
  }
  
  death_time=rsurv(the_survobj,n=1)
  return(death_time)
}


