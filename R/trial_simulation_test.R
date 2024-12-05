
############################################
# Load data
load("/Users/ahmadish/Desktop/OHSU project/Ishfaq folder/MCEDsim/R/lungoutfile.Rdata")
load("/Users/ahmadish/Desktop/OHSU project/Ishfaq folder/MCEDsim/R/liveroutfile.Rdata")

liverate <- liverout[[1]]$rate.matrix
lungrate <- lungout[[1]]$rate.matrix

early_sens_lung <- 0.3
late_sens_lung <- 0.9
early_sens_liver <- 0.2
late_sens_liver <- 0.8
specificity_lung <- 0.85
specificity_liver <- 0.9

############# Define emission matrices for lung and liver #############################
emission_matrix_lung <- create_emission_matrix(lungrate, early_sens_lung, late_sens_lung, specificity_lung)
emission_matrix_liver <- create_emission_matrix(liverate, early_sens_liver, late_sens_liver, specificity_liver)


############################# Test the function for an individual ####################################
#test1 <- sim_multiple_cancer_indiv(ID = 1,
#                                  cancer_sites = c("Lung", "Liver"),
#                                  rate_matrices = list(lungrate, liverate),
#                                  early_sensitivities = c(early_sens_lung, early_sens_liver),
#                                  late_sensitivities = c(late_sens_lung, late_sens_liver),
#                                  specificities = c(specificity_lung, specificity_liver),
#                                  obs.times = c(0, seq(50, 75, by = 2)),
#                                  end.time = 100,
#                                  start.time = 0,
#                                  start.state = 1)

########### Test the function for multiple individuals ####################
test2 <- sim_multiple_cancer_multiple_individuals(num_individuals = 10000,
                                                  cancer_sites = c("Lung", "Liver"),
                                                  rate_matrices = list(lungrate, liverate),
                                                  early_sensitivities = c(early_sens_lung, early_sens_liver),
                                                  late_sensitivities = c(late_sens_lung, late_sens_liver),
                                                  specificities = c(0,0),
                                                  obs.times = c(seq(50, 75, by = 2)),
                                                  end.time = 1000,
                                                  start.time = 0,
                                                  start.state = 1)




