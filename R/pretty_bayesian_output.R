#' Fancy Table Output Bayesian Linear and Logistic regression  
#' pretty_bayesian_output takes a stanreg fit object  ( family = binomial or guassian), and calculate estimates of 
#' regression coefficients or odds ratios. 
#'
#' @param fit  stanreg 
#' @param model_data data.frame or tibble  used to create model fits. Used for capturing variable labels, if they exist
#' 
#' @details 
#' 
#' Model type is determined by \code{fit} class, and also family if. If family is binomial,
#'  then the output is designed for a Logistic model (i.e. Odd Ratios), otherwise the output is designed for a linear model.
#'
#' @return
#' 
#' A tibble with:  \code{Variable}, \code{Level}, \code{Mean Est/OR },
#'  \code{mcse} (mcse monte carlo standard error), 
#'  \code{sd} (standard deviation), 
#'  \code{2.5} (percentile of posterior), 
#'  \code{50} (percentile of posterior), 
#'  \code{97.5} (percentile of posterior), 
#'  \code{n_eff} (effective number of simulation draws), 
#'  \code{Rhat} ( is essentially the ratio of between-chain variance to within-chain variance analogous to ANOVA.). 
#'   
#' @examples
#' 
#' # Basic linear model example
#' set.seed(542542522)
#' ybin <- sample(0:1, 100, replace = TRUE)
#' y <- rexp(100,.1)
#' x1 <- rnorm(100)
#' x2 <- y + rnorm(100)
#' x3 <- factor(sample(letters[1:4],100,replace = TRUE))
#' my_data <- data.frame(y, ybin, x1, x2, x3)
#' library(rstanarm)
#' lm_fit <- stan_glm(y ~ x1+ x2 + x3, data = my_data )
#' pretty_bayesian_output(fit = lm_fit, model_data = my_data)
#' 
#' library(dplyr)
#' # Logistic Regression
#' my_fit <- stan_glm(ybin ~ x1 + x2 + x3, data = my_data, family = binomial)
#' my_pretty_model_output <-pretty_bayesian_output(fit = my_fit, model_data = my_data)
#' 
#' # Printing of Fancy table in HTML
#' 
#' kableExtra::kable(my_pretty_model_output, 'html', caption = 'My Table') %>% 
#'    kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack')
#'   
#' # Real World Examples
#' data(Bladder_Cancer)
#' #surv_obj <- survival::Surv(Bladder_Cancer$Survival_Months, Bladder_Cancer$Vital_Status == 'Dead')  
#' #my_fit <- survival::coxph(surv_obj ~ Gender + PT0N0, data = Bladder_Cancer)
#' #my_output <- pretty_model_output(fit = my_fit, model_data = Bladder_Cancer)
#' #kableExtra::kable(my_output, 'html') %>% 
#' #    kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack')
#'   
#' @importFrom Hmisc label label<-
#' @importFrom utils globalVariables
#' 

#' @export


pretty_bayesian_output <- function(fit, model_data ) {
  
  if (any(class(fit) == 'stanreg') && fit$family$family %in% c('binomial')) {
    # Bayesian Logistic Regression 
    exp_output <- TRUE
    est_name <- 'OR'
  } else if (any(class(fit) == 'stanreg') && fit$family$family %in% c('???')) {
    # KEEPING THESE FOR NOW
    # Coxph Regression
    exp_output <- TRUE
    est_name <- 'HR'
  } else {
    # Not Logistic Regression or Coxph
    exp_output <- FALSE
    est_name <- 'Est'
  } 
  
  #Using Variable labels for output, if no label using variable name used later ####
  
  var_names <- all.vars(fit$terms)[-1] # remove outcome Y
  
  if (!all(var_names %in% colnames(model_data)))stop('All variables used in the "fit" must be in the "model_data" dataset')
  
  var_labels <- Hmisc::label(model_data)[var_names]
  
  if (any(var_labels == '')){ var_labels[var_labels == ''] <- gsub('_', ' ', var_names[var_labels == ''])}
  
  neat_fit = as.data.frame(summary(fit, digits = 2, prob=c(.025, .5, .975)))
  if (exp_output){  neat_fit[,c(1,3,4,5,6)] = exp(neat_fit[,c(1,3,4,5,6)]) }# exp to get the odds ratios 
  
  neat_fit = round(dplyr::slice(neat_fit, 1:(dplyr::n()-2)),3)#   remove bottom two rows
  colnames(neat_fit)[1] = "Mean OR"
  
  neat_fit <- neat_fit %>%  mutate(variable = row.names(neat_fit)) 
  
 
  if (length(fit$xlevels) > 0) {
    
    all_levels0 = fit$xlevels %>% tibble::enframe() %>% tidyr::unnest(.,cols = c("value")) %>%  
      dplyr::mutate(variable = paste0(name, value))
 
    continuousvar<-tibble(name = setdiff(var_names, all_levels0$name),
                          value = rep("",length(setdiff(var_names, all_levels0$name))),
                          variable = setdiff(var_names, all_levels0$name))
    all_levels <- dplyr::bind_rows(all_levels0, continuousvar)
 
    all_levels <- all_levels %>% dplyr::arrange(match(name,var_names))
  
    neat_fit <- dplyr::full_join(all_levels, neat_fit , by = "variable") %>%
      dplyr::mutate( name = ifelse(is.na(name), variable, name),
                     value = ifelse(is.na(value), '', value) )
  } else {
    neat_fit <- neat_fit %>%
      dplyr::mutate( name = variable,
                     value = ''  )
  }
  
  #  print("#1")
  # print(neat_fit) 
  neat_fit <- neat_fit  %>%
    dplyr::mutate(Variable2 = var_labels[match(name,var_names)] 
                  ,Variable2 = ifelse(is.na(Variable2), name, Variable2))  %>% 
    dplyr::select( -variable)  
  
  # Label the reference groups 
  var_names <- c("(Intercept)", all.vars(fit$terms))
  # print("the vars:")
  # print(var_names)
  neat_fit <- neat_fit %>%
    dplyr::mutate( est.label = ifelse(is.na(`Mean OR`), "1.0 (Reference)", `Mean OR` ),
                   oldname = name) %>% 
    dplyr::select( -`Mean OR`, -name) %>% #,-variable %>%  
    dplyr::select(name = Variable2 , Level = value, `Mean OR` = est.label,  everything()) %>% 
    dplyr::arrange(factor(name, levels = var_names))  %>% 
    dplyr::arrange(match(oldname, var_names))  %>%
    dplyr::select( -oldname) %>% 
    # rename estimate if OR or HR or lm model coefficient = Est eventually.    
    dplyr::rename(!!paste0("Mean ", est_name  ) := `Mean OR`)   
  
  
  neat_fit$name <- gsub("_"," ",neat_fit$name )
  neat_fit$name[neat_fit$name=="(Intercept)"] <- "Intercept";
  
  return(neat_fit)
  
} 

 

globalVariables(c("Mean OR","everything","n" ,"ROPE CI"))
 


#' Fancy Table of Tests of Bayesian Linear and Logistic regression 
#' 
#' pretty_bayesian_regression_tests takes a stanreg fit object  ( family = binomial or guassian), and conducts tests  
#'
#' @param fit  stanreg (currently only tested on linear model and  logistic stan_glm fit)
#' @param model_data data.frame or tibble  used to create model fits. Used for capturing variable labels, if they exist
#' 
#' @details 
#' 
#' Model type is determined by \code{fit} class, and also family if. If family is binomial,
#'  then the output is designed for a Logistic model (i.e. Odd Ratios), otherwise the output is designed for a linear model.
#'
#' @return
#' 
#' A tibble with: \code{Mean Est/OR/HR(CI)},
#'  \code{pd} (probabiity of direction), 
#'  \code{ROPE (CI)} (R.O.P.E and CI), 
#'  \code{2.5} (R.O.P.E. %),  
#'   
#' @examples
#' 
#' # Basic linear model example
#' set.seed(542542522)
#' ybin <- sample(0:1, 100, replace = TRUE)
#' y <- rexp(100,.1)
#' x1 <- rnorm(100)
#' x2 <- y + rnorm(100)
#' x3 <- factor(sample(letters[1:4],100,replace = TRUE))
#' my_data <- data.frame(y, ybin, x1, x2, x3)
#' library(rstanarm)
#' lm_fit <- stan_glm(y ~ x1+ x2 + x3, data = my_data )
#' pretty_bayesian_regression_tests(fit = lm_fit, model_data = my_data)
#' 
#' library(dplyr)
#' # Logistic Regression
#' my_fit <- stan_glm(ybin ~ x1 + x2 + x3, data = my_data, family = binomial)
#' my_pretty_model_output_tests <-pretty_bayesian_regression_tests(fit = my_fit, model_data = my_data)
#' 
#' # Printing of Fancy table in HTML
#' 
#' kableExtra::kable(my_pretty_model_output_tests, 'html', caption = 'My Table') %>% 
#'    kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack')
#'   
#'   
#' @importFrom Hmisc label label<-
#' @importFrom utils globalVariables
#' 

#' @export



pretty_bayesian_regression_tests <- function(fit, model_data){
  
  
  if (any(class(fit) == 'stanreg') && fit$family$family %in% c('binomial')) {
    # Bayesian Logistic Regression 
    exp_output <- TRUE
    est_name <- 'OR'
  } else if (any(class(fit) == 'stanreg') && fit$family$family %in% c('???')) {
    # KEEPING THESE FOR NOW
    # Coxph Regression
    exp_output <- TRUE
    est_name <- 'HR'
  } else {
    # Not Logistic Regression or Coxph
    exp_output <- FALSE
    est_name <- 'Est'
  }  
  
  
  PT0N0_bayesTest <- bayestestR::describe_posterior(fit)
  tempnames <- PT0N0_bayesTest[,1]
  # This will remove the first column since the rownames are already stored add back in below.
  PT0N0_bayesTest <- PT0N0_bayesTest[,-1]
  PT0N0_bayesTest <- round(PT0N0_bayesTest,3)
  
  neat_fit <- as.data.frame(PT0N0_bayesTest %>% mutate(variable = tempnames, 
                                                       MedianOR =   case_when(exp_output == TRUE ~ round(exp(Median),2),
                                                                              exp_output == FALSE ~ round(Median,2)),# round(exp(Median),2),
                                                       MORLL =   case_when(exp_output == TRUE ~ round(exp(CI_low),2),
                                                                           exp_output == FALSE ~ round(CI_low,2)),#round(exp(CI_low),3),
                                                       MORUL =   case_when(exp_output == TRUE ~ round(exp(CI_high),2),
                                                                           exp_output == FALSE ~ round(CI_high,2)),#round(exp(CI_high),3),
                                                       `Median OR (CI)` = paste(MedianOR," (",MORLL,", ",MORUL,")",sep=""),
                                                       `CI level` = CI,
                                                       `Median (CI)` = paste(Median," (",CI_low,", ",CI_high,")", sep=""),
                                                       `ROPE CI` = paste(ROPE_CI," (",ROPE_low,", ",ROPE_high,")", sep=""),
                                                       `Rope%` = ROPE_Percentage,
                                                       ESS = round(ESS,1)) %>%
                              select(variable, `Median OR (CI)` , 'Median (CI)' , pd, `ROPE CI`, `Rope%`))
  
  var_names <- all.vars(fit$terms)[-1] # remove outcome Y
  var_labels <- Hmisc::label(model_data)[var_names]
  
  
  # get all levels and combine 
  if (length(fit$xlevels) > 0) {
    
    all_levels0 = fit$xlevels %>% tibble::enframe() %>% tidyr::unnest(.,cols = c("value")) %>%  
      dplyr::mutate(variable = paste0(name, value))
    
    continuousvar<-tibble(name = setdiff(var_names, all_levels0$name),
                          value = rep("",length(setdiff(var_names, all_levels0$name))),
                          variable = setdiff(var_names, all_levels0$name))
    all_levels <- dplyr::bind_rows(all_levels0, continuousvar)
    
    all_levels <- all_levels %>%  dplyr::arrange(match(name,var_names))
    
    neat_fit <- dplyr::full_join(all_levels, neat_fit, by = "variable") %>%
      dplyr::mutate(  name = ifelse(is.na(name), variable, name),
                      value = ifelse(is.na(value), '', value)  )
  } else {
    neat_fit <- neat_fit %>%
      dplyr::mutate(  name = variable,
                      value = ''  )
  }
  
  neat_fit <- neat_fit  %>%
    dplyr::mutate(Variable2 = var_labels[match(name, var_names)] 
                  ,Variable2 = ifelse(is.na(Variable2), name, Variable2))  %>% 
    dplyr::select( -variable)  
  
  # Label the reference groups 
  var_names <- c("(Intercept)", all.vars(fit$terms)) 
  
  neat_fit <- neat_fit %>%
    dplyr::mutate( est.label = ifelse(is.na(`Median OR (CI)`), "1.0 (Reference)", `Median OR (CI)` ),
                   oldname = name) %>%
    dplyr::select( -`Median OR (CI)`, -name)  %>%
    dplyr::select( name = Variable2, Level = value,`Median OR (CI)` = est.label,  everything()) %>%
    dplyr::select( -`Median (CI)`) %>% 
    dplyr::arrange(factor(name, levels = var_names)) %>% 
    dplyr::arrange(match(oldname, var_names))  %>%
    dplyr::select( -oldname) %>% 
    # rename estimate if OR or HR or lm model coefficient = Est eventually.    
    dplyr::rename(!!paste0("Median ", est_name, " (CI)"  ) := `Median OR (CI)`)   
  
  neat_fit$name <- gsub("_"," ",neat_fit$name )
  neat_fit$name[neat_fit$name=="(Intercept)"] <- "Intercept";
  
  return(neat_fit)
}


globalVariables(c("CI", "CI_high", "CI_low", "ESS", "MORLL", "MORUL", "Median","Median (CI)","Median OR (CI)",
    "MedianOR", "ROPE", "CI", "ROPE_CI", "ROPE_Percentage", "ROPE_high", "ROPE_low", "Rope%",
    "Variable2", "bind_rows", "describe_posterior", "oldname", "pd" ))

#' Fancy Table for two group comparison
#' 
#' @param x two level variable
#' @param y variable to compare
#' @param data data.frame or tibble  used to create model fits. Used for capturing variable labels, if they exist
#' 
#' @details 
#' Two group comparison
#' 
#' @return
#' 
#' A tibble with:  \code{Variable}, \code{Parameter}, \code{Median CI},
#'  \code{pd} (probability of direction), 
#'  \code{ROPE CI} (R.O.P.E and CI), 
#'  \code{ROPE\%} (R.O.P.E precentage), 
#'  \code{BF} (bayes factor),  
#'  \code{Prior \[loc scale\]} (Prior distribution with location and scale parameters)
#'  
#' @examples
#' 
#' # Real World Examples
#' data(Bladder_Cancer)
#' library(bayestestR)
#' library(tidyverse)
#' Bladder_Cancer <- Bladder_Cancer %>% 
#' mutate(Cycles_cat = droplevels(Cycles_cat),
#'       Clinical_Stage_Model = recode_factor(Clinical_Stage_Grouped, 
#'                                            'Stage I/II (<=T2NxMx)' = 'Stage I/II (<=T2NxMx)',
#'                                            'Stage III (T3NxMx)' = 'Stage III/IV (T3/4NxMx)'
#'                                            ,'Stage IV (T4NxMx)' = 'Stage III/IV (T3/4NxMx)')
#')
#'
#' # For multiple variables use map_dfr function
#' vars_to_run <- c('PT0N0','Gender', 'Clinical_Stage_Model', 'Cycles_cat')
#' 
#' testmap <- purrr::map_dfr(
#'  vars_to_run, BayesFactorTest, y='Elix_Sum', data = Bladder_Cancer )
#' 
#' # kable(testmap, 'latex', booktabs = TRUE, linesep = '', 
#' # caption = 'Bayesian t-tests for Elix Sum')%>%
#' #  kable_styling(font_size = 7.0) 
#' 
#' # just one 
#' testone <- BayesFactorTest(x = 'Cycles_cat', y = 'Elix_Sum', data = Bladder_Cancer )
#'
#' # kable(testone, 'latex', booktabs = TRUE, linesep = '', 
#' # caption = 'Bayesian t-tests for Elix Sum')%>%
#' # kable_styling(font_size = 7.0) 
#'
#' @importFrom Hmisc label label<-
#' @importFrom utils globalVariables
#' 
#' @export


BayesFactorTest <- function(x, y, data){ 
  .x=x
  fo <- as.formula(  paste(eval(substitute(y), data) ,"~ ",eval(substitute(.x), data),sep=""))
  z <- BayesFactor::ttestBF( formula = fo, data = data)
  result0 <- bayestestR::describe_posterior(z) 
  results <- as.data.frame(result0 %>% dplyr::mutate_if(is.numeric, round, 3)  %>%
                             mutate(Variable =  gsub("_"," ",x ),
                                    `Median (CI)` = paste(Median," (",CI_low,", ",CI_high,")", sep=""),
                                    `ROPE CI` = paste(ROPE_CI," (",ROPE_low,", ",ROPE_high,")", sep=""),
                                    `ROPE%`= ROPE_Percentage,
                                    `Prior [loc, scale]` = paste0(Prior_Distribution," [",Prior_Location,", ",Prior_Scale,"]")) %>%
                             select("Variable","Parameter","Median (CI)","pd", 
                                    "ROPE CI",`ROPE%`,`Prior [loc, scale]`))
  return(results)
}


globalVariables(c("ROPE%", "Prior_Distribution", "Prior_Location",
                  "Prior_Scale", "Prior [loc, scale]"   ))

