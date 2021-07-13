#' Fancy Table Output of Poisson, Logistic, and Cox Models
#' 
#' pretty_model_output() takes a Poisson (glm family = poisson or quasipoisson), Logistic, and Cox model fit object and calculate estimates, odds ratios, or hazard ratios, respectively, with confidence intervals. P values are also produced. For categorical variables with 3+ levels overall Type 3 p values are calculated, in addition to p values comparing to the first level (reference).
#'
#' @param fit  glm, or coxph fit (currently only tested on logistic glm fit)
#' @param model_data data.frame or tibble  used to create model fits. Used for capturing variable labels, if they exist
#' @param title_name title to use (will be repeated in first column)
#' @param conf_level the confidence level required (default is 0.95).
#' @param overall_p_test_stat "Wald" (default) or "LR"; the test.statistic to pass through to the test.statistic param in car::Anova. Ignored for lm fits.
#' @param est_digits number of digits to round OR or HR to (default is 3)
#' @param p_digits number of digits to round p values (default is 4)
#' @param output_type output type, either NULL (default), "latex", or "html" (making special charaters latex friendly)
#' @param sig_alpha the defined significance level for highlighting. Default = 0.05 (Only used if output_type is not NULL)
#' @param background background color of significant values, or no highlighting if NULL. Default is "yellow" (Only used if output_type is not NULL)
#' @param ... other params to pass to \code{pretty_pvalues} (i.e. \code{bold} or \code{italic}) (Only used if output_type is not NULL)
#' 
#' @details 
#' 
#' Model type is determined by \code{fit} class, and also family if glm class. If the class is glm and  binomial or quasibinomial family, then the output is designed for a Logistic model (i.e. Odd Ratios), if the class is coxph the output is designed for a Cox model (i.e. Hazard Ratios), otherwise the output is designed for a Poisson model.
#'
#' @return
#' 
#' A tibble with: \code{Name} (if provided), \code{Variable}, \code{Level}, \code{Est/OR/HR (95\% CI)}, \code{P Value} (for categorical variables comparing to reference), \code{Overall P Value} (for categorical variables with 3+ levels). 
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
#' library(dplyr)
#' # Logistic Regression
#' my_fit <- glm(ybin ~ x1 + x2 + x3, data = my_data, family = binomial(link = "logit"))
#' pretty_model_output(fit = my_fit, model_data = my_data)
#' 
#' # Coxph Regression
#' my_fit <- survival::coxph(survival::Surv(y, ybin) ~ x1 + x2 + x3, data = my_data)
#' my_pretty_model_output <- pretty_model_output(fit = my_fit, model_data = my_data)
#' 
#' # Printing of Fancy table in HTML
#' 
#' kableExtra::kable(my_pretty_model_output, 'html', caption = 'My Table') %>% 
#'    kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack')
#'   
#' # Real World Examples
#' data(Bladder_Cancer)
#' surv_obj <- survival::Surv(Bladder_Cancer$Survival_Months, Bladder_Cancer$Vital_Status == 'Dead')  
#' my_fit <- survival::coxph(surv_obj ~ Gender + Clinical_Stage_Grouped + PT0N0, data = Bladder_Cancer)
#' my_output <- pretty_model_output(fit = my_fit, model_data = Bladder_Cancer)
#' kableExtra::kable(my_output, 'html') %>% 
#'     kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack')
#'   
#' @importFrom Hmisc label label<-
#' @importFrom utils globalVariables
#' 

#' @export
pretty_poisson_model_output <- function(fit, model_data, overall_p_test_stat = c('Wald', 'LR'), title_name = NULL, conf_level = 0.95, est_digits = 3, p_digits = 4, output_type = NULL, sig_alpha = 0.05, background = 'yellow', ...) {

   overall_p_test_stat <- match.arg(overall_p_test_stat)
  .check_numeric_input(est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(p_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(sig_alpha, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  .check_numeric_input(conf_level, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  
  if (any(class(fit) == 'glm') && fit$family$family %in% c('binomial', 'quasibinomial')) {
    # Logistic Regression
    
    #print("LOGISTIC ELSE")
    exp_output <- TRUE
    est_name <- 'OR'
  } else if (any(class(fit) == 'coxph')) {
    # Coxph Regression
    
    #print("COXPH ELSE")
    exp_output <- TRUE
    est_name <- 'HR'
  } else {  
    # print("POISSON ELSE")
    # quasipoisson
    any(class(fit) == 'glm') && fit$family$family %in% c('quasipoisson','poisson')
    exp_output <- TRUE
    est_name <- 'RR'
  }
  
  #Using Variable labels for output, is no label using variable name
  var_names <- all.vars(fit$terms)[-1]
  if (!all(var_names %in% colnames(model_data)))
    stop('All variables used in the "fit" must be in the "model_data" dataset')
  var_labels <- Hmisc::label(model_data)[var_names]
  if (any(var_labels == ''))
    var_labels[var_labels == ''] <- gsub('_', ' ', var_names[var_labels == ''])
  
  neat_fit = fit %>% broom::tidy(conf.int = TRUE, exponentiate = exp_output, conf.level = conf_level)
  
  if (!is.null(output_type)) {
    # P value highlighting if using for pdf output (latex)
    neat_fit$p.label = pretty_pvalues(neat_fit$p.value, digits = p_digits, trailing_zeros = TRUE, sig_alpha = sig_alpha,
                                      output_type = output_type, background = background, ...)
  } else {
    neat_fit$p.label = pretty_pvalues(neat_fit$p.value, digits = p_digits, trailing_zeros = TRUE, output_type = NULL)
  }
  
  neat_fit <- neat_fit %>%
    dplyr::select(variable = term, est = estimate, p.label, p.value, conf.low, conf.high) %>%
    dplyr::filter(variable != "(Intercept)")
  
  
  if (length(fit$xlevels) > 0) {
    all_levels = fit$xlevels %>% tibble::enframe() %>% tidyr::unnest() %>%
      dplyr::mutate(variable = paste0(name, value))
    
    neat_fit <- dplyr::full_join(all_levels, neat_fit, by = "variable") %>%
      dplyr::mutate(
        name = ifelse(is.na(name), variable, name),
        value = ifelse(is.na(value), '', value)
      )
  } else {
    neat_fit <- neat_fit %>%
      dplyr::mutate(
        name = variable,
        value = ''
      )
  }
  
  neat_fit <- neat_fit %>%
    dplyr::mutate(
      est.label = ifelse(is.na(est), "1.0 (Reference)",
                         stat_paste(est, conf.low, conf.high, digits = est_digits, trailing_zeros = TRUE)),
      p.label = ifelse(is.na(p.label), ifelse(!is.null(output_type) && output_type == 'latex', '---', '-'), p.label)
    ) %>%
    dplyr::select(name, Level = value, Est_CI = est.label, `P Value` = p.label) %>%
    dplyr::arrange(factor(name, levels = var_names)) %>%
    dplyr::rename(!!paste0(est_name, paste0(' (', round_away_0(conf_level, 2) * 100, ifelse(!is.null(output_type) && output_type == 'latex', '\\', '')), '% CI)') := Est_CI)
  
  # Dropping extra variable names (for overall p merging)
  neat_fit <- neat_fit %>%
    dplyr::mutate(name_sub = ifelse(duplicated(name), '', name),
                  Variable = var_labels[match(name,var_names)],
                  # Need to add in interaction names
                  Variable = ifelse(is.na(Variable), name, Variable))
  
  ## Type III Overall variable tests
  
  # Getting which vars we need overall tests for
  run_var<-as.vector(NA)
  for(i in 1:length(unique(neat_fit$Variable))){ 
    run_var[i] <- I(length(neat_fit$Variable[neat_fit$Variable == unique(neat_fit$Variable)[i]]) > 2 ) }
  
  overall_vars_needed <-  tibble(name = as.character(unique(neat_fit$Variable)), run_var)
  
  
  
  if (any(overall_vars_needed$run_var)) {
    type3_tests <- dplyr::full_join(broom::tidy(suppressWarnings(car::Anova(fit, type = 'III', test.statistic = overall_p_test_stat))),
                                    overall_vars_needed, by = c('term' = 'name'))
    
    if (!is.null(output_type)) {
      # P value highlighting if using for pdf output (latex)
      type3_tests$overall.p.label = pretty_pvalues(type3_tests$p.value, digits = p_digits,
                                                   trailing_zeros = TRUE, output_type = output_type,
                                                   sig_alpha = sig_alpha, background = background, ...)
    } else {
      type3_tests$overall.p.label = pretty_pvalues(type3_tests$p.value, digits = p_digits, trailing_zeros = TRUE, output_type = NULL)
    }
    
    type3_tests <- type3_tests %>%
      dplyr::filter(term != "(Intercept)" & run_var) %>%
      dplyr::select(variable = term, `Overall P Value` = overall.p.label)
    
    neat_fit <- dplyr::full_join(neat_fit, type3_tests, by = c("name_sub" = "variable")) %>%
      dplyr::mutate(`Overall P Value` = ifelse(is.na(`Overall P Value`), '', `Overall P Value`))
  } else {
    neat_fit <- neat_fit %>% dplyr::mutate(`Overall P Value` = '')
  }
  
  neat_fit <- neat_fit %>% dplyr::select(Variable, Level, dplyr::contains('CI'), dplyr::contains('P Value'))
  
  # Adding Title in front, if given
  if (!is.null(title_name)) dplyr::bind_cols(Name = rep(title_name,nrow(neat_fit)), neat_fit) else neat_fit
}   

globalVariables(c( "Overall P Value"))

#' Wrapper for Pretty Model Output
#' 
#' Wrapper for pretty_model_output(). This function takes a dataset, along with variables names for x (could be multiple), y, and possibly event status, for model fit.
#'
#' @param x_in name of x variables in model (can be vector of x names)
#' @param model_data data.frame or tibble that contains \code{x_in}, \code{time_in}, and \code{event_in} variables
#' @param y_in name of outcome measure for logistic and linear model, or name of time component in cox model
#' @param event_in name of event status variable. Shouled be left NULL for logistic and linear models. If \code{event_level} = NULL then this must be the name of a F/T or 0/1 variable, where F or 0 are considered the censored level, respectively.
#' @param event_level outcome variable event level for logistic model, and event status level for cox model.
#' @param title_name title to use (will be repeated in first column)
#' @param fail_if_warning Should program stop and give useful message if there is a warning message when running model (Default is TRUE)
#' @param conf_level the confidence level required (default is 0.95).
#' @param overall_p_test_stat "Wald" (default) or "LR"; the test.statistic to pass through to the test.statistic param in car::Anova. Ignored for lm fits.
#' @param est_digits number of digits to round OR or HR to (default is 3)
#' @param p_digits number of digits to round p values (default is 4)
#' @param output_type output type, either NULL (default), "latex", or "html" (making special charaters latex friendly)
#' @param sig_alpha the defined significance level for highlighting. Default = 0.05 (Only used if output_type not NULL)
#' @param background background color of significant values, or no highlighting if NULL. Default is "yellow" (Only used if output_type not NULL)
#' @param verbose a logical variable indicating if warnings and messages should be displayed. Default FALSE.
#' @param ... other params to pass to \code{pretty_pvalues} (i.e. \code{bold} or \code{italic})
#
#' 
#' @details 
#' \code{x_in} can be single variable name, or vector of variables to include in the model. All variables must be present in the \code{model_data} dataset.
#' 
#' \code{fail_if_warning} variable default to TRUE because most warnings should be addressed, such as the "Loglik converged before variable XX; beta may be infinite" warning.
#' 
#' @return
#' 
#' A tibble with: \code{Name} (if provided), \code{Variable}, \code{Level}, \code{Est/OR/HR (95\% CI)}, \code{P Value} (for categorical variables comparing to reference), \code{Overall P Value} (for categorical variables with 3+ levels), \code{n/n (event)}. 
#' 
#' @examples
#' 
#' # Basic linear model example
#' set.seed(542542522)
#' ybin <- sample(0:1, 100, replace = TRUE)
#' ybin2 <- sample(c('Male','Female'), 100, replace = TRUE)
#' ybin3 <- sample(c('Dead','Alive'), 100, replace = TRUE)
#' y <- rexp(100,.1)
#' x1 <- factor(sample(LETTERS[1:2],100,replace = TRUE))
#' x2 <- factor(sample(letters[1:4],100,replace = TRUE))
#' my_data <- data.frame(y, ybin, ybin2, ybin3, x1, x2)
#' Hmisc::label(my_data$x1) <- "X1 Variable"
#' library(dplyr)
#'  # Single runs 
#' run_pretty_model_output(x_in = 'x1', model_data = my_data, y_in = 'y', event_in = 'ybin')
#' run_pretty_model_output(x_in = 'x1', model_data = my_data, y_in = 'y', 
#'      event_in = 'ybin3', event_level = 'Dead')
#' run_pretty_model_output(x_in = c('x1','x2'), model_data = my_data, y_in = 'y', event_in = 'ybin')
#' run_pretty_model_output(x_in = 'x2', model_data = my_data, y_in = 'ybin',
#'  event_in = NULL, verbose = TRUE)
#' run_pretty_model_output(x_in = 'x2', model_data = my_data, y_in = 'y', event_in = NULL)
#' 
#' # Multiple runs for different variables
#'  
#' vars_to_run = c('x1', 'x2')
#' cox_models <- purrr::map_dfr(vars_to_run, run_pretty_model_output, model_data = my_data, 
#'      y_in = 'y', event_in = 'ybin')
#' 
#' kableExtra::kable(cox_models, 'html', caption = 'My Table') %>% 
#'   kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack', headers_to_remove = 1:2)
#' 
#' # Real World Example
#' data(Bladder_Cancer)
#' vars_to_run = c('Gender', 'Clinical_Stage_Grouped', 'PT0N0', 'Any_Downstaging')
#' 
#' univariate_output <- purrr::map_dfr(vars_to_run, run_pretty_model_output,
#'  model_data = Bladder_Cancer, 
#'       y_in = 'Survival_Months', event_in = 'Vital_Status', event_level = 'Dead')
#' kableExtra::kable(univariate_output, 'html') %>% 
#'       kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack', headers_to_remove = 1:2)
#' 
#' multivariable_output <- run_pretty_model_output(vars_to_run, model_data = Bladder_Cancer, 
#'       y_in = 'Survival_Months', event_in = 'Vital_Status', event_level = 'Dead')
#' kableExtra::kable(multivariable_output, 'html') %>% 
#'       kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack', headers_to_remove = 1:2)
#' 
#' @importFrom  Hmisc label
#' 
#' @export
#' 


run_pretty_poisson_model_output <- function(x_in, model_data, y_in, event_in = NULL, event_level = NULL, title_name = NULL, fail_if_warning = TRUE, conf_level = 0.95, overall_p_test_stat = c('Wald', 'LR'), est_digits = 3, p_digits = 4, output_type = NULL, sig_alpha = 0.05, background = 'yellow', verbose = FALSE, ...) {
  overall_p_test_stat <- match.arg(overall_p_test_stat)
  .check_numeric_input(est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(p_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(sig_alpha, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  .check_numeric_input(conf_level, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  if (!is.null(output_type) && !output_type %in% c('latex','html'))
    stop('"output_type" must be either NULL, "latex", or "html"')
  if (!all(x_in %in% colnames(model_data)))
    stop('All "x_in" (',paste0(x_in, collapse = ', '), ') must be in the "model_data" dataset')
  if (length(y_in) != 1) stop('"y_in" must be length of 1')
  if (all(y_in != colnames(model_data)))
    stop('"y_in" (',y_in, ') not in the "model_data" dataset')
  if (length(unique(na.omit(model_data[,y_in, drop = TRUE]))) <= 1)
    stop('"y_in" (',y_in, ') must have more than one unique value')
  
  x_in_paste = paste0(x_in, collapse = ' + ')
  
  if (is.null(event_in)) {
    tmp_formula <- as.formula(paste(y_in, " ~ ", x_in_paste))
    if (length(unique(na.omit(model_data[,y_in, drop = TRUE]))) == 2) {
      # Logistic Model
      # making y_in a factor
      if (!is.null(event_level)) {
        if (all(unique(model_data[, y_in, drop = TRUE]) != event_level))
          stop('"event_level" (',event_level, ') not present in "y_in" (',y_in, ')')
        model_data[, y_in] <- factor(model_data[, y_in, drop = TRUE] == event_level)
      } else {
        model_data[, y_in] <- factor(model_data[, y_in, drop = TRUE])
        if (verbose)
          message('Since no "event_level" specified setting "',levels(model_data[, y_in, drop = TRUE])[2], '" as outcome event level in logistic model')
      }
      if (nlevels(model_data[, y_in, drop = TRUE]) != 2)
        stop('"y_in" (',y_in, ') must have two levels for logistic model')
      tmp_fit <- tryCatch(expr =  glm(tmp_formula, data = model_data, family = binomial(link = "logit")),
                          error = function(c) stop('Logistic model with "',deparse(tmp_formula), '" formula has error(s) running'))
      if (fail_if_warning) {
        tmp_confint <- tryCatch(expr = suppressMessages(confint(tmp_fit)),
                                error = function(c) stop('Logistic model with "',deparse(tmp_formula), '" formula has error(s) calculating CI(s)'),
                                warning = function(c) stop('Logistic model with "',deparse(tmp_formula), '" formula has Inf CI(s); most likely a model error, most likely due to sparse counts or perfect seperation'))
      }
    } else {
      # Poisson Model
      tmp_fit <- tryCatch(expr =   glm(tmp_formula, data = model_data, family = "quasipoisson"),#lm(tmp_formula, data = model_data),
                          error = function(c) stop('(quasi)poisson model with "',deparse(tmp_formula), '" formula has error(s) calculating CI(s)'))
      if (fail_if_warning) {
        tmp_confint <- tryCatch(expr = suppressMessages(confint(tmp_fit)),
                                error = function(c) stop('(quasi)poisson  model with "',deparse(tmp_formula), '" formula has error(s) calculating CI(s)'),
                                warning = function(c) stop('(quasi)poisson  model with "',deparse(tmp_formula), '" formula has Inf CI(s); most likely a model error'))
      }
    } #end Poisson models...
    n_info <- paste0('n=',nrow(tmp_fit$model))
  } else {
    # Cox Model
    if (all(event_in != colnames(model_data)))
      stop('"event_in" (',event_in, ') not in the "model_data" dataset')
    if (length(unique(na.omit(model_data[,event_in, drop = TRUE]))) > 2)
      stop('"event_in" (',event_in, ') must have only two levels')
    
    if (!is.null(event_level)) {
      if (all(unique(model_data[, event_in, drop = TRUE]) != event_level))
        stop('"event_level" (',event_level, ') not present in "event_in" (',event_in, ')')
      model_data[,event_in] <- model_data[,event_in, drop = TRUE] == event_level
    }
    event_levels <- unique(model_data[,event_in, drop = TRUE])
    if (all(event_levels != TRUE))
      stop('"event_in" (',event_in, ') must have at least one event')
    
    
    tmp_formula <- as.formula(paste("survival::Surv(",y_in,",",event_in,") ~ ", x_in_paste))
    
    if (fail_if_warning) {
      tmp_fit <- tryCatch(expr = survival::coxph(tmp_formula, data = model_data),
                          error = function(c) stop('Cox model with "',deparse(tmp_formula), '" formula has error(s) running'),
                          # Need special error checking because coxph throws some warnings when warnPartialMatchArgs = TRUE
                          warning = function(c) if (any(grepl('converge', c)))
                            stop('Cox model with "',deparse(tmp_formula), '" formula has warnings(s) running')
                          else suppressWarnings(survival::coxph(tmp_formula, data = model_data)))
    } else {
      tmp_fit <- tryCatch(expr = survival::coxph(tmp_formula, data = model_data),
                          error = function(c) stop('Cox model with "',deparse(tmp_formula), '" formula has error(s) running'))
      
    }
    n_info <- paste0('n=',tmp_fit$n,' (',tmp_fit$nevent,')')
  }
  
  tmp_output <- pretty_poisson_model_output(fit = tmp_fit, model_data = model_data, title_name = title_name, conf_level = conf_level, overall_p_test_stat = overall_p_test_stat, est_digits = est_digits, p_digits = p_digits,  output_type =  output_type, sig_alpha = sig_alpha, background = background, ...)
  tmp_output <- dplyr::bind_cols(tmp_output, n =  c(n_info, rep("", nrow(tmp_output) - 1)))
  
  if (!is.null(event_in)) names(tmp_output)[names(tmp_output) == 'n'] <- 'n (events)'
  
  tmp_output
}


