#' Pasting Together Information for Two Groups
#'
#' Paste together information, often statistics, from two groups. There are two predefined combinations: mean(sd) and median[min,max], but user may also paste any single measure together.
#'
#'
#' @param data input dataset. User must use consistent naming throughout, \strong{with an underscore} to separate the group names from the measures (i.e. \code{Group1_mean} and \code{Group2_mean}). There also must be two columns with column names that exactly match the input for \code{first_name} and \code{second_name} (i.e. 'Group1' and 'Group2'), which are used to form the \code{Comparison} variable.
#' @param vars_to_paste vector of names of common measures to paste together. Can be the predefined 'median_min_max' or 'mean_sd', or any variable as long as they have matching columns for each group (i.e. Group1_MyMeasure and Group2_MyMeasure). Multiple measures can be requested. Default: "all" will run 'median_min_max' and 'mean_sd', as well as any pairs of columns in the proper format.
#' @param first_name name of first group (string before '_') . Default is 'Group1'.
#' @param second_name name of second group (string before '_'). Default is 'Group2'.
#' @param sep_val value to be pasted between the two measures. Default is ' vs. '.
#' @param na_str_out the character to replace missing values with.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". Will be used to determine the character to be pasted between the group names (\code{Comparison} variable).  Specifying "two.sided" will use the \code{sep_val} input.
#' @param digits integer indicating the number of decimal places to round to before pasting for numeric variables. Default is 0.
#' @param trailing_zeros logical indicating if trailing zeros should be included (i.e. 0.100 instead of 0.1). Note if set to TRUE output is a character vector.
#' @param keep_all logical indicating if all remaining, unpasted variables in \code{data} should be returned with the pasted variables. Default TRUE.
#' @param verbose a logical variable indicating if warnings and messages should be displayed. Default FALSE.
#' @details
#'
#' User must use consistant naming throughout, with a underscore to separate the group names from the measures (i.e. \code{Group1_mean} and \code{Group2_mean}). There also must be columns defining the group names (i.e. \code{Group1} and \code{Group2}), which are used to form the \code{Comparison} variable.
#'
#' \code{alternative} included as a parameter so the direction can easily be seen in one-sided test. If "two.sided" is selected the value to be pasted between the two group names will be set to \code{sep_val}, where "greater" will use " > " and "less" with use " < " as the pasting value.
#'
#'
#' @return data.frame with all the pasted values requested. Each name will have '_comparison' at the end of the names (i.e. mean_comparison, median_comparison, ...)
#' @examples
#'
#' # Same examples on data.table
#' library(data.table)
#' data(exampleData_BAMA)
#'
#' descriptive_stats_by_group <- exampleData_BAMA[, .(
#'      Group1 = unique(group[group == 1]), Group2 = unique(group[group == 2]),
#'      Group1_n = length(magnitude[group == 1]), Group2_n = length(magnitude[group == 2]),
#'      Group1_mean = mean(magnitude[group == 1]), Group2_mean = mean(magnitude[group == 2]),
#'      Group1_sd = sd(magnitude[group == 1]), Group2_sd = sd(magnitude[group == 2]),
#'      Group1_median = median(magnitude[group == 1]), Group2_median = median(magnitude[group == 2]),
#'      Group1_min = min(magnitude[group == 1]), Group2_min = min(magnitude[group == 2]),
#'      Group1_max = max(magnitude[group == 1]), Group2_max = max(magnitude[group == 2])
#' ), by = .(visitno,antigen)]
#'
#' paste_tbl_grp(data = descriptive_stats_by_group, vars_to_paste = 'all', 
#'    first_name = 'Group1', second_name = 'Group2', 
#'    sep_val = " vs. ", digits = 0, keep_all = TRUE)
#'
#' paste_tbl_grp(data = descriptive_stats_by_group, vars_to_paste = c("mean", "median_min_max"), 
#'    alternative= "less", keep_all = FALSE)
#'
#' paste_tbl_grp(data = descriptive_stats_by_group, vars_to_paste = 'all', 
#'    first_name = 'Group1', second_name = 'Group2', sep_val = " vs. ", 
#'    alternative = 'less', digits = 5, keep_all = FALSE)
#'
#'
#' # Same example with tidyverse (dplyr+tidyr) with some custom functions
#'
#' library(dplyr)
#' library(tidyr)
#'
#'q95_fun = function(x) quantile(x, 0.95)
#'N = function(x) length(x)
#'
#'exampleData_BAMA %>%
#'  mutate(group = paste0("Group", group)) %>%
#'  group_by(group, visitno, antigen) %>%
#'  summarise_at("magnitude", funs(N, mean, sd, median, min, max, q95_fun)) %>%
#'  gather(variable, value, -(group:antigen)) %>% # these three chains create a wide dataset
#'  unite(temp, group, variable) %>%
#'  spread(temp, value) %>%
#'  mutate(Group1 = "Group 1", Group2 = "Group 2") %>%
#'  paste_tbl_grp()
#'
#' @import data.table
#' @export


paste_tbl_grp <- function(data, vars_to_paste = 'all', first_name = 'Group1', second_name = 'Group2', sep_val = " vs. ", na_str_out = "---", alternative = c("two.sided", "less", "greater"), digits = 0, trailing_zeros = TRUE, keep_all = TRUE, verbose = FALSE){
  
  #####Checking variables being used
  
  data_here <- as.data.frame(data)
  alternative <- match.arg(alternative)
  .check_numeric_input(digits, lower_bound = 0, scalar = TRUE)
  if (first_name == second_name) stop('"first_name" and "second_name" must be different')
  if (sum(first_name == names(data_here)) != 1) stop('Expecting one column named "', first_name , '" in input dataset, but there are ', sum(first_name == names(data_here)), ' present')
  if (sum(second_name == names(data_here)) != 1) stop('Expecting one column named "', second_name , '" in input dataset, but there are ', sum(second_name == names(data_here)), ' present')
  if (length(vars_to_paste) != 1 & any(vars_to_paste == 'all')) {
    vars_to_paste = 'all'
    if (verbose) message('Since "all" specified other entries of vars_to_paste ignored.')
  }
  
  # Defining vars when vars_to_paste set to 'all'
  if (length(vars_to_paste) == 1 & any(vars_to_paste == 'all')) {
    #Need to address if one group name a subset of another
    temp_group1_names <- names(data_here)[(substr(names(data_here), 0, nchar(first_name)) == first_name) &
                                            (nchar(first_name) > nchar(second_name) | substr(names(data_here), 0, nchar(second_name)) != second_name)]
    temp_group2_names <- names(data_here)[(substr(names(data_here), 0, nchar(second_name)) == second_name) &
                                            (nchar(second_name) > nchar(first_name) | substr(names(data_here), 0, nchar(first_name)) != first_name)]
    temp_group1_measures <- gsub(paste0(first_name,'_'), '', temp_group1_names, fixed = TRUE)
    temp_group2_measures <- gsub(paste0(second_name,'_'), '', temp_group2_names, fixed = TRUE)
    vars_to_paste_here <- unique(intersect(temp_group1_measures, temp_group2_measures))
    
    # Adding speacial cases
    if (sum(vars_to_paste_here %in% c('median','min','max')) == 3) vars_to_paste_here <- c(vars_to_paste_here, 'median_min_max')
    if (sum(vars_to_paste_here %in% c('mean','sd')) == 2) vars_to_paste_here <- c(vars_to_paste_here, 'mean_sd')
    
    if (verbose) message('The following measures will be combined: ', paste0(vars_to_paste_here, collapse = ', '))
  } else {
    vars_to_paste_here <- unique(vars_to_paste)
  }
  # Giving a message if nothing to return
  if (any(vars_to_paste == 'all') & length(vars_to_paste_here) == 0) {
    if (verbose) message('"all" specified, but no matching columns to paste')
    if (keep_all) return(data) else return(NULL)
  }
  
  # Need to define which variables to check. Special considerations for the predefined values
  vars_to_check <- vars_to_paste_here[!vars_to_paste_here %in% c('median_min_max','mean_sd')]
  if (any(vars_to_paste_here == 'median_min_max')) vars_to_check <- unique(c(vars_to_check, 'median', 'min', 'max'))
  if (any(vars_to_paste_here == 'mean_sd')) vars_to_check <- unique(c(vars_to_check, 'mean', 'sd'))
  
  # Need to check the group1 and group2 version of each variable being pasted
  group1_vars_to_check <- paste0(first_name, '_', vars_to_check)
  group2_vars_to_check <- paste0(second_name, '_', vars_to_check)
  for (i in 1:length(vars_to_check)) {
    if (sum(group1_vars_to_check[i] == names(data_here)) != 1) stop('Expecting one column named "', group1_vars_to_check[i] , '" in input dataset, but there are ', sum(group1_vars_to_check[i] == names(data_here)), ' present')
    temp_group1_var <- data_here[, group1_vars_to_check[i] == names(data_here)]
    if (is.numeric(temp_group1_var) & any(!is.na(temp_group1_var))) .check_numeric_input(temp_group1_var)
    
    if (sum(group2_vars_to_check[i] == names(data_here)) != 1) stop('Expecting one column named "', group2_vars_to_check[i] , '" in input dataset, but there are ', sum(group2_vars_to_check[i] == names(data_here)), ' present')
    temp_group2_var <- data_here[, group2_vars_to_check[i] == names(data_here)]
    if (is.numeric(temp_group2_var) & any(!is.na(temp_group2_var))) .check_numeric_input(temp_group2_var)
  }
  
  ##### Pasting variables
  
  # Comparison variable
  comparison_sep <- switch(alternative,
                           two.sided = sep_val,
                           less = ' < ',
                           greater = ' > ')
  
  comparison_var <- paste0(data_here[, first_name == names(data_here)], comparison_sep, data_here[, second_name == names(data_here)])
  
  # Other variables
  pasted_results <- list()
  for (i in 1:length(vars_to_paste_here)) {
    if (vars_to_paste_here[i] == 'median_min_max') {
      pasted_results[[i]] <- paste0(
        stat_paste(stat1 = data_here[, paste0(first_name, '_median')],
                   stat2 = data_here[, paste0(first_name, '_min')],
                   stat3 = data_here[, paste0(first_name, '_max')],
                   digits = digits, bound_char = '[', sep = ', ', na_str_out = na_str_out, trailing_zeros = trailing_zeros
        ),
        sep_val,
        stat_paste(stat1 = data_here[, paste0(second_name, '_median')],
                   stat2 = data_here[, paste0(second_name, '_min')],
                   stat3 = data_here[, paste0(second_name, '_max')],
                   digits = digits, bound_char = '[', sep = ', ', na_str_out = na_str_out, trailing_zeros = trailing_zeros
        )
      )
    } else if (vars_to_paste_here[i] == 'mean_sd') {
      pasted_results[[i]] <-  paste0(
        stat_paste(stat1 = data_here[, paste0(first_name, '_mean')],
                   stat2 = data_here[, paste0(first_name, '_sd')],
                   digits = digits, bound_char = '(', na_str_out = na_str_out, trailing_zeros = trailing_zeros
        ),
        sep_val,
        stat_paste(stat1 = data_here[, paste0(second_name, '_mean')],
                   stat2 = data_here[, paste0(second_name, '_sd')],
                   digits = digits, bound_char = '(', na_str_out = na_str_out, trailing_zeros = trailing_zeros
        )
      )
    } else {
      first_var_here <- data_here[, paste0(first_name, '_', vars_to_paste_here[i])]
      second_var_here <- data_here[, paste0(second_name, '_', vars_to_paste_here[i])]
      both_var_here <- c(first_var_here, second_var_here)
      
      # Want to set digits to 0 if an integer
      if (is.numeric(both_var_here)) {
        if (any((both_var_here %% 1) != 0)) digits_here = digits else digits_here = 0
      } else digits_here = digits
      
      pasted_results[[i]] <- paste0(ifelse(is.na(first_var_here), na_str_out, .round_if_numeric(first_var_here, digits = digits_here, trailing_zeros = trailing_zeros)),
                                    sep_val,
                                    ifelse(is.na(second_var_here), na_str_out, .round_if_numeric(second_var_here, digits = digits_here, trailing_zeros = trailing_zeros))
      )
    }
  }
  names(pasted_results) <- paste0(vars_to_paste_here, '_comparison')
  
  pasted_results <- data.frame('Comparison' = comparison_var,pasted_results, stringsAsFactors = FALSE)
  
  # Returning all data if desired
  if (keep_all) {
    index_to_keep <- !names(data_here) %in% c(first_name, second_name, group1_vars_to_check, group2_vars_to_check)
    data.frame(data_here[, index_to_keep, drop = FALSE], pasted_results, stringsAsFactors = FALSE)
  } else {
    pasted_results
  }
  
}




#' Rounds and combines up to three numbers into table friendly presentation
#'
#' Takes in up to 3 numeric values, rounds each to a specified digit amount (if numeric), and then combines them accordingly.
#'
#' @param stat1 first statistic to be pasted.
#' @param stat2 second statistic to be pasted (optional).
#' @param stat3 third statistic to be pasted (optional).
#' @param digits positive integer of length 1 between 0 (default) and 14, giving the amount of digits to round to.
#' @param trailing_zeros logical indicating if trailing zeros should included (i.e. 0.100 instead of 0.1). Note is set to TRUE output is a character vector
#' @param bound_char the character to be used between stat1 and stat2/stat3. Available options are '(' (default), '[', '\{', and '|'.
#' @param sep the string to be used between stat2 and stat3. The default is ', '.
#' @param na_str_out the character to replace missing values with.
#' @return string of combined values
#'
#' @details
#'
#' One value provided - returns a rounded value or the missing character.
#' Two values - returns stat1 (stat2). e.g., mean (sd)
#' Three values - returns stat1 (stat2, stat3). e.g., estimate (95\% lower, 95\% upper) or median [min, max]
#'
#' Currently the method does work with string variables, but of course rounding would not be relevant for strings.
#'
#'
#' @examples
#'
#' stat_paste(5.109293)
#' stat_paste(NA)
#' stat_paste(5.109293, 2.145)
#' stat_paste(5.109293, 1.7645, 8.0345)
#' stat_paste(NA, NA, NA)
#' stat_paste(5.109, "p < 0.001", digits = 3)
#' stat_paste(c(rep(5,5),NA),c(1:5,NA),c(1,NA,2,NA,3,NA),bound_char = '[')
#'
#' library(data.table)
#' data(testData_BAMA)
#' testData_BAMA [!is.na(magnitude), .(
#'   median_min_max = stat_paste(
#'      median(magnitude), min(magnitude), max(magnitude)
#'      )), by = .(antigen, visit, group)]
#'
#' @export
stat_paste = function(stat1, stat2 = NULL, stat3 = NULL, digits = 0, trailing_zeros = TRUE, bound_char = c('(','[','{','|'), sep = ', ', na_str_out = "---"){
  bound_char <- match.arg(bound_char)
  end_bound_char <-   switch(bound_char,
                             `(` = ')',
                             `[` = ']',
                             `{` = '}',
                             `|` = '|'
  )
  
  stat1_pasted_obj <-  ifelse(is.na(stat1), na_str_out, as.character(.round_if_numeric(stat1, digits = digits, trailing_zeros = trailing_zeros)))
  if (is.null(stat2)) {
    pasted_output <- stat1_pasted_obj
  } else {
    stat2_pasted_obj <-  ifelse(is.na(stat2), na_str_out, as.character(.round_if_numeric(stat2, digits = digits, trailing_zeros = trailing_zeros)))
    if (is.null(stat3)) {
      pasted_output <- ifelse(is.na(stat1) & is.na(stat2), na_str_out, paste0(stat1_pasted_obj, " ", bound_char, stat2_pasted_obj, end_bound_char))
    } else {
      stat3_pasted_obj <-  ifelse(is.na(stat3), na_str_out, as.character(.round_if_numeric(stat3, digits = digits, trailing_zeros = trailing_zeros)))
      pasted_output <- ifelse(is.na(stat1) & is.na(stat2) & is.na(stat3), na_str_out, paste0(stat1_pasted_obj, " ", bound_char, stat2_pasted_obj, sep, stat3_pasted_obj, end_bound_char))
    }
  }
  pasted_output
}


#' Round and format a vector of p-values
#'
#' pretty_pvalues() takes a vector of p-values, rounds them to a specified digit amount,
#' allows options for emphasizing p-values < the defined significance level, and returns a character for missing.
#'
#' @param pvalues numeric vector of raw p-values to be formatted
#' @param digits number of digits to round to; values with zeros past this number of digits are truncated
#' @param missing_char character string that will replace missing values from the p-value vector. Default = "---"
#' @param include_p TRUE or FALSE: set to TRUE to print "p = " before each p-value
#' @param trailing_zeros TRUE or FALSE: default = TRUE, p-values are formatted with trailing zeros to the defined number of digits (i.e. 0.100 instead of 0.1 if digits= 3)
#' @param output_type output type, either NULL (default), "latex", or "html" (making special charaters latex friendly)
#' @param bold TRUE or FALSE: set to TRUE to bold p-values < the defined significance level
#' @param italic TRUE or FALSE: set to TRUE to italicize p-values < the defined significance level
#' @param background highlight color for p-values < the defined significance level. Default = NULL (no highlighting)
#' @param sig_alpha the defined significance level. Default = 0.05
#'
#' @return Vector of transformed p-values for table output
#'
#' @details
#'
#' With this function, there are two things to be noted:
#' Since the p-value vector formatting uses \code{cell_spec}, which generates raw HTML or LaTeX code, make sure you remember to put \code{escape = FALSE} into your kable code when generating your table. At the same time, you will need to escape special symbols manually.
#' Additionally, \code{cell_spec} needs a way to know whether you want HTML or LaTeX output. You can specify it locally in the function or globally using \code{options(knitr.table.format = "latex")}. If you don't provide anything, this function will output as HTML by default.
#'
#' @examples
#'
#' pvalue_example = c(1, 0.06, 0.0005, NA, 1e-6)
#'
#' pretty_pvalues(pvalue_example, background = "pink")
#'
#' pretty_pvalues(pvalue_example, digits = 4, missing_char = "missing", bold = TRUE)
#'
#' # How to use pretty_pvalues in reports
#' raw_pvals <- c(0.00000007, .05000001, NaN, NA, 0.783)
#' pretty_pvals <- pretty_pvalues(raw_pvals , digits = 3, 
#'    background = "green", italic = TRUE, bold = TRUE)
#' kableExtra::kable(pretty_pvals , format = "latex", escape = FALSE, col.names = c("P-values"))
#' @export
pretty_pvalues = function(pvalues, digits = 3, missing_char = '---', include_p = FALSE, trailing_zeros = TRUE, output_type = NULL, bold = FALSE, italic = FALSE, background = NULL, sig_alpha = 0.05){
  
  .check_numeric_input(pvalues, lower_bound = 0, upper_bound = 1)
  .check_numeric_input(sig_alpha, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  .check_numeric_input(digits, lower_bound = 1, upper_bound = 14, scalar = TRUE, whole_num = TRUE)
  if (!is.null(output_type) && !output_type %in% c('latex','html'))
    stop('"output_type" must be either NULL, "latex", or "html"')
  # Using Latex if user didn't specify 'output_type'
  if (is.null(output_type) & (bold == TRUE | italic == TRUE | !is.null(background)))
    output_type = 'latex'

  #Need to set options for no scientific notation, but set back to user preference on exit
  op <- options()
  options(scipen = 10)
  on.exit(options(op))
  
  lower_cutoff = 10^(-digits)
  
  ## relevant p-value indices for specific assignments
  missing_p = which(is.na(pvalues))
  below_cutoff_p = which(pvalues < lower_cutoff)
  sig_p = which(pvalues < sig_alpha)
  
  if (trailing_zeros) pvalues_new = round_away_0(pvalues, trailing_zeros = TRUE, digits = digits) else pvalues_new = as.character(round_away_0(pvalues, trailing_zeros = FALSE, digits))
  
  ## manipulate and assign pvalues as characters to output pvalue vector
  pvalues_new[missing_p] = missing_char
  pvalues_new[below_cutoff_p] = paste0("<", lower_cutoff)
  
  # the letter 'p' in front of values
  if (include_p) pvalues_new <- ifelse(pvalues_new < lower_cutoff, paste0('p',  pvalues_new), paste0('p=',  pvalues_new))
  
  # formatting
  if ((bold == TRUE | italic == TRUE | !is.null(background))) pvalues_new[sig_p] = kableExtra::cell_spec(pvalues_new[sig_p], format = output_type, bold = bold, italic = italic, background = background, escape = FALSE)
  
  pvalues_new
}


#' pretty_pvalues() wrapper for compareGroups table
#'
#' pretty_pvalues_compareGroups() takes a createTable object created from compareGroups::createTable(), and applies pretty_pvalues() to all p value columns.
#'
#' @param compareGroups_table_obj createTable object created from compareGroups::createTable() function, with at least one p value column
#' @param digits number of digits to round to; values with zeros past this number of digits are truncated
#' @param include_p TRUE or FALSE: set to TRUE to print "p = " before each p-value
#' @param trailing_zeros TRUE or FALSE: default = TRUE, p-values are formatted with trailing zeros to the defined number of digits (i.e. 0.100 instead of 0.1 if digits= 3)
#' @param output_type output type, either NULL (default), "latex", or "html" (making special charaters latex friendly)
#' @param bold TRUE or FALSE: set to TRUE to bold p-values < the defined significance level
#' @param italic TRUE or FALSE: set to TRUE to italicize p-values < the defined significance level
#' @param background highlight color for p-values < the defined significance level. Default = NULL (no highlighting)
#' @param sig_alpha the defined significance level. Default = 0.05
#'
#' @return createTable object with rounding, formating, and highlighting of any p value column.
#'
#' @details
#'
#' This function is design specifically for input from compareGroups tables. Use \code{pretty_pvalues} directly when working with normal data.frames and other objects.
#'
#' @examples
#'
#' data("Bladder_Cancer") 
#' library(compareGroups)
#' cycles_formula <- as.formula(Cycles_cat ~ Age_At_Diagnosis + Gender + Elix_Sum)
#' cycles_compare <- compareGroups::compareGroups(cycles_formula, data = Bladder_Cancer)
#' cycles_table <- compareGroups::createTable(cycles_compare, digits.p = 4, show.p.mul = TRUE)
#' cycles_table_fancy <- pretty_pvalues_compareGroups(cycles_table, background = 'yellow')
#' 
#' # Printing createTable object in report
#' compareGroups::export2latex(cycles_table_fancy, size = 'footnotesize'
#' , label = 'tab:Variable-of-Interest-Table'
#' , caption = 'Comparing Variables to Number of Cycles'
#' , header.labels = c('p.overall' = 'Overall P'), landscape = FALSE)
#' 
#' @export


pretty_pvalues_compareGroups <- function(compareGroups_table_obj, digits = 3, include_p = FALSE, trailing_zeros = TRUE, output_type = NULL, bold = FALSE, italic = FALSE, background = NULL, sig_alpha = 0.05){
  if (class(compareGroups_table_obj) != 'createTable') stop("'compareGroups_table_obj' must be a createTable class (output of compareGroups::createTable() function)")
  
  vars_to_adjust <- colnames(compareGroups_table_obj$descr)[grep('p\\.', colnames(compareGroups_table_obj$descr))]
  
  for (i in vars_to_adjust) {
    temp_values <- compareGroups_table_obj$descr[,i]
    temp_values[temp_values == '.'] <- ''
    numeric_values <- as.numeric(sub('<', '', temp_values))
    compareGroups_table_obj$descr[,i] <- pretty_pvalues(numeric_values, digits = digits, missing_char = '', include_p = include_p, trailing_zeros = trailing_zeros, output_type = output_type, bold = bold, italic = italic, background = background, sig_alpha = sig_alpha)
  }
  compareGroups_table_obj
}


#' Fancy Table Output of Linear, Logistic, and Cox Models
#' 
#' pretty_model_output() takes a Linear, Logistic, and Cox model fit object and calculate estimates, odds ratios, or hazard ratios, respectively, with confidence intervals. P values are also produced. For categorical variables with 3+ levels overall Type 3 p values are calculated, in addition to p values comparing to the first level (reference).
#'
#' @param fit lm, glm, or coxph fit (currently only tested on logistic glm fit)
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
#' Model type is determined by \code{fit} class, and also family if glm class. If the class is glm and  binomial or quasibinomial family, then the output is designed for a Logistic model (i.e. Odd Ratios), if the class is coxph the output is designed for a Cox model (i.e. Harzard Ratios), otherwise the output is designed for a linear model or other model where normal coefficient estimates are displayed.
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
#' # Linear Regression
#' my_fit <- lm(y ~ x1 + x2 + x3, data = my_data)
#' pretty_model_output(fit = my_fit, model_data = my_data)
#' 
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
#' @importFrom stats as.formula binomial chisq.test confint
#' cor.test fisher.test glm lm mcnemar.test
#' na.omit pchisq t.test
#' 
#' @export


pretty_model_output <- function(fit, model_data, overall_p_test_stat = c('Wald', 'LR'), title_name = NULL, conf_level = 0.95, est_digits = 3, p_digits = 4, output_type = NULL, sig_alpha = 0.05, background = 'yellow', ...) {
  overall_p_test_stat <- match.arg(overall_p_test_stat)
  .check_numeric_input(est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(p_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(sig_alpha, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  .check_numeric_input(conf_level, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  
  if (any(class(fit) == 'glm') && fit$family$family %in% c('binomial', 'quasibinomial')) {
    # Logistic Regression
    exp_output <- TRUE
    est_name <- 'OR'
  } else if (any(class(fit) == 'coxph')) {
    # Coxph Regression
    exp_output <- TRUE
    est_name <- 'HR'
  } else {
    # Not Logistic Regression or Coxph
    exp_output <- FALSE
    est_name <- 'Est'
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
   # all_levels = fit$xlevels %>% tibble::enframe() %>% tidyr::unnest() %>% 
     all_levels = fit$xlevels %>% tibble::enframe() %>% tidyr::unnest(.,cols = c("value")) %>% # <- added 07/20/2020
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
  # ADDED July 20th 2020 # 
      overall_vars_needed0 <-  tibble(name = as.character(unique(neat_fit$Variable)), run_var)
      fixnames <- data.frame(varnames = var_names, name = var_labels)
      overall_vars_needed <- dplyr::full_join(fixnames, overall_vars_needed0, by = "name") %>%  
        select(varnames, run_var) %>% dplyr::rename(name = varnames ) 
  # DONE with ADDED July 20th 2020 # 
  
 if (any( overall_vars_needed$run_var[is.na(overall_vars_needed$run_var)==FALSE])) {
 # if (any(overall_vars_needed$run_var)) {
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

  neat_fit <- neat_fit %>% dplyr::select(Variable, Level, dplyr::contains('CI'), dplyr::contains('P Value')) %>%
    dplyr::filter(is.na(Variable) == FALSE & is.na(Level) == FALSE )

  # Adding Title in front, if given
  if (!is.null(title_name)) dplyr::bind_cols(Name = rep(title_name,nrow(neat_fit)), neat_fit) else neat_fit
}  



#globalVariables(c("varnames","rename"))



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
#' run_pretty_model_output(x_in = 'x2', model_data = my_data, y_in = 'ybin'
#' , event_in = NULL, verbose = TRUE)
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
#' univariate_output <- purrr::map_dfr(vars_to_run,
#'  run_pretty_model_output, model_data = Bladder_Cancer, 
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
run_pretty_model_output <- function(x_in, model_data, y_in, event_in = NULL, event_level = NULL, title_name = NULL, fail_if_warning = TRUE, conf_level = 0.95, overall_p_test_stat = c('Wald', 'LR'), est_digits = 3, p_digits = 4, output_type = NULL, sig_alpha = 0.05, background = 'yellow', verbose = FALSE, ...) {
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
      # Linear Model
      tmp_fit <- tryCatch(expr =  lm(tmp_formula, data = model_data), 
                          error = function(c) stop('Linear model with "',deparse(tmp_formula), '" formula has error(s) calculating CI(s)'))
      if (fail_if_warning) {
        tmp_confint <- tryCatch(expr = suppressMessages(confint(tmp_fit)), 
                                error = function(c) stop('Linear model with "',deparse(tmp_formula), '" formula has error(s) calculating CI(s)'), 
                                warning = function(c) stop('Linear model with "',deparse(tmp_formula), '" formula has Inf CI(s); most likely a model error'))
      }
    }
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
                          # Need special error checking because coxph throws some wanrings when warnPartialMatchArgs = TRUE
                          warning = function(c) if (any(grepl('converge', c))) 
                            stop('Cox model with "',deparse(tmp_formula), '" formula has warnings(s) running')
                          else suppressWarnings(survival::coxph(tmp_formula, data = model_data)))
    } else {
      tmp_fit <- tryCatch(expr = survival::coxph(tmp_formula, data = model_data), 
                          error = function(c) stop('Cox model with "',deparse(tmp_formula), '" formula has error(s) running'))
      
    }
    n_info <- paste0('n=',tmp_fit$n,' (',tmp_fit$nevent,')')
  }
  
  tmp_output <- pretty_model_output(fit = tmp_fit, model_data = model_data, title_name = title_name, conf_level = conf_level, overall_p_test_stat = overall_p_test_stat, est_digits = est_digits, p_digits = p_digits,  output_type =  output_type, sig_alpha = sig_alpha, background = background, ...)
  tmp_output <- dplyr::bind_cols(tmp_output, n =  c(n_info, rep("", nrow(tmp_output) - 1)))
  
  if (!is.null(event_in)) names(tmp_output)[names(tmp_output) == 'n'] <- 'n (events)'
  tmp_output <- tmp_output %>% dplyr::filter(is.na(Variable) == FALSE & is.na(Level) == FALSE )
  tmp_output
}













#' Fancy Table Output of KM (survfit) Fit
#' 
#' This function takes a Kaplan-Meier model fit object (from survival::survfit) and calculate survival estimates at a specified time, and Median Survival Estimates. This can be performed on an overall KM fit or a fit including a categorical variable (strata).
#'
#' @param fit survfit object (with or without single strata variable)
#' @param time_est numerical vector of time estimates. If NULL (default) no time estimates are calculated
#' @param group_name strata variable name. If NULL and strata exists then using variable
#' @param title_name title to use
#' @param surv_est_prefix prefix to use in survival estimate names. Default is Time (i.e. Time:5, Time:10,...)
#' @param surv_est_digits number of digits to round p values for survival estimates for specified times
#' @param median_est_digits number of digits to round p values for Median Survival Estimates
#' @param output_type output type, either NULL (default), "latex", or "html" (making special charaters latex friendly)
#' 
#' @details 
#' Currently works with multiple strata in the fit (i.e. \code{survfit(Surv(time, event) ~ x1 + x2)}), although level and \code{Group} column names may be off.
#' 
#' @return
#' A tibble with: \code{Name} (if provided), \code{Group} (if strata variable in fit), \code{Level} (if strata variable in fit), \code{Median Estimate}, \code{Time:X} (Survival estimates for each time provided, if any). In no strata variable tibble is one row, otherwise nrows = number of strata levels.
#' 
#' @examples
#' 
#' # Basic linear model example
#' set.seed(542542522)
#' ybin <- sample(0:1, 100, replace = TRUE)
#' ybin2 <- sample(0:1, 100, replace = TRUE)
#' y <- rexp(100,.1)
#' x1 <- factor(sample(LETTERS[1:2],100,replace = TRUE))
#' x2 <- factor(sample(letters[1:4],100,replace = TRUE))
#' my_fit <- survival::survfit(survival::Surv(y, ybin) ~ 1)
#' my_fit2 <- survival::survfit(survival::Surv(y, ybin) ~ x1)
#' my_fit3 <- survival::survfit(survival::Surv(y, ybin) ~ x2)
#' my_fit_y2 <- survival::survfit(survival::Surv(y, ybin2) ~ 1)
#' 
#' pretty_km_output(fit = my_fit3, time_est = c(5,10), title_name = 'Overall Fit')
#' 
#' library(dplyr) 
#' km_info <- bind_rows(
#'   pretty_km_output(fit = my_fit, time_est = c(5,10), 
#'         group_name = 'Overall', title_name = 'Overall Survival---ybin'),
#'   pretty_km_output(fit = my_fit2, time_est = c(5,10), 
#'         group_name = NULL, title_name = 'Overall Survival---ybin'),
#'   pretty_km_output(fit = my_fit3, time_est = c(5,10), 
#'         group_name = 'x2', title_name = 'Overall Survival---ybin'),
#'   pretty_km_output(fit = my_fit_y2, time_est = c(5,10), 
#'         group_name = 'Overall', title_name = 'Overall Survival---ybin2'),
#' ) %>% select(Group, Level, everything())
#' 
#' kableExtra::kable(km_info, 'html', caption = 'Survival Percentage Estimates at 5 and 10 Years') %>% 
#'   kableExtra::collapse_rows(1:2, row_group_label_position = 'stack', headers_to_remove = 1:2)
#'   
#'   # Real World Examples
#'   data(Bladder_Cancer)
#' surv_obj <-survival::Surv(Bladder_Cancer$Survival_Months, Bladder_Cancer$Vital_Status == 'Dead')
#'     downstage_fit <- survival::survfit(surv_obj ~ PT0N0, data = Bladder_Cancer)
#'   
#'   pretty_km_output(fit = downstage_fit, time_est = c(24, 60), 
#'        surv_est_prefix = 'Month', surv_est_digits = 3)
#'   
#' 
#' @importFrom  dplyr %>%
#' @importFrom tibble tibble
#' @export
pretty_km_output <- function(fit, time_est = NULL, group_name = NULL, title_name = NULL, surv_est_prefix = 'Time', surv_est_digits = 2, median_est_digits = 1, output_type = NULL){

  
  # Input Checking
  if (!is.null(time_est)) {
    # Want to make sure time_est > 0
    .check_numeric_input(time_est, lower_bound = 1e-50)
  } else {
    # if no time_est requested, setting time_est to 0 and then removing later, to make coding much easier.
    time_est = 0
  }
  .check_numeric_input(surv_est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(median_est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  if (!is.null(output_type) && !output_type %in% c('latex','html'))
    stop('"output_type" must be either NULL, "latex", or "html"')
  
  # If group name not specified but using strata, will use var name
  if (is.null(group_name) && !is.null(fit$strata)) 
    group_name <- gsub('_', ' ', substr(names(fit$strata)[1], 1, regexpr('=', names(fit$strata)[1]) - 1))
  
  # Getting specific time_est estimates
  tmp_summary <- summary(fit, times = time_est, extend = TRUE)
  tmp_surv_est <-  stat_paste(tmp_summary$surv, tmp_summary$lower, tmp_summary$upper, 
                              digits = surv_est_digits, trailing_zeros = TRUE,  bound_char = '(', na_str_out = 'N.E.')
  if (!is.null(output_type) && output_type == 'latex') tmp_surv_est <- gsub('\\%','\\\\%', tmp_surv_est)
  
  if (!is.null(fit$strata)) {
    tmp_strata_levels <- substr(tmp_summary$strata, regexpr('=', tmp_summary$strata) + 1, nchar(as.vector(tmp_summary$strata)))
    tmp_all_levels <- tibble::tibble(Level = substr(names(fit$strata), regexpr('=',  names(fit$strata)) + 1, nchar( names(fit$strata))))
    tmp_med_info <- stat_paste(tmp_summary$table[,'median'], tmp_summary$table[,paste0(fit$conf.int,'LCL')], tmp_summary$table[,paste0(fit$conf.int,'UCL')], 
                               digits = median_est_digits, trailing_zeros = TRUE,  bound_char = '(', na_str_out = 'N.E.')
    n_events <- tmp_summary$table[,'events']
    # Getting max times for each level
    fit_for_max <- summary(fit, censored = TRUE)
    last_index <- length(fit_for_max$strata) - match(unique(fit_for_max$strata), rev(fit_for_max$strata)) + 1
    max_times <- dplyr::bind_cols(tmp_all_levels, max_times = fit_for_max$time[last_index])
  } else {
    tmp_strata_levels <- NULL
    tmp_med_info <-  stat_paste(tmp_summary$table['median'], tmp_summary$table[paste0(fit$conf.int,'LCL')], tmp_summary$table[paste0(fit$conf.int,'UCL')], 
                                digits = median_est_digits, trailing_zeros = TRUE,  bound_char = '(', na_str_out = 'N.E.')
    n_events <- tmp_summary$table['events']
    max_times <- max(summary(fit, censored = TRUE)$time)
  }
  tmp_surv_est_info_long <- dplyr::bind_cols(Level = tmp_strata_levels, 
                                             Time = tmp_summary$time, Est = tmp_surv_est)
  
  # Need to Replace times after last value
  if (!is.null(fit$strata)) {
    tmp_surv_est_info_long <- dplyr::full_join(tmp_surv_est_info_long, max_times, by = 'Level') %>% 
      dplyr::mutate(Est = ifelse(Time > max_times, 'N.E.', Est)) %>% dplyr::select(-max_times)
  } else {
    tmp_surv_est_info_long <- bind_cols(tmp_surv_est_info_long, max_times = rep(max_times, nrow(tmp_surv_est_info_long))) %>% 
      dplyr::mutate(Est = ifelse(Time > max_times, 'N.E.', Est)) %>% dplyr::select(-max_times)
  }
  names(tmp_surv_est_info_long)[ names(tmp_surv_est_info_long) == 'Time'] = surv_est_prefix
  
  tmp_surv_est_info <- tmp_surv_est_info_long %>% tidyr::spread_(surv_est_prefix, 'Est', sep = ':') 
  
  # Need to merge in time est with levels, in case entire strata missing. Then replace missing with N.E.
  if (!is.null(fit$strata))  tmp_surv_est_info <- dplyr::full_join(tmp_all_levels, tmp_surv_est_info, by = 'Level') 
  tmp_surv_est_info <- tmp_surv_est_info %>% replace(., is.na(.), 'N.E.')
  
  tmp_med_info <- gsub('NA', 'N.E.', tmp_med_info)
  
  # Stripping some unnecessary attributes
  attr(tmp_med_info, 'names') <- NULL
  attr(n_events, 'names') <- NULL
  
  tmp_output <- dplyr::bind_cols(Group = rep(group_name,length(tmp_med_info)),
                                 N = fit$n, `N Events` = n_events,
                                 tmp_surv_est_info,
                                 `Median Estimate` = tmp_med_info) %>% 
    # putting time variables at end
    dplyr::select(-dplyr::contains(paste0(surv_est_prefix, ':')), dplyr::everything(), dplyr::contains(paste0(surv_est_prefix, ':')))
  
  # Dropping if time = 0 (i.e. no time given)
  if (all(time_est == 0))  tmp_output <- tmp_output %>% dplyr::select(-dplyr::contains(surv_est_prefix))
  
  # Reorder if strata
  if (!is.null(fit$strata))  tmp_output <- tmp_output %>% dplyr::select(Group, Level, dplyr::everything())
  
  # Adding Title in front, if given
  if (!is.null(title_name)) dplyr::bind_cols(Name = rep(title_name,length(tmp_med_info)), tmp_output) else tmp_output
}


globalVariables(c("Est","Est_CI","Group","Level","Overall","P Value","Time","Variable", 
                  "conf.high","conf.low","est","est.label","estimate","name","overall.p.label",
                  "p.label","p.value","term","value","variable","bind_cols",'.'))




#' Wrapper for KM Model Output, with Log-Rank p value
#' 
#' This function takes a dataset, along with variables names for time and event status for KM fit, and possibly strata
#'
#' @param strata_in name of strata variable, or NA (default) if no strata desired
#' @param model_data dataset that contains \code{strata_in}, \code{time_in}, and \code{event_in} variables
#' @param time_in name of time variable component of outcome measure
#' @param event_in name of event status variable. If \code{event_level} = NULL then this must be the name of a F/T or 0/1 variable, where F or 0 are considered the censored level, respectively
#' @param event_level event level for event status variable.
#' @param time_est numerical vector of time estimates. If NULL (default) no time estimates are calculated
#' @param group_name strata variable name. If NULL and strata exists then using variable
#' @param title_name title to use
#' @param conf_level the confidence level required (default is 0.95).
#' @param surv_est_prefix prefix to use in survival estimate names. Default is Time (i.e. Time:5, Time:10,...)
#' @param surv_est_digits number of digits to round p values for survival estimates for specified times
#' @param median_est_digits number of digits to round p values for Median Survival Estimates
#' @param p_digits number of digits to round p values for Log-Rank p value
#' @param output_type output type, either NULL (default), "latex", or "html" (making special charaters latex friendly)
#' @param sig_alpha the defined significance level. Default = 0.05
#' @param background background color of significant values, or no highlighting if NULL. Default is "yellow"
#' @param ... other params to pass to \code{pretty_pvalues} (i.e. \code{bold} or \code{italic})
#
#'  
#' @return
#' A tibble with: \code{Name} (if provided), \code{Group} (if strata variable in fit), \code{Level} (if strata variable in fit), \code{Time:X} (Survival estimates for each time provided), \code{Median Estimate}. In no strata variable tibble is one row, otherwise nrows = number of strata levels.
#' 
#' @examples
#' 
#' # Basic survival model examples
#' set.seed(542542522)
#' ybin <- sample(0:1, 100, replace = TRUE)
#' ybin2 <- sample(0:1, 100, replace = TRUE)
#' ybin3 <- sample(c('Dead','Alive'), 100, replace = TRUE)
#' y <- rexp(100,.1)
#' x1 <- factor(sample(LETTERS[1:2],100,replace = TRUE))
#' x2 <- factor(sample(letters[1:4],100,replace = TRUE))
#' my_data <- data.frame(y, ybin, ybin2, ybin3, x1, x2)
#' Hmisc::label(my_data$x1) <- "X1 Variable"
#' 
#'  # Single runs 
#' run_pretty_km_output(strata_in = 'x1', model_data = my_data, 
#'      time_in = 'y', event_in = 'ybin', time_est = NULL)
#' run_pretty_km_output(strata_in = 'x1', model_data = my_data, 
#'      time_in = 'y', event_in = 'ybin', time_est = c(5,10))
#' run_pretty_km_output(strata_in = 'x2', model_data = my_data, 
#'      time_in = 'y', event_in = 'ybin3', event_level = 'Dead', time_est = c(5,10))
#' 
#' # Multiple runs for different variables
#' library(dplyr) 
#' vars_to_run = c(NA, 'x1', 'x2')
#' purrr::map_dfr(vars_to_run, run_pretty_km_output, model_data = my_data,
#'      time_in = 'y', event_in = 'ybin', event_level = '0', time_est = NULL) %>% 
#'    select(Group, Level, everything())
#'    
#' km_info <- purrr::map_dfr(vars_to_run, run_pretty_km_output, model_data = my_data, time_in = 'y', 
#'      event_in = 'ybin3', event_level = 'Dead', time_est = c(5,10), surv_est_prefix = 'Year', 
#'      title_name = 'Overall Survival') %>% 
#'    select(Group, Level, everything())
#'    
#' km_info2 <- purrr::map_dfr(vars_to_run, run_pretty_km_output, model_data = my_data, time_in = 'y', 
#'      event_in = 'ybin2', time_est = c(5,10), surv_est_prefix = 'Year', 
#'      title_name = 'Cancer Specific Survival') %>% 
#'    select(Group, Level, everything())
#' 
#' options(knitr.kable.NA = '')
#' kableExtra::kable(bind_rows(km_info, km_info2), escape = FALSE
#' , longtable = FALSE, booktabs = TRUE, linesep = '', 
#'      caption = 'Survival Percentage Estimates at 5 and 10 Years') %>%
#'   kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack'
#'   , headers_to_remove = 1:2) 
#' 
#' 
#'   # Real World Example
#'   data(Bladder_Cancer)
#'   
#'   vars_to_run = c(NA, 'Gender', 'Clinical_Stage_Grouped', 'PT0N0', 'Any_Downstaging')
#'   
#'   purrr::map_dfr(vars_to_run, run_pretty_km_output, model_data = Bladder_Cancer, 
#'        time_in = 'Survival_Months', event_in = 'Vital_Status', event_level = 'Dead', 
#'        time_est = c(24,60), surv_est_prefix = 'Month', p_digits=5) %>% 
#'    select(Group, Level, everything())
#'   
#' @importFrom  Hmisc label
#' 
#' @export
#' 
run_pretty_km_output <- function(strata_in = NA, model_data, time_in, event_in, event_level = NULL, time_est = NULL, group_name = NULL, title_name = NULL, conf_level = .95, surv_est_prefix = 'Time', surv_est_digits = 2, median_est_digits = 1, p_digits = 4, output_type = NULL, sig_alpha = 0.05, background = 'yellow', ...) {
  if (length(strata_in) != 1) stop('"strata_in" must be length of 1')
  .check_numeric_input(surv_est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(median_est_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(p_digits, lower_bound = 1, upper_bound = 14, whole_num = TRUE, scalar = TRUE)
  .check_numeric_input(sig_alpha, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  .check_numeric_input(conf_level, lower_bound = 0, upper_bound = 1, scalar = TRUE)
  if (!is.null(output_type) && !output_type %in% c('latex','html'))
    stop('"output_type" must be either NULL, "latex", or "html"')
  
  if (all(time_in != colnames(model_data)))
    stop('"time_in" must be in the "model_data" dataset')
  
  if (all(event_in != colnames(model_data)))
    stop('"event_in" must be in the "model_data" dataset')
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
  
  
  
  if (is.na(strata_in)) {
    if (is.null(group_name)) group_name <- 'Overall'
    tmp_formula <- as.formula(paste("survival::Surv(",time_in,",",event_in,") ~ 1"))
  } else {
    if (!all(na.omit(strata_in) %in% colnames(model_data)))
      stop('All variables used in the "strata_in" must be in the "model_data" dataset')
    if (is.null(group_name)) {
      #Using label if possible
      if (Hmisc::label(model_data[,strata_in]) != '') 
        group_name <- Hmisc::label(model_data[,strata_in]) 
      else group_name <- gsub('_', ' ', strata_in)
    }
    tmp_formula <- as.formula(paste("survival::Surv(",time_in,",",event_in,") ~ ", strata_in))
    tmp_pval <-  pchisq(survival::survdiff(formula = tmp_formula, data = model_data)$chisq, 
                        length(survival::survdiff(formula = tmp_formula, data = model_data)$n) - 1, 
                        lower.tail = FALSE)
    if (!is.null(output_type)) {
      tmp_p_info <- pretty_pvalues(tmp_pval, digits = p_digits, sig_alpha = sig_alpha, trailing_zeros = TRUE, 
                                   output_type = output_type, background = background, ...)
    } else {
      tmp_p_info <- pretty_pvalues(tmp_pval, digits = p_digits, trailing_zeros = TRUE, output_type = NULL)
    }
  }
  
  tmp_fit <- survival::survfit(tmp_formula, model_data, conf.int = conf_level)
  tmp_km_output <- pretty_km_output(fit = tmp_fit, time_est = time_est, group_name = group_name, title_name = title_name, surv_est_prefix = surv_est_prefix, surv_est_digits = surv_est_digits, median_est_digits = median_est_digits, output_type = output_type)
  
  if (is.na(strata_in)) 
    tmp_km_output else 
      dplyr::bind_cols(tmp_km_output, `Log-Rank P` = c(tmp_p_info, rep('', nrow(tmp_km_output) - 1)))
  
}
