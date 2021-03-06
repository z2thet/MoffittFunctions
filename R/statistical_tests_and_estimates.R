#' Rounding Using Round Away From 0 Method
#'
#' round_away_0 takes a numeric vector, rounds them to a specified digit amount using the round away from 0 method for ties (i.e. 1.5). This is the SAS method for rounding.
#'
#' @param x numeric vector (can include NA values).
#' @param digits positive integer of length 1 between 0 (default) and 14, giving the amount of digits to round to.
#' @param trailing_zeros logical indicating if trailing zeros should included (i.e. 0.100 instead of 0.1). Note is set to TRUE output is a character vector
#' @return if \code{trailing_zeros = TRUE} returns a character vector of rounded values with trailing zeros, otherwise returns a numeric vector of rounded values.
#' @details
#'
#' \code{round_away_0} is not designed for use at precision levels <= 1e-15
#'
#' @examples
#'
#' vals_to_round = c(NA,-3.5:3.5,NA)
#' # [1]   NA -3.5 -2.5 -1.5 -0.5  0.5  1.5  2.5  3.5   NA
#'
#' # round() will round to even numbers when tied at X.5
#' round(vals_to_round)
#' # [1] NA -4 -2 -2  0  0  2  2  4 NA
#'
#' # round_away_0() will round away from 0 when tied at X.5
#' round_away_0(vals_to_round)
#' # [1] NA -4 -3 -2 -1  1  2  3  4 NA
#'
#' # Can force trailing zeros (will output character vector)
#' round_away_0(vals_to_round, digits = 2, trailing_zeros = TRUE)
#'
#' @export

round_away_0 <- function(x, digits = 0, trailing_zeros = FALSE){
  .check_numeric_input(x, allow_NA = TRUE)
  .check_numeric_input(digits, lower_bound = 0, upper_bound = 14, scalar = TRUE, whole_num = TRUE)

  rounded_vals <- sign(x) * round(abs(x) + 1e-15, digits)

  if (trailing_zeros) {
    # Need to exclude NAs when doing formatting
    rounded_vals[!is.na(rounded_vals)] <- formatC(rounded_vals[!is.na(rounded_vals)], digits, format = 'f')

    #Need to change -0.00... to 0.00...
    neg_to_change <- paste0('-0.',paste0(rep(0,digits), collapse = ''))
    if (any(rounded_vals == neg_to_change, na.rm = TRUE)) rounded_vals[rounded_vals == neg_to_change] <- substr(neg_to_change, 2, nchar(neg_to_change))
  }
  rounded_vals
}

#' Wrapper for round_away_0 to account for non-numeric values
#'
#' Internal wrapper for round_away_0
#' 
#' @noRd
#' 
#' @param x vector (can include NA values).
#' @param digits positive integer of length 1 between 0 and 14, giving the amount of digits to round to.
#' @param trailing_zeros logical indicating if trailing zeros should included (i.e. 0.100 instead of 0.1). Note is set to TRUE output is a character vector
#' @return if \code{x} non-numeric vector returns x, else if \code{trailing_zeros = TRUE} returns a character vector of rounded values with trailing zeros, otherwise returns a numeric vector of rounded values.
#' @details
#'
#' \code{round_away_0} is not designed for use at precision levels <= 1e-15
#'
#' @examples
#'
#' .round_if_numeric(c(NA,-3.5:3.5,NA))
#' .round_if_numeric(c(NA,letters[1:5],NA))
#'
#'

.round_if_numeric <- function(x, digits = 0, trailing_zeros = FALSE){
  if (is.numeric(x)) round_away_0(x, digits = digits, trailing_zeros = trailing_zeros) else x
}


#' Continuous Variable Compared to Binary Variable Test
#'
#' Either Wilcox or T-Test Performed, for unpaired or paired data
#'
#' @param x numeric vector (can include NA values).
#' @param y vector with only 2 levels (can include NA values unless \code{paired = TRUE}).
#' @param method what test to run (wilcox or t-test).
#' @param paired a logical indicating whether you want a paired test.
#' @param verbose a logical variable indicating if warnings and messages should be displayed.
#' @param ... parameters to pass to wilcox_test or t.test functions. For example the testing direction (\code{alternative}) in either call or the \code{var.equal} in the t.test function.
#' @return p-value for comparing x at the different levels of y.
#' @details
#'
#' Runs \code{wilcox_test()} in the coin package, with "exact" distribution and mid-ranks ties method.
#'
#' For one sided tests if \code{y} is a factor variable the level order is respected, otherwise the levels will set to alphabetical order (i.e. if \code{alternative = less} then testing a < b ).
#'
#' If \code{paired = TRUE} assumes the first observations of the first group matches the first observation of the second group, and so on. Also if \code{paired = TRUE} then \code{y} must have the same number of samples for each level.
#'
#' @examples
#'
#' set.seed(5432322)
#' x <- c(rnorm(10,0,3), rnorm(10,3,3))
#' y <- c(rep('a', 10), rep('b', 10))
#' two_samp_cont_test(x = x, y = y, method = 'wilcox', paired = FALSE)
#' two_samp_cont_test(x = x, y = y, method = 'wilcox', paired = TRUE)
#' two_samp_cont_test(x = x, y = y, method = 't', paired = FALSE)
#' two_samp_cont_test(x = x, y = y, method = 't', paired = TRUE)
#'
#' @export


two_samp_cont_test <- function(x, y, method = c('wilcox', 't.test'), paired = FALSE, verbose = FALSE, ...){
  # Input checking
  method <- match.arg(method)
  .check_numeric_input(x)
  .check_binary_input(y, paired = paired)
  y <- droplevels(factor(y))

  # Removing cases where x and y are both NA and returning p-value where no complete cases or only one distinct value
  rm_na_and_check_output <- .rm_na_and_check(x, y, y_type = 'binary', verbose = verbose)
  if (is.data.frame(rm_na_and_check_output)) data_here <- rm_na_and_check_output else return(rm_na_and_check_output)

  if (method == 'wilcox') {
    if (paired) {
      as.double(coin::pvalue(coin::wilcoxsign_test(data_here$x[data_here$y == levels(data_here$y)[1]] ~ data_here$x[data_here$y == levels(data_here$y)[2]], distribution = "exact", zero.method = "Pratt", ...)))
    } else {
      as.double(coin::pvalue(coin::wilcox_test(x ~ y, data = data_here, distribution = "exact", ties.method = "mid-ranks", ...)))
    }
  } else {
    #If both groups have only one distinct value t.test will throw error
    if (any(by(data_here$x[!is.na(data_here$y)], data_here$y[!is.na(data_here$y)], function(xx) {length(unique(xx[!is.na(xx)])) > 1}))) {
      as.double(t.test(data_here$x[data_here$y == levels(data_here$y)[1]], data_here$x[data_here$y == levels(data_here$y)[2]], data = data_here, paired = paired, ...)$p.value)
    } else {
      if (verbose) message('t.test can not run when both levels of "y" have only 1 unique "x" value, so p=NA returned')
      NA
    }
  }
}





#' Binary (Response) Variable Compared to Binary (Group) Variable Test
#'
#' Either Barnard, Fisher's, or Chi-sq test performed for unpaired data and McNemar's test for paired data
#'
#' @param x vector with only 2 levels (can include NA values).
#' @param y vector with only 2 levels (can include NA values unless \code{method = 'mcnemar'}).
#' @param method what test to run ("barnard", "fisher" ,"chi.sq" , "mcnemar"). No default so user must enter one of the four selections
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter. Only "two.sided" available for \code{method = 'chi.sq' or 'mcnemar'}
#' @param verbose a logical variable indicating if warnings and messages should be displayed.
#' @param ... parameters to pass to wilcox_test or t.test functions. For example the testing direction (\code{alternative}) in either call or the \code{var.equal} in the t.test function.
#' @return p-value for comparing x at the different levels of y.
#' @details
#'
#'
#' For one sided tests if \code{y} is a factor variable the level order is respected, otherwise the levels will set to alphabetical order (i.e. if \code{alternative = less} then testing a < b ).
#'
#' If \code{method = 'mcnemar'} assumes the first observations of the first group matches the first observation of the second group, and so on. Also if \code{method = 'mcnemar'} then \code{y} must have the same number of samples for each level.
#'
#' @examples
#'
#' set.seed(5432322)
#' x <- c(sample(0:1,10,replace = TRUE, prob = c(.75,.25)), 
#'    sample(0:1,10,replace = TRUE, prob = c(.25,.75)))
#' y <- c(rep('a', 10), rep('b', 10))
#' two_samp_bin_test(x,y, method = 'barnard')
#' two_samp_bin_test(x,y, 'fisher')
#' two_samp_bin_test(x,y, 'chi.sq')
#' two_samp_bin_test(x,y, 'mcnemar')
#'
#' @export


two_samp_bin_test <- function(x, y, method = NA, alternative = c("two.sided", "less", "greater"), verbose = FALSE, ...){
  # Input checking
  if (!method %in% c('barnard', 'fisher' ,'chi.sq' , 'mcnemar')) stop('"method" must be one of these choices: "barnard", "fisher", "chi.sq", "mcnemar"')
  alternative <- match.arg(alternative)
  if (method == 'chi.sq' & alternative != 'two.sided') stop('When "method" is chi.sq then "alternative" must be two.sided')
  if (method == 'mcnemar' & alternative != 'two.sided') stop('When "method" is mcnemar then "alternative" must be two.sided')
  .check_binary_input(x)
  .check_binary_input(y, paired = ifelse(method == 'mcnemar', TRUE, FALSE))
  y <- droplevels(factor(y))

  # Removing cases where x and y are both NA and returning p-value where no complete cases or only one distinct value
  rm_na_and_check_output <- .rm_na_and_check(x, y, x_type = ifelse(method == 'barnard', 'fixed_binary', 'binary'), y_type = 'binary', verbose = verbose)
  if (is.data.frame(rm_na_and_check_output)) data_here <- rm_na_and_check_output else return(rm_na_and_check_output)

  if (method == 'barnard') {
    pval_out <- as.double(Exact::exact.test(table(data_here), method = 'Z-pooled', to.plot = FALSE, alternative = alternative)$p.value)
  }
  if (method == 'fisher') {
    pval_out <- as.double(fisher.test(data_here$x, data_here$y, alternative = alternative)$p.value)
  }
  if (method == 'chi.sq') {
    pval_out <- as.double(chisq.test(data_here$x, data_here$y)$p.value)
  }
  if (method == 'mcnemar') {
    pval_out <- as.double(mcnemar.test(data_here$x[data_here$y == levels(data_here$y)[1]], data_here$x[data_here$y == levels(data_here$y)[2]])$p.value)
  }
  pval_out
}




#' Correlation Test for Two Continuous Variables
#'
#' A wrapper for cor.test function, except if spearman selected and ties in at least one variable, in which case this is a wrapper for coin::spreaman_test in with approximate method.
#'
#'
#' @param x numeric vector (can include NA values)
#' @param y numeric vector (can include NA values)
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated
#' @param seed seed (only used if \code{method = "spearman"})
#' @param B number of reps (only used if \code{method = "spearman"})
#' @param exact Should exact method be used. Ingorned it \code{method = "spearman"} and ties present
#' @param verbose a logical variable indicating if warnings and messages should be displayed
#' @return spearman_test pvalue
#' @details
#'
#' To always get reproducible results when using approximate method we need to set seed inside of the call, and order the data
#'
#' @examples
#'
#' set.seed(5432322)
#' x <- rnorm(20,0,3)
#' y <- x + rnorm(20,0,5)
#' cor_test(x,y, method = 'pearson')
#' cor_test(x,y, method = 'kendall')
#' cor_test(x,y, method = 'spearman')
#'
#' @export


cor_test <- function(x, y, method = c("pearson", "kendall", "spearman"), seed = 68954857, B = 10000, exact = TRUE, verbose = FALSE){
  method <- match.arg(method)
  .check_numeric_input(x)
  .check_numeric_input(y)
  if (method == "spearman") {
    .check_numeric_input(seed, lower_bound = -2^30, upper_bound = 2^30, scalar = TRUE, whole_num = TRUE)
    .check_numeric_input(B, lower_bound = 1, upper_bound = 2^20, scalar = TRUE, whole_num = TRUE)
  }

  rm_na_and_check_output <- .rm_na_and_check(x, y, x_type = 'continuous', y_type = 'continuous', verbose = verbose)
  if (is.data.frame(rm_na_and_check_output)) data_here <- rm_na_and_check_output else return(rm_na_and_check_output)

  if (method == "spearman") {
    #if no ties calling cor.test, otherwise spearman test
    if (any(duplicated(data_here$x)) | any(duplicated(data_here$y))) {
      set.seed(seed)
      as.double(coin::pvalue(coin::spearman_test(x~y, data = data_here, distribution = coin::approximate(B))))
    } else {
      as.double(cor.test(data_here$x,data_here$y, method = 'spearman', exact = exact)$p.value)
    }
  } else {
    as.double(cor.test(data_here$x,data_here$y, method = method, exact = exact)$p.value)
  }
}


