#' Fancy Table Output of Mixed effects model from glmer
#' 
#' pretty_mixed_effects_glmer() takes a mixed effects model from glmer() fit object and calculates estimates, odds ratios with confidence intervals. P values are also produced. For categorical variables with 3+ levels overall Type 3 p values are calculated, in addition to p values comparing to the first level (reference).
#'
#' @param fit  glmer fit
#' @param df is data.frame or tibble  used to create model fits. Used for capturing variable labels, if they exist
#' @param rdintterm  random term in the model
#' 
#' @return
#' 
#' A data.frame with: \code{Variable}, \code{Level}, \code{ OR (95\% CI)}, \code{P Value} (for categorical variables comparing to reference), \code{Overall P Value} (for categorical variables with 3+ levels). 
#'   
#' @importFrom  dplyr %>%
#' @importFrom  dplyr case_when mutate group_by select ungroup 
#' @importFrom stats anova
#' @importFrom broom.mixed tidy
#' 
#' @examples
#' 

#'# glmer fit ####
#' library(magrittr)
#' library(ggeffects)
#' library(sjmisc)
#' library(lme4)
#' library(splines)
#' library(coxme)
#' library(car) 
#' library(survival)
#' library(broom.mixed)
#' library(rlang)
#' library(kableExtra)
#' data(PDACdata)
#' gmlmerfit <- glmer(I(PDACdata$mortality90day)=="Dead" ~   caci3
#'                   + surgery_volume_recode
#'                   + (1|puf_facility_id)
#'                   ,family = binomial(link = "logit"), 
#'                   data = PDACdata )
#'  
#' # CALL function include fit name, data frame name and random intercept term name in quotes ####
#' outputtable <- pretty_mixed_effects_glmer(gmlmerfit, PDACdata, "puf_facility_id")
#' 
#' kableExtra::kable(outputtable, 'html', caption = 'My Table') %>% 
#'    kableExtra::collapse_rows(c(1:2), row_group_label_position = 'stack')
#' 
#' # PDF output
#' #kable(outputtable, 'latex', escape = FALSE, longtable = T, booktabs = TRUE, 
#' # linesep = '', caption = 'Multivariable Mixed Effects Logistic Regression 
#' # Model Results for 90 day mortality with CACI, Volume, and Facility type')%>%
#' # footnote(number = c('OR are odds of death within 90 days
#' #  after the most definitive primary site surgery'))
#'
#' 
#' @export

pretty_mixed_effects_glmer <- function(fit,df,rdintterm){# fit is glmer object, and df = data frame is the data used.
  #rdintterm = is the name of the variable used for random intercept in quotes
  df<-df;
  fitresults<-fit;
  confidenceintervalsexp2<-tidy(fit,conf.int=TRUE, exponentiate=TRUE, effects="fixed");  
  
  car::Anova(fit,  test = "Chisq") -> LRfitanova;  
  LRanovastats <- as.data.frame(cbind(Variable = rownames(LRfitanova), 
                                      'Overall P' = MoffittFunctions::pretty_pvalues(LRfitanova$'Pr(>Chisq)') ));
  
  reppedvec <- function(fitresults,df){
    facid<-match("puf_facility_id",names(attributes(fit)$frame));
    thevars <-names(attributes(fit)$frame)[c(-1,-facid)]#attributes(LRmodel_rint0$terms)$term.labels;
    
    thevars <- gsub("`","", thevars)
    nvars <- length(thevars);
    therepeats <- rep(NA,nvars)
    for (i in 1:nvars){ therepeats[i] <- ifelse("numeric" %in% class(df[,thevars[i]]), 1, length(unique(df[,thevars[i]]))- 1)  }
    return(therepeats)
  }
  
  
  # NEED match names to get rid of facilty id ####
  facid<- match(rdintterm, names(attributes(fit)$frame))
  
  Variable <- rep(names(attributes(fit)$frame)[c(-1,-facid)], times = reppedvec(fit,df ) )
  
  confidenceintervalsexp2$Variable <- c("Intercept",Variable);
  
  lrmetable <- as.data.frame(  dplyr::left_join(confidenceintervalsexp2, LRanovastats, by = "Variable") %>% 
                                 mutate( 'P Value' = MoffittFunctions::pretty_pvalues(p.value)) %>%
                                 group_by(Variable) %>%
                                 mutate(first = 1:dplyr::n(),
                                        total = dplyr::n()) %>%
                                 mutate(
                                   'Overall P' = case_when(first == 1 & (first != total) ~ as.character(`Overall P`) ,
                                                           first == total ~ "",
                                                           first != total ~ "" ))%>%
                                 ungroup() %>%
                                 mutate(Variable = as.character(Variable),
                                        Level = as.character(term)  )
  )#end as data .frame
  
  # remove variable names from levels ####
  test   <- rep(NA,NROW(lrmetable))
  for(i in 1:NROW(lrmetable)){
    test[i] <- ifelse((lrmetable$first[i] == 1 & (lrmetable$first[i] == lrmetable$total[i])) ,
                      " ",stringr::str_replace(lrmetable[i,"Level"], lrmetable[i,"Variable"], "")  ) } 
  lrmetable$test <- test
  
  lrmetable <- lrmetable  %>%
    mutate( Level = test) %>%
    mutate( Variable = case_when(first == 1 ~ as.character(Variable),
                                 first != 1 ~ " ") ) %>%
    select(-first, -total, -test)
  # ta da! ####     
  
  confidenceintervalsexp2 <- lrmetable  
  confidenceintervalsexp2 <- confidenceintervalsexp2 %>%
    mutate('OR 95\\% CI' = paste(round(estimate,3)," (",round(conf.low,3),", ",round(conf.high,3),")",sep=""),
           'P value' = MoffittFunctions::pretty_pvalues(p.value,  background = "yellow")) %>% 
    mutate(numpvalue = as.numeric(gsub("<","",`Overall P`)),
           testOverallP = ifelse(numpvalue < 0.0500001 & is.na(numpvalue) == FALSE, paste("\\cellcolor{yellow}{",`Overall P`,"}",sep=""), `Overall P`),
           `Overall P` = testOverallP)  %>% 
    select(Variable, Level, 'OR 95\\% CI','P value','Overall P')   %>%
    mutate(Variable = gsub("_"," ", Variable))
  return(confidenceintervalsexp2)
} #end o'function

globalVariables(c("Overall P","numpvalue","total","testOverallP" ))


#' Fancy Table Output of Mixed effects COX PH model
#' 
#' pretty_mixed_effects_coxme()  Cox mixed effects model using the coxme() function.  Hazard ratios are produced with confidence intervals. P values are also produced. For categorical variables with 3+ levels overall Type 3 p values are calculated, in addition to p values comparing to the first level (reference).
#'
#' @param fit  glmer fit
#' @param df is data.frame or tibble  used to create model fits. Used for capturing variable labels, if they exist
#' @return
#' 
#' A data.frame  with:    \code{Variable}, \code{Level}, \code{HR (95\% CI)}, \code{P Value} (for categorical variables comparing to reference), \code{Overall P Value} (for categorical variables with 3+ levels). 
#'   
#' @importFrom lme4 fixef
#' @examples
#' # Mixed effect cox ph model example - random intercepts
#' library(magrittr)
#' library(ggeffects)
#' library(sjmisc)
#' library(lme4)
#' library(splines)
#' library(coxme)
#' library(car)
#' library(survival)
#' library(broom.mixed)
#' library(rlang)
#' library(kableExtra)
#' data(PDACdata)
#' set.seed(542542522)
#' #USE THIS ONE for 3 months to 2 years (24 months)
#'surv_months_obj <- survival::Surv(time = PDACdata$dx_lastcontact_death_months,
#'                                  event = PDACdata$puf_vital_status == 0)
#'    cactfit  <- coxme(surv_months_obj ~   
#'      caci3 + surgery_volume_recode +  (1|puf_facility_id),
#'      data = PDACdata,
#'                     control=coxme.control(eps = 1e-04, toler.chol = .Machine$double.eps^0.75,
#'                             iter.max = 200, inner.iter = Quote(max(4, fit0$iter+1)),
#'                             sparse.calc = NULL,
#'                             optpar = list(method = "BFGS", control=list(reltol = 1e-3)),
#'                             refine.df=4, refine.detail=FALSE, refine.method="control",
#'                             sparse=c(50, .02),
#'                             varinit=c(.02, .1, .4, .8)^2, corinit = c(0, .3))  )
#'  cactfitcmetable <-  pretty_mixed_effects_coxme(cactfit, PDACdata)
#' # for PDF output 
#' #  kable(cactfitcmetable , 'latex', escape = FALSE, longtable = T, booktabs = TRUE, linesep = '', 
#' #  caption = 'Full Mulitvariable Mixed effects Cox Proportional-Hazards Regression Model with 
#' #  interaction for Overall Survival2 to 10 years')
#' 
#' @export

pretty_mixed_effects_coxme <- function(fit, df){# fit is glmer object, and data frame is the data used. 
  
  #code cme table ####
  
  mod<-fit;
  extract_coxme_table <- function (mod){
    beta <- fixef(mod)
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar;
     #  print("the class of mod$var is: ")
     # print(class(mod$var))
    # print(diag(mod$var)[1:3])
    # print(as.numeric(diag(mod$var)[1:3]))
    # print(class(diag(mod$var)[nfrail + 1:nvar]))
    #  
    # print(slotNames(mod$var))
    # print(class(as.matrix(mod$var )))
    se <- sqrt(diag(as.matrix(mod$var )))[nfrail + 1:nvar] #bdsmatrix
    z<- round(beta/se, 2)
    p<- signif(1 - pchisq((beta/se)^2, 1), 2)
    tableout=data.frame(cbind(beta,se,z,p))
    return(tableout)
  }
  reppedvec <- function(fitresults, df){
    thevars <- attributes(fitresults$terms)$term.labels;
    nvars <- length(attributes(fitresults$terms)$term.labels); 
    thevars <- gsub("`","", thevars)
    therepeats <- rep(NA,nvars)
    for (i in 1:nvars){#print(i);
      therepeats[i] <- ifelse("numeric" %in% class(df[,thevars[i]]), 1, length(unique(df[,thevars[i]]))- 1)
    }
    return(therepeats)
  }
  Variable <- rep(attributes(fit$terms)$term.labels, times = reppedvec(fit, df) )
  testout <- cbind(round(cbind(exp( fixef(fit)), exp(confint(fit, level = 0.95))), 3), 'P Value'= extract_coxme_table(fit)$p)
  
  coxmeoutput00 <- as.data.frame(testout) %>% mutate('HR 95\\%' = paste(testout[, 1]," (",testout[, 2],", ",testout[, 3],")",sep=""))
  coxmeoutput0 <- cbind(Variable, Level = rownames(testout), coxmeoutput00[, c(5,4)] )
  
  anova(fit) -> fitanova
  anovastats <- as.data.frame(cbind(Variable = rownames(fitanova)[-1],
                                    'Overall P' = MoffittFunctions::pretty_pvalues(fitanova$"Pr(>|Chi|)"[-1]) ))
  
  #MERGE in overall pvalue from anova function
  
  cmetable <- as.data.frame(  dplyr::left_join(coxmeoutput0, anovastats, by = "Variable")  %>%
                                mutate( 'P Value' = MoffittFunctions::pretty_pvalues(`P Value`, background = "yellow")) %>%
                                group_by(Variable) %>%
                                mutate(first = 1:dplyr::n(),
                                       total = dplyr::n()) ) %>%
    mutate('Overall P' = case_when(first == 1 & (first != total) ~ as.character(`Overall P`) ,
                                   first == total ~ "",
                                   first != total ~ ""  ) )  %>%
    ungroup() %>%
    mutate(Variable = as.character(Variable),
           Level = as.character(Level)  )
  
  #remove variable names from levels
  test   <- rep(NA,NROW(cmetable))
  for(i in 1:NROW(cmetable)){
    test[i] <- ifelse((cmetable$first[i] == 1 & (cmetable$first[i] == cmetable$total[i])) ,
                      " ", stringr::str_replace(cmetable[i,2], cmetable[i,1], "")  ) }
  cmetable$test <- test
  
  cmetable <- cmetable  %>%
    mutate( Level = test) %>%
    mutate( Variable = case_when(first == 1 ~ as.character(Variable),
                                 first != 1 ~ " ") ) %>% 
    mutate(numpvalue = as.numeric(gsub("<","",`Overall P`)),
           testOverallP = ifelse(numpvalue < 0.0500001 & is.na(numpvalue) == FALSE, 
                                 paste("\\cellcolor{yellow}{",`Overall P`,"}",sep=""), `Overall P`),
           `Overall P` = testOverallP)  %>% 
    select(  -first, -total, -test, -numpvalue, -testOverallP) %>%
    mutate(Variable = gsub("_"," ", Variable))
  
  return(cmetable)
} # End pretty coxme o'function
