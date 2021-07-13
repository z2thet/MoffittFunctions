\name{Bladder_Cancer}
\alias{Bladder_Cancer}
\docType{data}
\title{
Bladder_Cancer
}
\description{
  A concise  description of the dataset. 
}
\usage{data("Bladder_Cancer")}
\format{
  A data frame with 166 observations on the following 37 variables.
  \describe{
    \item{\code{PTID}}{a numeric vector}
    \item{\code{Age_At_Diagnosis}}{a labelled}
    \item{\code{Gender}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{Race}}{a factor with levels \code{White} \code{Black} \code{Other}}
    \item{\code{Ethnicity}}{a factor with levels \code{Hispanic/Latino} \code{Non-Hispanic/Latino}}
    \item{\code{Marital_Status}}{a factor with levels \code{Married/living together} \code{Separated/divorced} \code{Single (never married)} \code{Widowed}}
    \item{\code{Education_Status}}{a factor with levels \code{High School or Less} \code{College or Some College} \code{Graduate/Professional Degree} \code{Unknown}}
    \item{\code{Primary_Insurance}}{a factor with levels \code{Private} \code{Medicare} \code{Medicaid} \code{Self-paying/uninsured/other}}
    \item{\code{Primary_Insurance_More_Cats}}{a factor with levels \code{Private} \code{Medicare} \code{Medicaid} \code{Self-paying/other} \code{Uninsured}}
    \item{\code{Surgery_Year}}{a labelled}
    \item{\code{Elix_Sum}}{a labelled}
    \item{\code{Comorbidities}}{a labelled}
    \item{\code{Histology}}{a factor with levels \code{Urothelial Neoplasia} \code{UC w/Squamous Differentiation} \code{UC w/Glandular Differentiation} \code{Micropapillary Carcinoma} \code{Nested Variant} \code{Plasmacytoid} \code{Sarcomatoid Carcinoma/Carcinosarcom} \code{Neuroendocrine Carcinoma}}
    \item{\code{Histology_Grouped}}{a factor with levels \code{Pure Urothelial Neoplasia} \code{Mixed Tumors} \code{Variant Histology}}
    \item{\code{Clinical_Stage}}{a factor with levels \code{Stage I (<=T1NxMx)} \code{Stage II (T2NxMx)} \code{Stage III (T3NxMx)} \code{Stage IV (T4NxMx)}}
    \item{\code{Clinical_Stage_Grouped}}{a factor with levels \code{Stage I/II (<=T2NxMx)} \code{Stage III (T3NxMx)} \code{Stage IV (T4NxMx)}}
    \item{\code{Pathologic_Stage}}{a factor with levels \code{Stage 0 (T0/Ta/isN0M0)} \code{Stage I (T1N0M0)} \code{Stage II (T2N0M0)} \code{Stage III (T3N0M0)} \code{Stage IV (T4N0-3M0-1)}}
    \item{\code{Neoadjuvant_Chemo}}{a factor with levels \code{Yes}}
    \item{\code{Chemo_Cat}}{a factor with levels \code{GEM/CIS}}
    \item{\code{Adjuvant_Chemo}}{a labelled}
    \item{\code{Cycles}}{a labelled}
    \item{\code{Cycles_cat}}{a factor with levels \code{<3} \code{3} \code{4+}}
    \item{\code{Cycles_cat_4group}}{a factor with levels \code{<3} \code{3} \code{4} \code{5+}}
    \item{\code{Urinary_Diversion}}{a factor with levels \code{Ileal Conduit} \code{Pouch} \code{Neobladder} \code{Cutaneous Ureterostomy} \code{Sigmond}}
    \item{\code{Urinary_Diversion_Grouped}}{a factor with levels \code{Ileal Conduit} \code{Pouch} \code{Neobladder} \code{Other}}
    \item{\code{PT0N0}}{a factor with levels \code{No Completed Response} \code{Complete Response}}
    \item{\code{Any_Downstaging}}{a factor with levels \code{No Downstaging} \code{Downstaging}}
    \item{\code{Path_N_Stage}}{a factor with levels \code{pNX} \code{pN0} \code{pN1} \code{pN2} \code{pN3} \code{Missing/Other}}
    \item{\code{Lympho_invasion}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Peri_invasion}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Carc_in_situ}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Positive_Margins}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Vital_Status}}{a labelled}
    \item{\code{Cancer_Specific_Vital_Status}}{a labelled}
    \item{\code{Survival_Days}}{a labelled}
    \item{\code{Survival_Months}}{a labelled}
    \item{\code{Survival_Years}}{a labelled}
  }
}
\details{  more details than the description above 
}
\source{ Not sure
}
\references{ Don't know
}
\examples{
data(Bladder_Cancer)
## maybe str(Bladder_Cancer) ; plot(Bladder_Cancer) ...
}
\keyword{datasets}
