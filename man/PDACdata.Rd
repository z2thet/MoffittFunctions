\name{PDACdata}
\alias{PDACdata}
\docType{data}
\title{
PDAC data set
}
\description{
Pancreatic cancer data set
}
\usage{data("PDACdata")}
\format{
  A data frame with 29571 observations on the following 184 variables.
  \describe{
    \item{\code{puf_case_id}}{a character vector}
    \item{\code{puf_facility_id}}{a character vector}
    \item{\code{facility_type_cd}}{a numeric vector}
    \item{\code{facility_location_cd}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{sex}}{a numeric vector}
    \item{\code{race}}{a numeric vector}
    \item{\code{spanish_hispanic_origin}}{a numeric vector}
    \item{\code{insurance_status}}{a numeric vector}
    \item{\code{med_inc_quar_00}}{a numeric vector}
    \item{\code{no_hsd_quar_00}}{a numeric vector}
    \item{\code{ur_cd_03}}{a numeric vector}
    \item{\code{med_inc_quar_12}}{a numeric vector}
    \item{\code{no_hsd_quar_12}}{a numeric vector}
    \item{\code{ur_cd_13}}{a numeric vector}
    \item{\code{crowfly}}{a numeric vector}
    \item{\code{charlson_recode}}{a numeric vector}
    \item{\code{sequence_number}}{a character vector}
    \item{\code{class_of_case}}{a numeric vector}
    \item{\code{year_of_diagnosis}}{a numeric vector}
    \item{\code{primary_site}}{a character vector}
    \item{\code{laterality}}{a numeric vector}
    \item{\code{histology}}{a numeric vector}
    \item{\code{behavior}}{a numeric vector}
    \item{\code{grade}}{a numeric vector}
    \item{\code{diagnostic_confirmation}}{a numeric vector}
    \item{\code{tumor_size}}{a numeric vector}
    \item{\code{regional_nodes_positive}}{a numeric vector}
    \item{\code{regional_nodes_examined}}{a numeric vector}
    \item{\code{dx_staging_proc_days}}{a numeric vector}
    \item{\code{rx_summ_dxstg_proc}}{a numeric vector}
    \item{\code{tnm_clin_t}}{a character vector}
    \item{\code{tnm_clin_n}}{a character vector}
    \item{\code{tnm_clin_m}}{a character vector}
    \item{\code{tnm_clin_stage_group}}{a character vector}
    \item{\code{tnm_path_t}}{a character vector}
    \item{\code{tnm_path_n}}{a character vector}
    \item{\code{tnm_path_m}}{a character vector}
    \item{\code{tnm_path_stage_group}}{a character vector}
    \item{\code{tnm_edition_number}}{a numeric vector}
    \item{\code{analytic_stage_group}}{a numeric vector}
    \item{\code{cs_mets_at_dx}}{a character vector}
    \item{\code{cs_mets_eval}}{a character vector}
    \item{\code{cs_extension}}{a character vector}
    \item{\code{cs_tumor_sizeext_eval}}{a character vector}
    \item{\code{cs_mets_dx_bone}}{a numeric vector}
    \item{\code{cs_mets_dx_brain}}{a numeric vector}
    \item{\code{cs_mets_dx_liver}}{a numeric vector}
    \item{\code{cs_mets_dx_lung}}{a numeric vector}
    \item{\code{lymph_vascular_invasion}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_1}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_2}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_3}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_4}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_5}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_6}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_7}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_8}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_9}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_10}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_11}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_12}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_13}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_14}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_15}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_16}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_17}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_18}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_19}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_20}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_21}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_22}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_23}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_24}}{a numeric vector}
    \item{\code{cs_sitespecific_factor_25}}{a numeric vector}
    \item{\code{cs_version_latest}}{a numeric vector}
    \item{\code{dx_rx_started_days}}{a numeric vector}
    \item{\code{dx_surg_started_days}}{a numeric vector}
    \item{\code{dx_defsurg_started_days}}{a numeric vector}
    \item{\code{rx_summ_surg_prim_site}}{a numeric vector}
    \item{\code{rx_hosp_surg_appr_2010}}{a numeric vector}
    \item{\code{rx_summ_surgical_margins}}{a numeric vector}
    \item{\code{rx_summ_scope_reg_ln_sur}}{a numeric vector}
    \item{\code{rx_summ_surg_oth_regdis}}{a numeric vector}
    \item{\code{surg_discharge_days}}{a numeric vector}
    \item{\code{readm_hosp_30_days}}{a numeric vector}
    \item{\code{reason_for_no_surgery}}{a numeric vector}
    \item{\code{dx_rad_started_days}}{a numeric vector}
    \item{\code{rx_summ_radiation}}{a numeric vector}
    \item{\code{rad_location_of_rx}}{a numeric vector}
    \item{\code{rad_treat_vol}}{a numeric vector}
    \item{\code{rad_regional_rx_modality}}{a numeric vector}
    \item{\code{rad_regional_dose_cgy}}{a numeric vector}
    \item{\code{rad_boost_rx_modality}}{a numeric vector}
    \item{\code{rad_boost_dose_cgy}}{a numeric vector}
    \item{\code{rad_num_treat_vol}}{a numeric vector}
    \item{\code{rx_summ_surgrad_seq}}{a numeric vector}
    \item{\code{rad_elapsed_rx_days}}{a numeric vector}
    \item{\code{reason_for_no_radiation}}{a numeric vector}
    \item{\code{dx_systemic_started_days}}{a numeric vector}
    \item{\code{dx_chemo_started_days}}{a numeric vector}
    \item{\code{rx_summ_chemo}}{a numeric vector}
    \item{\code{dx_hormone_started_days}}{a numeric vector}
    \item{\code{rx_summ_hormone}}{a numeric vector}
    \item{\code{dx_immuno_started_days}}{a numeric vector}
    \item{\code{rx_summ_immunotherapy}}{a numeric vector}
    \item{\code{rx_summ_trnsplnt_endo}}{a numeric vector}
    \item{\code{rx_summ_systemic_sur_seq}}{a numeric vector}
    \item{\code{dx_other_started_days}}{a numeric vector}
    \item{\code{rx_summ_other}}{a numeric vector}
    \item{\code{palliative_care}}{a numeric vector}
    \item{\code{rx_summ_treatment_status}}{a numeric vector}
    \item{\code{puf_30_day_mort_cd}}{a numeric vector}
    \item{\code{puf_90_day_mort_cd}}{a numeric vector}
    \item{\code{dx_lastcontact_death_months}}{a numeric vector}
    \item{\code{puf_vital_status}}{a numeric vector}
    \item{\code{rx_hosp_surg_prim_site}}{a numeric vector}
    \item{\code{rx_hosp_chemo}}{a numeric vector}
    \item{\code{rx_hosp_immunotherapy}}{a numeric vector}
    \item{\code{rx_hosp_hormone}}{a numeric vector}
    \item{\code{rx_hosp_other}}{a numeric vector}
    \item{\code{puf_mult_source}}{a numeric vector}
    \item{\code{reference_date_flag}}{a numeric vector}
    \item{\code{rx_summ_scope_reg_ln_2012}}{a numeric vector}
    \item{\code{rx_hosp_dxstg_proc}}{a numeric vector}
    \item{\code{palliative_care_hosp}}{a numeric vector}
    \item{\code{tumor_size_summary}}{a numeric vector}
    \item{\code{mets_at_dx_other}}{a logical vector}
    \item{\code{mets_at_dx_distant_ln}}{a logical vector}
    \item{\code{mets_at_dx_bone}}{a logical vector}
    \item{\code{mets_at_dx_brain}}{a logical vector}
    \item{\code{mets_at_dx_liver}}{a logical vector}
    \item{\code{mets_at_dx_lung}}{a logical vector}
    \item{\code{no_hsd_quar_16}}{a numeric vector}
    \item{\code{med_inc_quar_16}}{a numeric vector}
    \item{\code{medicaid_expn_code}}{a numeric vector}
    \item{\code{age70recode}}{a character vector}
    \item{\code{sex_recode}}{a labelled}
    \item{\code{hispanic_recode}}{a labelled}
    \item{\code{race_recode}}{a labelled}
    \item{\code{facility_type_recode}}{a labelled}
    \item{\code{facility_academic_binary}}{a character vector}
    \item{\code{sequence_recode}}{a character vector}
    \item{\code{vital_status}}{a labelled}
    \item{\code{surgery_recode}}{a labelled}
    \item{\code{grade_recode}}{a labelled}
    \item{\code{histology_recode}}{a labelled}
    \item{\code{clinical_stage_recode}}{a labelled}
    \item{\code{hospital_count}}{a numeric vector}
    \item{\code{hospital_volume}}{a numeric vector}
    \item{\code{surgery_year}}{a numeric vector}
    \item{\code{volume_quartiles}}{a numeric vector}
    \item{\code{insurance_recode}}{a labelled}
    \item{\code{surgery_volume_recode}}{a labelled}
    \item{\code{surgery_volume_median}}{a character vector}
    \item{\code{lvi}}{a labelled}
    \item{\code{path_n_8ed}}{a labelled}
    \item{\code{path_t}}{a character vector}
    \item{\code{path_t_recode}}{a character vector}
    \item{\code{path_t_8ed}}{a labelled}
    \item{\code{margin_recode}}{a labelled}
    \item{\code{radiation_seq_recode}}{a character vector}
    \item{\code{chemo_seq_recode}}{a character vector}
    \item{\code{chemo_neoadjuvant_recode}}{a labelled}
    \item{\code{chemo_adjuvant_recode}}{a labelled}
    \item{\code{age_adjustment}}{a character vector}
    \item{\code{caci}}{a numeric vector}
    \item{\code{caci_recode}}{a character vector}
    \item{\code{caci_categorical}}{a labelled}
    \item{\code{mortality90day}}{a labelled}
    \item{\code{volume_high_low}}{a labelled}
    \item{\code{facility_type_cd0}}{a numeric vector}
    \item{\code{num_obs_by_fac}}{a labelled}
    \item{\code{surgery_year12}}{a numeric vector}
    \item{\code{volumequantiles1}}{a factor with levels \code{(0,2.42]} \code{(2.42,7.42]} \code{(7.42,17.4]} \code{(17.4,79.7]}}
    \item{\code{volumequantiles}}{a labelled}
    \item{\code{ageadjust}}{a numeric vector}
    \item{\code{CACIcheck}}{a numeric vector}
    \item{\code{caci3}}{a labelled}
    \item{\code{chemo_adjuvant_recode0}}{a character vector}
    \item{\code{chemo_neoadjuvant_recode0}}{a character vector}
    \item{\code{lvi0}}{a character vector}
    \item{\code{dx_lastcontact_death_months0}}{a numeric vector}
    \item{\code{puf_vital_status0}}{a numeric vector}
  }
}
 
\examples{
data(PDACdata)
## maybe str(PDACdata) ; plot(PDACdata) ...
}
\keyword{datasets}
