#' @title vecmatch: Vector Matching for Generalized Propensity Scores
#'
#' @description The `vecmatch` package implements vector matching methods based
#' on generalized propensity scores (GPS) for balancing multiple treatment
#' groups, along with diagnostic tools and visualization functions.
#'
#' @section Life cycle:
#' `vecmatch` is currently in a *stable* stage. The core matching and
#' diagnostic functionality is implemented and covered by tests, and the
#' package is actively maintained. Minor API changes may still occur in
#' response to user feedback and further methodological development, but
#' major breaking changes are not anticipated. Future work will focus on
#' extending the set of matching algorithms, adding further diagnostics,
#' and improving performance and documentation.
#'
#'
#' @keywords internal
"_PACKAGE"

#' @title Patients with Colorectal Cancer and Adenoma
#'
#' @description This is a synthetically generated dataset containing metadata
#'   for healthy individuals and patients diagnosed with colorectal cancer or
#'   adenomas. The primary purpose of this dataset in the context of matching is
#'   to balance the `status` groups across various covariates and achieve
#'   optimal matching quality.
#'
#' @format A data frame (`cancer`) with 1,224 rows and 5 columns:
#' \describe{
#'   \item{status}{Patient's health status, which can be one of the following:
#'   \code{healthy}, \code{adenoma}, \code{crc_benign}
#'   (benign colorectal carcinoma), or \code{crc_malignant}
#'   (malignant colorectal carcinoma).}
#'   \item{sex}{Patient's biological sex, recorded as either \code{M} (male) or
#'   \code{F} (female).}
#'   \item{age}{Patient's age, represented as a continuous numeric variable.}
#'   \item{bmi}{Patient's Body Mass Index (BMI), represented as a continuous
#'   numeric variable.}
#'   \item{smoker}{Smoking status of the patient,
#'   recorded as \code{yes} or \code{no}.}
#' }
#' @name cancer
#' @docType data
#' @usage data(cancer)
#' @source Data generated artificially
"cancer"
