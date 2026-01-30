#' GLOW11M dataset
#'
#' This dataset was originally shared by Hosmer et al. (2013).
#'
#' @format A data frame with 238 rows and 16 variables:
#' \describe{
#'   \item{pair}{Matched pair identifier (each case-control pair shares the same value)}
#'   \item{sub_id}{Identification Code}
#'   \item{site_id}{Study Site (1–6)}
#'   \item{phy_id}{Physician ID code (128 unique codes)}
#'   \item{priorfrac}{History of prior fracture (1 = Yes, 0 = No)}
#'   \item{age}{Years}
#'   \item{weight}{Kilograms}
#'   \item{height}{Centimeters}
#'   \item{bmi}{Body mass index (kg/m\eqn{^2})}
#'   \item{premeno}{Menopause before age 45 (1 = Yes, 0 = No)}
#'   \item{momfrac}{Mother had hip fracture (1 = Yes, 0 = No)}
#'   \item{armassist}{Arms are needed to stand from a chair (1 = Yes, 0 = No)}
#'   \item{smoke}{Former or current smoker (1 = Yes, 0 = No)}
#'   \item{raterisk}{Self-reported risk of fracture (1 = Less than others of the same age, 2 = Same as others of the same age, 3 = Greater than others of the same age)}
#'   \item{fracscore}{Composite risk score}
#'   \item{fracture}{Any fracture in first year (1 = Yes, 0 = No)}
#' }
#' @source Dataset shared by Hosmer et al. (2013) with permission. Also appears in the \pkg{aplore3} package.
"glow11m"
