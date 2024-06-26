#' Log transformed cutaneous delivery of econazole (ECZ) from bioequivalent products on porcine skin
#'
#' @docType data
#' @usage data(skin)
#'
#' @description Original data were collected in the same way as described in Quariter et. al. (2019),
#' and represents cutaneous delivery of econazole nitrate (ECZ in ng/cm^2) on porcine skin from a reference medicinal
#' product and an approved bioequivalent product. The dataset contains 17 pairs of comparable porcine skin samples on which
#' measurement of ECZ deposition was gathered, and log transformed, using both creams.
#'
#' @format A `data.frame` with 16 rows and 2 columns:
#' \describe{
#'   \item{Reference}{Econazole nitrate delivery for the reference product.}
#'   \item{Generic}{Econazole nitrate delivery for the generic bioequivalent product.}
#'   \item{Obs.}{The observation corresponds to a given skin on which the log ECZ delivery was collected for each of the reference and generic cream.}
#' }
#' @references Quartier, Julie, et al. "Cutaneous Biodistribution: A High-Resolution Methodology
#' to Assess Bioequivalence in Topical Skin Delivery", Pharmaceutics, (2019).
#' Boulaguiem, Younes, et al. "Finite Sample Adjustments for Average Equivalence Testing", Statistics in Medicine, (2023)
#'
#' @examples
#' data(skin)
#' theta <- diff(apply(skin,2,mean))
#' nu <- nrow(skin)-1
#' sigma_nu <- sd(apply(skin,1,diff))/sqrt(nu)
"skin"

#' ????????
#'
#' @docType data
#' @usage data(ticlopidine)
#'
#' @description ??????
#'
#' @format ????
#' @references ????
#'
#' @examples
#' data(ticlopidine)
"ticlopidine"
