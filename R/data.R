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

#' Ticlopidine Bioequivalence Data
#'
#' @docType data
#' @usage data(ticlopidine)
#'
#' @description Bioequivalence data for ticlopidine. [PLACEHOLDER: Source and detailed description needed]
#'
#' @format A data frame. [PLACEHOLDER: Specify dimensions and variable descriptions]
#' @references [PLACEHOLDER: Add source reference]
#'
#' @examples
#' data(ticlopidine)
"ticlopidine"

#' Multivariate Skin Bioequivalence Data
#'
#' @docType data
#' @usage data(skin_mvt)
#'
#' @description Multivariate skin bioequivalence data measuring drug delivery across multiple
#' skin layers. Contains measurements from stratum corneum, viable epidermis, upper dermis,
#' and lower dermis. [PLACEHOLDER: Source and detailed description needed]
#'
#' @format A `data.frame` with 12 rows and 4 columns:
#' \describe{
#'   \item{stratum corneum}{Log-transformed drug delivery measurement in stratum corneum layer.}
#'   \item{viable epidermis}{Log-transformed drug delivery measurement in viable epidermis layer.}
#'   \item{upper dermis}{Log-transformed drug delivery measurement in upper dermis layer.}
#'   \item{lower dermis}{Log-transformed drug delivery measurement in lower dermis layer.}
#' }
#' @references [PLACEHOLDER: Add source reference]
#'
#' @examples
#' data(skin_mvt)
#' head(skin_mvt)
"skin_mvt"

#' Biodistribution Data
#'
#' @docType data
#' @usage data(biodistribution)
#'
#' @description Biodistribution data for bioequivalence analysis.
#' [PLACEHOLDER: Source and detailed description needed]
#'
#' @format A `data.frame` with 6 rows and 2 columns:
#' \describe{
#'   \item{X}{Measurement for treatment X.}
#'   \item{Y}{Measurement for treatment Y.}
#' }
#' @references [PLACEHOLDER: Add source reference]
#'
#' @examples
#' data(biodistribution)
"biodistribution"
