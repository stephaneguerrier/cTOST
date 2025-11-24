#' Log-transformed cutaneous delivery of econazole (ECZ) from bioequivalent products on porcine skin
#'
#' @docType data
#' @usage data(skin)
#'
#' @description Original data were collected using the cutaneous biodistribution method described in Quartier et. al. (2019),
#' and represents cutaneous delivery of econazole nitrate (ECZ in ng/cm^2) on porcine skin from a reference medicinal
#' product and an approved bioequivalent product. The dataset contains 17 pairs of comparable porcine skin samples on which
#' measurement of ECZ deposition was gathered, and log transformed, using both creams.
#'
#' @format A `data.frame` with 17 rows and 2 columns:
#' \describe{
#'   \item{Reference}{Log-transformed econazole nitrate delivery for the reference product.}
#'   \item{Generic}{Log-transformed econazole nitrate delivery for the generic bioequivalent product.}
#' }
#' @references
#' Quartier J. et al. "Cutaneous Biodistribution: A High-Resolution Methodology
#' to Assess Bioequivalence in Topical Skin Delivery", Pharmaceutics, 2019.
#'
#' Boulaguiem Y. et al. "Finite Sample Adjustments for Average Equivalence Testing", Statistics in Medicine, 2024.
#'
#' @examples
#' data(skin)
#' theta <- diff(apply(skin,2,mean))
#' n <- nrow(skin)
#' nu <- n - 1
#' sigma <- var(apply(skin,1,diff)) / n
"skin"

#' Ticlopidine hydrochloride data from a crossover design
#'
#' @docType data
#' @usage data(ticlopidine)
#'
#' @description Original data were collected in Marzo et. al. (2002) to assess the
#' bioequivalence of a new formulation of ticlopidine hydrochloride with the formulation
#' that was marketed at that time. The original study involved 24 healthy male
#' volunteers who received both formulations in the form of a tablet containing 250 mg
#' of active ingredient in a 2 by 2 by 2 crossover design. The dataset contains differences between
#' the two formulations for each pair of pharmacokinetic outcomes after applying the
#' logarithmic transformation, where evident outliers were removed bringing the available sample size to 20.
#'
#' @format A `data.frame` with 20 rows and 4 columns:
#' \describe{
#'   \item{t_half}{Differences for the log-transformed elimination half-life.}
#'   \item{AUC}{Differences for the log-transformed area under the concentration-time curve from time zero to the last measurable concentration.}
#'   \item{AUC_inf}{Differences for the log-transformed area under the concentration-time curve from time zero to infinity.}
#'   \item{C_max}{Differences for the log-transformed maximum plasma concentration.}
#' }
#' @references
#' Marzo A. et al. "Bioequivalence of ticlopidine hydrochloride administered in single
#' dose to healthy volunteers", Pharmacological Research, 2002.
#'
#' Pallmann P. and Jaki T. "Simultaneous confidence regions for multivariate bioequivalence",
#' Statistics in Medicine, 2017.
#'
#' Boulaguiem Y. et al. "Multivariate adjustments for average equivalence testing",
#' Statistics in Medicine, 2025.
#'
#' @examples
#' data(ticlopidine)
#' n <- nrow(ticlopidine)
#' nu <- n - 1
#' theta <- apply(ticlopidine,2,mean)
#' Sigma <- cov(ticlopidine)/n
"ticlopidine"

#' Multivariate measurements for log-transformed cutaneous delivery of econazole (ECZ) from bioequivalent products on porcine skin
#'
#' @docType data
#' @usage data(skin_mvt)
#'
#' @description Original data were collected using the cutaneous biodistribution method described in Quartier et. al. (2019),
#' and represents cutaneous delivery of econazole nitrate (ECZ in ng/cm^2) on porcine skin from a reference medicinal
#' product and an approved bioequivalent product across multiple skin layers (from the surface to a depth of approximately 800 μm).
#' The dataset contains differences of log-transformed measurements of ECZ deposition for the two creams,
#' across 12 comparable porcine skin samples and four anatomical regions.
#'
#' @format A `data.frame` with 12 rows and 4 columns:
#' \describe{
#'   \item{`stratum corneum`}{Differences of log-transformed econazole nitrate delivery for the two creams across the stratum corneum (0-20 μm).}
#'   \item{`viable epidermis`}{Differences of log-transformed econazole nitrate delivery for the two creams across the viable epidermis (20-160 μm).}
#'   \item{`upper dermis`}{Differences of log-transformed econazole nitrate delivery for the two creams across the upper dermis (160-400 μm).}
#'   \item{`lower dermis`}{Differences of log-transformed econazole nitrate delivery for the two creams across the lower dermis (400-800 μm).}
#' }
#' @references
#' Quartier J. et al. "Cutaneous Biodistribution: A High-Resolution Methodology
#' to Assess Bioequivalence in Topical Skin Delivery", Pharmaceutics, 2019.
#'
#' Insolia L. et al. "Bioequivalence Assessment for Locally Acting Drugs: A Framework for Feasible
#' and Efficient Evaluation", arXiv, 2025.
#'
#' @examples
#' data(skin_mvt)
#' n <- nrow(skin_mvt)
#' nu <- n - 1
#' theta <- apply(skin_mvt,2,mean)
#' Sigma <- cov(skin_mvt) / n
"skin_mvt"

#' Log-transformed cutaneous delivery of "Molecule X" on human skin as measured by two operators
#'
#' @docType data
#' @usage data(biodistribution)
#'
#' @description Original data were collected using the cutaneous biodistribution method described in Quartier et. al. (2019),
#' and represents cutaneous delivery of "Molecule X" (in ng/cm^2) on human abdominal skin.
#' The study in Wu et al. (2025) compared the reproducibility of an identical experimental protocol performed by two different operators.
#' The dataset contains measurements from 6 independent skin samples per operator, after a log-transformation.
#'
#' @format A `data.frame` with 6 rows and 2 columns:
#' \describe{
#'   \item{X}{Log-transformed measurements from the "reference" operator.}
#'   \item{Y}{Log-transformed measurements from the "target" operator.}
#' }
#' @references
#' Quartier J. et al. "Cutaneous Biodistribution: A High-Resolution Methodology
#' to Assess Bioequivalence in Topical Skin Delivery", Pharmaceutics, 2019.
#'
#' Wu J. et al. "Bridging the gap between experimental burden and statistical power for quantiles
#' equivalence testing", arXiv, 2025.
#'
#' @examples
#' data(biodistribution)
#' x_bar <- mean(biodistribution$X)
#' x_sd <- sd(biodistribution$X)
#' y_bar <- mean(biodistribution$Y)
#' y_sd <- sd(biodistribution$Y)
"biodistribution"
