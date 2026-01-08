#' Example Function for Test Association Between Alpha Diversity and a Covariate
#'
#' This function fits a simple linear regression model to test the association
#' between an alpha diversity value (\code{Y}) and a given covariate (\code{X}).
#' It is typically used to assess whether alpha diversity differs across groups
#' (e.g., case vs control).
#'
#' @param X A numeric or binary covariate vector (e.g., 0/1 group indicator or
#' continuous variable).
#' @param Y A numeric response vector representing alpha diversity values (e.g.,
#' Shannon index).
#'
#' @return A named numeric vector with the following components:
#' \describe{
#'   \item{Estimate}{The estimated regression coefficient for \code{X}.}
#'   \item{Std.Error}{The standard error of the coefficient.}
#'   \item{Test.Stat}{The t-statistic for testing the null hypothesis of no
#'                    association.}
#'   \item{P.val}{The p-value corresponding to the test statistic.}
#' }
#'
#' @details This function performs a univariate test of association by fitting
#' the linear model \code{Y ~ X}. It is often used as a simple association
#' screen between diversity metrics and clinical or experimental covariates.
#' @importFrom stats lm
#' @examples
#' # Simulate binary covariate and diversity values
#' set.seed(42)
#' X <- rbinom(50, 1, 0.5)  # binary group (0/1)
#' Y <- 5 + 0.8 * X + rnorm(50)  # alpha diversity with group effect
#' test.func.bin(X, Y)
#'
#' @export
test.func.bin = function(X, Y){

  lmfit.bin = summary(lm(Y ~ X))
  b = lmfit.bin$coefficients[2,"Estimate"]
  v = lmfit.bin$coefficients[2,"Std. Error"]
  test.stat = lmfit.bin$coefficients[2,"t value"]
  pval = lmfit.bin$coefficients[2,"Pr(>|t|)"]

  res = c(b, v, test.stat, pval)
  names(res) = c("Estimate", "Std.Error", "Test.Stat", "P.val")

  return(res)
}
