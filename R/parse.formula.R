#' Parse the formula "response ~ covariates | groups"
#' to lists: response, covariates, groups
#'
#' syntax var1:var2 will add 2 variables to a corresponding list
#'
#' @noRd
# @export
parse.formula <- function(formula) {
  vars <- terms(as.formula(formula))
  y <- if(attr(vars, "response"))
    nlme::getResponseFormula(formula)
  x <- nlme::getCovariateFormula(formula)
  z <- nlme::getGroupsFormula(formula)
  list(response = all.vars(y),
       covariates = all.vars(x),
       groups = all.vars(z))
}
