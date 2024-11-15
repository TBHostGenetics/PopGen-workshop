runLaaaModelSummary <- function(model.frame) {
  return (summary(glm(Class ~  Age + Gender + SanAncestry + AfricanAncestry + EuropeanAncestry + SouthEastAsianAncestry  + allele_dose + san_dose + allele_san_dose, data=model.frame, family = "binomial")))
}

runNullModel <- function(model.frame) {
  return (glm(Class ~ Age + Gender + SanAncestry + AfricanAncestry + EuropeanAncestry + SouthEastAsianAncestry , data=model.frame, family = "binomial"))
}

runLaaaModel <- function(model.frame) {
  return (glm(Class ~ Age + Gender + SanAncestry + AfricanAncestry + EuropeanAncestry + SouthEastAsianAncestry  + allele_dose + san_dose + allele_san_dose, data=model.frame, family = "binomial"))
}

runAlleleModelSummary <- function(model.frame) {
  return (summary(glm(Class ~ Age + Gender + SanAncestry + AfricanAncestry + EuropeanAncestry + SouthEastAsianAncestry  + allele_dose, data=model.frame, family = "binomial")))
}


