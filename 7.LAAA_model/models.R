runLaaaModelSummary <- function(model.frame) {
  return (summary(glm(Class ~  Age + Sex + AFR + SAN + allele_dose + bantu_dose + allele_bantu_dose, data=model.frame, family = "binomial")))
}

runNullModel <- function(model.frame) {
  return (glm(Class ~ Age + Sex + AFR + SAN, data=model.frame, family = "binomial"))
}

runLaaaModel <- function(model.frame) {
  return (glm(Class ~ Age + Sex + SAN + AFR + allele_dose + bantu_dose + allele_bantu_dose, data=model.frame, family = "binomial"))
}

runAlleleModelSummary <- function(model.frame) {
  return (summary(glm(Class ~ Age + Sex + SAN + AFR + allele_dose, data=model.frame, family = "binomial")))
}


