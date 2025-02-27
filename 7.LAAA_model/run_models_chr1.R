source("models.R")
out.file.name <- "/mnt/lustre/users/yswart/Dosages_popgen/BANTU/LAAA_bantu_chr1.txt"
pheno.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/Phenotype.txt", stringsAsFactors = F)
allele.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/BANTU/chr1_allele_dose.txt", stringsAsFactors = F)
bantu.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/BANTU/chr1_bantu_dose.txt", stringsAsFactors = F)
allele.bantu.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/BANTU/chr1_allele_bantu_dose.txt", stringsAsFactors = F)
cat("position\tref\talt\talt_frq\tunadj_allele_dose_beta\tunadj_allele_dose_beta_se\tunadj_allele_dose_p\tallele_dose_beta\tallele_dose_beta_se\tallele_dose_p\tbantu_dose_beta\tbantu_dose_beta_se\tbantu_dose_p\tallele_bantu_dose_beta\tallele_bantu_dose_beta_se\tallele_bantu_dose_p\tanova_p\n", sep="", file=out.file.name, append=F)

for (position in allele.frame$position) {
  #Merge in the allele dose
  allele.col.frame <- 
    data.frame(sample_id=names(allele.frame)[-c(1:3)],
               allele_dose=t(allele.frame[allele.frame$position == position,-c(1:3)]))
  names(allele.col.frame)[2] <- "allele_dose"
  model.frame <- merge(pheno.frame, allele.col.frame)
  
  #Merge in the bantu dose
  bantu.col.frame <- 
    data.frame(sample_id=names(bantu.frame)[-c(1:3)],
               bantu_dose=t(bantu.frame[bantu.frame$position == position,-c(1:3)]))
  names(bantu.col.frame)[2] <- "bantu_dose"
  model.frame <- merge(model.frame, bantu.col.frame)
  
  #Merge in the allele.bantu dose
  allele.bantu.col.frame <- 
    data.frame(sample_id=names(allele.bantu.frame)[-c(1:3)],
               allele.bantu_dose=t(allele.bantu.frame[allele.bantu.frame$position == position,-c(1:3)]))
  names(allele.bantu.col.frame)[2] <- "allele_bantu_dose"
  model.frame <- merge(model.frame, allele.bantu.col.frame)
  
  #Fit the model
  model <- runLaaaModelSummary(model.frame)
  m.null <- runNullModel(model.frame)
  m.laaa <- runLaaaModel(model.frame)
  anova_p <- anova(m.null, m.laaa, test="Chisq")[2,5]
  allele.model <- runAlleleModelSummary(model.frame)
  
  #Get the ref and alt alleles
  ref <- allele.frame$ref[allele.frame$position == position]
  alt <- allele.frame$alt[allele.frame$position == position]
  
  #Estimate the alternate allele frequency
  frq <- sum(model.frame$allele_dose)/(dim(model.frame)[1]*2)

  #Get the model output
  if ("allele_dose" %in% rownames(model$coefficients)) {
    allele_dose_beta <- model$coefficients["allele_dose", "Estimate"]
    allele_dose_beta_se <- model$coefficients["allele_dose", "Std. Error"]
    allele_dose_p <- model$coefficients["allele_dose", "Pr(>|z|)"]      
  } else {
    allele_dose_beta <- NA
    allele_dose_beta_se <- NA
    allele_dose_p <- NA  
  }
  if ("bantu_dose" %in% rownames(model$coefficients)) {
    bantu_dose_beta <- model$coefficients["bantu_dose", "Estimate"]
    bantu_dose_beta_se <- model$coefficients["bantu_dose", "Std. Error"]
    bantu_dose_p <- model$coefficients["bantu_dose", "Pr(>|z|)"]      
  } else {
    bantu_dose_beta <- NA
    bantu_dose_beta_se <- NA
    bantu_dose_p <- NA  
  }
  if ("allele_bantu_dose" %in% rownames(model$coefficients)) {
    allele_bantu_dose_beta <- model$coefficients["allele_bantu_dose", "Estimate"]
    allele_bantu_dose_beta_se <- model$coefficients["allele_bantu_dose", "Std. Error"]
    allele_bantu_dose_p <- model$coefficients["allele_bantu_dose", "Pr(>|z|)"]      
  } else {
    allele_bantu_dose_beta <- NA
    allele_bantu_dose_beta_se <- NA
    allele_bantu_dose_p <- NA  
  }
  if ("allele_dose" %in% rownames(allele.model$coefficients)) {
    unadj_allele_dose_beta <- allele.model$coefficients["allele_dose", "Estimate"]
    unadj_allele_dose_beta_se <- allele.model$coefficients["allele_dose", "Std. Error"]
    unadj_allele_dose_p <- allele.model$coefficients["allele_dose", "Pr(>|z|)"]      
  } else {
    unadj_allele_dose_beta <- NA
    unadj_allele_dose_beta_se <- NA
    unadj_allele_dose_p <- NA  
  }
  
  #Write the output
  cat(position,"\t", ref, "\t", alt, "\t", frq, "\t",
      unadj_allele_dose_beta,"\t",
      unadj_allele_dose_beta_se,"\t",
      unadj_allele_dose_p,"\t",
      allele_dose_beta,"\t",
      allele_dose_beta_se,"\t",
      allele_dose_p,"\t",
      bantu_dose_beta,"\t",
      bantu_dose_beta_se,"\t",
      bantu_dose_p,"\t",
      allele_bantu_dose_beta,"\t",
      allele_bantu_dose_beta_se,"\t",
      allele_bantu_dose_p, "\t",
      anova_p, "\n", file=out.file.name, append=T, sep="")
  
}

