source("models.R")
out.file.name <- "/mnt/lustre/users/yswart/Dosages_popgen/BANTU/LAAA_bantu_chr1.txt"
pheno.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/Phenotype.txt", stringsAsFactors = F)
allele.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/BANTU/chr1_allele_dose.txt", stringsAsFactors = F)
bantu.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/BANTU/chr1_bantu_dose.txt", stringsAsFactors = F)
allele.san.frame <- read.delim("/mnt/lustre/users/yswart/Dosages_popgen/BANTU/chr1_allele_bantu_dose.txt", stringsAsFactors = F)
cat("position\tref\talt\talt_frq\tunadj_allele_dose_beta\tunadj_allele_dose_beta_se\tunadj_allele_dose_p\tallele_dose_beta\tallele_dose_beta_se\tallele_dose_p\tbantu_dose_beta\tbantu_dose_beta_se\tbantu_dose_p\tallele_bantu_dose_beta\tallele_bantu_dose_beta_se\tallele_bantu_dose_p\tanova_p\n", sep="", file=out.file.name, append=F)

for (position in allele.frame$position) {
  #Merge in the allele dose
  allele.col.frame <- 
    data.frame(sample_id=names(allele.frame)[-c(1:3)],
               allele_dose=t(allele.frame[allele.frame$position == position,-c(1:3)]))
  names(allele.col.frame)[2] <- "allele_dose"
  model.frame <- merge(pheno.frame, allele.col.frame)
  
  #Merge in the san dose
  san.col.frame <- 
    data.frame(sample_id=names(san.frame)[-c(1:3)],
               san_dose=t(san.frame[san.frame$position == position,-c(1:3)]))
  names(san.col.frame)[2] <- "san_dose"
  model.frame <- merge(model.frame, san.col.frame)
  
  #Merge in the allele.san dose
  allele.san.col.frame <- 
    data.frame(sample_id=names(allele.san.frame)[-c(1:3)],
               allele.san_dose=t(allele.san.frame[allele.san.frame$position == position,-c(1:3)]))
  names(allele.san.col.frame)[2] <- "allele_san_dose"
  model.frame <- merge(model.frame, allele.san.col.frame)
  
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
  if ("san_dose" %in% rownames(model$coefficients)) {
    san_dose_beta <- model$coefficients["san_dose", "Estimate"]
    san_dose_beta_se <- model$coefficients["san_dose", "Std. Error"]
    san_dose_p <- model$coefficients["san_dose", "Pr(>|z|)"]      
  } else {
    san_dose_beta <- NA
    san_dose_beta_se <- NA
    san_dose_p <- NA  
  }
  if ("allele_san_dose" %in% rownames(model$coefficients)) {
    allele_san_dose_beta <- model$coefficients["allele_san_dose", "Estimate"]
    allele_san_dose_beta_se <- model$coefficients["allele_san_dose", "Std. Error"]
    allele_san_dose_p <- model$coefficients["allele_san_dose", "Pr(>|z|)"]      
  } else {
    allele_san_dose_beta <- NA
    allele_san_dose_beta_se <- NA
    allele_san_dose_p <- NA  
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
      san_dose_beta,"\t",
      san_dose_beta_se,"\t",
      san_dose_p,"\t",
      allele_san_dose_beta,"\t",
      allele_san_dose_beta_se,"\t",
      allele_san_dose_p, "\t",
      anova_p, "\n", file=out.file.name, append=T, sep="")
  
}

