## Plot a Manhattan and QQ-plot

# Run these command on RStudio
install.packages("qqman")
library(qqman)
head(gwasResults)

manhattan(gwasResults)
qq(gwasResults$P)
