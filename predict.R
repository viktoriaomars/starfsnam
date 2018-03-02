library("e1071")

info <- read.delim("output.txt")
info$GT = factor(info$GT, levels=c(0,1,2), labels=c("HOMO_REF", "HET", "HOMO_ALT"))

svmfit <- svm(GT ~ REF_AD + ALT_AD + ALT_AD/(ALT_AD + REF_AD), data = info, kernel = "linear", cost = 0.1, scale = FALSE)

table(svmfit$fitted, info$GT)

plot(svmfit, info, GT ~ REF_AD + ALT_AD)

