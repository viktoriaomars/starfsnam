library("e1071")

info <- read.delim("output.txt")
info$GT = as.factor(info$GT)

svmfit <- svm(GT ~ REF_AD + ALT_AD, data = info, kernel = "linear", cost = 0.1, scale = FALSE)

table(svmfit$fitted, info$GT)

plot(svmfit, info, GT ~ REF_AD + ALT_AD)

