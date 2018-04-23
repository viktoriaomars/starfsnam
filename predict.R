library("e1071")
library("kernlab")

ad <- read.delim("output_adtest.txt")
ad$GT <- factor(ad$GT == 0)
#factor(ad$GT, levels=c(0,1,2), labels=c("HOMO_REF", "HET", "HOMO_ALT"))

#zero <- as.numeric(ad$GT == "HOMO_REF")
one <- as.numeric(ad$GT == "HET")
two <- as.numeric(ad$GT == "HOMO_ALT")

for (i in 1:length(one))
{
  if (two[i] == 1)
  {
    one[i] <- 2
  }
}

info <- read.delim("output_info.txt")
total <- merge(ad,info,by="Variant")
#svmfit <- svm(GT ~ ., data = total, kernel = "linear", cost = 0.1, scale = FALSE)
svmfit <- ksvm(GT ~ REF_AD + ALT_AD, data = total, type ="C-svc", kernel = "vanilladot")
#total$GT <- one

table(svmfit$fitted, total$GT)


#plot(svmfit, data = total, GT ~ REF_AD + ALT_AD)
plot(svmfit, data = total)

