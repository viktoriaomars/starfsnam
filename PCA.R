train1 <- read.delim("output_adtrain.txt")
test1 <- read.delim("output_adtest.txt")
info <- read.delim("output_info.txt")
train <- merge(train1,info,by="Variant")
test <- merge(test1,info,by="Variant")

test$GT <- 1

combi <- rbind(train, test)

#combi$REF_AD <- ifelse(combi$REF_AD == 0, median(combi$REF_AD), combi$REF_AD)
#combi$ALT_AD <- ifelse(combi$ALT_AD == 0, median(combi$ALT_AD), combi$ALT_AD)

#Find columns with no variance, wan to remove that
no_var <- which(apply(combi, 2, var)==0)

my_data <- subset(combi, select = -c(GT, Sample, Variant, no_var)) #Hafa GT hér eğa ekki, virğist ekki breyta neinu

colnames(my_data)
str(my_data)

#divide data
pca.train <- my_data[1:nrow(train),]
pca.test <- my_data[-(1:nrow(train)),]

#principal component analysis
prin_comp <- prcomp(pca.train, scale. = T)
names(prin_comp)

#outputs the mean of variables
prin_comp$center

#outputs the standard deviation of variables
prin_comp$scale

#outputs prinicipal component loading
prin_comp$rotation

#plot resultant principal components
biplot(prin_comp, scale = 0)

#compute standard deviation of each principal component
std_dev <- prin_comp$sdev

#compute variance
pr_var <- std_dev^2

#check variance of first 10 components
pr_var

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex

#scree plot
plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained",type = "b")

#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", type = "b")



#add a training set with principal components
train.data <- data.frame(GT = train$GT, prin_comp$x)

#we are interested in first n PCAs, the ones that account for 95% of the variance
n = 0
cumulative <- cumsum(prop_varex)
for(val in cumulative) {
  if(val <= 0.95){
    n = n + 1
  }
}
train.data <- train.data[,1:n]

#run a decision tree
pca.train$GT <- combi$GT[1:nrow(train)] #Hafa GT meğ svona? ss ekki í PCA sjálfu??? virğist ekki breyta neinu
pca.train$GT <- as.factor(pca.train$GT)

library(tree)
tree.model <- tree(GT ~ .,data = pca.train)
tree.model <- prune.tree(tree.model, best=3)
tree.model

plot(tree.model)
text(tree.model)

#transform test into PCA
#test.data <- predict(prin_comp, newdata = pca.test)
#test.data <- as.data.frame(test.data)

#select the first n components, they account for 95% of the variance
#test.data <- test.data[,1:n]

#make prediction on test data
tree.prediction <- predict(tree.model, pca.test)

test$GT <- tree.prediction

write.csv(test, "pca.csv", row.names = F)

