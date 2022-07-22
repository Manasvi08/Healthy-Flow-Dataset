hfd<- read.csv("hfd_sub.csv")

gates<- hfd[,5]

hfd[,1:4]<-scale(hfd[,1:4])  

size.train<- floor(0.5*nrow(hfd))

size.test<-floor(0.25*nrow(hfd))

train_lab<- hfd[1:size.train,5]

test_lab<- hfd[(size.train+1):(size.train+size.test),5]

valid_lab<- hfd[(size.train+size.test+1):3864,5]

train_data<-hfd[1:size.train,1:4]

test_data<- hfd[(size.train+1):(size.train+size.test),1:4]

valid_data<- hfd[(size.train+size.test+1):3864,1:4]

library(class)
result<- knn(train_data,test_data,train_lab,4)
result


class_agree<-table(result,test_lab)
class_agree
sum_agree<-sum(diag(class_agree))
sum_agree
(nrow(test_data)-sum_agree)/nrow(test_data)


kmax <- 50
k <- 1:kmax
p <- rep(0, kmax)
ntest <- nrow(test_data)
k_summary <- cbind(k, p)
colnames(k_summary) <- c("k","% misclassified")
for(i in 1:kmax){
  result <- knn(train_data, test_data, cl =train_lab, k = i)
  class_agree <- table(result, test_lab)
  sum_agree <- sum(diag(class_agree))
  k_summary[i, 2] <- (ntest - sum_agree) / ntest
}
k_summary[1:10, ]

kmin<-min(k_summary)#k=32

plot(k_summary, type = "l",main="Missclassification Plot")

#validation 
result_validation<- knn(train_data,valid_data,train_lab,32)
class_agree1 <- table(result_validation, valid_lab)
class_agree1
sum_agree1 <- sum(diag(class_agree1))
sum_agree1
missfit <- (ntest - sum_agree1) / ntest
missfit*100



library(MASS)
plot(hfd[,-5], col = as.factor(hfd[,5]))


size.train.qda<- floor(0.8*nrow(hfd.norm))
size.test.qda<-floor(0.2*nrow(hfd.norm))
train_lab.qda<- hfd[1:size.train.qda,5]
test_lab.qda<- hfd[(size.train.qda+1):(size.train.qda+size.test.qda),5]
train_data.qda<-hfd[1:size.train.qda,1:4]
test_data.qda<- hfd[(size.train.qda+1):(size.train.qda+size.test.qda),1:4]

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

hfd.norm <- as.data.frame(lapply(hfd, min_max_norm))

qsol <- qda(train_data.qda, grouping = train_lab.qda)
qsol$prior
qsol$means
qsol$scaling
predict(qsol, test_data.qda)
qsol

      


