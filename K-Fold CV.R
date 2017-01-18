###K-fold cross validation loop that was written, however not used in the
###Final project due to time constraints

#SVM-KNN Hybrid
library(e1071)
library(kknn)
library(plyr)
library(splines)

#Classifies expressions as 0 or 1
pred <- read.csv("D:/pred_fin.csv")
exp <- read.csv("D:/expression_fin.csv")
pred.exp <- merge(pred,exp, by="X", all.y = T)
pred.exp$V1 <- as.factor(ifelse(pred.exp$V1 > 0, "1", "0"))
write.csv(pred.exp,file="D:/predexp.csv")

#Creates random subsets
k<-3
pred.exp$id<-sample(1:k,nrow(pred.exp),replace=TRUE)
groups<-1:k
tk<-train.kknn(V1~.,pred.exp,ks=c(1:15),distance=2,scale=FALSE)
best.k<-tk$best.parameters$k

#Euclidean Distance
means<-vector()
accs<-vector()

pb<-create_progress_bar("text")
pb$init(15)

for (i in 1:15){
  for (j in 1:3){
    set.seed(0226)
    train.temp<-subset(pred.exp,id %in% groups[-j])
    test.temp<-subset(pred.exp,id %in% c(j))
    
    train<-train.temp[,-4874]
    test<-test.temp[,-4874]
    
    ge.svm<-svm(V1~.,data=train,scale=FALSE)
    
    ge.knn.train<-train[ge.svm$index,]
    
    preds<-kknn(V1~.,train,test,k=i,kernel="rectangular",distance=2)
    meansk<-mean(preds$fitted.values==test$V1)
    means[j]<-meansk
    acc<-mean(means)
  }
  accs[i]<-acc
  pb$step()
}
#Plots the accuracy of each K after 3 fold cross validation
plot(1:15,accs,xlab="K",ylab="Accuracy",main="Accuracy for Varying Levels of K")
lines(1:15,accs,type="l")

#Chooses the best K
best.acc.pos<-which(accs==max(accs))
best.k<-max(best.acc.pos)

#Predicts the entire data set for the chromosome using the best K found
svm.data<-pred.exp[,-4874]
model.svm<-svm(V1~.,data=svm.data,scale=FALSE)
final.train<-pred.exp[model.svm$index,-4874]
final.test<-pred.exp[-4874]
model.knn<-kknn(V1~.,final.train,final.test,k=best.k,kernel="rectangular",dist=2)
mean(model.knn$fitted.values==pred.exp$V1)

