data <- read.csv("fetal_health.csv")
data <- data[,-6]
n <- floor(0.75 * nrow(data))
set.seed(123)
rows <- sort(sample(seq_len(nrow(data)), size = n))
train.data <- data[rows, ]
test.data <- data[-rows, ]
likelihood <- function(mean, sd, X){
  exponent <- exp(-((X-mean)**2 / (2 * sd**2 )))
  return(log((1 / (sqrt(2 * pi) * sd)) * exponent))
}
gaussian.classifier <- function(train.data, test.data){
  prior = c()
  for(i in 1:3 ){
    prior[i] = sum(train.data$fetal_health==i)/nrow(train.data)
  }
  X1 <- subset(train.data, fetal_health==1)
  X2 <- subset(train.data, fetal_health==2)
  X3 <- subset(train.data, fetal_health==3)
  preds = c()
  for(j in 1:nrow(test.data)){
    prob1 = 0
    prob2 = 0
    prob3 = 0
    for(i in 1:20){
      x <- test.data[j,i]
      prob1 <- prob1 + likelihood(mean(X1[,i]), sd(X1[,i]), x)
      prob2 <- prob2 + likelihood(mean(X2[,i]), sd(X2[,i]), x)
      prob3 <- prob3 + likelihood(mean(X3[,i]), sd(X3[,i]), x)
    }
    probX1 <- log(prior[1]) + prob1
    probX2 <- log(prior[2]) + prob2
    probX3 <- log(prior[3]) + prob3
    temp <- max(probX1, probX2, probX3)
    if(temp == probX1){
      preds <- append(preds, 1)
    }else if(temp == probX2){
      preds <- append(preds, 2)
    }else{
      preds <- append(preds, 3)
    }
  }
  return(preds)
}
accuracy <- function(pred, y){
  correct <- 0
  for(i in 1:length(y)){
    if(pred[i]==y[i]){
      correct  <- correct + 1
    }
  }
  return(correct/length(y))
}
k.fold.cross.vaidation <- function(dataset, folds){
  start <- 1
  end <- folds
  acc <- 0
  count <- 0
  while(end <= nrow(dataset)){
    test.data <- dataset[start:end,]
    train.data <- dataset[-c(start:end),]
    preds <- gaussian.classifier(test.data = test.data, train.data = train.data)
    acc <- acc + accuracy(preds, test.data[,21])
    count <- count + 1
    start <- end + 1
    end <- end + folds
  }
  return(acc/count)
}
best.accuracy <- 0
acc.scores = 0
for(i in seq(from=200, to=600, by=50)){
  acc <- k.fold.cross.vaidation(train.data, i)
  acc.scores = acc.scores + acc
  best.accuracy <- max(best.accuracy, acc)
  print(paste("For k: ",i," accuracy: ", acc*100,"%"))
}
print(paste("Best accuracy score: ", best.accuracy*100))
preds <- gaussian.classifier(train.data, test.data)
print(accuracy(preds, test.data[,21])*100)

for(i in 1:20){
  x1<- subset(train.data, fetal_health==1)
  x2 <- subset(train.data, fetal_health==2)
  x3 <- subset(train.data, fetal_health==3)
  
  dist1 <- dnorm(x1[,i], mean=mean(x1[,i]), sd=sd(x1[,i]))
  dist2 <- dnorm(x2[,i], mean=mean(x2[,i]), sd=sd(x2[,i]))
  dist3 <- dnorm(x3[,i], mean=mean(x3[,i]), sd=sd(x3[,i]))
  
  par(mfrow=c(1, 3))
  
  plot(x1[,i], dist1, col="red")
  plot(x2[,i], dist2, col="blue")
  plot(x3[,i], dist3, col="green")
}


tab <- table(test.data[,21], preds)
true.positive <- function(class){
  return(tab[class, class])
}
true.negative <- function(class){
  temp <- 0
  for(i in 1:3){
    for(j in 1:3){
      if(i!=class && j!=class){
        temp <- temp + tab[i,j]
      }
    }
  }
  return(temp)
}
false.positive <- function(class){
  temp <- 0
  for(i in 1:3){
    if(i != class){
      temp <- temp + tab[i, class]
    }
  }
  return(temp)
}
false.negative <- function(class){
  temp <- 0
  for(i in 1:3){
    if(i != class){
      temp <- temp + tab[class, i]
    }
  }
  return(temp)
}
mat <- matrix(0, nrow = 3, ncol = 4)
for(i in 1:3){
  tp <- true.positive(i)
  tn <- true.negative(i)
  fp <- false.positive(i)
  fn <- false.negative(i)
  mat[i, ] <- c(tp,tn,fp,fn)
}
#accuracy
accFun <- function(class){
  return((mat[class,1]+mat[class,2])/(mat[class,1]+mat[class,2]+mat[class,3]+mat[class,4]))
}
#positive predictive value
precision <- function(class){
  return(mat[class,1]/(mat[class,1]+mat[class,3]))
}

#True positive rate
recall <- function(class){
  return(mat[class,1]/(mat[class,1]+mat[class,4]))
}

f1score <- function(class){
  beta <- 1
  p <- precision(class)
  r <- recall(class)
  return((1+beta^2)*p*r/((beta^2 * p)+r))
}

#false positive rate
fpr <- function(class){
  return(mat[class,3]/(mat[class,3]+mat[class,2]))
}

#False negative rate
fnr <- function(class){
  return(mat[class,4]/(mat[class,1]+mat[class,4]))
}

#True negative rate
specificity <- function(class){
  return(mat[class,2]/(mat[class,2]+mat[class,3]))
}

#Detection Rate
detectionRate <- function(class){
  return(mat[class,1]/(mat[class,1]+mat[class,2]+mat[class,3]+mat[class,4]))
}

#Detection Prevalence
detectionPrevalence <- function(class){
  return((mat[class,1]+mat[class,3])/(mat[class,1]+mat[class,2]+mat[class,3]+mat[class,4]))
}
#balanced accuracy
balanceAccuracy <- function(class)
{
  return((specificity(class)+recall(class))/2)
}
# Result Matrix
diffmat <- matrix(0,nrow = 9,ncol = 3,dimnames = list( c("Precision","Recall","Specificity","False Positive Rate","False Negative rate","f1-score","Detection Rate","Detection Prevalence","Balanced Accuracy"),c("class: 1","class: 2","class: 3")))

for(i in 1:3){
  diffmat[,i] <- c(precision(i),recall(i),specificity(i),fpr(i),fnr(i),f1score(i),detectionRate(i),detectionPrevalence(i),balanceAccuracy(i))
}



