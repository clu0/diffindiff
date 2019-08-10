library(Matrix)
library(foreach)
library(magrittr)
library(glmnet)
library(diffindiff)
library(EQL)
library(rpart)
library(foreign)

profits = read.dta("./data/drugsdata/Table6/profits.dta")

profits.mat = as.matrix(profits)
# last 16 are group indicators
X = matrix(NA, dim(profits)[1], 9+16)
cov_index = c(1,3,35,36,79)
X[,1:5] = profits.mat[,cov_index]

for (i in 1:39){
  profits.mat[,i+37] = profits.mat[,i+37]*(20+i)
}

X[,6] = rowSums(profits.mat[,38:76])
X[X[,6]==0,6] = NA

# 26, 28 for 1993 and 1995, 31, 33 for 1998 and 2000
# X_real is data in 93 vs 95, whereas X_control is 98 vs 00
X_real = X
X_real[profits.mat[,26]==1,7] = 0
X_real[profits.mat[,28]==1,7] = 1
growing = c(8,23,24)
nongrowing = c(7,9,10,12,13,15,16,17,18,19,21,22,25)
X_real[rowSums(profits.mat[,growing])==1,8] = 1
X_real[rowSums(profits.mat[,nongrowing])==1,8] = 0
X_real[,10:12] = profits.mat[,growing]
X_real[,13:25] = profits.mat[,nongrowing]
X_real[,9] = profits.mat[,78]

# throw away null values because otherwise file is too large
X_real = X_real[complete.cases(X_real),]

X_real = X_real[sample(1:dim(X_real)[1],dim(X_real)[1]),]

write.csv(X_real,"./experiments/ak07data_93_95.csv")


X_control = X
X_control[profits.mat[,31]==1,7] = 0
X_control[profits.mat[,33]==1,7] = 1
X_control[rowSums(profits.mat[,growing])==1,8] = 1
X_control[rowSums(profits.mat[,nongrowing])==1,8] = 0
X_control[,10:12] = profits.mat[,growing]
X_control[,13:25] = profits.mat[,nongrowing]
X_control[,9] = profits.mat[,78]

X_control = X_control[complete.cases(X_control),]

X_control = X_control[sample(1:dim(X_control)[1],dim(X_control)[1]),]

write.csv(X_control,"./experiments/ak07data_98_00.csv")

# also have files for testing urban vs rural in 93 vs 95
X_urban = X_real[(X_real[,5]==0),]
X_rural = X_real[(X_real[,5]==1),]
write.csv(X_urban,"./experiments/ak07data_urban.csv")
write.csv(X_rural,"./experiments/ak07data_rural.csv")
