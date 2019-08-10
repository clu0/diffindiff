rm(list = ls())

library(xtable)
library(data.table)


filenames = list.files("./experiments/results/results_non_const", pattern="*", full.names=TRUE)

param.names = c("alg", "setup", "n", "p", "constant_eff", "NREP")

setup.values = c('B','D','E','F')

raw = c()

for (i in 1:length(filenames)) {

  #output2 = (read.csv(fnm)[,-1])
  output = t(matrix(unlist(read.csv(filenames[i])[,-1]), ncol=100, byrow=FALSE))
  params = strsplit(filenames[i], "-")[[1]][2:7]

  err = as.numeric(output[,1])
  mse.mean = sqrt(mean(err^2))
  if (mse.mean >= 100){
    mse.mean = 9999
  }
  bias = mean(-err)

  n = as.numeric(params[3])
  Vhat = as.numeric(output[,2])
  sd.err = sqrt(Vhat/n)
  if (params[1] == "OLS"){
    sd.err = output[,2]
  }
  covg = mean(abs(err)/as.numeric(sd.err) < qnorm(0.975))
  if (is.na(covg)){
    covg = 9999
  }

  mse_row = c(paste(params[1],"mse",sep='.'), params[-1],
              value=sprintf("%.8f", round(mse.mean, 2)))
  bias_row = c(paste(params[1],"bias",sep='.'), params[-1],
              value=sprintf("%.8f", round(bias, 2)))
  covg_row = c(paste(params[1],"covg",sep='.'), params[-1],
              value=sprintf("%.8f", round(covg, 2)))

  raw = rbind(raw,rbind(rbind(mse_row,bias_row),covg_row))
}
raw = data.frame(raw)

rownames(raw) = 1:nrow(raw)
names(raw) = c(param.names,
               "mean")

#raw<-raw[raw$learner==learner & (raw$setup %in% setup.values), ] # only look at boost or lasso results

options(stringsAsFactors = FALSE)
raw = data.frame(apply(raw, 1:2, as.character))
#
# raw = raw[order(as.numeric(raw$p)),]
# raw = raw[order(as.numeric(raw$n)),]
# raw = raw[order(as.character(raw$setup)),]
# rownames(raw) = 1:nrow(raw)
raw = dcast(setDT(raw), setup + n + p + constant_eff ~ alg, value.var=c("mean"))

raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.character(raw$setup)),]
#raw = raw[order(as.character(raw$use_spline)),]
raw = raw[order(as.character(raw$constant_eff)),]
rownames(raw) = 1:nrow(raw)

raw.round = raw
for (col in 6:7){
  raw.round[,col] <- round(as.numeric(unlist(raw[,..col])),2)
}
raw = data.frame(apply(raw, 1:2, as.character))
raw.round = data.frame(apply(raw.round, 1:2, as.character))

# write raw csv output file
write.csv(raw, file=paste0("output_non_const_full_cov", ".csv"))

# write latex tables
raw <- read.csv(file="output_non_const_full_cov_old.csv", header=TRUE, sep=",")
raw = raw[-1]
raw = raw[-1]
#raw = round(raw,3)
tab.all = cbind("", raw)
rmse.idx = c(6, 9, 12,15)
for(iter in 1:nrow(tab.all)) {
  best.idx = rmse.idx[which(as.numeric(tab.all[iter,c(6,9,12,15)]) == min(as.numeric(tab.all[iter,c(6,9,12,15)])))]
  tab.all[iter,best.idx] = paste("\\bf", tab.all[iter,best.idx])
}

xtab.all = xtable(tab.all, align=c("c","c", "|", "|", rep("c",5), "|", rep("c", 3), "|", "|", rep("c", 3), "|", rep("c", 3), "|"))
print(xtab.all, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity,
      hline.after = c(-1, -1, 0, 8, 16, 24, 32, 32), file = "simulation_results.tex")

