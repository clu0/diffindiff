rm(list = ls())

library(xtable)
library(data.table)


filenames = list.files("./results/results_const", pattern="*", full.names=TRUE)

param.names = c("alg", "setup", "n", "p", "constant_eff", "NREP")
setup.values = c('A', 'B', 'C','D')

raw = data.frame(t(sapply(filenames, function(fnm) {

  #output2 = (read.csv(fnm)[,-1])
  output = t(matrix(unlist(read.csv(fnm)[,-1]), ncol=100, byrow=FALSE))
  params = strsplit(fnm, "-")[[1]][2:7]

  err = as.numeric(output[,1])
  mse.mean = mean(err^2)
  bias = mean(-err)

  n = as.numeric(params[3])
  Vhat = as.numeric(output[,2])
  sd.err = sqrt(Vhat/n)
  covg = mean(abs(err)/as.numeric(sd.err) < qnorm(0.975))


  #covg = mean(abs(err)/as.numeric(output[,2]) < qnorm(0.975))


  c(params,
    mse=sprintf("%.8f", round(mse.mean, 8)))
})))


rownames(raw) = 1:nrow(raw)
# names(raw) = c(param.names,
#                "mse","bias","covg")
names(raw) = c(param.names,"mean")

#raw<-raw[raw$learner==learner & (raw$setup %in% setup.values), ] # only look at boost or lasso results

options(stringsAsFactors = FALSE)
raw = data.frame(apply(raw, 1:2, as.character))

raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(as.character(raw$setup)),]
rownames(raw) = 1:nrow(raw)
raw = dcast(setDT(raw), setup + n + p + constant_eff  ~ alg, value.var=c("mean"))

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
write.csv(raw, file=paste0("output_const_full", ".csv"))

# get a dataframe for each setup
#raw.by.setup = lapply(c(setup.values), function(x) raw.round[raw.round$setup==x, ])

# write to latex tables
# for (i in seq_along(setup.values)){
#   tab.setup = cbind("", raw.by.setup[[i+1]][,-1])
#   mse.idx = 1 + c(5:8)
#   for(iter in 1:nrow(tab.setup)) {
#     best.idx = mse.idx[which(as.numeric(tab.setup[iter,mse.idx]) == min(as.numeric(tab.setup[iter,mse.idx])))]
#     for (j in 1:length(best.idx)) {
#       best.idx.j = best.idx[j]
#       tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,best.idx.j])
#     }
#   }
#   tab.setup = tab.setup[,-1]
#   print(setup.values[i])
#   print(tab.setup)
#   # if (learner == "boost") {
#   #   xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean-squared error running \\texttt{boosting} from Setup ", setup.values[i], ". Results are averaged across 200 runs, rounded to two decimal places, and reported on an independent test set of size $n$."), align="ccccccccccc", label=paste0("table:setup",i))
#   # } else if  (learner == "lasso"){
#   #   xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean-squared error running \\texttt{lasso} from Setup ", setup.values[i], ". Results are averaged across 500 runs, rounded to two decimal places, and reported on an independent test set of size $n$."), align="ccccccccccc", label=paste0("table:setup",i))
#   # }
#   xtab.setup = xtable(tab.setup, caption = paste0("\\tt Mean-squared error running \\texttt{lasso} from Setup ", setup.values[i], ". Results are averaged across 500 runs, rounded to two decimal places, and reported on an independent test set of size $n$."), align="ccccccccc", label=paste0("table:setup",i))
#   names(xtab.setup) <- c(c('n','d','$\\sigma$'), algs.tex)
#   print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste("tables", "/simulation_results_setup_", setup.values[i], "_", ".tex", sep=""))
# }
#
#
