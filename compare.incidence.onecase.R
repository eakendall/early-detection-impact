setwd("/Users/emilykendall/Library/CloudStorage/GoogleDrive-emily.kendall@gmail.com/My Drive/DALY impact of ACF 2024/step_change_TB")

load('nointervention_results_emily.rda', verbose = T)
noint.result.base <- noint.result

load('preventedcase_results_emily.rda', verbose = T)
noint.result.onecase <- noint.result

load('preventedcase_symptomatic_results_emily.rda', verbose = T)
noint.result.onecase.symptomatic <- noint.result


z.base <- lapply(noint.result.base,function(x) x$y[,-c(1,26,27,28,29,30,31,32,33,34,35,36,37)])
y.base <- lapply(noint.result.base,function(x) x$y[,])
z.onecase <- lapply(noint.result.onecase,function(x) x$y[,-c(1,26,27,28,29,30,31,32,33,34,35,36,37)])
y.onecase <- lapply(noint.result.onecase,function(x) x$y[,])
z.onecase.symptomatic <- lapply(noint.result.onecase.symptomatic,function(x) x$y[,-c(1,26,27,28,29,30,31,32,33,34,35,36,37)])
y.onecase.symptomatic <- lapply(noint.result.onecase.symptomatic,function(x) x$y[,])

cum.cases.base <- sapply(y.base, function(x) (x[-(1:10),28]+ x[-(1:10),29])/rowSums(x[511,-c(1,26,27,28,29,30,31,32,33,34,35,36,37)])*1000000)
cum.cases.onecase <- sapply(y.onecase, function(x) (x[-(1:10),28]+ x[-(1:10),29])/rowSums(x[511,-c(1,26,27,28,29,30,31,32,33,34,35,36,37)])*1000000)
cum.cases.onecase.symptomatic <- sapply(y.onecase.symptomatic, function(x) (x[-(1:10),28]+ x[-(1:10),29])/rowSums(x[511,-c(1,26,27,28,29,30,31,32,33,34,35,36,37)])*1000000)

cum.cases.difference <- cum.cases.base - cum.cases.onecase
cum.cases.difference.symptomatic <- cum.cases.base - cum.cases.onecase.symptomatic

plot(1:501,cum.cases.difference[,1], type = 'l', ylim=c(-0.1,2))
for (i in 2:1170) lines(1:501,cum.cases.difference[,i])

plot(1:501,cum.cases.difference.symptomatic[,1], type = 'l', ylim=c(-0.1,2))
for (i in 2:1170) lines(1:501,cum.cases.difference.symptomatic[,i])

summary(cum.cases.difference[501,])
quantile(cum.cases.difference[501,], c(0.025, 0.975))

summary(cum.cases.difference.symptomatic[501,])
quantile(cum.cases.difference.symptomatic[501,], c(0.025, 0.975))

# timing of averted cases
summary(apply(cum.cases.difference, 2, function(x) (which.min(x < 0.5*x[501]) - 200)/10))
summary(apply(cum.cases.difference.symptomatic, 2, function(x) (which.min(x < 0.5*x[501]) - 200)/10))
