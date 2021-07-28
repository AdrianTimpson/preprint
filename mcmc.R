#---------------------------------------------------------------------------------------------------------------------
# preferred prior beliefs require the argument frequency.col='alexandra'
# additional sensitivity tests performed with the argument frequency.col='emily'
source('functions.R')
data <- read.csv('data with frequency priors.csv',na.strings='')
draw.priors <- read.csv('draw priors.csv',row.names=1)
N <- 100 # N <- 100000
#---------------------------------------------------------------------------------------------------------------------
# null hypothesis
#---------------------------------------------------------------------------------------------------------------------
null.EMP <- mcmc(data, start.pars=c(70, 250, 0.0032), jump=c(13, 10, 0.0018), objective.function, proposal.function, type='null', period='EMP', frequency.col='alexandra', island=NULL, draw.priors, N)
null.LMP <- mcmc(data, start.pars=c(108, 28, 0.008) , jump=c(22, 5, 0.0056) , objective.function, proposal.function, type='null', period='LMP', frequency.col='alexandra', island=NULL, draw.priors, N)
null.LP <- mcmc(data, start.pars=c(53, 119, 0.0136) , jump=c(7.64, 2.51, 0.0058) , objective.function, proposal.function, type='null', period='LP', frequency.col='alexandra', island=NULL, draw.priors, N)
#---------------------------------------------------------------------------------------------------------------------
# independent hypothesis
#---------------------------------------------------------------------------------------------------------------------
null.EMP <- mcmc(data, start.pars=c(70, 250, 0.0032), jump=c(13, 10, 0.0018), objective.function, proposal.function, type='null', period='EMP', frequency.col='alexandra', island=NULL, draw.priors, N)

#---------------------------------------------------------------------------------------------------------------------

citation('arrangements')
citation('PreciseSums')
citation('extraDistr')