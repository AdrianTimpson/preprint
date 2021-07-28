#---------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
#---------------------------------------------------------------------------------------------------------------------------------------------
# the following packages need to be installed:
citation('arrangements')
citation('PreciseSums')
citation('extraDistr')
#---------------------------------------------------------------------------------------------------------------------------------------------
presence.probability <- function(draws,probs){

	# Calculates the overall probability of S observed species for presence-only data using fast algorithm (method.2)
	# Draws: integer number of draws from the urn
	# probs: vector of S species frequencies and missingness (last value). Must sum to 1.
	# example 1: 4 draws from 2 species of equal frequency, with no missingness: presence.probability(4,c(0.5,0.5,0))
	# example 2: 4 draws from 2 species with a proportion of 5:1, with 10% missingness: presence.probability(4,c(0.75,0.15,0.1))
	
	require(arrangements)
	require(PreciseSums)

	# floating point issue: ensure probs sum to 1
	while(sum(probs)!=1)probs <- probs/sum(probs)

	# any extreme conditions
	if(draws==0)return(0)
	if(is.infinite(draws))return(1)

	# expression containing N sets of values, such that the 1st will be included, the 2nd excluded, etc
	N <- length(probs)
	expr <- list()
	if(N==1){
		k1 <- 1 - probs 
		res <- kahanSum(k1^draws) # more precise version of sum()
		}
	if(N>1){
		for(n in 1:N){
			k1 <- 1 - rowSums(cbind(combinations(v=probs[1:(N-1)],k=n-1),probs[N]))
			expr[[n]] <- (k1^draws)
			}

		# floating point issue: to prevent floating point catastrophic cancellation,
		# pair up all inclusions and exclusions by similar order of magnitudes
		include <- sort(unlist(expr[seq(1,N,by=2)]),decreasing=T)
		exclude <- sort(unlist(expr[seq(2,N,by=2)]),decreasing=T)
		remaining <- include-exclude
		res <- kahanSum(remaining)
		}

	# safeguard incase of exceptional floating point problem
	if(res<0)res <- 0	

return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------
presence.probability.slow <- function(draws,probs){
	
	# performs the same job as presence.probability() using the intuitive algorithm (method.1) that considers all compositions
	# example 1: 4 draws from 2 species of equal frequency, with no missingness: presence.probability.slow(4,c(0.5,0.5,0))
	# example 2: 4 draws from 2 species with a proportion of 5:1, with 10% missingness: presence.probability.slow(4,c(0.75,0.15,0.1))

	require(arrangements)

	# all compositions that satisfy the presence-only observations
	x <- cbind(compositions(draws,length(probs)-1),0)

	# probability of each possible set, summed
	res <- sum(apply(x,1,dmultinom, prob=probs))

return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------
presence.lik <- function(draws,prior.ratios,presence,missing){

	# Wrapper function
	# Can handle presence-only data (all presence values = T)
	# Can also handle presence-absence data (some presence values = F), which are become additional missingness
	# example 1: 4 draws from presence-only data of 2 species of equal frequency, with no missingness: presence.lik(4,c(1,1),c(T,T),0)
	# example 2: 4 draws from presence-only data of 2 species with a proportion of 5:1, with 10% missingness: presence.lik(4,c(1,5),c(T,T),0.1)
	# example 3: 4 draws from presence-absence data of 2 species of equal frequency, one present and one absent, with no missingness: presence.lik(4,c(1,1),c(T,F),0)
	# example 4: 4 draws from presence-absence data of 2 species with a proportion of 5:1, first present and second absent, with 10% missingness: presence.lik(4,c(5,1),c(T,F),0.1)

	if(draws<sum(presence))return(0)
	if(draws==0)return(1)
	if(draws>0 & sum(presence)==0)return(0)

	tmp <- c((prior.ratios[presence]/sum(prior.ratios)),sum(prior.ratios[!presence])/sum(prior.ratios))

	# vector of probabilities, including the missing
	probs <- tmp*(1-missing)
	probs[length(probs)] <- probs[length(probs)] + missing

	res <- presence.probability(draws=draws,probs=probs)
return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------
mcmc <- function(data, start.pars, jump, objective.function, proposal.function, type, period, frequency.col, island, draw.priors, N){ 

	all.pars <- matrix(,N,length(start.pars))
	all.probs <- numeric(N)

	# mcmc loop
	for(n in 1:N){
		if(n==1){
			pars <- start.pars
			lprob <- -objective.function(pars, data, type, period, frequency.col, island, draw.priors)
			}

		all.pars[n,] <- pars
		all.probs[n] <- lprob

		prop.pars <- proposal.function(pars, jump)
		prop.lprob <- -objective.function(prop.pars, data, type, period, frequency.col, island, draw.priors)

		ratio <- min(exp(prop.lprob-lprob),1)

		if(is.infinite(lprob))ratio <- 1 # always move if the current pars are crap
		if(is.infinite(prop.lprob))ratio <- 0 # never move if the proposed pars are crap
		move <- sample(c(T,F),size=1,prob=c(ratio,1-ratio))
		if(move){
			pars <- prop.pars
			lprob <- prop.lprob
			}
		if(n%in%seq(0,N,length.out=101))print(c(n,'iterations of',N, pars, lprob ))
		}

	# acceptance ratio
	tmp <- rowSums(apply(all.pars, 2, diff))
	ar <- 1-sum(tmp==0)/length(tmp)

return(list(all.pars=all.pars, acceptance.ratio=ar, all.probs=all.probs))}
#---------------------------------------------------------------------------------------------------------------------------------------------
proposal.function <- function(pars, jump){

	# exclude the impossible or ridiculous
	min <- c(rep(1,length(pars)-1),0)
	max <- c(rep(50000,length(pars)-1),1)

	new.pars <- rnorm(length(pars),pars,jump)
	i <- new.pars>1
	new.pars[i] <- round(new.pars[i])
	cond <- sum(new.pars<min) + sum(new.pars>max)
	if(cond>0)new.pars <- pars
return(new.pars)}
#---------------------------------------------------------------------------------------------------------------------------------------------
objective.function <- function(pars, data, type, period, frequency.col, island, draw.priors){
	
	require(extraDistr)
	P <- format.pars(pars,type)
	D <- format.data(data,type,period,frequency.col,island)
	W <- format.draw.priors(draw.priors,type,period,island)

	if(grepl('null',type)){
		lik.malta <- presence.lik(draws=P$malta$draws, prior.ratios=D$malta$ratio, presence=D$malta$presence, missing=P$missing)
		lik.sicily <- presence.lik(draws=P$sicily$draws, prior.ratios=D$sicily$ratio, presence=D$sicily$presence, missing=P$missing)
		loglik <-  log(lik.malta) + log(lik.sicily)
		prior.log.prob.malta <- ddgamma(P$malta$draws, W$malta$shape, W$malta$rate,log=TRUE)
		prior.log.prob.sicily <- ddgamma(P$sicily$draws, W$sicily$shape, W$sicily$rate,log=TRUE)
		prior.log.prob <- prior.log.prob.malta + prior.log.prob.sicily
		}
	if(grepl('inde',type)){
		lik <- presence.lik(draws=P$draws, prior.ratios=D$ratio, presence=D$presence, missing=P$missing)
		loglik <- log(lik)
		prior.log.prob <- ddgamma(P$draws, W$shape, W$rate,log=TRUE)
		}

	res <- (prior.log.prob + loglik) * (-1)
return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------
extract.priors <- function(data,period,frequency.col){

	cols <- grep(period,names(data))
	rows <- rowSums(is.na(data[,cols,drop=F]))!=2
	PA <- !is.na(data[rows,cols])
	colnames(PA) <- unlist(strsplit(colnames(PA),split='.',fixed=T))[c(2,4)]
	freq <- data[rows,frequency.col,drop=F]
	colnames(freq) <- 'freq'
	sub <- cbind(freq,PA)
return(sub)}
#---------------------------------------------------------------------------------------------------------------------------------------------
format.data <- function(data,type,period,frequency.col,island){
	
	e <- extract.priors(data,period,frequency.col)
	if(grepl('null',type)){
		malta <- data.frame(ratio=e$freq,presence=e$malta)
		sicily <- data.frame(ratio=e$freq,presence=e$sicily)
		res <- list(malta=malta,sicily=sicily)
		}
	if(grepl('inde',type)){
		if(is.null(island))stop('if type = inde, island must be sicily or malta')
		ratio <- e$freq[e[,island]]
		presence <- rep(TRUE,length(ratio))
		res <- list(ratio=ratio,presence=presence)
		}
return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------
format.pars <- function(pars,type){

	if(grepl('null',type)){
		malta <- data.frame(draws =pars[1])
		sicily <- data.frame(draws = pars[2])
		missing <- pars[3]
		res <- list(malta=malta,sicily=sicily, missing=missing)
		}

	if(grepl('inde',type)){	
		draws <- pars[1]
		missing <- pars[2]
		res <- list(draws=draws, missing=missing)
		}

return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------
format.draw.priors <- function(draw.priors,type,period,island){

	if(grepl('null',type)){
		malta <- as.data.frame(t(draw.priors[,paste(period,'malta',sep='.'),drop=F]))
		sicily <- as.data.frame(t(draw.priors[,paste(period,'sicily',sep='.'),drop=F]))
		res <- list(malta=malta,sicily=sicily)
		}
	if(grepl('inde',type)){
		if(is.null(island))stop('if type = inde, island must be sicily or malta')
		res <- as.data.frame(t(draw.priors[,paste(period,island,sep='.'),drop=F]))
		}
return(res)}
#---------------------------------------------------------------------------------------------------------------------------------------------


