#' Simulates data from Barker capture-recapture model 
#'
#' Simulates capture-recapture data from Burnham (both)  model of releases of multiple cohorts.
#' Allows general model for s, p, r and F or constant or time-specific model
#' depending on input arguments. Optionally outputs file to outfile with .inp extension for input to program MARK.
#'  
#' The design.data dataframe contains a record for each animal for each time period it is contained in the study.  The default
#' fieldnames are id(unique # in cohort), cohort (cohort number from 1:num.cohorts), time (i:num.cohorts,
#' where i is cohort #) and age (assumed to increment by 1 each occasion and defaults with 0 initial age)
#' The function create.design.data can be used to create an initial dataframe which can be supplemented 
#' with other covariates. If a formula is given but no design.data is not provided this function
#' 
#' The design.data dataframe contains a record for each animal for each time period it is contained in the study.  The default
#' fieldnames are id(unique # in cohort), cohort (cohort number from 1:num.cohorts), time (i:num.cohorts,
#' where i is cohort #) and age (assumed to increment by 1 each occasion and defaults with 0 initial age)
#' The function create.design.data can be used to create an initial dataframe which can be supplemented 
#' with other covariates. If a formula is given but no design.data is not provided this function
#' calls create.design.data to construct a default dataframe for the specified problem.
#'
#' For designation of models for s or p, there are 3 options here:
#'         constant model: s=list(par=value) 
#'         time model:     s=list(par=c(val1,val2,...valk)) where k=num.cohorts-1 is number of survival intervals
#'         general model:  s=list(par=c(val1,val2,...valk),formula=~yourformula)) k is number of cols in model matrix
#' See \code{\link{simcjs}} for more details.
#' 
#' @param num.cohorts  number of cohorts; design is square with same number of c-r eventsas num.cohorts; number of recapture events is num.cohorts-1
#' @param cohort.sizes a scalar giving constant size of each cohort or a vector of sizes of length num.cohorts
#' @param s a list defining the survival model with the following elements (see details)
#'                    par      - a vector of parameter values
#'                    formula  - a formula to use with design.data to construct model
#'                    link     - link function used with model to create probabilites
#' @param p a list defining the capture probability model (same structure as s)
#' @param r a list defining the recovery probability model (same structure as s)
#' @param F a list defining the model of fidelity to capture region (same structure as s)
#' @param design.data a dataframe with design data that allows model construction for probabilities (see details).
#' @param outfile prefix name of the output file for the ch data. extension .inp is always added for MARK
#' @export simburnham
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @return dataframe with ch (capture history) as character
#' @examples
#'do.burnham <- function(s,p,r,F,reps,output=FALSE)
#'#
#'#  do.burnham -  a simple example piece of code to show simulation/analysis loop
#'#                with a Burnham model
#'#
#'#  Arguments:
#'#
#'#  s   -  parameter list for s
#'#  p   -  parameter list for p
#'#  r   -  parameter list for r
#'#  F   -  parameter list for F
#'#  reps-  number of simulation reps
#'#
#'#  Value:
#'# 
#'#  results - for this simple example it will be a matrix of the real parameter estimates
#'#
#'#  Functions used: simburnham, mark
#'#
#'#
#'{
#'	results=NULL
#'	for(i in 1:reps)
#'	{
#'		cat("rep ",i,"\n")
#'		xx=simburnham(3,500,s=s,p=p,r=r,F=F)
#'		mod<-mark(xx,title="Burnham test",model="Burnham",output=output)
#'		results=rbind(results,mod$results$real$estimate)
#'	}
#'	return(results)
#'}
#' 
simburnham <- function(num.cohorts,cohort.sizes,s=list(),p=list(),r=list(),F=list(),design.data=NULL,outfile=NULL)
{
#
# Setup default values for link and allow simple parameter vectors
#
	if(!is.list(s))
		s=list(par=s,link="identity")
	if(!is.list(r))
		r=list(par=r,link="identity")
	if(!is.list(p))
		p=list(par=p,link="identity")
	if(!is.list(F))
		F=list(par=F,link="identity")
	if(is.null(s$link))s$link="logit"
	if(is.null(p$link))p$link="logit"
	if(is.null(r$link))r$link="logit"
	if(is.null(F$link))F$link="logit"
#
# Check to make sure that s$par and p$par are defined
#
	if(is.null(s$par)) stop("Survival parameters are missing\n")
	if(is.null(p$par)) stop("Capture parameters are missing\n")
	if(is.null(r$par)) stop("Recovery parameters are missing\n")
	if(is.null(F$par)) stop("Fidelity parameters are missing\n")
#       
# Setup cohorts and cohort sizes
#
	ch=NULL
	if(num.cohorts>1 & length(cohort.sizes)==1)
		cohort.sizes=rep(cohort.sizes,num.cohorts)
	else
	if(num.cohorts !=length(cohort.sizes))
		stop("Number of cohorts doesn't match vector of cohort sizes\n")
#
# Check if p,s are using formula and then make sure design.data exists;
# if it doesn't exist create default design data
#
	if(!is.null(s$formula) | !is.null(p$formula) | !is.null(r$formula) | !is.null(F$formula))
		if(is.null(design.data))
			design.data=create.design.data(num.cohorts,cohort.sizes,initial.age=0)
#
# Set up survival parameters
#
	if(is.null(s$formula))
	{
		sformula=FALSE
		if(length(s$par)==1)
			sconstant=TRUE
		else
		if(length(s$par)==num.cohorts)
			sconstant=FALSE
		else
			stop("Incorrect number of survival parameters; doesn't match number of occasions\n")
		if(s$link=="identity")   
			s=s$par
		else
			s=exp(s$par)/(1+exp(s$par))
	}
	else
	{
		if(s$link!="logit")cat("Note: Logit link will be used with formula\n")
		sformula=TRUE
		smat=model.matrix(s$formula,design.data)
		if(dim(smat)[2]!=length(s$par))
			stop(paste("Dimension of model matrix",dim(smat)[2],"is not consistent with length of survival parameter vector",length(s$par),"\n"))
	}
#
# Set up capture parameters
#
	if(is.null(p$formula))
	{
		pformula=FALSE
		if(length(p$par)==1)
			pconstant=TRUE
		else
		if(length(p$par)==(num.cohorts-1))
			pconstant=FALSE
		else
			stop("Incorrect number of capture parameters; doesn't match number of occasions\n")
		if(p$link=="identity")   
			p=p$par
		else
			p=exp(p$par)/(1+exp(p$par))
	}
	else
	{
		if(p$link!="logit")cat("Note: Logit link will be used with formula\n")
		pformula=TRUE
		pmat=model.matrix(p$formula,design.data[design.data$cohort<num.cohorts&design.data$time>design.data$cohort,])
		if(dim(pmat)[2]!=length(p$par))
			stop(paste("Dimension of model matrix",dim(pmat)[2],"is not consistent with length of capture parameter vector",length(p$par),"\n"))
	}
#
# Set up recovery parameters
#
	if(is.null(r$formula))
	{
		rformula=FALSE
		if(length(r$par)==1)
			rconstant=TRUE
		else
		if(length(r$par)==num.cohorts)
			rconstant=FALSE
		else
			stop("Incorrect number of recovery parameters; doesn't match number of occasions\n")
		if(r$link=="identity")   
			r=r$par
		else
			r=exp(r$par)/(1+exp(r$par))
	}
	else
	{
		if(r$link!="logit")cat("Note: Logit link will be used with formula\n")
		rformula=TRUE
		rmat=model.matrix(r$formula,design.data)
		if(dim(rmat)[2]!=length(r$par))
			stop(paste("Dimension of model matrix",dim(rmat)[2],"is not consistent with length of recovery parameter vector",length(r$par),"\n"))
	}
#
# Set up fidelity parameters
#
	if(is.null(F$formula))
	{
		Fformula=FALSE
		if(length(F$par)==1)
			Fconstant=TRUE
		else
		if(length(F$par)==(num.cohorts-1))
			Fconstant=FALSE
		else
			stop("Incorrect number of fidelity parameters; doesn't match number of occasions\n")
		if(F$link=="identity")   
			F=F$par
		else
			F=exp(F$par)/(1+exp(F$par))
	}
	else
	{
		if(F$link!="logit")cat("Note: Logit link will be used with formula\n")
		Fformula=TRUE
		Fmat=model.matrix(F$formula,design.data[design.data$cohort<num.cohorts&design.data$time<num.cohorts,])
		if(dim(Fmat)[2]!=length(F$par))
			stop(paste("Dimension of model matrix",dim(Fmat)[2],"is not consistent with length of fidelity parameter vector",length(F$par),"\n"))
	}
#
# Loop over cohorts and generate simulated data for each cohort
#
	for (i in 1:(num.cohorts-1))
	{
#
#  Setup cohort specific survival parameters
#
		if(!sformula)
		{
			if(!sconstant)  
				svalues=s[i:num.cohorts]
			else
				svalues=s
		} 
		else
		{
			svalues=smat[design.data$cohort==i,] %*% s$par
			svalues=exp(svalues)/(1+exp(svalues))
			select=design.data$cohort==i&design.data$time>=i
			svalues=tapply(as.vector(svalues),list(as.factor(design.data$id[select]),as.factor(design.data$time[select])),sum)
			if(dim(svalues)[1]!=cohort.sizes[i]) 
				stop("Number of rows in s model matrix not consistent with number in cohort\n")
		}
#
#  Setup cohort specific capture parameters
#
		if(!pformula)
		{
			if(!pconstant)
				pvalues=p[i:(num.cohorts-1)]       
			else
				pvalues=p
		} 
		else
		{
			pvalues=pmat[design.data$cohort[design.data$time>design.data$cohort&design.data$cohort<num.cohorts]==i,] %*% p$par
			pvalues=exp(pvalues)/(1+exp(pvalues))
			select=design.data$cohort==i&design.data$time>i
			pvalues=tapply(as.vector(pvalues),list(as.factor(design.data$id[select]),as.factor(design.data$time[select])),sum)
			if(dim(pvalues)[1]!=cohort.sizes[i]) 
				stop("Number of rows in p model matrix not consistent with number in cohort\n")
		}
#
#  Setup cohort specific fidelity parameters
#
		if(!Fformula)
		{
			if(!Fconstant)  
				Fvalues=F[i:(num.cohorts-1)]
			else
				Fvalues=F
		} 
		else
		{
			Fvalues=Fmat[design.data$cohort[design.data$time<num.cohorts]==i,] %*% F$par
			Fvalues=exp(Fvalues)/(1+exp(Fvalues))
			select=design.data$cohort==i&design.data$time>=i&design.data$time<num.cohorts
			Fvalues=tapply(as.vector(Fvalues),list(as.factor(design.data$id[select]),as.factor(design.data$time[select])),sum)
			if(dim(Fvalues)[1]!=cohort.sizes[i]) 
				stop("Number of rows in F model matrix not consistent with number in cohort\n")
		}
#
#  Setup cohort specific recovery parameters
#
		if(!rformula)
		{
			if(!rconstant)  
				rvalues=r[i:num.cohorts]
			else
				rvalues=r
		} 
		else
		{
			rvalues=rmat[design.data$cohort==i,] %*% r$par
			rvalues=exp(rvalues)/(1+exp(rvalues))
			select=design.data$cohort==i
			rvalues=tapply(as.vector(rvalues),list(as.factor(design.data$id[select]),as.factor(design.data$time[select])),sum)
			if(dim(rvalues)[1]!=cohort.sizes[i]) 
				stop("Number of rows in r model matrix not consistent with number in cohort\n")
		}
#
#  Create simulated capture histories for the cohort and append to ch; if this is 
#  not the first cohort, the appropriate number of 0's are used as prefix to the
#  cohort capture history.
#
		ch=c(ch,paste(rep(paste(rep("00",i-1),collapse=""),cohort.sizes[i]),simburnham.cohort(cohort.sizes[i],svalues,pvalues,rvalues,Fvalues,num.cohorts-i+1),sep=""))   
	}
#
#  If a filename specified output ch to file for input to MARK
#
	if(!is.null(outfile))
		write.table(paste(ch,rep(" 1 ;",length(ch))),paste(outfile,".inp",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
#
#  Return ch vector as value 
#
	return(data=data.frame(ch=I(ch)))
}
"simburnham.cohort" <-
		function(N,s,p,r,F,nocc)
#
# simburnham - simulates capture histories for a single cohort of N release animals 
#             for nocc occasions using the Burnham model 
#
# Arguments: 
#
#  N        - size of cohort
#  s        - survival probability (can be either scalar(constant model), 
#             vector (time model - length = nocc-1 or heterogeneity model - length=N)
#             or matrix (N by nocc-1)
#  p        - capture probability (same types as survival)
#  r        - recovery probability (same types as survival)
#  F        - probability of at risk animal staying at risk (same types as survival)
#  nocc     - number of occasions (ch length is 2*nocc)
#
# Value:
#  ch       - character vector of capture histories 
#
{
#
# Setup parameter matrix depending on what was input; it can be a scalar (constant model), a vector with
# nocc-1 values (time model) or N values (heterogeneity) or a matrix (N rows by nocc-1 cols) for a completely 
# general specification.  In each case it is transformed to a matrix.
#
	s=create.parmat(s,nocc+1,N)
	r=create.parmat(r,nocc+1,N)
	p=create.parmat(p,nocc,N)
	F=create.parmat(F,nocc,N)
#
# Set up as release; ch is the vector of capture histories
#
	ch=rep("1",N)
#
# alive keeps track of whether animals are alive = 1; at release all are alive
#
	alive=rep(1,N)
#
# atrisk keeps track of whether animals are in capture area (at risk of capture); at release
# all are at risk
#
	atrisk=rep(1,N)
#
# Loop over each occasion. First release is first occasion
#
	for (i in 1:nocc)
	{
#
#  Generate Bernoulli rv for survival (sx) and recovery (rx)
#
		sx=rbinom(N,1,s[,i])
		rx=rbinom(N,1,r[,i])
#
#  Compute value of history between capture occasions for recovery of dead animal
#
		ch=paste(ch,as.character(alive*(1-sx)*rx),sep="")        
#
#  Only continue on if this is not the last occasion
#
		if(i!=nocc)
		{
# 
#     Determine which animals are at risk of capture and those captured
#       Fx      - bernoulli rv for animal at risk to stay at risk
#       px      - bernoulli rv for capture
#
			Fx=rbinom(N,1,F[,i])
			px=rbinom(N,1,p[,i])
#
#     Compute those currently at risk 
#
			atrisk=atrisk*Fx  
#
#     Update who is alive 
#
			alive=alive*sx
#
#     Update ch with those captured on the occasion    
#
			ch=paste(ch,as.character(alive*px*atrisk),sep="")        
		}
	}
#
#    Return a vector of capture histories with a 1 ; appended for input to MARK;
#    Vector can be written to file with command
#    write.table(ch,filename,row.names=FALSE,col.names=FALSE,quote=FALSE)
#
	return(ch)
}

