\name{simburnham}
\alias{simburnham}
\title{Simulates data from Barker capture-recapture model}
\usage{
  simburnham(num.cohorts, cohort.sizes, s = list(),
    p = list(), r = list(), F = list(), design.data = NULL,
    outfile = NULL)
}
\arguments{
  \item{num.cohorts}{number of cohorts; design is square
  with same number of c-r eventsas num.cohorts; number of
  recapture events is num.cohorts-1}

  \item{cohort.sizes}{a scalar giving constant size of each
  cohort or a vector of sizes of length num.cohorts}

  \item{s}{a list defining the survival model with the
  following elements (see details) par - a vector of
  parameter values formula - a formula to use with
  design.data to construct model link - link function used
  with model to create probabilites (not used at present)}

  \item{p}{a list defining the capture probability model
  (same structure as s)}

  \item{r}{a list defining the recovery probability model
  (same structure as s)}

  \item{F}{a list defining the model of fidelity to capture
  region (same structure as s)}

  \item{design.data}{a dataframe with design data that
  allows model construction for probabilities (see
  details).}

  \item{outfile}{prefix name of the output file for the ch
  data. extension .inp is always added for MARK}
}
\value{
  dataframe with ch (capture history) as character
}
\description{
  Simulates capture-recapture data from Burnham (both)
  model of releases of multiple cohorts. Allows general
  model for s, p, r and F or constant or time-specific
  model depending on input arguments. Optionally outputs
  file to outfile with .inp extension for input to program
  MARK.
}
\details{
  The design.data dataframe contains a record for each
  animal for each time period it is contained in the study.
  The default fieldnames are id(unique # in cohort), cohort
  (cohort number from 1:num.cohorts), time (i:num.cohorts,
  where i is cohort #) and age (assumed to increment by 1
  each occasion and defaults with 0 initial age) The
  function create.simdesign.data can be used to create an
  initial dataframe which can be supplemented with other
  covariates. If a formula is given but no design.data is
  not provided this function

  The design.data dataframe contains a record for each
  animal for each time period it is contained in the study.
  The default fieldnames are id(unique # in cohort), cohort
  (cohort number from 1:num.cohorts), time (i:num.cohorts,
  where i is cohort #) and age (assumed to increment by 1
  each occasion and defaults with 0 initial age) The
  function create.simdesign.data can be used to create an
  initial dataframe which can be supplemented with other
  covariates. If a formula is given but no design.data is
  not provided this function calls create.simdesign.data to
  construct a default dataframe for the specified problem.

  For designation of models for s or p, there are 3 options
  here: constant model: s=list(par=value) time model:
  s=list(par=c(val1,val2,...valk)) where k=num.cohorts-1 is
  number of survival intervals general model:
  s=list(par=c(val1,val2,...valk),formula=~yourformula)) k
  is number of cols in model matrix See
  \code{\link{simcjs}} for more details.
}
\examples{
do.burnham <- function(s,p,r,F,reps,output=FALSE)
#
#  do.burnham -  a simple example piece of code to show simulation/analysis loop
#                with a Burnham model
#
#  Arguments:
#
#  s   -  parameter list for s
#  p   -  parameter list for p
#  r   -  parameter list for r
#  F   -  parameter list for F
#  reps-  number of simulation reps
#
#  Value:
#
#  results - for this simple example it will be a matrix of the real parameter estimates
#
#  Functions used: simburnham, mark
#
#
{
	results=NULL
	for(i in 1:reps)
	{
		cat("rep ",i,"\\n")
		xx=simburnham(3,500,s=s,p=p,r=r,F=F)
		mod<-mark(xx,title="Burnham test",model="Burnham",output=output)
		results=rbind(results,mod$results$real$estimate)
	}
	return(results)
}
}
\author{
  Jeff Laake <jeff.laake@noaa.gov>
}

