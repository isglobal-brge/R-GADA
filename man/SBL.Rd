\name{SBL}
\alias{SBL}
\alias{print.SBL}

\title{Sparse Bayesian Learning (SBL) Segmentation Algorithm}

\description{
  Fits the a sparse Bayesian learning (SBL) model on a single sample setupGADA object 
}

\usage{
SBL(x, sigma2, aAlpha = 0.2, 
    estim.sigma2 = FALSE, maxit = 10000, tol = 1e-08, 
    verbose = FALSE, saveInfo = TRUE)
}

\arguments{
  \item{x}{an object of class 'setupGADA' prepared using
  'setupGADA'}
  \item{sigma2}{the array noise level. See details}
  \item{aAlpha}{sparseness hyperparameter. See details}
  \item{estim.sigma2}{Should the 'sigma2' be estimated? the default is FALSE}
  \item{maxit}{maximum number of iterations in the SBL algorithm}
  \item{tol}{tolerance criteria to stop the SBL algorithm}
  \item{verbose}{print verbose information,  the default is FALSE (useful to debug errors)}
  \item{saveInfo}{TRUE if the annotation data is transfered to the
  retuned SBL object. The default is TRUE}
}


\details{

  This function fits a SBL model on a single DNA array observation. 
  The underlying copy number is assumed assumed to be piecewise constant
  (PWC) with a sparse number of breakpoints but the observed array is
  degraded by noise.

  The array noise level can be provided by the 'sigma2' parameter or if
  it is unkown we can use 'estim.sigma2=TRUE' to estimate the
  hybridization noise level. 
  
  The SBL model places an hierarchical Bayesian prior over the
  breakpoints delimiting the probes that fall into a copy number a
  piece-wise constant (PWC) vector. This Bayesian prior is uninformative
  about the magnitude and position of the breakpoints but enforces
  sparseness (i.e., assumes that only very few breakpoints are true
  positives). The hyperparameter 'aAlpha' is used to control the
  sparseness level. Instead of adjusting we recomment to use
  a high sensitivity value (e.g. 'aAlpha=0.2', default value) and adjust
  the False Discovery Rate (FDR) using the
  'BackwardElimination' procedure. The 'aAlpha=0.2' can be
  increased to obtain a faster result which may be interesting
  in case of very high density arrays, but we may not be able
  to recover some of the segments that would be otherwise
  obtained with a high sensitivity setting with the BackwardElimination
  procedure.

  The SBL model is fit using an expectation maximization (EM)
  algorithm. The 'tol' parameter sets the maximum allowed change
  on the model parameters to consider that the algorithm has
  converged. The 'maxit' parameter establishes a maximum number of
  iterations. If the algorithm does not converges before 'maxit' number
  of iterations have been computed a warning is returned with the last
  change in magnitude. 
   
}


\value{

  An object of class 'SBL' to be used with the 'BackwardElimination'
  algorithm. 

  'print' returns the number of discontinuities or segments by
  chromosome and convergence information about the SBL algorithm.

}


\references{ 

Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples". BMC Bioinformatics, 2010;11:380.

Pique-Regi R, Monso-Varona J,Ortega A, Seeger RC, Triche TJ, Asgharzadeh S. "Sparse representation and Bayesian detection of the genome copy number alterations from microarray data", Bioinformatics , 2008;24(3):309-18.

}


\seealso{\code{\link{setupGADA}}, \code{\link{BackwardElimination}}, \code{\link{parSBL}}, \code{\link{parBE}}}


\examples{
 \dontrun{
  See the vignette
 } 
}

\keyword{models}
