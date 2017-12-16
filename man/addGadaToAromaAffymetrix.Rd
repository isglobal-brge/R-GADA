\name{addGadaToAromaAffymetrix}
\alias{addGadaToAromaAffymetrix}

\title{Connect gada and aroma.affymetrix packages}

\description{
  Defines the GadaModel object and methods for calling GADA from within the aroma.affymetrix package framework. 
}


\usage{
addGadaToAromaAffymetrix()
}

\details{
  This function loads the Aroma.Affymetrix package if not loaded
  already, and creates the GadaModel object and methods such that we
  can use the gada package from within the aroma.affymetrix package.

  See the example in the gada user manual for more information, and the
  aroma.affymetrix package documentation.
}

\value{
  No return value, just prepares gada to be used under aroma.affymetrix package
}


\references{ 
Pique-Regi R, Caceres A, Gonzalez JR. "R-Gada: a package for fast detection and visualization of copy number alterations on multiple samples",  BMC Bioinformatics , Submitted Nov 2009
}


\keyword{utilities}

