# R-GADA
R implementation of GADA (Genetic Alteration Detector Algorithm) - used to detect CNVs from aCGH and intensity array SNP data

Package is installed by:

`devtools::install_github("isglobal-brge/R-GADA")`

and it can be loaded into R by:

`library(gada)`

GADA requires C and Fortran compilation. Therefore, gfortran is required to be installed for those researchers
using Windows. For those of you having problems with this, you can install the compiled .zip file 
