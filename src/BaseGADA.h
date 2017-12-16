/*=================================================================
* BaseGenomeBreaks.h 
*
*=================================================================*/
/*

This File is part of GADA
GADA v1.0 Genome Alteration Detection Algorithm 
Copyright (C) 2008  Childrens Hospital of Los Angeles

GADA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

GADA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GADA.  If not, see <http://www.gnu.org/licenses/>.

Author: 
Roger Pique-Regi    piquereg@usc.edu
*/

#ifndef _BaseGADA_H_
#define _BaseGADA_H_

//#include "matlabdefines.h"

void reconstruct (double *wr,int M,double *aux_vec);
void intBubbleSort (int *I,int L);
void doubleBubbleSort (double *D,int *I,int L);
void TrisolveREG(double *t0,double *tu,double *tl,double *coef,double *sol,int sizeh0);
// DiagOfTriXTri  Diagonal of L*R where L,R tridiagonal (MexDiagOfTriXTri)
void   
DiagOfTriXTri(
			  const double *ll, //(In) ll Lower diagonal
			  const double *l0, //(In) l0 Central diagonal
			  const double *lu, //(In) lu Upper diagonal
			  const double *rl, //(In) rl Lower diagonal
			  const double *r0, //(In) r0 Central diagonal
			  const double *ru, //(In) ru Upper diagonal
			  double *d,  //(Out) d=DiagOfTriXTri(ll,l0,lu,rl,r0,ru)
			  int N  //(In) Number of variables or length of l0, 
			  );

void tridiagofinverse(double *t0,double *tl,double *itl,double *it0,double *itu,int N,double *d,double *e);

void ForwardElimination(double *A,int N);

void BackSubstitution(double *A,int N);

void BackwardElimination(double *A,int N);

void TriSolveINV(double *AA,int M, int N, double *x,double *d,double *e);

void ComputeH(double *h0, double *h1, int M);

void ComputeFdualXb(int M, double *b);

// 20080119 REMOVED void ComputeHs(int *s,double *a,int M,int Ms,double *h0,double *h1);

void ComputeHs(int *s,int M,int Ms,double *h0,double *h1);

void TridSymGaxpy(double *t0, double *t1, double *x, int M, double *y);

void ComputeT(double *h0,double *h1,int M,double *alfa,double sigma,double *t0,double *tl,double *tu);

int findminus(double *alpha,int Ms,double maxalpha,int *sel);

int simpletresholding(double *inputvector,int N,double thres,double *disc);

void computesegmentmeans(double *inputvector,int N,double *disc,int numdisc,double *amp);

void reconstructoutput(double *rec,int N,double *disc,int numdisc,double *amp);

int BEthresh( //To eliminate...   

			 double *Scores,

			 int Nscores,

			 double *wr,

			 int *indsel,

			 int *pointNumRem,

			 double *pointTau

			 );


int SBLandBE( //Returns breakpoint list lenght.

			 double *tn,

			 int M,  //length of the noisy signal tn

			 double *sigma2, //If sigma2 < 0, compute sigma2 (Input/Output)

			 double a,      // SBL parameter

			 double T,      // Threshold to prune

			 int MinSegLen, //Minimum length of the segment.

			 int **pI,   //Returns breakpoint positions

			 double **pw //Returns breakpoint weights.

			 //int *pK    

			 );





void Project(

			 double *y,

			 int M,

			 int *I,

			 int L,

			 double *xI,

			 double *wI

			 );


/// Implements the SBL algorithm for single sample using a the normalized PWC representation
int
Single_SBL_PWC_norm(
    double *y 	///[in] The observed data (one single sample and chromosome
    ,int *Iext 	///[in/out] Initial (Final) candidate breakpoints  [0,i_1,..,i_K,M]
    ,double *alpha  ///[in/out] Estimated brkp hyperparameters  [alpha_1,...,alpha_K]
    ,double *Wext   ///[in/out] Estimated brkp posterior mean [MeanY,MeanW_1,MeanW_2,...,MeanW_K]
    ,double *Wsigma2 ///[in/out] Estimated brkp posterior variance [SigmaW_1,SigmaW_2,...,SigmaW_K]
    ,int *pK ///[in/out] Number of breakpoints 
    ,double *sigma2 /// [in/out] Noise estimation
		    /// Algorithm parameters
    ,double aAlpha /// a -- hyperparameter alpha_i ~ Gamma(a,b)
    ,double bAlpha /// b -- hyperparameter alpha_i ~ Gamma(a,b) (usually 0 -- flat prior)
    //	,double aSigma /// a -- hyperparameter 1/sigma2 ~ Gamma(a,b)
    //	,double bSigma /// b -- hyperparameter 1/sigma2 ~ Gamma(a,b) 
    //  ,int toggleSigma2est /// 1 with noise estimation 0 without noise estimation
    ,double maxAlpha /// maximum setting for the alpha parameter before removal
    ,int isRemovingOK /// 0 -- no basis elements are removed, 1 -- elements whose alpha > maxlpha are removed
    ,int maxit /// Maximum number of EM iterations
    ,double *tol /// Threshold on the minimum change on the alpha parameters to consider that the algorithm has converged
    ,int debug /// verbose mode if higher than 0, also asseses numerical stability.
    );




void IextToSegLen(

				  int *Iext, // Enters Iext

				  int *SegLen,		// Outputs SegLen.. can be the same as Iext?

				  int K			// Length Iext - 1

				  );

void IextWextToSegAmp(

					  int *Iext, 

					  double *Wext,

					  double *SegAmp,

					  int K 

					  );

void CompZ(// computes z=F'y for entire possilbe breakpoint positions (normalized PWC)

		   double *y,

		   double *z,

		   int M

		   );

void ComputeHsIext(

				   //input variables:

				   int *Iext,     // Indices of selection,

				   int K,     // Length of the indices, 

				   double *h0, // Returning diagonal of H,

				   double *h1  // Returning upper diagonal of H

				   );

void ProjectCoeff ( //IextYobs2Wext

				   double *y,

				   int M,

				   int *Iext,

				   int K,

				   double *Wext

				   );

void 

CollapseAmpTtest(//Uses a T test to decide which segments collapse to neutral
				 double *SegAmp, //Segment Amplitudes (input output)
				 int *SegLen, //Segment Lengths
				 int L, //Number of segments.
				 double BaseAmp, //Reference amplitude to compare
				 double sigma2,  //Reference noise 
				 double T		//Critical value that decides when to colapse
				 );

double 
CompBaseAmpMedianMethod( 
    int *SegLen    /// Lengths corresponding to the segment amplitudes
    ,double *SegAmp /// Amplitudes of the segments.
    ,int K /// Number of segments. 
    );

/// Classifies Segment amplitudes to baseline reference level
void 
ClassifySegments(
    double *SegAmp ///Segment Amplitudes (input output)
    ,double *SegState ///Segment State, Loss (-1),Gain(1),Neutral(0)
    ,int *SegLen ///Segment Lengths
    ,int K ///Number of segments.
    ,double BaseAmp ///Reference amplitude to compare (see CompBaseAmpMedianMethod)
    ,double sigma2  ///Estimate of the noise variance. 
    ,double T	    ///Critical value of the T test to decide Gain/Loss against Neutral as null hyphotesis
    );

void 
EstimateSegmentDiscreteCN(
    double *SegAmp 
    ,int *SegLen,
    double *SegState,
    int K,
    double BaseAmp,
    double ploidy,
    double sigma2,  //Reference noise 
    double T		//Critical value that decides when to colapse
    );





void
ComputeTScores(

			   const double *Wext,

			   const int *Iext,

			   double *Scores,

			   int K,

			   int start,

			   int end

			   );



int BEwTscore(

			  double *Wext,  //IO Breakpoint weights extended notation...

			  int *Iext,     //IO Breakpoint positions in extended notation...

			  double *tscore,

			  int *pK,       //IO Number breakpoint positions remaining.

			  double T      //IP  Threshold to prune

			  );



int BEwTandMinLen( //Returns breakpoint list lenght. with T and MinSegLen

				  double *Wext,  //IO Breakpoint weights extended notation...

				  int *Iext,     //IO Breakpoint positions in extended notation...

				  int *pK,       //IO Number breakpoint positions remaining.

				  double sigma2, //IP If sigma2 

				  double T,      //IP  Threshold to prune,  T=T*sqrt(sigma2);

				  int MinSegLen  //IP Minimum length of the segment.

				  );

int

RemoveBreakpoint(

				 double *Wext,

				 int *Iext,

				 int K,

				 int jrem

				 );

void ReconstructFromIextWext(
		int *Iext ///(in) Breakpoint set defining the segmentation.
		,double *Wext ///(in) Breakpoint coefficients defining the segmentation.
		,double *xRec ///(out) Full reconstruction 
		,int K ///(in) Number of breakpoints
		);

#ifdef _DEPRECATED_

int SBL(

		double *y, //I -- 1D array with the input signal

		int *I, //IO -- 1D array with the initial (final) candidate breakpoints

		double *alpha, //I -- 1D array with the initial (final) hyperparameter inv. varainces.

		double *w, //O -- 1D array containing the breakpoint weigths or posterior mean. 

		double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha and H

		int M, //Initial size of the array in y

		int *K, //Size of the I alpha w



		//Algorithm parameters:

		double sigma2, //Noise estimated 

		double a,      //

		double b,       

		double maxalpha,  //Basis reduction parameter 

		int    maxit,     //Max number of iterations

		double tol,       //Tolerance for convergence

		int debug       //verbosity... set equal to 1 to see messages  0 to not see them

		);


#endif //_DEPRECATED_


#endif //_BaseGADA_H_





