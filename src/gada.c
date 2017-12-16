/*=================================================================
 * gada.c 
 *
 * R calling
 *
 *=================================================================*/

#include <R.h>
#include "BaseGADA.h"
//#include <math.h>
#include "RGADAdefines.h" // We would put specific defines to handle R calling in our code.

/// Implements the R calling to the Single_SBL_PWC_norm
void Rcall_Single_SBL_PWC_norm(
			       double *y 	///[in] The observed data (one single sample and chromosome
			       ,int *Iext 	///[in/out] Initial (Final) candidate breakpoints  [0,i_1,..,i_K,M]
			       ,double *alpha  ///[in/out] Estimated brkp hyperparameters  [alpha_1,...,alpha_K]
			       ,double *Wext   ///[in/out] Estimated brkp posterior mean [MeanY,MeanW_1,MeanW_2,...,MeanW_K]
			       ,double *Wsigma2 ///[in/out] Estimated brkp posterior variance [SigmaW_1,SigmaW_2,...,SigmaW_K]
			       ,int *K ///[in/out] Number of breakpoints 
			       ,double *sigma2 /// [in/out] Noise estimation
			       /// Algorithm parameters
			       ,double *aAlpha /// a -- hyperparameter alpha_i ~ Gamma(a,b)
			       ,double *bAlpha /// b -- hyperparameter alpha_i ~ Gamma(a,b) (usually 0 -- flat prior)
			       //	,double aSigma /// a -- hyperparameter 1/sigma2 ~ Gamma(a,b)
			       //	,double bSigma /// b -- hyperparameter 1/sigma2 ~ Gamma(a,b) 
			       //  ,int toggleSigma2est /// 1 with noise estimation 0 without noise estimation 
			       ,double *maxAlpha /// maximum setting for the alpha parameter before removal (1E8)
			       ,int *isRemovingOK /// 0 -- no basis elements are removed, 1 -- elements whose alpha > maxlpha are removed (1)
			       ,int *maxit /// Maximum number of EM iterations (10000)
			       ,double *tol /// Threshold on the minimum change on the alpha parameters to consider that the algorithm has converged (1E-10)
			       ,int *pdebug /// verbose mode if higher than 0, also asseses numerical stability assessment. 
			       )
{

  int k,M,m;
  double myaux;

  ///Vectors may be shortened by recalloc after running the function
  ///All vectors are assumed to be properlly allocated in memory
	
  M=Iext[*K+1];
	
  if(*pdebug>0)
    Rprintf("#(Rcall_Single_SBL_PWC_norm) K=%d,M=%d,sigma2=%g,a=%g,b=%g,maxAlpha=%g,RemovingFlag=%d,maxit=%d,tol=%g,debug=%d\n",K,M,*sigma2,*aAlpha,*bAlpha,*maxAlpha,*isRemovingOK,*maxit,*tol,*pdebug);

	    
  //If sigma2 < 0, compute sigma2
  if(*sigma2<0){
    *sigma2=0;
    for(m=1;m<M;m++){
      myaux=y[m]-y[m-1];
      *sigma2+=(0.5*myaux*myaux);
    }
    *sigma2=*sigma2/(M-1);
    if(*pdebug>0)
      Rprintf("(Rcall_Single_SBL_PWC_norm) Estimated sigma2=%g\n");
  }

  //Mean computation and removal
  myaux=0;
  for(m=0;m<M;m++)
    myaux+=y[m];
  myaux=myaux/M;
  for(m=0;m<M;m++)
    y[m]=y[m]-myaux;
  Wext[0]=myaux;
  
  //Call to SBL
  if(*pdebug>0)
    Rprintf("#(Rcall_Single_SBL_PWC_norm) Calling Single_SBL_PWC_norm\n");  
  *maxit=Single_SBL_PWC_norm(y,Iext,alpha,Wext,Wsigma2,K,sigma2,*aAlpha,*bAlpha,*maxAlpha,*isRemovingOK,*maxit,tol,*pdebug);
  if(*pdebug>0)
    Rprintf("#(Rcall_Single_SBL_PWC_norm) After SBL K=%d,M=%d,sigma2=%g,a=%g,b=%g,maxAlpha=%g,RemovingFlag=%d,maxit=%d,tol=%g,debug=%d\n",K,M,*sigma2,*aAlpha,*bAlpha,*maxAlpha,*isRemovingOK,*maxit,*tol,*pdebug);
  
   
  if(*pdebug>200){
    //Debugging
    Rprintf("#I:");
    for (k=0;k<=*K+1;++k)
      Rprintf("%d ",Iext[k]);
    Rprintf("\n");
  }
	
	
  if(*pdebug>200){
    //Debugging
    Rprintf("#W:");
    for (k=0;k<=*K;++k)
      Rprintf("%f ",Wext[k]);
    Rprintf("\n");
  }
}

void RcallSBLandBE(
	double *y,
	int *M,
	double *sigma2, 
	double *a,
	double *T, 
	int *MinSegLen, 
	int *pK, 
	int *Iext, 
	double *Wext,
	int *pVerboseFlag)
{
	int k,K;
	int *auxIext;
	double *auxWext;


	K=SBLandBE(y,*M,sigma2,*a,*T,*MinSegLen,&auxIext,&auxWext); 

	if(*pVerboseFlag>0)
	  Rprintf("# After SBLandBE -- %d discontinuities \n",K);

	*pK=K;

	for (k=0;k<=K+1;++k)
		Iext[k]=auxIext[k];
	

	if(*pVerboseFlag>2){
	//Debugging
	  Rprintf("#I:");
	  for (k=0;k<=K+1;++k)
	    Rprintf("%d ",Iext[k]);
	  Rprintf("\n");
	}
	
	for (k=0;k<=K;++k)
		Wext[k]=auxWext[k];
	
	if(*pVerboseFlag>2){
	//Debugging
	  Rprintf("#W:");
	  for (k=0;k<=K;++k)
	    Rprintf("%f ",Wext[k]);
	  Rprintf("\n");
	}


	myFree(auxIext);
	myFree(auxWext);
}

void RcallBEwTandMinLen
(
	double *Wext, //Input and output breakpoint weights
	int *Iext,    //Input and output breakpoint locations
	int *K,       // Input and output numof breakpoints
	const double *sigma2, //Input, noise estimation
	const double *T,      //Input, Backward Elimination Critical value
	const int *MinSegLen //Minimum number of probes in an alteration
)
{

	BEwTandMinLen(Wext,Iext,K,*sigma2,*T,*MinSegLen); 
	//Rprintf("# After BEwTandMinLen -- %d discontinuities \n",*K);

}

void RcallWextIextToSegments
(
	double *Wext, //Input and output breakpoint weights
	int *Iext,    //Input and output breakpoint locations
	int *K,       // Input and output numof breakpoints
	double *SegAmp, // Segment amplitudes
	int *SegLen // Segment lengths
)
{
//	Rprintf("#  WextIextToSegments -- %d discontinuities \n",*K);
	IextToSegLen(Iext,SegLen,*K);
	IextWextToSegAmp(Iext,Wext,SegAmp,*K);
}

///Returns the median intensity of all segments to estimate the base-line hybridization intensity associated with the ploidy copy number
void RcallCompAmpMedianMethod(
			      int *SegLen /// Lengths corresponding to the segment amplitudes
			      ,double *SegAmp /// Amplitudes of the segments.
			      ,int *K ///Number of Segments
			      ,double *RefLevel ///Reference base-line hybridization intensity
			      )
{
  //Rprintf("# RcallCompAmpMedianMethod_ K=%d\n",*K);
  *RefLevel=CompBaseAmpMedianMethod(SegLen,SegAmp,*K);
  //Rprintf("# RcallCompAmpMedianMethod_ BaseAmp=%g\n",*RefLevel);
}

/// Classifies Segment amplitudes to baseline reference level
void RcallClassifySegments(
    int *SegLen ///Segment Lengths
    ,double *SegAmp ///Segment Amplitudes (input output)
    ,double *SegState ///Segment State, Loss (-1),Gain(1),Neutral(0)
    ,int *K ///Number of segments.
    ,double *BaseAmp ///Reference amplitude to compare (see CompBaseAmpMedianMethod)
    ,double *sigma2  ///Estimate of the noise variance. 
    ,double *T	    ///Critical value of the T test to decide Gain/Loss against Neutral as null hyphotesis
    )
{
  //Rprintf("# RcallClassifySegments_ K=%d BaseAmp=%g sigma2=%g T=%g \n",*K,*BaseAmp,*sigma2,*T);
  ClassifySegments(SegAmp,SegState,SegLen,*K,*BaseAmp,*sigma2,*T);
}
    			  


