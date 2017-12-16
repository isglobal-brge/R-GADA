/*=================================================================
 * BaseGADA.c 
 * All basic functions for SBL, BE, etc. 
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

  Authors: 
  Roger Pique-Regi    piquereg@usc.edu
  Jordi Monso-Varona

*/
#include "BaseGADA.h"
#include <math.h>

#ifdef _MATLAB_
#include "matlabdefines.h"
#endif //_MATLAB_  memory allocation and messages for Matlab

#ifdef _RGADA_
#include "RGADAdefines.h"
#endif //_RGADA_ memory allocation and messages for R

#ifdef _CONSOLE_   //
#include "consoledefines.h"
#endif //_CONSOLE_ memory allocation and messages for console application

#define min(x,y) x<y?x:y

/* General Notation definitions
   K		= number of breakpoints, num of segments is K+1 
   M		= number of probes in chromosome unit of analysis (can be entire chromosome or arm)
   Ir Ired	= Breakpoint set, double array I[0] contains first breakpoint I[K-1] contains Kth breakpoint length 0 if no breakpoints
   Ie Iext	= Breakpoint set, extended notation I[0] contains 0 I[K] contains Kth breakpoint I[K+1]=M
   I        = Can refer to Ir or Ie types of vectors
   We Wext	= Aligned with Ie We[0] contains the overall mean and We[K+1] nothing
   SegAmp   = Segment amplitudes, lenght K+1
   SegLen   = Segment longitudes in probes, lenght K+1 and should sum up to M?
*/

/// integer list Bubble Sort
void intBubbleSort (
    int *I  ///[in,out] Input vector to sort
    ,int L  ///[in]  Size of the input vector
    )
{
    int i=0;
    int aux=0;
    
    for (i=0;i<L-1;i++) {
	if (I[i]>I[i+1]){
	    aux=I[i+1];
	    I[i+1]=I[i];
	    I[i]=aux;
	    i=-1;
	}
    }
}

/// Double list Bubble Sort 
void doubleBubbleSort (
    double *D ///[in,out] Input vector to sort
    ,int *I   ///[in,out] The same permutation is applied to this vector as in D, useful to find the indices, if we want to sort another vector with same order.
    ,int L	 ///[in] Size of the list to sort
    )
{
    int i=0;
    double Daux=0;
    int Iaux=0;

    for (i=0;i<L-1;i++) {
	if (D[i]>D[i+1]){
	    Daux=D[i+1];
	    D[i+1]=D[i];
	    D[i]=Daux;
			
	    if(I!=NULL)
	    {
		Iaux=I[i+1];
		I[i+1]=I[i];
		I[i]=Iaux;
	    }

	    i=-1;
	}
    }
}

/// Caluculates the norm infinity of the difference between two vectors
double ///Returns the norm
normInfOfDiff( 
    double *x1  /// First vector
    ,double *x2 /// Second vector
    ,int N /// Size of x1 and x2
    )
{   
    int n;
    double norm,aux;

    norm=0;
    for (n=0;n<N;n++)
    {
	aux=(x1[n]-x2[n]);
	if (aux<0)aux=-aux;
	if (aux>norm)norm=aux;          
    }             

    return norm;
}

/// Caluculates the norm 1 of the difference between two vectors
double ///Returns the norm
norm1OfDiff( 
    double *x1  /// First vector
    ,double *x2 /// Second vector
    ,int N /// Size of x1 and x2
    )
{   
    int n;
    double norm,aux;

    norm=0;
    for (n=0;n<N;n++)
    {
	aux=(x1[n]-x2[n]);
	if (aux<0)aux=-aux;
	norm+=aux;          
    }             

    return norm;
}
    
/// Calculates the norm 2 of the difference between two vectors
double ///returns the norm
norm2OfDiff( 
    double *x1  /// First vector
    ,double *x2 /// Second vector
    ,int N /// Size of x1 and x2
    )
{
    int n;
    double norm;

    norm=0;
    for (n=0;n<N;n++)
	norm+=(x1[n]-x2[n])*(x1[n]-x2[n]);

    return sqrt(norm);
}


/************************************************************************/
/*************** TRIDIAGONAL MATRIX OPERATION FUNCTIONS *****************/
          
// DiagOfTridXTrid Diagonal of L*R where L,R tridiagonal (MexDiagOfTridXTrid)
/// DiagOfTridXTrid computes \f${\bf d}=diag({\bf L R} )\f$ where \f$\bf L\f$ and \f$\bf R\f$  are tridiagonal matrices
void   
DiagOfTridXTrid(
    const double *ll ///(in) ll Lower diagonal left matrix L.
    ,const double *l0 ///(in) l0 Central diagonal left matrix L.
    ,const double *lu ///(in) lu Upper diagonalleft matrix L.
    ,const double *rl ///(in) rl Lower diagonal right matrix R.
    ,const double *r0 ///(in) r0 Central diagonal rigth matrix R.
    ,const double *ru ///(in) ru Upper diagonal rigth matrix R.
    ,double *d  ///(out) d diagonal of the product L R. Can replace r0 or l0.
                //d=DiagOfTridXTrid(ll,l0,lu,rl,r0,ru)
    ,int N  ///(in) Number of variables or length of l0, 
    )
{
    int i;
    
    d[0]=l0[0]*r0[0]+lu[0]*rl[0];
    for (i=1;i<N-1;i++){ //It can be run in parallel
        d[i]=ll[i-1]*ru[i-1]+l0[i]*r0[i]+lu[i]*rl[i];
    }
    d[N-1]=ll[N-2]*ru[N-2]+l0[N-1]*r0[N-1];
    
}

/// Computes the tridiagonal part of inv(T) given the LDU and UDL decompositions
void   
TridOfInvTrid(
    double *Fu ///(in) Upper diagonal of U in LDU .
    ,double *Fd ///(in) Diagonal of D in LDU.
    ,double *Bl ///(in) Lower diagonal of L in UDL.
    ,double *itl ///(out) itl First lower diagonal of inv(T).
    ,double *it0 ///(out) it0 Central diagonal of inv(T).
    ,double *itu ///(out) itu First upper diagonal of inv(T).
    ,int N  ///(in) N Number of variables.
    )
{   
    int i;    

    

    // Compute main diagonal terms of inv(T)
    for(i=0;i<N-1;i++)
	it0[i]=1/(Fd[i] * (1-Fu[i]*Bl[i]));
    it0[N-1]=1/Fd[N-1];
   
    
    //Computes off-diagonal terms
    for (i=0;i<(N-1);i++)
	itl[i]=-Bl[i]*it0[i];
    for (i=0;i<(N-1);i++)
	itu[i]=-Fu[i]*it0[i+1];

}


///  LDU factorization of a tridiagonal matrix T using Forward Gaussian Elimination
void
LDUofTrid(
    double *tu  ///[in] 1st Upper Diagonal T matrix 
    ,double *tc ///[in] Central Diagonal T matrix   
    ,double *tl ///[in] 1st Lower Diagonal T matrix 
    ,int N ///[in] Number of variables
    ,double *u ///[out] Upper Diagonal of the U matrix, can overwrite tu
    ,double *d ///[out] Diagonal of the D matrix, can overwrite t0
    ,double *l ///[out] Lower Diagonal of the L matrix, can overwrite tl
    )
{ 
    int i;

    ///\note 
    ///- If the output parameter u or l are Null pointers then, they will not be computed.
    ///- It takes 3N-2 flops to compute d, and N-1 flops to compute u and l.  
 

    d[0]=tc[0];
    for(i=0;i<N-1;i++) 
	d[i+1]=tc[i+1]-tl[i]*tu[i]/d[i];
  
    if(l!=NULL)
	for(i=0;i<N-1;i++)
	    l[i]=tl[i]/d[i];

    if(u!=NULL)
	for(i=0;i<N-1;i++)
	    u[i]=tu[i]/d[i];  
   
}


/// Solve the tridiagonal system Tx=b using existing T=LDU (Forward Elimination / Back Substitution)
void 
TridSolveUsingLDU(
    double *u ///[in] Upper Diagonal of the U matrix, can overwrite tu
    ,double *d ///[in] Diagonal of the D matrix, can overwrite t0
    ,double *l ///[in] Lower Diagonal of the L matrix, can overwrite tl
    ,int N ///[in] Number of variables       
    ,double *b ///[in] Right side of the equation. 
    ,double *x ///[out] Vector with the solutions can safely replace b
    )
{
    int i;

    ///\note: 5N flops are used.
  
    //Forward Elimination
    x[0]=b[0];
    for(i=1;i<N;i++)
	x[i]=b[i]-l[i-1]*x[i-1];
    for(i=0;i<N;i++)
	x[i]=x[i]/d[i];

    //Back Substitution
    //x[N-1]=x[N-1];
    for(i=N-2;i>=0;i--)
	x[i]=x[i]-u[i]*x[i+1];  
}

/// Solve the tridiagonal system Tx=b by first computing  T=LDU (Forward Elimination / Back Substitution)
void 
ForwardTridSolve(
    double *tu  ///[in] 1st Upper Diagonal T matrix 
    ,double *tc ///[in] Central Diagonal T matrix   
    ,double *tl ///[in] 1st Lower Diagonal T matrix 
    ,int N ///[in] Number of variables
    ,double *u ///[out] Upper Diagonal of the U matrix, can overwrite tu
    ,double *d ///[out] Diagonal of the D matrix, can overwrite tc
    ,double *l ///[out] Lower Diagonal of the L matrix, can overwrite tl
    ,double *b ///[in] Right side of the equation.
    ,double *x ///[out] Vector with the solutions can safely replace b
    )
{ 
    LDUofTrid(tu,tc,tl,N,u,d,l);
    TridSolveUsingLDU(u,d,l,N,b,x);
}


///  UDL factorization of a tridiagonal matrix T using Backward Gaussian Elimination
void
UDLofTrid(
    double *tu  ///[in] 1st Upper Diagonal T matrix 
    ,double *tc ///[in] Central Diagonal T matrix   
    ,double *tl ///[in] 1st Lower Diagonal T matrix 
    ,int N ///[in] Number of variables
    ,double *u ///[out] Upper Diagonal of the U matrix, can overwrite tu
    ,double *d ///[out] Diagonal of the D matrix, can overwrite t0
    ,double *l ///[out] Lower Diagonal of the L matrix, can overwrite tl
    )
{ 
    int i;


    ///\note 
    ///- If the output parameter u or l are Null pointers then, they will not be computed.
    ///- It takes 3N-2 flops to compute d, and N-1 flops to compute u and l.  
 

    d[N-1]=tc[N-1];
    for(i=N-2;i>=0;i--) 
	d[i]=tc[i]-tl[i]*tu[i]/d[i+1];
  
    if(l!=NULL)
	for( i = 0 ; i < N-1 ; i++ )
	    l[i]=tl[i]/d[i+1];

    if(u!=NULL)
	for(i=0;i<N-1;i++)
	    u[i]=tu[i]/d[i+1];  
   
}

/// Solve the tridiagonal system Tx=b using existing T=UDL (Backward Elimination / Forward Substitution)
void 
TridSolveUsingUDL(
    double *u ///[in] Upper Diagonal of the U matrix, can overwrite tu
    ,double *d ///[in] Diagonal of the D matrix, can overwrite t0
    ,double *l ///[in] Lower Diagonal of the L matrix, can overwrite tl
    ,int N ///[in] Number of variables       
    ,double *b ///[in] Right side of the equation. 
    ,double *x ///[out] Vector with the solutions can safely replace b
    )
{
    int i;

    ///\note: 5N flops are used.
  
    //Backward Elimination
    x[N-1]=b[N-1];
    for(i=N-2;i>=0;i--)
	x[i]=b[i]-u[i]*x[i+1];
    for(i=0;i<N;i++)
	x[i]=x[i]/d[i];

    //Forward Substitution
    //x[0]=x[0];
    for(i=1;i<N;i++)
	x[i]=x[i]-x[i-1]*l[i-1];  
}

/// Solve the tridiagonal system Tx=b by first computing  T=UDL (Backward Elimination / Forward Substitution)
void 
BackwardTridSolve(
    double *tu  ///[in] 1st Upper Diagonal T matrix 
    ,double *tc ///[in] Central Diagonal T matrix   
    ,double *tl ///[in] 1st Lower Diagonal T matrix 
    ,int N ///[in] Number of variables
    ,double *u ///[out] Upper Diagonal of the U matrix, can overwrite tu
    ,double *d ///[out] Diagonal of the D matrix, can overwrite tc
    ,double *l ///[out] Lower Diagonal of the L matrix, can overwrite tl
    ,double *b ///[in] Right side of the equation.
    ,double *x ///[out] Vector with the solutions can safely replace b
    )
{
    UDLofTrid(tu,tc,tl,N,u,d,l);
    TridSolveUsingUDL(u,d,l,N,b,x);
}


/// Tridiagonal matrix symetric T Gaxpy  y <-- Tx + y
void 
TridSymGaxpy(
    //input variables:
    double *t0,
    double *t1,
    double *x,
    int M,
    double *y
    
    )
{
    int i;
    
    if (M==1){
        y[0]=t0[0]*x[0];
    }
    else
    {        
        y[0]=t0[0]*x[0]+t1[0]*x[1];        
        for (i=1;i<M-1;i++){
            y[i]=t0[i]*x[i]+t1[i]*x[i+1]+t1[i-1]*x[i-1];
        }
        y[M-1]=t0[M-1]*x[M-1]+t1[M-2]*x[M-2];
    }
}

/********************* End Tridiagonal operations *****************************/
/******************************************************************************/

/******************************************************************************/
/* PWC specific functions, considering specific  normalization of F   *********/
/******************************************************************************/

// Alternative formulations with differen normalization constants may lead to different implementations of F
// which may or may not have different properties.
// Here we use that the columns of F have zero mean and unit norm.

/// Extract the segment amplitudes from segmentation (Iext,Wext)
void IextWextToSegAmp(
    int *Iext ///(in) Breakpoint set defining the segmentation.
    ,double *Wext ///(in) Breakpoint coefficients defining the segmentation.
    ,double *SegAmp ///(out) Amplitude of the segments
    ,int K ///(in) Number of breakpoints
    )
{
    int k;
    double M,TotalMean,AuxMean;
    M=Iext[K+1];

    TotalMean=Wext[0];
    SegAmp[0]=0;
    for(k=1;k<=K;k++)
	SegAmp[k]=Wext[k]/sqrt((double)(M-Iext[k])*(double)Iext[k]/M)+SegAmp[k-1];
    AuxMean=0;
    for(k=0;k<=K;k++)
	AuxMean=AuxMean+SegAmp[k]*(double)(Iext[k+1]-Iext[k]);
    AuxMean=AuxMean/M;
    for(k=0;k<=K;k++)
	SegAmp[k]=SegAmp[k]-AuxMean+TotalMean;

}

/// Reconstructs the output from the segmentation 
void ReconstructFromIextWext(
    int *Iext ///(in) Breakpoint set defining the segmentation.
    ,double *Wext ///(in) Breakpoint coefficients defining the segmentation.
    ,double *xRec ///(out) Full reconstruction 
    ,int K ///(in) Number of breakpoints
    )
{
    int m,k;
    double *SegAmp;

    SegAmp=myCalloc(K+1,sizeof(double));

    IextWextToSegAmp(Iext,Wext,SegAmp,K);

    for(k=0;k<=K;k++)
	for(m=Iext[k];m<Iext[k+1];m++)
	    xRec[m]=SegAmp[k];

    myFree(SegAmp);

}

/// Extract projection coefficients onto breakpoint set aka IextYobs2Wext
void ProjectCoeff ( 
    double *y       ///[in]  Observed vector.
    ,int M		   ///[in]  Size of the observed vector, M=Iext[K+1]
    ,int *Iext      ///[in]  Breakpoint positions
    ,int K		   ///[in]  Number of breakpoints
    ,double *Wext   ///[out] Breakpoint coefficients
    )
{
    // Intern variables declaration
    double ymean=0;
    int i=0;
    double *h0;
    double *h1;
    double *z;
    
    // Variables inizialitation
    h0=myCalloc(K,sizeof(double));
    h1=myCalloc(K-1,sizeof(double));
    z=myCalloc(M-1,sizeof(double));
    
    // Compute the mean
    ymean=0;    
    for (i=0;i<M;i++)
	ymean=ymean+y[i];
    ymean=ymean/M;

    // Remove y's mean, Not necessary
    //    for (i=0;i<M;i++)
    //        y[i]=y[i]-ymean;	

    if (K>0)
    {
	//I sort -> assumed that I is already ordered
	// intBubbleSort(Iext,K+2);
	
	CompZ(y,z,M);
	//		myPrintf("\n CHECKING: ymean=%g,M=%d,K=%d\n",ymean,M,K);//
	//		myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M-3]=%g,z[M-2]=%g\n",z[0],z[1],z[2],z[M-3],z[M-2]);

	for (i=1;i<=K;i++)
	    z[i-1]=z[Iext[i]-1];
	//    	myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[K-2]=%g,z[K-1]=%g\n",z[0],z[1],z[2],z[K-2],z[K-1]);   

	ComputeHsIext(Iext,K,h0,h1);
	//myPrintf("\n h0 CHECKING: h0[0]=%g,h0[1]=%g,h0[2]=%g,h0[K-2]=%g,h0[K-1]=%g,h0[K]=%g\n",h0[0],h0[1],h0[2],h0[K-2],h0[K-1]);   
	//myPrintf("\n h1 CHECKING: h1[0]=%g,h1[1]=%g,h1[2]=%g,h1[K-2]=%g,h1[K-1]=%g\n",h1[0],h1[1],h1[2],h1[K-2]);   
		
	for (i=0;i<K+1;i++)
	    Wext[i]=0.0;

	TridSymGaxpy(h0,h1,z,K,Wext+1);
	//   	myPrintf("\n w CHECKING: w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]);   
    }
    Wext[0]=ymean;
    
    myFree(h0);
    myFree(h1);
    myFree(z);
   
}

/// Computes b <- inv(F'F)Fb When F comprises the entire data. (normalized PWC)
void 
ComputeFdualXb(
    int M,
    double *b
    )
{
    int i=0;
    double a=0;
    
    /* diff */
    for (i=0;i<M-1;i++){
        b[i]=b[i+1]-b[i];
        
    }
    b[M-1]=0;
    for (i=0;i<M-1;i++){
        a=(double)(M-1-i)*(double)(i+1)/(double)M;
        b[i]=b[i]*sqrt(a);
    }
}

/// Computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
void CompZ(
    double *y, /// Input with M elements
    double *z, /// Output with M-1 elements already allocated in memory.
    int M      
    )
{
    double LeftSum;
    double RightSum;
    double *aux;
    int m;

    aux=myCalloc(M-1,sizeof(double));
	
    LeftSum=0;
    RightSum=0;
    for(m=1;m<M;m++)
    {
	LeftSum=LeftSum+y[m-1];
	aux[m-1]=(-1)*sqrt( (double)(M-m) / ((double)M*(double)m) )*LeftSum;	
    }
    for(m=M-1;m>=1;m--)
    {
	RightSum=RightSum+y[m];
	z[m-1]=aux[m-1]+sqrt( (double)m / ((double)(M-m)*(double)(M)) )*RightSum;	
    }
//	myPrintf("\n CHECKING: Rsum=%g,Lsum=%g,M=%d\n",RightSum,LeftSum,M);
//	myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M-3]=%g,z[M-2]=%g\n",z[0],z[1],z[2],z[M-3],z[M-2]);
//	myPrintf("\n z CHECKING: z[0]=%g,z[1]=%g,z[2]=%g,z[M-1]=%g,z[M]=%g\n",aux[0],aux[1],aux[2],aux[M-1],aux[M]);

    myFree(aux);
}



/// Deprecated -- Compute H   fuction [h0,h1]=CompH(dim) 
void   
ComputeH(
    //Input variables:
    double *h0,
    double *h1,
    int M  //Number of variables
    )
{
    int i;

    for (i=0;i<M-1;i++){
        h0[i]=((double)(M-1-i)*(double)(i+1))/(double)M;
    }
    for (i=0;i<M-2;i++){
        h1[i]=-sqrt((h0[i+1]*h0[i]));
    }
    for (i=0;i<M-1;i++){
        h0[i]=2*h0[i];
    }      
}


/// deprecated Computes the H at the vector of indices...
void 
ComputeHs(
    //input variables:
    int *s,     // Indices of selection,
    int M,      // Length of the chormosome, 
    int K,     // Length of the indices, 
    double *h0, // Returning diagonal of H,
    double *h1  // Returning upper diagonal of H
    )
{
    int i;
    double iC,iL,iR;
    //double M;
    //M=(double)MM;
    
    i=0;
    
    if (K==1)
    {   
        iL=0;iC=(double)(s[i]+1);iR=M;    
        h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
        //h1[i]=0;   
    }
    else
    {
        iL=0;iC=(double)(s[i]+1);iR=(double)(s[i+1]+1);
        h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
        h1[i]=-sqrt((M-iC)*iC*(M-iR)*iR)/(M*(iR-iC));
        for (i=1;i<K-1;i++){
            //        h0[i]=((double)((M-s[i]-1)*(s[i]+1))/(double)M)*(double)(s[i+1]-s[i-1])/(double)((s[i+1]-s[i])*((s[i]-s[i-1])));
            //        h1[i]=sqrt((double)(((M-s[i]-1)*(s[i]+1))*(M-s[i+1]-1)*(s[i+1]+1)))/(double)(M*(s[i+1]-s[i]));
            iL=iC;iC=iR;iR=(double)(s[i+1]+1);
            h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
            h1[i]=-sqrt((M-iC)*iC*(M-iR)*iR)/(M*(iR-iC));
        }
        //i=K-1
        iL=iC;iC=iR;iR=M;
        h0[i]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
        //    h0[i]=((double)((M-s[i]-1)*(s[i]+1))/(double)M)*(double)(M-s[i-1]-1)/(double)((M-s[i]-1)*((s[i]-s[i-1]))); //s[K]=M;
    }
       
}


void 
ComputeHsIext(
    //input variables:
    int *Iext,     // Indices of selection,
    int K,     // Length of the indices, 
    double *h0, // Returning diagonal of H,
    double *h1  // Returning upper diagonal of H
    )
{
    int i,M;
    double iC,iL,iR;
    //double M;
    //M=(double)MM;
    
    M=Iext[K+1];

    //iL=0;iC=(double)Iext[1];iR=(double)s[2];
    for (i=1;i<K;i++){
	//iL=iC;iC=iR;iR=(double)(Iext[i+1]);
	iL=Iext[i-1];iC=Iext[i];iR=(double)Iext[i+1];
	h0[i-1]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));
	h1[i-1]=-sqrt((M-iC)*iC*(M-iR)*iR)/(M*(iR-iC));
    }
    //i=K
    //iL=iC;iC=iR;iR=M;
    iL=Iext[i-1];iC=Iext[i];iR=(double)Iext[i+1];
    h0[i-1]=((M-iC)*iC/M)*(iR-iL)/((iR-iC)*((iC-iL)));       
}



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/// Extract the segment lengths from breakpoint set Iext
void IextToSegLen(
    int *Iext /// (In) Extended breakpoint set
    ,int *SegLen /// (Out) Segment lengths. Can be the same as Iext?
    ,int K      /// Number of Breakpoints. Length of Iext - 1
    )
{
    int k;
    for(k=1;k<=K+1;k++)
	SegLen[k-1]=Iext[k]-Iext[k-1];
}

        
/*************************************************************************/

/// Collapses all the segment amplitudes close to baseline reference level
void 
CollapseAmpTtest(
    double *SegAmp ///Segment Amplitudes (input output)
    ,int *SegLen ///Segment Lengths
    ,int L ///Number of segments.
    ,double BaseAmp ///Reference amplitude to compare (see CompBaseAmpMedianMethod)
    ,double sigma2  ///Estimate of the noise variance. 
    ,double T	    ///Critical value of the T test to decide Gain/Loss against Neutral as null hyphotesis
    )
{
    int k;

//	if(L==1) // If only one segment...
//	{
//		if( fabs(SegAmp[0]-BaseAmp)/sqrt(sigma2/(double)SegLen[0]) < T)
//			SegAmp[0]=BaseAmp;
//	}
    for(k=0;k<L;k++)
	if( fabs(SegAmp[k]-BaseAmp)/sqrt(sigma2/(double)SegLen[k]) < T)
	{
	    //but do it only if one of the neigboring ones have been collapsed
	    if((k>0)&&(fabs(SegAmp[k-1]-BaseAmp)/sqrt(sigma2/(double)SegLen[k-1]) < T))
		SegAmp[k]=BaseAmp;
	    if((k<L-1)&&(fabs(SegAmp[k+1]-BaseAmp)/sqrt(sigma2/(double)SegLen[k+1]) < T))
		SegAmp[k]=BaseAmp;
	    //or it is the initial segment or the final segment of the unit, since we assume that the unseen neighbors where in collapsed state
	    if((k==0)||(k==L-1))
		SegAmp[k]=BaseAmp;
            //or it has larger size than the neighboring segments
	    if((k>0)&&(SegLen[k]>SegLen[k-1]))
		SegAmp[k]=BaseAmp;
	    if((k<L-1)&&(SegLen[k]>SegLen[k+1]))
		SegAmp[k]=BaseAmp;

	}

/* 	
//Classify amplitudes into Gain +1, Loss -1, or Neutral=0
for(k=0;k++;k<K)
{
if(SegAmp[k]>BaseAmp)
SegAmp=+1; //Gain
else if(SegAmp[k]<BaseAmp)
SegAmp=-1;
else
SegAmp=0;		
}
*/
}


/*************************************************************************/


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
    )
{
    int k;

    for(k=0;k<K;k++)
      SegState[k]=SegAmp[k];

    CollapseAmpTtest(SegAmp,SegLen,K,BaseAmp,sigma2,T);
    
    for(k=0;k<K;k++)
      if(SegAmp[k]>BaseAmp)
	SegState[k]=1;
      else if(SegAmp[k]<BaseAmp)
	SegState[k]=-1;
      else
	SegState[k]=0;
       
}




/*************************************************************************/

/// Classify segment amplitudes into a discrete copy number state
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
    )
{
    int k;
    double c;
    double aux;
	
    for(k=0;k<=K;k++)
    {
	if( fabs(SegAmp[k]-BaseAmp)/sqrt(sigma2/(double)SegLen[k]) < T)
	{
	    SegState[k]=ploidy;
	}
	else if((SegAmp[k]-BaseAmp)>0)
	{
	    for(c=ploidy;c<100;c++)
	    {
		aux=(SegAmp[k]-BaseAmp-log2(c)+log2(ploidy))/sqrt(sigma2/(double)SegLen[k]);
		if((SegAmp[k]-BaseAmp-log2(c)+log2(ploidy))/sqrt(sigma2/(double)SegLen[k]) < T )
		{
		    //					c=c-1;
		    break;
		}
	    }
	    SegState[k]=c;
	}
	else if((SegAmp[k]-BaseAmp)<0)
	{
	    for(c=ploidy;c>0;c--)
	    {
		if((SegAmp[k]-BaseAmp-log2(c)+log2(ploidy))/sqrt(sigma2/(double)SegLen[k]) > -T )
		{
		    break;
		}
	    }
	    SegState[k]=c;
	}
    }
}



///Returns the median intensity of all segments to estimate the base-line hybridization intensity associated with the ploidy copy number
double 
CompBaseAmpMedianMethod( 
    int *SegLen    /// Lengths corresponding to the segment amplitudes
    ,double *SegAmp /// Amplitudes of the segments.
    ,int K /// Number of segments. 
    )
{	
    int M,k,RunLen;
    double BaseAmp=0;


    ///\retruns Returns BaseAmp corresponding to the baseline level.
	
    //If they need to be sorted use the following...
    double *D;
    int *I;
    D=myCalloc(K,sizeof(double));
    I=myCalloc(K,sizeof(double));
    for(k=0;k<K;k++)
    {
	D[k]=SegAmp[k];
	I[k]=SegLen[k];
    }
    doubleBubbleSort(D,I,K); //I need indexes of the sort
    SegAmp=D;
    SegLen=I;

    M=0;
    for(k=0;k<K;k++)
	M=M+SegLen[k];
#ifdef _DebugCompBaseAmpMedianMethod_
    myPrintf("_DebugCompBaseAmpMedianMethod_: M%d K%d M/2%d\n",M,K,M/2);
#endif
	
    RunLen=0;
    k=0;
    while(RunLen<M/2)RunLen=RunLen+SegLen[k++];
#ifdef _DebugCompBaseAmpMedianMethod_
    myPrintf("_DebugCompBaseAmpMedianMethod_: k%d RunLen=%d SegAmp[k-1]%g SegAmp[k]%g\n",k,RunLen,SegAmp[k-1],SegAmp[k]);
#endif

    BaseAmp=SegAmp[k-1];

#ifdef _DebugCompBaseAmpMedianMethod_
    printf("_DebugCompBaseAmpMedianMethod_: BaseAmp=%g\n",BaseAmp);
#endif

    return BaseAmp;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/************************************************/
/// Backward Elimination Procedure with T score and Segment Length control
int /// Returns breakpoint list lenght. with T and MinSegLen
BEwTandMinLen( 
    double *Wext   ///IO Breakpoint weights extended notation...
    ,int *Iext     ///IO Breakpoint positions in extended notation...
    ,int *pK       ///IO Number breakpoint positions remaining.
    ,double sigma2 ///IP Sigma2 
    ,double T      ///IP Threshold to prune,  T=T*sqrt(sigma2);
    ,int MinSegLen ///IP Minimum length of the segment.
    )
{
    int i,K,imin,M; 
    //    double myaux;
    double vmin;
    double *tscore; //Statistical Scores,  
    int *L; //Vector with the smallest of the two neighboring segments of each breakpoint.
    int smallinside; //Variable that indicates that there are still small segments to eliminate.

    K=*pK;       //Number of breakpoints
    M=Iext[K+1]; //Total length		
    T=T*sqrt(sigma2); //Adjusting T to the noise power
    
    if(MinSegLen>0)
	smallinside=1; 
    else
	smallinside=0;
        
    L=myCalloc(K+1,sizeof(int)); //Bring from outside?
    tscore=myCalloc(K+1,sizeof(double)); //
    
    //Computing scores
#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("_BEwTandMinLen_ Computing scores\n");
#endif    
    for (i=0;i<K+1;i++)
	tscore[i]=0;    

    ComputeTScores(Wext,Iext,tscore,K,1,K);

#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("\n_BEwTandMinLen_:K%d w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",K,Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]); 
    myPrintf("\n_BEwTandMinLen_:K%d t[0]=%g,t[1]=%g,t[2]=%g,t[K-1]=%g,t[K]=%g\n",K,tscore[0],tscore[1],tscore[2],tscore[K-1],tscore[K]); 
#endif
#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("_BEwTandMinLen_ Backward Elimination sigma2=%g T=%g MinSegLen %d\n",sigma2,T,MinSegLen);
#endif

    BEwTscore(Wext,Iext,tscore,&K,T);

#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("\n_BEwTandMinLen_:K%d w[0]=%g,w[1]=%g,w[2]=%g,w[K-1]=%g,w[K]=%g\n",K,Wext[0],Wext[1],Wext[2],Wext[K-1],Wext[K]); 
    myPrintf("\n_BEwTandMinLen_:K%d t[0]=%g,t[1]=%g,t[2]=%g,t[K-1]=%g,t[K]=%g\n",K,tscore[0],tscore[1],tscore[2],tscore[K-1],tscore[K]); 
#endif


#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("_BEwTandMinLen_ Small segment prunning ML=%d SI=%d \n",MinSegLen,smallinside);
#endif    
    while(smallinside)
    {
	//Compute segment lengths of flanking breakpoints
	for(i=1;i<K+1;i++)
	    L[i]=min(Iext[i]-Iext[i-1],Iext[i+1]-Iext[i]);
        
	//Find smallest segment with lowest score to remove
	imin=-1;
	vmin=1E100;
	for(i=1;i<K+1;i++)
	    if((tscore[i]<vmin)&&(L[i]<MinSegLen)){
			vmin=tscore[i];
			imin=i;
	    }

	if(imin<0)
	    smallinside=0;
	else
        {
#ifdef _DEBUG_BEwTandMinLen_
	    myPrintf("_BEwTandMinLen_ Removing %d %d %g Numrem(old)=%d",Iext[imin],imin,vmin,K);
#endif
	    tscore[imin]=-1;
	    BEwTscore(Wext,Iext,tscore,&K,T);

#ifdef _DEBUG_BEwTandMinLen_
	    myPrintf("_BEwTandMinLen_ K = %d \n",K);
#endif
        }
	if(K<1)smallinside=0;
    }    
#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("T=%g K=%d\n",T,K);
    myPrintf("_BEwTandMinLen_ Small segment prunning ends imin=%d vmin=%g\n",imin,vmin);
#endif    

#ifdef _DEBUG_BEwTandMinLen_
    myPrintf("_BEwTandMinLen_ After BE K=%d \n",K);
#endif
    
    //Freeing memory
    myFree(tscore);
    myFree(L);

    *pK=K;
    return K;    
}


/******************************************************/
/// BEwTscore  Backward Breakpoint Elimination procedure with T score control
int BEwTscore(
    double *Wext,   // IO Breakpoint weights extended notation...
    int *Iext,      // IO Breakpoint positions in extended notation...
    double *tscore, // IO T scores of the breakpoints.
    int *pK,        // IO Number breakpoint positions remaining.
    double T        // IP  Critical value of the minimum Tscore that sets when to stop the BE procedure.
    )
{
    int i,K,jmin,M; 
    double vmin;

    K=*pK;       //Number of breakpoints
    M=Iext[K+1]; //Total length		

#ifdef _DebugBEwTscore_
    myPrintf("_DebugBEwTscore_ BE starts K=%d M=%d T=%g\n",K,M,T);
#endif    


    vmin=0;
    while((vmin<T)&&(K>0))
    {
	//Find breakpoint with lowest score to remove
	jmin=-1;
	vmin=1E100;
	for(i=1;i<K+1;i++)
	    if(tscore[i]<vmin){
		vmin=tscore[i];
		jmin=i;
	    }
#ifdef _DebugBEwTscore_
	myPrintf("_DebugBEwTscore_ Smallest breakpoint imin=%d vmin=%g T=%g\n",jmin,vmin,T);
#endif    
	if(vmin<T) //Remove breakpoint at imin
	{
#ifdef _DebugBEwTscore_
	    myPrintf("_DebugBEwTscore_ Removing imin=%d vmin=%g T=%g\n",jmin,vmin,T);
#endif    

#ifdef _DebugBEwTscore_
	    myPrintf("_DebugBEwTscore_ Removing imin=%d vmin=%g T=%g\n",jmin,vmin,T);
#endif    

#ifdef _DebugBEwTscore_
	    myPrintf("_DebugBEwTscore_ %d, W[jmin-2]=%g W[jmin-1]=%g W[jmin]=%g W[jmin+1]=%g W[jmin+2]=%g T=%g\n",jmin,Wext[jmin-2],Wext[jmin-1],Wext[jmin],Wext[jmin+1],Wext[jmin+2],T);
	    myPrintf("_DebugBEwTscore_ %d, I[jmin-2]=%d I[jmin-1]=%d I[jmin]=%d I[jmin+1]=%d I[jmin+2]=%d T=%g\n",jmin,Iext[jmin-2],Iext[jmin-1],Iext[jmin],Iext[jmin+1],Iext[jmin+2],T);
	    myPrintf("_DebugBEwTscore_ %d, t[jmin-2]=%g t[jmin-1]=%g t[jmin]=%g I[jmin+1]=%g t[jmin+2]=%g T=%g\n",jmin,tscore[jmin-2],tscore[jmin-1],tscore[jmin],tscore[jmin+1],tscore[jmin+2],T);
#endif    
	    for(i=jmin;i<K;i++)
		tscore[i]=tscore[i+1];

	    K=RemoveBreakpoint(Wext,Iext,K,jmin);
	    ComputeTScores(Wext,Iext,tscore,K,jmin-1,jmin);

#ifdef _DebugBEwTscore_
	    myPrintf("_DebugBEwTscore_ %d, W[jmin-2]=%g W[jmin-1]=%g W[jmin]=%g W[jmin+1]=%g W[jmin+2]=%g T=%g\n",jmin,Wext[jmin-2],Wext[jmin-1],Wext[jmin],Wext[jmin+1],Wext[jmin+2],T);
	    myPrintf("_DebugBEwTscore_ %d, I[jmin-2]=%d I[jmin-1]=%d I[jmin]=%d I[jmin+1]=%d I[jmin+2]=%d T=%g\n",jmin,Iext[jmin-2],Iext[jmin-1],Iext[jmin],Iext[jmin+1],Iext[jmin+2],T);
	    myPrintf("_DebugBEwTscore_ %d, t[jmin-2]=%g t[jmin-1]=%g t[jmin]=%g I[jmin+1]=%g t[jmin+2]=%g T=%g\n",jmin,tscore[jmin-2],tscore[jmin-1],tscore[jmin],tscore[jmin+1],tscore[jmin+2],T);
#endif    

	}
    }
#ifdef _DebugBEwTscore_
    myPrintf("_DebugBEwTscore_ BE ends imin=%d vmin=%g T=%g K=%d\n",jmin,vmin,T,K);
#endif  

    *pK=K;
    return K;
}

int
RemoveBreakpoint(
    double *Wext,
    int *Iext,
    int K,
    int jrem
    )
{
    int j;
    double iC,iL,iR,M;

    M=(double)Iext[K+1];
    iL=(double)Iext[jrem-1];
    iC=(double)Iext[jrem];
    iR=(double)Iext[jrem+1];

    //Change coefficients
    if(jrem>1)
	Wext[jrem-1] = Wext[jrem-1] + sqrt((M-iL)/(M-iC)*iL/iC)*(iR-iC)/(iR-iL) * Wext[jrem]; 
    if(jrem<K)
	Wext[jrem+1] = Wext[jrem+1] + sqrt((M-iR)/(M-iC)*iR/iC)*(iC-iL)/(iR-iL) * Wext[jrem]; 
    Wext[jrem]=0; //

    //Shorten list
    for(j=jrem;j<K;j++)
	Wext[j]=Wext[j+1];
    for(j=jrem;j<K+1;j++)
	Iext[j]=Iext[j+1];
    return K-1;
}

void
ComputeTScores(
    const double *Wext,
    const int *Iext,
    double *Scores,
    int K,
    int start,
    int end
    )
{
    int j;
    double h0,M;
	
    M=(double)Iext[K+1];

    for(j=start;j<=end;j++)
    {
	h0=(double)(M-Iext[j]) * (double)Iext[j] / M * (double)(Iext[j+1]-Iext[j-1]) / (double)(Iext[j+1]-Iext[j]) / (double)(Iext[j]-Iext[j-1]);
	Scores[j]=fabs(Wext[j])/sqrt(h0);
    }
}/********* Compute TScores *********/




/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/// Computes the T matrix of the SBL algorithm
void ComputeT(
    double *h0,
    double *h1,
    int M,
    double *alfa, 
    double sigma,
    double *t0,
    double *tl,
    double *tu   
    )
{
    int i=0;
    
    for (i=0;i<M;i++){
        t0[i]=(h0[i]*alfa[i]*sigma)+1;
        if (i<M-1){
            tu[i]=h1[i]*alfa[i+1]*sigma;
            tl[i]=h1[i]*alfa[i]*sigma;
        }
    }
}        


/************************************************/

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
    )
{
    int n,i,j,K;
    int M; // Size of the observed vector y

    double mynorm=0,error=0;

    double *WrB,*WrF,*Wr0;
    double *Wold;

    double *z;  // z=F'y
    double *h0; // H=inv(F(I)'F(I)) matrix diagonals
    double *h1; // 

    double *t0; // T matrix (T=I+sHA)
    double *tl;
    double *tu;
    double *Fd;
    double *Fl;
    double *Fu;
    double *Bd;
    double *Bl;
    double *Bu;

    ///\returns Returns the number of EM iterations


    K=*pK;
    M=Iext[K+1];
    WrF=Wext+1;
    
    if(debug>1)
	myPrintf("## Input K=%d,M=%d\n",K,M);

    // Non discontinuities case
    if (K==0){
	if(debug>0)
	    myPrintf("#(Single_SBL_PWC_norm) Nothing to be done K=0\n");
	return 0; // Nothing needs to be done.
    }


    // Memory needed to perform operations (Needs to be freed before exiting the function)
    z = myCalloc(M-1,sizeof(double));    
    h0=myCalloc(K,sizeof(double));
    h1=myCalloc(K-1,sizeof(double));

    Wold=myCalloc(K,sizeof(double));
    WrB=myCalloc(K,sizeof(double));
    Wr0=myCalloc(K,sizeof(double));

    t0=myCalloc(K,sizeof(double));
    tl=myCalloc(K-1,sizeof(double));
    tu=myCalloc(K-1,sizeof(double));
    Fd=myCalloc(K,sizeof(double));
    Fl=myCalloc(K-1,sizeof(double));
    Fu=myCalloc(K-1,sizeof(double));
    Bd=myCalloc(K,sizeof(double));
    Bl=myCalloc(K-1,sizeof(double));
    Bu=myCalloc(K-1,sizeof(double));


    CompZ(y,z,M); 

    //  myPrintf("\nCOMPUTE H\n");
    ComputeHsIext(Iext,K,h0,h1);
    //myPrintf("H0\n h0[0]:%g\nh0[1]:%g\nh0[2]:%g\nh0[3]:%g\nh0[%d]:%g\n",h0[0],h0[1],h0[2],h0[3],K-1,h0[K]);
    //myPrintf("H1\n h1[0]:%g\nh1[1]:%g\nh1[2]:%g\nh1[3]:%g\nh1[%d]:%g\n",h1[0],h1[1],h1[2],h1[3],K-2,h1[K-2]); 

    TridSymGaxpy(h0,h1,z,K,Wr0);


/********************************/
/******* EM loop         ********/
/********************************/    
    
    for (n=0;n<maxit;n++)
    {
	if(debug>2)
	    myPrintf(".");

	for(j=0;j<K;j++)
	    Wold[j]=WrF[j];
    
	//*** E-STEP ***
	ComputeT(h0,h1,K,alpha,*sigma2,t0,tl,tu);
    
	if (K==1){
	    WrF[0]=Wr0[0]/t0[0]; 
	    Wsigma2[0]=h0[0]*(*sigma2)/(t0[0]);
	}
	else{
	    BackwardTridSolve(tu,t0,tl,K,Bu,Bd,Bl,Wr0,WrB);
	    ForwardTridSolve(tu,t0,tl,K,Fu,Fd,Fl,Wr0,WrF);
	    if(debug>0)  //Error assesment between backward and forward.
	    {
		error=normInfOfDiff(WrB,WrF,K);
		if(error>1E-6)
		    myPrintf("#(Single_SBL_PWC_norm) EM iteration %d, !!!WARNING: The difference between Forward and Backward solution is %g\n",n,error);
	    }

	    TridOfInvTrid(Fu,Fd,Bl,tl,t0,tu,K); //Inverse stored in tl t0 tu
	    DiagOfTridXTrid(tl,t0,tu,h1,h0,h1,Wsigma2,K);
	    for(j=0;j<K;j++)
		Wsigma2[j]=(*sigma2)*Wsigma2[j];
	}
	//*** end E-STEP ***

	//*** M-STEP for the alphas
	for (j=0;j<K;j++)
	    alpha[j]=(1+2*aAlpha)/(WrF[j]*WrF[j]+Wsigma2[j]+2*bAlpha);


	//*** Test if the algorithm has converged
	//euclidean norm of the change in the weights
	//mynorm=sqrt(norm2OfDiff(Wold,WrF,K));
	//mynorm=sqrt(mynorm);  
	mynorm = normInfOfDiff(Wold,WrF,K);
       
	if (mynorm<*tol){ // Converged !!
	    if (debug>0)
		myPrintf("#(Single_SBL_PWC_norm): Converged after %d iterations within tolerance %g, M=%d \n",n,*tol,K);	    
	    break; 
	}


	//***  [optional] -- Elimination of bases     
	i=0;j=0;

	if(isRemovingOK){
	    j=0;    
	    for(i=0;i<K;i++){
		if (alpha[i]<maxAlpha){	    
		    alpha[j]=alpha[i];
		    Iext[j+1]=Iext[i+1];
		    j++;
		}else{
		    //nRem++;
		    if(debug>3)
			myPrintf("##(Single_SBL_PWC_norm) NumIt=%d  Removing %d -- %d breakpoint w=%g alpha=%g\n",n,i+1,Iext[i+1],Wext[i+1],alpha[i]);
		}	  
	    }
	    //      myPrintf("## i %d j %d\n",i,j);
	    Iext[j+1]=Iext[K+1];
	    // j++;
	}else{ // Do not eliminate but prevent overflow of alpha

	    for(i=0;i<K;i++)
		if (alpha[i]>maxAlpha)
		    alpha[i]=maxAlpha;

	    j=i; //Means nothing has been removed
	}

	if(j<i){ //Something has been removed -- So I need to rebuild H and W0
	    K=j;
	    if(K>0){
		//Recompute H,W0	
		if(debug>1)
		    myPrintf("##(Single_SBL_PWC_norm) NumIt=%d  Removed %d bases, new K=%d\n",n,i-j,K);

		ComputeHsIext(Iext,K,h0,h1);

		for(j=1;j<=K;j++)
		    t0[j-1]=z[Iext[j]-1];

		TridSymGaxpy(h0,h1,t0,K,Wr0);

		for(j=0;j<K;j++)
		    Wold[j]=Wr0[j];

	    }else{
		myPrintf("#(Single_SBL_PWC_norm): After %d iterations, %d breakpoints found M=%d \n",n,K,M);      
        break;
	    }
	}

	//*** End Optional elimination of bases       

    }/*** End EM loop ***/
    /*********************/
 
    if (debug>0)
	if (n>=maxit)              
	    myPrintf("#      SBL: Converged??? Stopped after %d iterations with change %g, K=%d, M=%d\n",n,mynorm,K,M);       
       
    if (debug>0)
	myPrintf("#(Single_SBL_PWC_norm) Reffiting K=%d\n",K);
       
    // Projection onto Breakpoint set... (REFITTING)
    if ((K>0)&&(isRemovingOK)) //NOTE: Add another flag? to select this step
	for (i=0;i<K;i++)
	    WrF[i]=Wr0[i];
    
    //Memory freeing....
    myFree(z);

    myFree(h0);
    myFree(h1);

    myFree(Wold);
    myFree(WrB);
    myFree(Wr0);

    myFree(t0);
    myFree(tl);
    myFree(tu);
    myFree(Fd);
    myFree(Fl);
    myFree(Fu);
    myFree(Bd);
    myFree(Bl);
    myFree(Bu);
    
    
    *pK=K;
    
    *tol=mynorm;
    return n;
		
}//Single_SBL_PWC_norm
/*****************************************/

int SBLandBE( //Returns breakpoint list lenght.
    double *y, 
    int M,  //length of the noisy signal tn
    double *Psigma2, //If sigma2 < 0, compute sigma2
    double a,      // SBL parameter
    double T,      // Threshold to prune
    int MinSegLen, //Minimum length of the segment.
    int **pIext,   //Returns breakpoint positions
    double **pWext //Returns breakpoint weights.
    //int *pK    
    )
{
    int i,K;
    int debug=0; //I should convert that into a parameter
    double myaux;
    double sigma2;
    double ymean;
    double tol,maxalpha,b;
    double *alpha,*Wext,*tn,*aux;
    int *Iext;
    int maxit;
    int NumEMsteps;

    sigma2=*Psigma2;
		
    tn=myCalloc(M,sizeof(double));
    for(i=0;i<M;i++)
	tn[i]=y[i];
    
    //If sigma2 < 0, compute sigma2
    if(sigma2<0){
        sigma2=0;
        for(i=1;i<M;i++){
            myaux=tn[i]-tn[i-1];
            sigma2+=(0.5*myaux*myaux);
        }
        sigma2=sigma2/(M-1);
    }      
    
    //Mean removal
    ymean=0;
    for(i=0;i<M;++i)
        ymean+=tn[i];
    ymean=ymean/M;
    for(i=0;i<M;++i)
        tn[i]=tn[i]-ymean;
    
    // SBL optimization parameters
    tol=1E-10;     //1E-10 or 1E-8 seems to work well for this parameter. -- => ++ conv time
    maxalpha=1E8;  //1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
    maxit=10000;   //Maximum number of iterations to reach convergence...
    b=1E-20;       //
    
    //Call to SBL 
#ifdef _DEBUG_SBLBE_
    myPrintf("#(SBLandBE) Memory initialization\n");
#endif

    K=M-1;
    Iext=myCalloc(K+2,sizeof(int)); //
    Wext=myCalloc(K+1,sizeof(double)); //
    alpha=myCalloc(K,sizeof(double));   
    aux=myCalloc(K,sizeof(double));

    
#ifdef _DEBUG_SBLBE_
    myPrintf("#(SBLandBE) Breakpoint Initialization\n");
#endif
    //Initialize breakpoints
    for (i=0;i<K;i++)
        alpha[i]=0.0;
    for (i=0;i<=K;i++)
        Wext[i]=0.0;
    for (i=0;i<=K+1;i++)
        Iext[i]=i;
    
#ifdef _DEBUG_SBLBE_
    myPrintf("#(SBLandBE) SBL starts\n");
#endif
    //    NumEMsteps=SBL(tn,Iext,alpha,Wext+1,aux,M,&K,sigma2,a,b,maxalpha,maxit,tol,1);
    NumEMsteps=Single_SBL_PWC_norm(tn,Iext,alpha,Wext,aux,&K,&sigma2,a,b,maxalpha,1,maxit,&tol,debug);

    //Freeing memory
    myFree(alpha);
    myFree(tn);
    myFree(aux);

	//mxAssert(Iext[K+1]==M,"after SBL");

    //Convert Iext and Wext to the extended notation.
    Iext=myRealloc(Iext,(K+2)*sizeof(int));
    //for(i=(K+1);i>0;i--)
    //	Iext[i]=Iext[i-1]+1;
    //Iext[0]=0;
    //Iext[K+1]=M;

    Wext=myRealloc(Wext,(K+1)*sizeof(double));
    //for(i=K;i>0;i++)
    //	Wext[i]=Wext[i-1];
    Wext[0]=ymean;
    
    if(debug>0)
	myPrintf("#(SBLandBE) Backward Elimination T=%g MinLen=%d\n",T,MinSegLen);
    BEwTandMinLen(Wext,Iext,&K,sigma2,T,MinSegLen); 
   //mxAssert(Iext[K+1]==M,"after BE");


    Iext=myRealloc(Iext,(K+2)*sizeof(int));
    Wext=myRealloc(Wext,(K+1)*sizeof(double));

    if(debug>0)      
        myPrintf("#(SBLandBE) After BE K=%d \n",K);

    
    
    *pIext=Iext;
    *pWext=Wext;
    *Psigma2=sigma2;

    return K;    
}//SBLandBE


/*****************************************/

#ifdef _DEPRECATED_

/** deprecated simpletresholding 
 * Applies a simple tresholding algorithm to find the discontinuities
 */
int   //Returns the number of discontinuities
simpletresholding(
    //Input variables:
    double *inputvector, //1D Array containing the input vector
    int N,  //Vector length
    double thres, //Treshold value
    //Output variables:
    double *disc //1D empty array, with memory already allocated for finding up to N discontinuities positions
    )
{
    int i=0;
    int numdisc=0; //Number of discontinuities
    int state=0; //To keep up the state 0(<tresh),1(>tresh)
    
    //First step: Find discontinuities
    if(inputvector[0]>=thres) state=1;    
    for(i=0;i<N;i++){
        if (state==0)
        {
            if (inputvector[i]>=thres) 
            {
                disc[numdisc]=(double)(i+1); //Matlab uses 1..N indices instead of 0..(N-1)
                numdisc++;
                state=1;
//                myPrintf("i%d numdisc%d State change state%d.\n",i,numdisc,state);
            }
        }
        else
        {
            if (inputvector[i]<thres) 
            {
                disc[numdisc]=(double)(i+1); //Matlab uses 1..N indices instead of 0..(N-1)
                numdisc++;
                state=0;
		//              myPrintf("i%d numdisc%d State change state%d.\n",i,numdisc,state);
            }            
        }
    }    
    
    return numdisc;
}


/** Deprecated computesegmentmeans, Convert this function to IextYobsToSegAmp to obtain projection in a different way. 
    Using 2N flops.
    *  Computed segment means.
    */

void
computesegmentmeans(
    //Input variables:
    double *inputvector, //1D Array containing the input vector
    int N,  //Vector length
    double *disc, //1D empty array, with memory already allocated for finding up to N discontinuities positions
    int numdisc, //Length of the discontinuity vector
    //Output variables:    
    double *amp  //1D empty array, containing the amplitudes of the segments between discotinuities, 
    )
{
    int i,j;
    int prevdisc=0;

    for(i=0;i<numdisc;i++)
    {
        amp[i]=0;
        for(j=prevdisc;j<(int)(disc[i]-1);j++)
            amp[i]+=inputvector[j];
        amp[i]=amp[i]/((double)(disc[i]-1-prevdisc));
        prevdisc=(int)(disc[i]-1);
    }
    amp[i]=0;
    for(j=prevdisc;j<N;j++)
        amp[i]+=inputvector[j];
    amp[i]=amp[i]/(double)(N-prevdisc);   
}

/* reconstructoutput 
 * Reconstruct signal from the discontinuity location, and segment 
 * amplitudes
 */
void
reconstructoutput(
    //Output variables:
    double *rec, //1D Array returning the reconstructed vector
    //Input variable:
    int N,  //Vector length
    double *disc, //1D empty array, with memory already allocated for finding up to N discontinuities positions
    int numdisc, //Length of the discontinuity vector
    //Output variables:    
    double *amp  //1D empty array, containing the amplitudes of the segments between discotinuities, 
    )
{
    int i,j;
    
    j=0;
    for(i=0;i<numdisc;i++)
        for(;j<(int)(disc[i]-1);j++)
            rec[j]=amp[i];
    //If no discontinuities or from the last to the end 
    for(;j<N;j++)
        rec[j]=amp[i];    
}

#endif //_DEPRECATED_
