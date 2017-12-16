/*=================================================================
 * consoledefines.h 
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
#ifndef _RGADADEFINES_H_
#define _RGADADEFINES_H_

//Matlab printing and memory allocation 
#include <R.h>
//#include <malloc.h>
//#include <RS.h>  //R memory allocation via Calloc and Free
//#include <Rmath.h>

#define myPrintf Rprintf //Syntax of myPrintf similar to that of fprintf('asfd',var,var)
//#define myCalloc(n,size) R_alloc(n, size)
//#define myFree(ptr) ptr=ptr //do nothing, does not free memory, leave to Garbage_collection???
#define myCalloc(n,size) calloc(n , size)
#define myRealloc(n,size) realloc(n,size)
#define myFree(ptr) free(ptr) 

#define log2(x) log(x)/log(2)


//#define _DEBUG_BEwTandMinLen_
//#define _DebugBEwTscore_

#endif //_RGADADEFINES_H_
