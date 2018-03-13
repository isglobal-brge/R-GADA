#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Rcall_Single_SBL_PWC_norm(double *
                                      ,int *
                                        ,double *
                                        ,double *
                                        ,double *
                                        ,int *
                                        ,double *
                                        ,double *
                                        ,double *
                                        ,double *
                                        ,int *
                                        ,int *
                                        ,double *
                                        ,int *);
extern void RcallBEwTandMinLen(double *,
                               int *,
                               int *,
                               const double *,
                               const double *,
                               const int *);
extern void RcallClassifySegments(int *
                                    ,double *
                                    ,double *
                                    ,int *
                                    ,double *
                                    ,double *
                                    ,double *);
extern void RcallCompAmpMedianMethod(int *
                                       ,double *
                                       ,int *
                                       ,double *);
extern void RcallWextIextToSegments(double *,
                                    int *,
                                    int *,
                                    double *,
                                    int *);

/* .Fortran calls */
extern void F77_NAME(floc)(double *, double *, double *, double *, 
                     int *, int *, int *, int *, double *);

static const R_CMethodDef CEntries[] = {
  {"Rcall_Single_SBL_PWC_norm", (DL_FUNC) &Rcall_Single_SBL_PWC_norm, 14},
  {"RcallBEwTandMinLen",        (DL_FUNC) &RcallBEwTandMinLen,         6},
  {"RcallClassifySegments",     (DL_FUNC) &RcallClassifySegments,      7},
  {"RcallCompAmpMedianMethod",  (DL_FUNC) &RcallCompAmpMedianMethod,   4},
  {"RcallWextIextToSegments",   (DL_FUNC) &RcallWextIextToSegments,    5},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {"floc", (DL_FUNC) &F77_NAME(floc), 9},
  {NULL, NULL, 0}
};

void R_init_gada(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}