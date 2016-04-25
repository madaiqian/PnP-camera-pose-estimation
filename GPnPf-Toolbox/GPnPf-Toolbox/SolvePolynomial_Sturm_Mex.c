/*
Using Sturm Sequences to Bracket Real Roots of Polynomial Equations
by D.G. Hook and P.R. McAree
from "Graphics Gems", Academic Press, 1990
*/

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "solve.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
double *Srcc;
poly	sseq[MAX_ORDER];
double 	min, max, roots[MAX_ORDER];
int		i, j, order, nroots, nchanges, np, atmin, atmax;
mwSize *q1dim;
int dims[2];
double *sL, *sU;

if (nrhs!=3)
{
   mexErrMsgTxt("Input one coefficient vector and lower/upper bounds\n");
}

q1dim = mxGetDimensions(prhs[0]);
order = q1dim[0]-1;
if (order < q1dim[1]-1)
{
    order = q1dim[1]-1;
}
   
if (order >= MAX_ORDER)
{
   mexPrintf("Exceeds the maximum order!\n");
   return;
}

Srcc = mxGetPr(prhs[0]);
sL = mxGetPr(prhs[1]); 
sU = mxGetPr(prhs[2]);

for (i = order; i >= 0; i--) 
{
    sseq[0].coef[i] = Srcc[order-i];
}

/*
 * build the Sturm sequence
 */
np = buildsturm(order, sseq);

/* 
 * get the number of real roots
 */
nroots = numroots(np, sseq, &atmin, &atmax);

/*
 * calculate the bracket that the roots live in
 */
min = sL[0];
nchanges = numchanges(np, sseq, min);

if (nchanges != atmin) {
        atmin = nchanges;
}

max = sU[0];
nchanges = numchanges(np, sseq, max);

if (nchanges != atmax) {
        atmax = nchanges;
}

nroots = atmin - atmax;

/*
 * perform the bisection.
 */
sbisect(np, sseq, min, max, atmin, atmax, roots);

dims[0] = nroots;
dims[1] = 1; 
plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
memcpy (mxGetData(plhs[0]), roots, nroots*sizeof(double));
}
