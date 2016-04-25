#include "mex.h";
#include <math.h>

void Solvecubic(double *Srcc, double *Dsts);

int Solvef_f3t4_Transt_1_Mex(int num, double *sL, double *sU, double *Srct, double *Srca, double *Srcb, double  *outf, double *outt)
{
    int i,j; 
    int numf;
    double t, t2, t3, t4;
    double eqaf[4], eqbf[4];
    double f[3];
    double minv, tempv, tempf;
    
    numf = 0;
    for (i=0;i<num;i++)
    {
        t = Srct[i];
        t2 = pow(t,2.0);
        t3 = pow(t,3.0);
        t4 = pow(t,4.0);
        eqaf[0] = Srca[0]*t4; 
        eqaf[1] = Srca[1]*t4 + Srca[4]*t3 + Srca[7]*t2;
        eqaf[2] = Srca[2]*t4 + Srca[5]*t3 + Srca[8]*t2 + Srca[10]*t + Srca[12];
        eqaf[3] = Srca[3]*t4 + Srca[6]*t3 + Srca[9]*t2 + Srca[11]*t + Srca[13];
        
        eqbf[0] = Srcb[0]*t4; 
        eqbf[1] = Srcb[1]*t4 + Srcb[4]*t3 + Srcb[7]*t2;
        eqbf[2] = Srcb[2]*t4 + Srcb[5]*t3 + Srcb[8]*t2 + Srcb[10]*t + Srcb[12];
        eqbf[3] = Srcb[3]*t4 + Srcb[6]*t3 + Srcb[9]*t2 + Srcb[11]*t + Srcb[13];
        
        Solvecubic(eqaf, f);
        minv = 1e20;
        for (j=0;j<3;j++)
        {
            tempv = fabs(eqbf[0]*pow(f[j],3)+eqbf[1]*pow(f[j],2)+eqbf[2]*f[j]+eqbf[3]);
            if (tempv < minv)
            {
                tempf = f[j];
                minv = tempv;
            }
        }
        
        if (tempf >= sL[0] && tempf <= sU[0])
        {
            outf[numf] = tempf;
            outt[numf] = t+1; //due to translation
            numf = numf+1;
        }
 
    }
    return numf;
}

void Solvecubic(double *Srcc, double *Dsts)
{
int	i, num;
double	sub,
	A, B, C,
	sq_A, p, q,
	cb_p, D;
double u, phi, t, sqrt_D;

double M_PI = 3.14159265358979323846;

// normalize the equation:x ^ 3 + Ax ^ 2 + Bx  + C = 0
A = Srcc[1] / Srcc[0];
B = Srcc[2] / Srcc[0];
C = Srcc[3] / Srcc[0];

// substitute x = y - A / 3 to eliminate the quadric term: x^3 + px + q = 0

sq_A = A * A;
p = 1.0/3.0 * (-1.0/3.0 * sq_A + B);
q = 1.0/2.0 * (2.0/27.0 * A *sq_A - 1.0/3.0 * A * B + C);

// use Cardano's formula

cb_p = p * p * p;
D = q * q + cb_p;

if (isZero(D))
    {
    if (isZero(q))
	{
	// one triple solution
	Dsts[0] = 0.0;
	num = 1;
	}
    else
	{
	// one single and one double solution
	u = pow(-q,1.0/3.0);
	Dsts[0] = 2.0 * u;
	Dsts[1] = - u;
	num = 2;
	}
    }
else
    if (D < 0.0)
	{
	// casus irreductibilis: three real solutions
	phi = 1.0/3.0 * acos(-q / sqrt(-cb_p));
    t = 2.0 * sqrt(-p);
	Dsts[0] = t * cos(phi);
	Dsts[1] = -t * cos(phi + M_PI / 3.0);
	Dsts[2] = -t * cos(phi - M_PI / 3.0);
	num = 3;
	}
    else
	{
	// one real solution
	sqrt_D = sqrt(D);
	u = pow(sqrt_D + fabs(q),1.0/3.0);
	if (q > 0.0)
	    Dsts[0] = - u + p / u ;
	else
	    Dsts[0] = u - p / u;
	num = 1;
	}

// resubstitute
sub = 1.0 / 3.0 * A;
for (i = 0; i < num; i++)
    Dsts[i] -= sub;
}

static int isZero(double x)
{
return x > -1e-10 && x < 1e-10;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

double *Srct;
double *Srca;
double *Srcb;
double *sL;
double *sU;
double *Dstf;
double *Dstt;
int num, numf; 
mwSize *q1dim;
int dims[2];
double f[14];
double t[14];

if (nrhs!=5)
{
   mexErrMsgTxt("Five input vectors.");
}

q1dim = mxGetDimensions(prhs[0]);
num = q1dim[0];
if (num < q1dim[1])
{
    num = q1dim[1];
}

Srct= mxGetPr(prhs[0]);
Srca= mxGetPr(prhs[1]);
Srcb= mxGetPr(prhs[2]);
sL= mxGetPr(prhs[3]);
sU= mxGetPr(prhs[4]);

numf = Solvef_f3t4_Transt_1_Mex(num, sL, sU, Srct, Srca, Srcb, f, t);

dims[0] = numf;
dims[1] = 1;

plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

memcpy(mxGetData(plhs[0]), f, numf*sizeof(double));
memcpy(mxGetData(plhs[1]), t, numf*sizeof(double));
}

