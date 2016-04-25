#include "mex.h";
#include <math.h>

void SolveCubic_Mex(double *Srcc, double  *Dsts)
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

double *Srcc;
double *Dsts;

if (nrhs!=1)
{
   mexErrMsgTxt("One coefficient vector required.");
}

Srcc= mxGetPr(prhs[0]);

plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);

Dsts = mxGetPr(plhs[0]);

SolveCubic_Mex(Srcc, Dsts);

}

