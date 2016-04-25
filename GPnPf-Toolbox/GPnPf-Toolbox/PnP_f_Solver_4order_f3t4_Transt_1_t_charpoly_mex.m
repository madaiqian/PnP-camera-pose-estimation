%A characteristic polynomial method to solve the p4p+f problem
%1. closed-form coefficients; 2. sturm-sequence for eigenvalue (t); 
%3. closed-form solution for f
function [ftemp ttemp] = PnP_f_Solver_4order_f3t4_Transt_1_t_charpoly_mex(a,b)
%calculate characteristic polynomial coefficient
coeff = CalculateCoefficient_f3t4_Transt_1_Mex(a,b);

%solve t1 using the sturm sequence (-0.99<t1<0.99)
ttemp = SolvePolynomial_Sturm_Mex(coeff,-0.99, 0.99);

%solve f via close-form solution (0.1<f<100, f is normalized)
[ftemp ttemp] = Solvef_f3t4_Transt_1_Mex(ttemp,a,b,0.1,100.0);

end
