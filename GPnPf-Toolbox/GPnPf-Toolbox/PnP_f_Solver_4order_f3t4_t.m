% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [ftemp ttemp] = PnP_f_Solver_4order_f3t4_t(a,b)
  %nn = [ f^3*t1^4, f^2*t1^4, f*t1^4, t1^4, f^3*t1^3, f^2*t1^3, f*t1^3, t1^3, f^3*t1^2, f^2*t1^2, f*t1^2, t1^2, f^3*t1, f^2*t1, f*t1, t1, f^3, f^2, f, 1];
    a1 =  a(1);  a2 =  a(2);  a3 =  a(3);  a4 =  a(4);  a5 =  a(5);  a6 =  a(6);  a7 =  a(7);  a8 =  a(8);  a9 =  a(9);  a10 = a(10);
    a11 = a(11); a12 = a(12); a13 = a(13); a14 = a(14); a15 = a(15); a16 = a(16); a17 = a(17); a18 = a(18); a19 = a(19); a20 = a(20);
       
    b1 =  b(1);  b2 =  b(2);  b3 =  b(3);  b4 =  b(4);  b5 =  b(5);  b6 =  b(6);  b7 =  b(7);  b8 =  b(8);  b9 =  b(9);  b10 = b(10);
    b11 = b(11); b12 = b(12); b13 = b(13); b14 = b(14); b15 = b(15); b16 = b(16); b17 = b(17); b18 = b(18); b19 = b(19); b20 = b(20);
 
Mat0 = [[ a17, a18, a19, a20,   0,   0]
[   0, a17, a18, a19, a20,   0]
[   0,   0, a17, a18, a19, a20]
[ b17, b18, b19, b20,   0,   0]
[   0, b17, b18, b19, b20,   0]
[   0,   0, b17, b18, b19, b20]];
 
 
Mat1 = [[ a13, a14, a15, a16,   0,   0]
[   0, a13, a14, a15, a16,   0]
[   0,   0, a13, a14, a15, a16]
[ b13, b14, b15, b16,   0,   0]
[   0, b13, b14, b15, b16,   0]
[   0,   0, b13, b14, b15, b16]];
 
 
Mat2 = [[ a9, a10, a11, a12,   0,   0]
[  0,  a9, a10, a11, a12,   0]
[  0,   0,  a9, a10, a11, a12]
[ b9, b10, b11, b12,   0,   0]
[  0,  b9, b10, b11, b12,   0]
[  0,   0,  b9, b10, b11, b12]];
 
 
Mat3 = [[ a5, a6, a7, a8,  0,  0]
[  0, a5, a6, a7, a8,  0]
[  0,  0, a5, a6, a7, a8]
[ b5, b6, b7, b8,  0,  0]
[  0, b5, b6, b7, b8,  0]
[  0,  0, b5, b6, b7, b8]];
 
 
Mat4 = [[ a1, a2, a3, a4,  0,  0]
[  0, a1, a2, a3, a4,  0]
[  0,  0, a1, a2, a3, a4]
[ b1, b2, b3, b4,  0,  0]
[  0, b1, b2, b3, b4,  0]
[  0,  0, b1, b2, b3, b4] ];
 
[X E] = polyeig(Mat0,Mat1,Mat2,Mat3,Mat4);

I = find(abs(imag(E))<1e-6);
ttemp = real(E(I));
X = real(X(:,I));
X = X(1:end-1,:)./X(2:end,:);
ftemp = median(X,1).';

end
