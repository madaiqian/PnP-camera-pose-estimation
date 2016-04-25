function [ftemp ttemp] = PnP_f_Solver_4order_f3t4_Transt_1_t(a,b)
  %nn = [ f^3*t1^4, f^2*t1^4, f*t1^4, t1^4, f^3*t1^3, f^2*t1^3, f*t1^3, t1^3, f^3*t1^2, f^2*t1^2, f*t1^2, t1^2, f^3*t1, f^2*t1, f*t1, t1, f^3, f^2, f, 1];
    a1 =  a(1);  a2 =  a(2);  a3 =  a(3);  a4 =  a(4);  a5 =  a(5);  a6 =  a(6);  a7 =  a(7);  a8 =  a(8);  a9 =  a(9);  a10 = a(10);
    a11 = a(11); a12 = a(12); a13 = a(13); a14 = a(14); 
       
    b1 =  b(1);  b2 =  b(2);  b3 =  b(3);  b4 =  b(4);  b5 =  b(5);  b6 =  b(6);  b7 =  b(7);  b8 =  b(8);  b9 =  b(9);  b10 = b(10);
    b11 = b(11); b12 = b(12); b13 = b(13); b14 = b(14); 
 
Mat0 = [[ 0, 0, a13, a14,   0,   0]
[ 0, 0,   0, a13, a14,   0]
[ 0, 0,   0,   0, a13, a14]
[ 0, 0, b13, b14,   0,   0]
[ 0, 0,   0, b13, b14,   0]
[ 0, 0,   0,   0, b13, b14]];
 
 
Mat1 = [[ 0, 0, a11, a12,   0,   0]
[ 0, 0,   0, a11, a12,   0]
[ 0, 0,   0,   0, a11, a12]
[ 0, 0, b11, b12,   0,   0]
[ 0, 0,   0, b11, b12,   0]
[ 0, 0,   0,   0, b11, b12]];
 
 
Mat2 = [[ 0, a8, a9, a10,   0,   0]
[ 0,  0, a8,  a9, a10,   0]
[ 0,  0,  0,  a8,  a9, a10]
[ 0, b8, b9, b10,   0,   0]
[ 0,  0, b8,  b9, b10,   0]
[ 0,  0,  0,  b8,  b9, b10]];
 
 
Mat3 = [[ 0, a5, a6, a7,  0,  0]
[ 0,  0, a5, a6, a7,  0]
[ 0,  0,  0, a5, a6, a7]
[ 0, b5, b6, b7,  0,  0]
[ 0,  0, b5, b6, b7,  0]
[ 0,  0,  0, b5, b6, b7]];
 
 
Mat4 = [[ a1, a2, a3, a4,  0,  0]
[  0, a1, a2, a3, a4,  0]
[  0,  0, a1, a2, a3, a4]
[ b1, b2, b3, b4,  0,  0]
[  0, b1, b2, b3, b4,  0]
[  0,  0, b1, b2, b3, b4]];
 
[X E] = polyeig(Mat0,Mat1,Mat2,Mat3,Mat4);
% Mat{1} = Mat1; Mat{2} = Mat2; Mat{3} = Mat3; Mat{4} = Mat4;
% A = eye(6*4);
% A(1:6,1:6) = Mat0;
% B = diag(ones(6*(4-1),1),-6);
% j = 1:6;
% for k = 1:4
%     B(1:6,j) = - Mat{k};
%     j = j+6;
% end
% [YYY XXX] = eig(A,B);
I = find(abs(imag(E))<1e-6);
ttemp = real(E(I));
X = real(X(:,I));
X = X(1:end-1,:)./X(2:end,:);
ftemp = median(X,1).';

end
