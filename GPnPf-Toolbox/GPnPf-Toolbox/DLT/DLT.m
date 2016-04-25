function [f R t] = DLT(XX,xx)
%direct linear transformation for projection matrix P
P = CalibNormDLT(xx,XX);

%decompose P into K, R, t
[K,R,t] = DecompPMat(P);

if det(R) < 0;
    R = -R;
end

%focal length
f = (K(1,1)+K(2,2))/2;
end