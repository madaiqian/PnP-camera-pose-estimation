function Alph=compute_alphas(Xw,Cw)

n=size(Xw,1); %number of 3d points

%generate auxiliar matrix to compute alphas
C=[Cw';ones(1,4)];
X=[Xw';ones(1,n)];
Alph_=inv(C)*X;


Alph=Alph_';
