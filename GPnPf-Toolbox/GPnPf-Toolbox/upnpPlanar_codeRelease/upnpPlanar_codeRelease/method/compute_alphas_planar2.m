function Alph=compute_alphas_planar2(Xw,Cw)
%clear all; close all; load compute_alphas_planar2;

n=size(Xw,1); %number of 3d points

C=[Cw';ones(1,3)];
X=[Xw';ones(1,n)];
Alph_=inv(C'*C)*C'*X;

Alph=Alph_';





