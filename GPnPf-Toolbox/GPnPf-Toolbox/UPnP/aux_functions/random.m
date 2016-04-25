function [x]=random(xini,xend,n)

a=xend-xini;
b=xini;

xu=rand(n,1);
x=xu*a+b;



