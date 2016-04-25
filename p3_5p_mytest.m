addpath('GPnPf-Toolbox/GPnPf-Toolbox');
addpath('p4pf');
A=load('1.in');
B=load('2.in');
fgt=1436.694355;
f=fgt;

error=zeros(6,1);
ferror=[1e-5,1e-4,1e-3,1e-2,1e-1,1,10000];
dx=[-5,-4,-3,-2,-1,0];

[m,n]=size(A);
tol=0;
for i=1:m
    pp=0;
    i
    t=0;
    for j=A(i,1):A(i,2)
        pp=pp+1;            %每四个点为一组
        m(1,pp)=B(j,1);
        m(2,pp)=B(j,2);
        M(1,pp)=B(j,3);
        M(2,pp)=B(j,4);
        M(3,pp)=B(j,5);
        if pp==4
            [min_p35p]=p3_5p(M,m,fgt);
            [min_Bujnak]=p4pf_func(M,m,fgt);
            [f_R, f_A, f_AT ,f_ATc] = RPnP_f_Test(M,m);
            tol=tol+1;     %组号
            ans(tol,1)=min_p35p;
            ans(tol,2)=min_Bujnak;  
            if length(f_R)>0 && length(f_A)>0 && length(f_AT)>0 && length(f_ATc)>0
               ans(tol,3) = min(abs(f_R-f)/f);
               ans(tol,4) = min(abs(f_A-f)/f);
               ans(tol,5) = min(abs(f_AT-f)/f);
               ans(tol,6) = min(abs(f_ATc-f)/f);
            end
            pp=0;
            t=t+1;
        end
        if t==1
            break;
        end
    end
end