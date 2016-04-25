clc;
clear;
A=load('1.in');
B=load('2.in');
fgt=1436.694355;

error=zeros(6,1);
ferror=[1e-5,1e-4,1e-3,1e-2,1e-1,1,10000];
dx=[-5,-4,-3,-2,-1,0];

[m,n]=size(A);
tol=0;
for i=1:m
    pp=0;
    t=0;
    for j=A(i,1):A(i,2)
        pp=pp+1;
        m(1,pp)=B(j,1);
        m(2,pp)=B(j,2);
        M(1,pp)=B(j,3);
        M(2,pp)=B(j,4);
        M(3,pp)=B(j,5);
        if pp==4
            [min]=p4pf_func(M,m,fgt);
            tol=tol+1;
            ans(tol)=min;
            pp=0;
            for k=1:6
                if min>=ferror(k)&&min<=ferror(k+1)
                    error(k)=error(k)+1;
                end
            end
            t=t+1;
        end
        if t==10
            break;
        end
    end
end

dy=error/1410;

plot(dx,dy,'-*b');