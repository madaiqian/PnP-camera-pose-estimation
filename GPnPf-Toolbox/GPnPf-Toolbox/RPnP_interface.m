function [f R t] = RPnP_interface(XX,xx,K)
xx = K\[xx;ones(1,size(xx,2))];
xx = xx(1:2,:)./repmat(xx(3,:),2,1);

[R t]= RPnP_revise(XX,xx);

f = K(1,1)*ones(1,size(t,2));
end