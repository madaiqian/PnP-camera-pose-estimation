function [f R t] = upnp_planar_interface(XX,xx,K)

if nargin < 3 
    K = diag([0,0,1]);
end

[R,t,Xc,f] = upnp_planar(XX.',xx.',K);
end