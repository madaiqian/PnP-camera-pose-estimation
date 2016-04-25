function [f R t] = upnp_GN_interface(XX,xx,K)

if nargin < 3 
    K = diag([0,0,1]);
end

Number_eq_sets = 500;   % Number of equations sets to explore
Explore_N3 = true;      % Set this to 'false' to avoid exploration of Kernel dimensionality N=3

[R,t,Xc,f] = upnp_GN(XX.',xx.',K,Number_eq_sets,Explore_N3);
end