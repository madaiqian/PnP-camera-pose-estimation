function [f] = Best_Ratio_HVR(M3D,m2D)

% shift 3D data so that variance = sqrt(2), mean = 0
mean3d = (sum(M3D,2) / 4);
M3D = M3D - repmat(mean3d, 1, 4);

% variance (isotropic)
var = (sum( sqrt(sum( M3D.^2 ) ) ) / 4);
M3D = (1/var)*M3D;

% scale 2D data
var2d = (sum( sqrt(sum( m2D.^2 ) ) ) / 4);
m2D = (1/var2d)*m2D;

%call the solver
[f] = pnpf_hv_peig(m2D, M3D);

%denormalize
f=var2d*f;
end