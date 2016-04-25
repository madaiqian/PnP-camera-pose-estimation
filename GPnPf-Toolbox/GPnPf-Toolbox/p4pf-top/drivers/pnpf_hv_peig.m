% calc focal length, rotation and translation such that
% lambda*m = [f 0 0;0 f 0;0 0 1]*[R t]*M
%
% by Martin Bujnak, nov2007
%
% m - 2*N vector (at most first 7 points are used)
% M - 3*N vector with 3D correspondence to 'm'

function [f] = pnpf_hv_peig(m, M)

    [f10 f11 f12 f13 f14] = BuildP3Pfpoly_coefs(m(:,1), m(:,2), m(:,3), M(:,1), M(:,2), M(:,3));
    [f20 f21 f22 f23 f24] = BuildP3Pfpoly_coefs(m(:,1), m(:,4), m(:,3), M(:,1), M(:,4), M(:,3));

    % 4pt
    F0 = [0 0 0 f10
          0 0 0 f20
          0 0 f10 0 
          0 0 f20 0 
          0 f10 0 0
          0 f20 0 0
          f10 0 0 0
          f20 0 0 0];

    F1 = [0 0 0 f11
          0 0 0 f21
          0 0 f11 0 
          0 0 f21 0 
          0 f11 0 0
          0 f21 0 0
          f11 0 0 0
          f21 0 0 0];

    F2 = [0 0 0 f12
          0 0 0 f22
          0 0 f12 0 
          0 0 f22 0 
          0 f12 0 0
          0 f22 0 0
          f12 0 0 0
          f22 0 0 0];

    F3 = [0 0 0 f13
          0 0 0 f23
          0 0 f13 0 
          0 0 f23 0 
          0 f13 0 0
          0 f23 0 0
          f13 0 0 0
          f23 0 0 0];

    F4 = [0 0 0 f14
          0 0 0 f24
          0 0 f14 0 
          0 0 f24 0 
          0 f14 0 0
          0 f24 0 0
          f14 0 0 0
          f24 0 0 0];
      
    f = polyeig(F0, F1, F2, F3, F4);
    f = sqrt(f(find((f > 0) & abs(imag(f)) < 1e-5)));
end