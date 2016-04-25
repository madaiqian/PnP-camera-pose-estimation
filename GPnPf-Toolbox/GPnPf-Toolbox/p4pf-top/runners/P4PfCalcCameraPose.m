% create camera matrix based on depths and 2D-3D observations
%
% by Martin Bujnak (c) may 2012

function [f R t] = P4PfCalcCameraPose(za, zb, zc, zd, f, m2D, M3D)

    % recover camera rotation and translation
    % create p3d points in a camera coordinate system (using depths)

    p3dc(:, 1) = za * [m2D(:, 1); f];
    p3dc(:, 2) = zb * [m2D(:, 2); f];
    p3dc(:, 3) = zc * [m2D(:, 3); f];
    p3dc(:, 4) = zd * [m2D(:, 4); f];

    % calc camera
    [R t] = GetRigidTransform2(M3D, p3dc, [], false);
end