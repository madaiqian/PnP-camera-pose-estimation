function [f R t] = Best_Distance_GB(M3D,m2D)

% shift 3D data so that variance = sqrt(2), mean = 0
mean3d = (sum(M3D,2) / 4);
M3D = M3D - repmat(mean3d, 1, 4);

% variance (isotropic)
var = (sum( sqrt(sum( M3D.^2 ) ) ) / 4);
M3D = (1/var)*M3D;

% scale 2D data
var2d = (sum( sqrt(sum( m2D.^2 ) ) ) / 4);
m2D = (1/var2d)*m2D;

%coefficients of 5 reduced polynomials in these monomials mon
glab = (sum((M3D(:,1)-M3D(:,2)).^2));
glac = (sum((M3D(:,1)-M3D(:,3)).^2));
glad = (sum((M3D(:,1)-M3D(:,4)).^2));
glbc = (sum((M3D(:,2)-M3D(:,3)).^2));
glbd = (sum((M3D(:,2)-M3D(:,4)).^2));
glcd = (sum((M3D(:,3)-M3D(:,4)).^2));

%call the best distance solver (561x581, 20 solutions)
[f za zb zc zd] = p4pfd_1_179([glab,glac,glad,glbc,glbd,glcd], m2D(:,1), m2D(:,2), m2D(:,3), m2D(:,4));

if nargout == 1
    %denormalize
    f=var2d*sqrt(f);
%recover the rotation and translation
else
    % recover camera rotation and translation
    lcnt = length(f);
    R = zeros(3,3,lcnt);
    t = zeros(3,lcnt);
    for i=1:lcnt
        %
        f(i) = sqrt(f(i));
        
        % create p3d points in a camera coordinate system (using depths)
        p3dc(:, 1) = za(i) * [m2D(:, 1); f(i)];
        p3dc(:, 2) = zb(i) * [m2D(:, 2); f(i)];
        p3dc(:, 3) = zc(i) * [m2D(:, 3); f(i)];
        p3dc(:, 4) = zd(i) * [m2D(:, 4); f(i)];

        % fix scale (recover 'za')
        d(1) = sqrt(glab / (sum((p3dc(:,1)-p3dc(:,2)).^2)));
        d(2) = sqrt(glac / (sum((p3dc(:,1)-p3dc(:,3)).^2)));
        d(3) = sqrt(glad / (sum((p3dc(:,1)-p3dc(:,4)).^2)));
        d(4) = sqrt(glbc / (sum((p3dc(:,2)-p3dc(:,3)).^2)));
        d(5) = sqrt(glbd / (sum((p3dc(:,2)-p3dc(:,4)).^2)));
        d(6) = sqrt(glcd / (sum((p3dc(:,3)-p3dc(:,4)).^2)));

        % all d(i) should be equal...but who knows ;)
        %gta = median(d);
        gta = sum(d) ./ 6;

        p3dc = gta * p3dc;

        % calc camera
        [Rr tt] = GetRigidTransform2(M3D, p3dc, [],false);

        R(:,:,i)=Rr;
        t(:,i)=var*tt - Rr*mean3d;
        f(i)=var2d*f(i);
    end
end
end