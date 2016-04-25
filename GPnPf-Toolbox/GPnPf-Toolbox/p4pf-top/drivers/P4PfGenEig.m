% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008
%
% unknowns
%	f zb zc zd
% known
%	glab glac glad glbc glbd glcd a1 a2 b1 b2 c1 c2 d1 d2

function [f R t] = P4PfGenEig(m2D, M3D)

%%
    tol = 2.2204e-10;

    %normalize 2D, 3D
    
    % shift 3D data so that variance = sqrt(2), mean = 0
    mean3d = (sum(M3D') / 4);
    M3D = M3D - repmat(mean3d', 1, 4);

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
    
    if glab*glac*glad*glbc*glbd*glcd < tol
        
        % initial solution degeneracy - invalid input
        R = [];
        t = [];
        f = [];
        return;
    end

%%    
	[M1 M2] = CalcP4PfCoefs(glab, glac, glad, glbc, glbd, glcd, m2D(1,1), m2D(2,1), m2D(1,2), m2D(2,2), m2D(1,3), m2D(2,3), m2D(1,4), m2D(2,4));

    [v, w] = polyeig(M2,M1);
    w

    [v w] = eig(-M2, M1);
    ff=w(find(abs(imag(w(:))) < eps & abs(w(:)<5) & abs(w(:)>0)))
    sqrt(ff)
    sqrt(1/ff)
    1/ff
    
    %                                                                                                15 16 17 18 
    % zc*zb^2 zc^2*zb zb^2*zd zd*zc*zb zd*zc^2 zd^2*zb zd^2*zc zd^3 zb^2 zc*zb zc^2 zb*zd zd*zc zd^2 zb zc zd 1 

    R = [];
    t = [];
    f = [];

    for i=1:length(w)

        if (isfinite(w(i)) && isreal(w(i)) && w(i) > 0)

            za = v(18, i);
            zd = v(17, i);
            zc = v(16, i);
            zb = v(15, i);
            
            % create p3d points in a camera coordinate system (using depths)

            fi = sqrt(1/w(i));
            
            p3dc(:, 1) = za * [m2D(:, 1); fi];
            p3dc(:, 2) = zb * [m2D(:, 2); fi];
            p3dc(:, 3) = zc * [m2D(:, 3); fi];
            p3dc(:, 4) = zd * [m2D(:, 4); fi];

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
            [Rr tt] = GetRigidTransform2(M3D, p3dc, [], true);

            R(:,:,i)=Rr;
            t(:,:,i)=var*tt - Rr*mean3d';
            f(i)=var2d*fi;
        end
    end
end
