%
% calculate poses for all p4pf solvers

function [D] = GetAllPoses(filter, m, M, priors)

    solvers_ratio = {@p4pfr_1_116 @p4pfr_2_287 @p4pfr_3_393 @p4pfr_4_139};
    solvers_dist = {@p4pfd_1_179 @p4pfd_2_10 @p4pfd_3_126 @p4pfd_4_7678za};
    solvers_DIAC = {@mbp4pf_np_1_7036 @mbp4pf_np_2_7015 @mbp4pf_np_3_2036};
    solvers_quat1 = {@p4pf_quat2 @p4pf_quat2_wdcb @p4pf_quat3_wdcb};
    solvers_quat2 = {@p4pf_quat101};
    
    allAlgs = 19;
    
    if nargin == 1
        
        % create fake 2D-3D coordinates in order to be able to read
        % metadata
        [Pgt M m] = GenerateScene(4, [40 20], 1, 80, 250, 0, 0, diag([1 1 1]), 'randombox', [], [], [], true);
        m = m{1};
        m(3,:) = 1;
        M(4,:) = 1;
        
        %    Solid line, Dash-dot line, Dashed line, Dotted line
        cfgstyles = {'-' '-.*' '--' ':'};
        
        % generate abs algs info
        D.algsNames = {'p3p'};
        D.algsFullNames = {'p3p'};
        D.algsColors = [1 0 0];
        D.algsStyle = {'-'};

        % read params from metadata
        for i=1:length(solvers_ratio)
            
            [fs R t meta Mr] = EvalP4PfSolverUni(solvers_ratio{i}, m(1:2,:), M(1:3,:));
            r = size(Mr, 1);
            c = size(Mr, 2);
            rc = ['ratio-' int2str(r) '\times' int2str(c)];
            D.algsNames = [D.algsNames rc];
            D.algsFullNames = [D.algsFullNames ['p4p+focal+' rc]];
            D.algsColors = [D.algsColors; [0 1 0]];
            D.algsStyle = [D.algsStyle cfgstyles{i}];
        end

        for i=1:length(solvers_ratio)
            
            [fs R t meta Mr] = EvalP4PfSolverUni(solvers_dist{i}, m(1:2,:), M(1:3,:));
            r = size(Mr, 1);
            c = size(Mr, 2);
            rc = [int2str(r) '\times' int2str(c)];
            D.algsNames = [D.algsNames ['dist-' rc]];
            D.algsFullNames = [D.algsFullNames ['p4p+focal+distance-' rc]];
            D.algsColors = [D.algsColors; [0 0 1]];
            D.algsStyle = [D.algsStyle cfgstyles{i}];
        end

        for i=1:length(solvers_quat1)
            
            [f Ps meta Mr] = P4PfQuat(solvers_quat1{i}, m, M);
            r = size(Mr, 1);
            c = size(Mr, 2);
            rc = [int2str(r) '\times' int2str(c)];
            D.algsNames = [D.algsNames ['quat-' rc]];
            D.algsFullNames = [D.algsFullNames ['p4p+focal+quaternion-' rc]];
            D.algsColors = [D.algsColors; [0 1 1]];
            D.algsStyle = [D.algsStyle cfgstyles{i}];
        end
        
        for i=1:length(solvers_quat2)
            
            [f Ps meta Mr] = P4PfQuat2(solvers_quat2{i}, m, M);
            r = size(Mr, 1);
            c = size(Mr, 2);
            rc = [int2str(r) '\times' int2str(c)];
            D.algsNames = [D.algsNames ['quat-2-' rc]];
            D.algsFullNames = [D.algsFullNames ['p4p+focal+quaternion2-' rc]];
            D.algsColors = [D.algsColors; [1 1 0]];
            D.algsStyle = [D.algsStyle cfgstyles{i}];
        end

        for i=1:length(solvers_DIAC)
            
            [f R t meta Mr] = EvalP4PfNonPlanarSolver(solvers_DIAC{i}, m, M);
            r = size(Mr, 1);
            c = size(Mr, 2);
            rc = ['DIAC-' int2str(r) '\times' int2str(c)];
            D.algsNames = [D.algsNames rc];
            D.algsFullNames = [D.algsFullNames ['p4p+focal+' rc]];
            D.algsColors = [D.algsColors; [1 0 1]];
            D.algsStyle = [D.algsStyle cfgstyles{i}];
        end
        
        D.algsNames = [D.algsNames 'hiddenvar'];
        D.algsFullNames = [D.algsFullNames 'hidden variable'];
        D.algsColors = [D.algsColors; [0 0 0]];
        D.algsStyle = [D.algsStyle cfgstyles{1}];

        D.algsNames = [D.algsNames 'cvpr08'];
        D.algsFullNames = [D.algsFullNames 'cvpr08 original'];
        D.algsColors = [D.algsColors; [0 0 0]];
        D.algsStyle = [D.algsStyle cfgstyles{2}];

        D.algsNames = [D.algsNames 'planar'];
        D.algsFullNames = [D.algsFullNames 'planar (H)'];
        D.algsColors = [D.algsColors; [0.7 0.7 0.7]];
        D.algsStyle = [D.algsStyle cfgstyles{1}];

        D.algsCount = length(D.algsNames);

        
        D.measurementSize = 3;
        D.measurementNames = {'Camera' 'Focal length' 'Radial distortion'};
        D.min_sample = 4;
                                      
        % select algorithms
        if ~isempty(filter)
            
            selection = find(filter == 1);
            D.algsNames = D.algsNames(selection);
            D.algsFullNames = D.algsFullNames(selection);
            D.algsColors = D.algsColors(selection);
            D.algsStyle = D.algsStyle(selection);
        end
        
        return;
    end

    m = m(:,1:4);
    M = M(:,1:4);
    m(3,:) = 1;
    M(4,:) = 1;
    
    if (isempty(filter))
        filter = ones(1,allAlgs);
    end
    
    algsCnt = sum(filter);
    D = cell(1,algsCnt);
    a = 1;

    Kgt = priors.K;
    lgt = priors.l;

    algi = 1;
   
    %
    % calibrated p3p
    if (filter(algi))
        
        [P] = CalcP3P(Kgt, m(1:3,1:3), M(1:3,1:3));
        cnt = size(P, 1)/3;
        for i=1:cnt

            Pi = P((3*i-2):(3*i), :);
            D{a}.P{i} = Kgt*Pi;
            D{a}.f{i} = Kgt(1);
            D{a}.l{i} = lgt;
        end
        D{a}.sols = cnt;
        a = a + 1;
    end
    algi = algi + 1;
   
    
    %
    % ratio solvers
    for s=solvers_ratio

        if (filter(algi))
           
            try
                [fs Rs ts] = EvalP4PfSolverUni(s{1}, m(1:2,:), M(1:3,:));
                
                cnt = size(fs, 2);
                for i=1:cnt

                    f = fs(i);
                    R = Rs(:,:,i);
                    t = ts(:,i);

                    D{a}.P{i} = diag([f f 1])*[R t];
                    D{a}.f{i} = f;
                    D{a}.l{i} = lgt;
                end
                D{a}.sols = cnt;
            catch e
                D{a}.sols = 0;
            end
            a = a + 1;
        end       
        algi = algi + 1;
    end
    
    
    %
    % distance solvers
    for s=solvers_dist

        if (filter(algi))
            
            try
                [fs Rs ts] = EvalP4PfSolverUni(s{1}, m(1:2,:), M(1:3,:));
                cnt = size(fs, 2);
                for i=1:cnt

                    f = fs(i);
                    R = Rs(:,:,i);
                    t = ts(:,i);

                    D{a}.P{i} = diag([f f 1])*[R t];
                    D{a}.f{i} = f;
                    D{a}.l{i} = lgt;
                end
                D{a}.sols = cnt;
            catch e
                D{a}.sols = 0;
            end
            a = a + 1;
        end       
        algi = algi + 1;
    end

    
    %
    % quaternion solvers
    for s=solvers_quat1

        if (filter(algi))
            
            try
                [fs Ps] = P4PfQuat(s{1}, m, M);
                cnt = size(fs, 2);
                for i=1:cnt

                    f = fs(i);
                    P = Ps(:,:,i);

                    D{a}.P{i} = P;
                    D{a}.f{i} = f;
                    D{a}.l{i} = lgt;
                end
                D{a}.sols = cnt;
            catch e
                D{a}.sols = 0;
            end
            a = a + 1;
        end       
        algi = algi + 1;
    end
    
    for s=solvers_quat2

        if (filter(algi))
            
            try
                [fs Ps] = P4PfQuat2(s{1}, m, M);
                cnt = size(fs, 2);
                for i=1:cnt

                    f = fs(i);
                    P = Ps(:,:,i);

                    D{a}.P{i} = P;
                    D{a}.f{i} = f;
                    D{a}.l{i} = lgt;
                end
                D{a}.sols = cnt;
            catch e
                D{a}.sols = 0;
            end
            a = a + 1;
        end       
        algi = algi + 1;
    end
    
    
    %
    % non-planar diac solvers
    for s=solvers_DIAC

        if (filter(algi))
            
            try
                [fs Rs ts] = EvalP4PfNonPlanarSolver(s{1}, m, M);
                cnt = size(fs, 2);
                for i=1:cnt

                    f = fs(i);
                    R = Rs(:,:,i);
                    t = ts(:,i);

                    D{a}.P{i} = diag([f f 1])*[R t];
                    D{a}.f{i} = f;
                    D{a}.l{i} = lgt;
                end
                D{a}.sols = cnt;
            catch
                D{a}.sols = 0;
            end
            a = a + 1;
        end       
        algi = algi + 1;
    end   

    
    %
    % hidden variable solver
    if (filter(algi))
        
        D{a}.sols = 0; 
        try
            [fs] = pnpf_hv_peig(m(1:2,1:4), M(1:3,1:4));
            cnt = size(fs, 1);
            if (size(fs, 1) == 0) 
                cnt = 0;
            end

            ij = 0;
            for i=1:cnt

                f = fs(i);
                if (f > 0 && f < 10000)
                
                    % techically, it is p3p once focal length is found                
                    % calculare rotation-translation 

                    K = diag([f f 1]);
                    [P] = CalcP3P(K, m(1:3,1:3), M(1:3,1:3));
                    subcnt = size(P, 1)/3;
                    for j=1:subcnt

                        Pi = P((3*j-2):(3*j), :);
                        ij = ij + 1;
                        D{a}.P{ij} = K*Pi;
                        D{a}.f{ij} = f;
                        D{a}.l{ij} = lgt;
                    end
                    D{a}.sols = D{a}.sols + subcnt;
                end
            end
        catch
        end
        a = a + 1;
    end
    algi = algi + 1;
    
    
    %
    % original CVPR 2008 solver
    if (filter(algi))
        
        try
            [fs,Rs,ts] = P4PfocalGBr(m(1:2,1:4), M(1:3,1:4));
            cnt = size(fs, 2);
            for i=1:cnt

                f = fs(i);
                R = Rs(:,:,i);
                t = ts(:,i);

                D{a}.P{i} = diag([f f 1])*[R t];
                D{a}.f{i} = f;
                D{a}.l{i} = lgt;
            end
            D{a}.sols = cnt;
        catch
            D{a}.sols = 0;
        end
        a = a + 1;
    end
    algi = algi + 1;
    
    
    %
    % planar solver
    if (filter(algi))
        
        D{a}.sols = 0;
        try
            [fs RTs] = CalcP4PfPlanar(M(1:3,1:4), m(1:2,1:4));

            if (isreal(fs))
        
                D{a}.P{1} = diag([fs fs 1])*RTs;
                D{a}.f{1} = fs;
                D{a}.l{1} = lgt;
                D{a}.sols = 1;
            end
        catch
        end
        a = a + 1;
    end
    algi = algi + 1;    
end