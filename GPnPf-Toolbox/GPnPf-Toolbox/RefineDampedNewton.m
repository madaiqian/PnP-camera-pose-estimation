function [f_r, q_r, obj_pre] = RefineDampedNewton(f0,q0,MTM,MTb,Q)
%q=[a b c d];
%vec = [1 a^2 ab ac ad b^2 bc bd c^2 cd d^2];
%objective function: vec.'*Q*vec
%refine the solution by using damped Newton method

%maximum allowed iterations (one-step in our implementation)
maxItr = 10; 
%damped factor
lambda = 1e-8;
%max lambda
maxLambda = 1e2;

%min lambda
minLambda = 1e-8;

%flag 
flag = 0; 

fq = [f0; q0];
f = fq(1);
a = fq(2); b = fq(3); c = fq(4); d = fq(5);
vec = [ a^2*f, a*b*f, a*c*f, a*d*f, b^2*f, b*c*f, b*d*f, c^2*f, c*d*f, d^2*f, a^2, a*b, a*c, b^2, b*d, c^2, c*d, d^2];
obj_pre = [vec 1]*Q*[vec 1].';

itr = 1; 
%iteration 
while itr <= maxItr 
        
    f = fq(1); a = fq(2); b = fq(3); c = fq(4); d = fq(5);
    vec = [a^2*f, a*b*f, a*c*f, a*d*f, b^2*f, b*c*f, b*d*f, c^2*f, c*d*f, d^2*f, a^2, a*b, a*c, b^2, b*d, c^2, c*d, d^2];
    
    jac = [   a^2, a*b, a*c, a*d,   b^2, b*c, b*d,   c^2, c*d,   d^2,   0, 0, 0,   0, 0,   0, 0,   0 
            2*a*f, b*f, c*f, d*f,     0,   0,   0,     0,   0,     0, 2*a, b, c,   0, 0,   0, 0,   0
                0, a*f,   0,   0, 2*b*f, c*f, d*f,     0,   0,     0,   0, a, 0, 2*b, d,   0, 0,   0 
                0,   0, a*f,   0,     0, b*f,   0, 2*c*f, d*f,     0,   0, 0, a,   0, 0, 2*c, d,   0 
                0,   0,   0, a*f,     0,   0, b*f,     0, c*f, 2*d*f,   0, 0, 0,   0, b,   0, c, 2*d ];
    g = jac*(MTb - MTM*vec.');
 
    fq_temp = fq;
    while (lambda < maxLambda)
        %increment
        delta = (jac*MTM*jac.'+lambda*eye(5))\g;

        %update parameter
        fq = fq_temp + delta;
         
        f = fq(1); a = fq(2); b = fq(3); c = fq(4); d = fq(5);
        vec = [a^2*f, a*b*f, a*c*f, a*d*f, b^2*f, b*c*f, b*d*f, c^2*f, c*d*f, d^2*f, a^2, a*b, a*c, b^2, b*d, c^2, c*d, d^2];
        
        %evaluate the sampson error 
        obj_cur = [vec 1]*Q*[vec 1].';

        %check convergence
        if obj_cur >= obj_pre
            lambda = 10*lambda;
            continue;
        else
            obj_pre = obj_cur; 
            lambda = 0.1*lambda;
            break;
        end        
    end
    if lambda >= maxLambda
        fq = fq_temp;
        break;
    end
    
    if lambda <= minLambda
        lambda = minLambda;
    end
    itr = itr + 1;
end

%account for the sign ambiguity
f_r = fq(1); q_r = fq(2:end); 
if f_r < 0
    f_r = -f_r;
    q_r = [q_r(4); q_r(3); -q_r(2); -q_r(1)];
end
end

