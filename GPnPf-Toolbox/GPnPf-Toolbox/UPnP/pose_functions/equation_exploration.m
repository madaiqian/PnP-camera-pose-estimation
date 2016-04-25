function [beta1,beta2,beta3,f] = equation_exploration(ordered_equations,Cw,Alph,Xw,U,A_UPnP,lambdas_,Kd,Kmf1,Kmf2,Kmf3,max_iter)
%clear all; close all; load ransac_epnp; intwarning('off')

lambda_size=6; % Number of lambda equations to use
lambda_window=8; % max 15
betaf_size=5; % Number of beta equations to use
betaf_window=7; % max 12
min_error_consensus=20;
prosac_threshold=10;

%%%%%%%%%%%%%%
lambda_counter= 0;
betaf_counter= 0;

failed= 0;
for iterator=1:max_iter
    
    while true
        try
            lambda_pts=randsample(lambda_window,lambda_size);
            betaf_pts=randsample(betaf_window,betaf_size);

            lambda_eq = ordered_equations.lambda_probability(lambda_pts,:);
            betaf_eq = ordered_equations.betaf_probability(betaf_pts,:);

            betas_and_f= solve_for_equation_combination( sort(lambda_eq(:,2))', sort(betaf_eq(:,2))', lambdas_, Kd );
            
            % Xci = beta*Vi1
            % Yci = beta*Vi2
            % Zci = f*beta*Vi3
            X3=betas_and_f(1)*Kmf1+betas_and_f(2)*Kmf2+betas_and_f(3)*Kmf3;
            X3_aux=X3;
            X3_aux(3:3:size(Kmf1,1))=X3(3:3:size(Kmf1,1))*betas_and_f(4);
            
            [Cc,Xc]=compute_norm_sign_scaling_factor(X3_aux,Cw,Alph,Xw);
            
            % Get new A matrix
            A_UPnP(1,1) = betas_and_f(4);
            A_UPnP(2,2) = betas_and_f(4);
            
            [R,T]=getrotT(Xw,Xc);  % Solve exterior orientation
            error_aux=reprojection_error_usingRT(Xw,U,R,T,A_UPnP);
            
            if ((iterator==1) || (error_aux < smallest_error))
                beta1= betas_and_f(1);
                beta2= betas_and_f(2);
                beta3= betas_and_f(3);
                f= betas_and_f(4);
                smallest_error=error_aux;
            end
            
            break;
        catch exception
            % Expand the initial windows
            if (lambda_window < 15)
                lambda_counter = lambda_counter + 1;
                if (lambda_counter >= prosac_threshold)
                    lambda_window = lambda_window + 1;
                    lambda_counter= 0;
                end
            end
            if (betaf_window < 15)
                betaf_counter = betaf_counter + 1;
                if (betaf_counter >= prosac_threshold)
                    betaf_window = betaf_window + 1;
                    betaf_counter= 0;
                end
            end
            
            failed = failed + 1;
        end
    end
    
    if (smallest_error < min_error_consensus)
        break;
    else
        % Expand the initial windows
        if (lambda_window < 15)
            lambda_counter = lambda_counter + 1;
            if (lambda_counter >= prosac_threshold)
                lambda_window = lambda_window + 1;
                lambda_counter= 0;
            end
        end
        if (betaf_window < 15)
            betaf_counter = betaf_counter + 1;
            if (betaf_counter >= prosac_threshold)
                betaf_window = betaf_window + 1;
                betaf_counter= 0;
            end
        end
    end
end


