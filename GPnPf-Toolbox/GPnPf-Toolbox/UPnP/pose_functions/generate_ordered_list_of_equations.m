function ordered_equations= generate_ordered_list_of_equations(n_samples, load_order_data)

try
    filename = 'data/upnp_equation_order_data.mat';
    
    % If data has already been generated, we load it
    if load_order_data
        load(filename,'ordered_equations');
        return;
    end
    
    %% HOW TO CREATE TRAINING SAMPLES
    % We need the data that has been previously calculated by
    % "training_upnp.m" to get the probabilities for each combination of
    % equations. Save the output in this file format:
    %       [ 'equationOrderData/upnp_data' 1:n_samples '.mat' ];
    %
    % Our data is not included beacuse it takes up to 2.2 GB.
    %%
    
    lambda_probability = zeros(1,15);
    betaf_probability = zeros(1,12);
    ordered_equations.lambda_probability = zeros(1,15);
    ordered_equations.betaf_probability = zeros(1,12);
    
    for iterator=1:n_samples
        filename2 = [ 'equationOrderData/upnp_data' num2str(iterator) '.mat' ];
        load(filename2,'sol');
        
        if(iterator==1)
            sum_sol= sol(3);
        else
            sum_sol.case_error_Rec(:)= sum_sol.case_error_Rep(:) + sol(3).case_error_Rep(:);
        end
        
        nelements = size(sum_sol.case_error_Rep,1);
        
        for iterator2=1:nelements
            lambda_probability(sol(3).lambda_eq(iterator2,1)) = lambda_probability(sol(3).lambda_eq(iterator2,1)) + sol(3).case_error_Rep(iterator2);
            lambda_probability(sol(3).lambda_eq(iterator2,2)) = lambda_probability(sol(3).lambda_eq(iterator2,2)) + sol(3).case_error_Rep(iterator2);
            lambda_probability(sol(3).lambda_eq(iterator2,3)) = lambda_probability(sol(3).lambda_eq(iterator2,3)) + sol(3).case_error_Rep(iterator2);
            lambda_probability(sol(3).lambda_eq(iterator2,4)) = lambda_probability(sol(3).lambda_eq(iterator2,4)) + sol(3).case_error_Rep(iterator2);
            lambda_probability(sol(3).lambda_eq(iterator2,5)) = lambda_probability(sol(3).lambda_eq(iterator2,5)) + sol(3).case_error_Rep(iterator2);
            
            betaf_probability(sol(3).beta_eq(iterator2,1)) = betaf_probability(sol(3).beta_eq(iterator2,1)) + sol(3).case_error_Rep(iterator2);
            betaf_probability(sol(3).beta_eq(iterator2,2)) = betaf_probability(sol(3).beta_eq(iterator2,2)) + sol(3).case_error_Rep(iterator2);
            betaf_probability(sol(3).beta_eq(iterator2,3)) = betaf_probability(sol(3).beta_eq(iterator2,3)) + sol(3).case_error_Rep(iterator2);
            betaf_probability(sol(3).beta_eq(iterator2,4)) = betaf_probability(sol(3).beta_eq(iterator2,4)) + sol(3).case_error_Rep(iterator2);
        end
        
        ordered_equations.lambda_probability = ordered_equations.lambda_probability + lambda_probability/max(lambda_probability);
        ordered_equations.betaf_probability = ordered_equations.betaf_probability + betaf_probability/max(betaf_probability);
    end
    
    ordered_equations.lambda_probability = [ordered_equations.lambda_probability/n_samples; 1:15];
    ordered_equations.betaf_probability = [ordered_equations.betaf_probability/n_samples; 1:12];
    
    ordered_equations.lambda_probability = sort_matrix(ordered_equations.lambda_probability',1,'ascend');
    ordered_equations.betaf_probability = sort_matrix(ordered_equations.betaf_probability',1,'ascend');
    
    save(filename,'ordered_equations','-v7.3');
catch exception
    disp(exception);
end