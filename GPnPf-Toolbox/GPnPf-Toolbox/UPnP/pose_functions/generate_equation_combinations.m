function combinations=generate_equation_combinations(number_of_equations, number_of_unknowns)

% basic case
if (number_of_unknowns <= 1)
    combinations= zeros(number_of_equations, 1);
    for i=1:number_of_equations
        combinations(i,1)= i;
    end
else
    % recursive case
    aux_combinations= generate_equation_combinations(number_of_equations, number_of_unknowns-1);
    previous_combination_size= size(aux_combinations,1);
    
    combinations= zeros(1,number_of_unknowns);
    index= 1;
    for i=1:previous_combination_size
        % We start from 1 value plus de depth of the recursive call
        for j=number_of_unknowns:number_of_equations
            % we don't consider trivial solutions
            aux_value= find(aux_combinations(i,:)==j);
            [nrows,ncolumns]= size(aux_value);
            aux_tam= nrows*ncolumns;
            if (aux_tam == 0) && (aux_combinations(i,number_of_unknowns-1)<j)
                combinations(index,:)= [aux_combinations(i,:), j];
                index= index+1;
            end
        end
    end
end
