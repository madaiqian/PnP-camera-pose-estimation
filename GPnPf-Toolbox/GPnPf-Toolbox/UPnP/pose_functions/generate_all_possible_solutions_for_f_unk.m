% Calculates all possible solutions for f unknown cases...
%
% We return the solutions in the rows
function [solutions, equation_combination]=generate_all_possible_solutions_for_f_unk(Betas_fij, number_of_Bi)

aux_number_of_equations= cumsum(1:number_of_Bi);
number_of_equations= aux_number_of_equations(number_of_Bi)*2;
number_of_unknowns= number_of_Bi+1; % we add f to the unknowns
number_of_Bij= round(number_of_equations/2);

% We could also use nchoosek
combinations= generate_equation_combinations(number_of_equations, number_of_unknowns);

number_of_solutions= size(combinations,1);

if (combinations(1,1)~=0)
    
    % We introduce
    linearization_coordinates= generate_linearization_coordinates(number_of_equations/2, number_of_Bi);
    
    solutions= zeros(1,number_of_unknowns);
    equation_combination= zeros(1,number_of_unknowns);
    index= 1;
    for i=1:number_of_solutions
        % We will apply logarithms on both sides in order to calculate the
        % absolute value
        matrix_to_resolve= zeros(number_of_unknowns);
        independent_term= zeros(1,number_of_unknowns);
        for iterator=1:number_of_unknowns
            % solution = [B1,B2, .. ,f]
            % If the equation is one of the Bfij
            if (combinations(i,iterator) > number_of_Bij)
                coordinates= linearization_coordinates(combinations(i,iterator)-number_of_Bij,:);
                % We introduce the f value for the vector
                matrix_to_resolve(iterator, number_of_unknowns)= 2;
            else
                coordinates= linearization_coordinates(combinations(i,iterator),:);
            end
            
            if (coordinates(1)==coordinates(2))
                matrix_to_resolve(iterator, coordinates(1))= 2;
            else
                matrix_to_resolve(iterator, coordinates(1))= 1;
                matrix_to_resolve(iterator, coordinates(2))= 1;
            end
            
            independent_term(iterator)= log(abs(Betas_fij(combinations(i,iterator))));
        end
        
        if (rank(matrix_to_resolve) == number_of_unknowns)
            solutions(index,:)= exp(inv(matrix_to_resolve)*independent_term');
            
            for iterator=2:(number_of_unknowns-1)
                solutions(index,iterator)= solutions(index,iterator)*sign(Betas_fij(iterator));
            end
            solutions(index,number_of_unknowns)= abs(solutions(index,number_of_unknowns));
            equation_combination(index,:)= combinations(i,:);
            
            index= index+1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We create the vector to trasform the coordinates of the vector to the
% coordinates of the linearization
function linearization_coordinates=generate_linearization_coordinates(number_of_equations, number_of_Bi)

coordinate_matrix= zeros(number_of_Bi);
counter= 1;
for i=1:number_of_Bi
    for j=i:number_of_Bi
        coordinate_matrix(i,j)= counter;
        counter= counter+1;
    end
end

linearization_coordinates= zeros(number_of_equations,2);
for iterator=1:number_of_equations
    [row,column]= find(coordinate_matrix==iterator);
    linearization_coordinates(iterator,:)= [row, column];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
