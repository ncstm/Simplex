function simplex = simplex(constrains)
%assuming constrains are bigger than variables

%coeff_var_in_constrain constrains from left of constrain equation 
%type_of_constrain[1,:] is the coefficients from right of the constrain equation
%type_of_constrain has rows corresponding to equations of constrains 0 means equal, -1 less (or equal) than and 1 more (or equal) than
%is_min is true or false. Assumed Max if false
%maximize is the coefficients of the equation we are trying to solve

%Improevments could be made by not calculating the columns of the
%artificial variables and just return the variables in the base with their
%corespondign value with 2 arrays.

%Notes: For system of equations with not a large amount of coefficients use matrix type for mostly empty matrix for performance, can add check
%Assuis_ming number of variables >= number of equations

	num_constrain_equ = length(constrains(:,1)) - 1;
	num_var = length(constrains(1,:)) - 2*num_constrain_equ - 1;
	%Makes sure we are solving a maximization problem
    %Restrain equations. 	
	
	%Has positions of variables we are working with. e.g. For 3 var x_1->1, s_1->4 
	%current_solution = zeros(num_constrain_equ, 1);

    %first solution
    %m_space = zeros(num_constrain_equ + 1, num_constrain_equ + 2*num_var);
    %not necessarly in basis
    %first_artificial = 2*num_constrain_equ + num_var + 2;
    m_c_vector = -1*[zeros(1,num_var + num_constrain_equ + 1) constrains(1,(end-num_constrain_equ + 1):end)]';
    % tic
    artificial_zone = (1 + 1 + num_var + num_constrain_equ):(1 + num_var + 2*num_constrain_equ);
    artificial_zone(constrains(1, artificial_zone)==0) = artificial_zone(constrains(1, artificial_zone)==0) - num_constrain_equ;
    current_solution = artificial_zone;
    constrains(1, (end-num_constrain_equ + 1):end) = zeros(num_constrain_equ,1)';

    change_to = 0;

    while true
        zi_m_version = sum(constrains(2:(num_constrain_equ+1), 2:(end - num_constrain_equ)).*m_c_vector(current_solution));
        %not needed because mc_vector will be 0 for non artificials and
        %since they are not considered for reeintroduction to the basis 

        min_val = min(zi_m_version);
        if(min_val <= 0)
            ind = find(zi_m_version == min_val) + 1;
            if (length(ind) ~= 1) || (min_val == 0)
                zi_normal = sum(constrains(2:(num_constrain_equ+1), ind).*constrains(1, current_solution)') - (constrains(1, ind));
                [temp , change_to] = min(zi_normal);
                if temp >= 0
                    if(max(current_solution) > 2*num_constrain_equ + 1)
                        disp("Infesible")
                        simplex = [[0 0 1:(2*num_constrain_equ + num_var )];[0; (current_solution-1)'] constrains];
                        return 
                    end
                    simplex = constrains;
                    return
                end
                change_to = ind(change_to) - 1;
            else 
                [~ , change_to] = min(zi_m_version);
            end
%change_to is what column has a change lowest z


            div_col = constrains(2:end, 1)./constrains(2:end, change_to + 1);
            %what about constrains(2:end, 1) < 0
            div_col(constrains(2:end, change_to + 1) <= 0) = Inf;
            [~, new_ind] = min(div_col);
            current_solution(new_ind) = change_to + 1;
            
            constrains(new_ind + 1, :) = constrains(new_ind + 1,:)./constrains(new_ind + 1, change_to + 1);


            temp_elimination  = [false (constrains(2:end, change_to + 1) ~= 0)'];
            temp_elimination(new_ind + 1) = false;
            constrains(temp_elimination,:) = constrains(temp_elimination,:) - constrains(new_ind + 1, :).*constrains(temp_elimination,change_to+1);

        elseif min_val > 0
            %Because we can't have a z with M's without an artificial in
            %the basis
            %Have to check that the normals are not gone
            disp("Infesible")
            simplex = [current_solution ; constrains(2:end, 1)];
            return 
        end
    end
    %simplex = [current_solution ; constrains(2:end, 1)];
end

% 