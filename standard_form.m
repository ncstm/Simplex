function standard_form = standard_form(constrains, type_of_constrain, is_min)
%right side in the first column
%constrains first column has 
%type_of_constrain is a horizontal vector with entries corresponding to equations of constrains 0 means equal, -1 less (or equal) than and 1 more (or equal) than
%is_min is true or false. Assumed Max if false

%Notes: For system of equations with not a large amount of coefficients use matrix type for mostly empty matrix for performance, can add check
%Assuis_ming number of variables >= number of equations

%We put to max and equalice with slack, surplus and artificial, with ss at
%the on same columns and artifials after
    num_constr_equ = length( constrains(:, 1) ) - 1;
    %num_var = length(constrains(1,:)) - 1;

    if(is_min)
		constrains(1,:) = -1*constrains(1,:);
    end

    %For efficiency maybe less computation iff we fill the extended and
    %then do the check for less than 0 mambo jambo
    %Now be careful for the artifical if ennacted on the change
	%Makes sure we have all constrains coefficients more than zero to the right of the equation
    temp0 = constrains(2:end,1) < 0;
    type_of_constrain(temp0) = type_of_constrain(temp0)*-1;
    temp = [false ; temp0];
    constrains(temp, :) = constrains(temp, :)*-1;


    standard_form = [constrains [zeros(1, num_constr_equ) ; diag(type_of_constrain*-1)] [(type_of_constrain > -1)' ; diag(type_of_constrain > -1)]]; 

    %since surplus is substracted from equation so biger than is
    %substracted


end