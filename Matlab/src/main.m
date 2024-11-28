function main(polynomials,vars)
arguments
    polynomials (:,1);
    vars (1,3) = [sym("x","real"),sym("y","real"),sym("z","real")]
end
pos_matrix = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,10];
pos_matrix(:,:,2) = [0,0,0,0;0,0,0,0;0,0,0,6;0,0,9,16];
pos_matrix(:,:,3) =    [0,0,0,0;0,0,0,3;0,0,5,13;0,8,15,19];
pos_matrix(:,:,4) =    [0,0,0,1;0,0,2,11;0,4,12,17;8,14,18,20];
n = numel(polynomials);
C = zeros(n,20);
for i=1:n
    p = polynomials(i);
    % for s =
    [c,~] = coeffs(p,vars,"All");
    for x_ind=min(size(c,1),4)-1:-1:0
        for y_ind=min(size(c,2),4)-1:-1:0
            for z_ind = min(size(c,3),4)-1:-1:0
                if pos_matrix(4-x_ind,4-y_ind,4-z_ind) > 0
                    C(i,pos_matrix(4-x_ind,4-y_ind,4-z_ind)) = c(size(c,1)-x_ind,size(c,2)-y_ind,size(c,3)-z_ind);
                end
            end
        end
    end
end
    C
    if C(:,1:10) == zeros(n,10)
        if n==3
            E3Q3(C(:,11:20))
        else
            error("Invalid number of equations for the E3Q3 solver.")
        end
    else
        switch n
            case 3
                E3C3(C)
            case 4
                E4C3(C)
            case 5
                E5C3(C)
            otherwise
                error("Invalid number of equations for the solver.")
        end
    end

end