function [lin_vars, p_var, P2] = split_matrices_5C3(C, x, y, z, options)
%SPLIT_MATRICES_5C3 Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_5C3(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices Q and P for all variable choices.If Q is invertable,
%   the matrix and the combinations of variables for the approach are returned.
arguments
    C (5,20) {mustBeReal};
    x (1,1) sym = sym("x","real"); % Symbolic variable for x-direction
    y (1,1) sym = sym("y","real"); % Symbolic variable for y-direction
    z (1,1) sym = sym("z","real"); % Symbolic variable for z-direction
    options.verbose (1,1) {mustBeInteger, mustBeInRange(options.verbose,0,2)} = 0;
    options.error (1,1) {mustBeNumericOrLogical} = true;
end

assert(~any(C(:,5)))
A = -[C(:,7),C(:,8),C(:,9),C(:,10),C(:,5)*x+C(:,15)];
B = -[C(:,1),C(:,3),C(:,6),C(:,10),C(:,5)*y+C(:,13)];
C = -[C(:,1),C(:,2),C(:,4),C(:,7),C(:,5)*z+C(:,12)];
i =   [ 4,14,6,16,2,12,18,3,13,19,1,11,17;...
    2,11,9,16,4,12,17,8,15,19,7,14,18;...
    3,11,8,14,6,13,17,9,15,18,10,16,19];

Q_list = {A,B,C};
v = [x;y;z];

if any([rank(A),rank(B),rank(C)]==5)
    conds = real([vpa(cond(A)),vpa(cond(B)),vpa(cond(C))]);
    [~,I] = min(conds,[],"omitnan");
    if options.verbose > 0
        fprintf("Using P(%s)\n",string(v(I)))
    end
    Q = Q_list{I};
    p_var = v(I);
    
    P  = [C(:,i(I,1)).*p_var+C(:,i(I,2)),...
        C(:,i(I,3)).*p_var+C(:,i(I,4)),...
        C(:,i(I,5)).*p_var^2+C(:,i(I,6)).*p_var+C(:,i(I,7)),...
        C(:,i(I,8)).*p_var^2+C(:,i(I,9)).*p_var+C(:,i(I,10)),...
        C(:,i(I,11)).*p_var^3+C(:,i(I,12)).*p_var^2+C(:,i(I,13)).*p_var+C(:,20)];
    lin_vars = [v(setdiff(1:3,I)).^2;v(setdiff(1:3,I));1];
    P2 = Q\P;
else
    if options.error
        error("All Matrices are singular!")
    else
        warning("All Matrices are singular!")
    end
end