function [lin_vars, p_var, P2] = split_matrices_5C3(c, x, y, z, options)
%SPLIT_MATRICES_5C3 Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_5C3(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices Q and P for all variable choices.If Q is invertable,
%   the matrix and the combinations of variables for the approach are returned.
arguments
    c (5,20) {mustBeReal};
    x (1,1) sym = sym("x","real"); % Symbolic variable for x-direction
    y (1,1) sym = sym("y","real"); % Symbolic variable for y-direction
    z (1,1) sym = sym("z","real"); % Symbolic variable for z-direction
    options.verbose (1,1) {mustBeInteger, mustBeInRange(options.verbose,0,2)} = 0;
    options.error (1,1) {mustBeNumericOrLogical} = true;
end

assert(~any(c(:,5)))
A = -[c(:,7),c(:,8),c(:,9),c(:,10),c(:,5)*x+c(:,15)];
B = -[c(:,1),c(:,3),c(:,6),c(:,10),c(:,5)*y+c(:,13)];
C = -[c(:,1),c(:,2),c(:,4),c(:,7),c(:,5)*z+c(:,12)];
i =   [ 4,14,6,16,2,12,18,3,13,19,1,11,17;...
    2,11,9,16,4,12,17,8,15,19,7,14,18;...
    3,11,8,14,6,13,17,9,15,18,10,16,19];

Q_list = {A,B,C};
v = [x;y;z];

if any([rank(A),rank(B),rank(C)]==5)
    conds = real([vpa(cond(A)),vpa(cond(B)),vpa(cond(C))]);
    [M,I] = min(conds,[],"omitnan");
else
    M = Inf;
end


if ~isinf(M)
    if options.verbose > 0
        fprintf("Using P(%s)\n",string(v(I)))
    end
    Q = Q_list{I};
    p_var = v(I);
    
    P  = [c(:,i(I,1)).*p_var+c(:,i(I,2)),...
        c(:,i(I,3)).*p_var+c(:,i(I,4)),...
        c(:,i(I,5)).*p_var^2+c(:,i(I,6)).*p_var+c(:,i(I,7)),...
        c(:,i(I,8)).*p_var^2+c(:,i(I,9)).*p_var+c(:,i(I,10)),...
        c(:,i(I,11)).*p_var^3+c(:,i(I,12)).*p_var^2+c(:,i(I,13)).*p_var+c(:,20)];
    lin_vars = [v(setdiff(1:3,I)).^2;v(setdiff(1:3,I));1];
    P2 = Q\P;
    if options.verbose >= 2
        fprintf("Q:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(Q))
        fprintf("P:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(P))
        fprintf("P2:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(P2))
    end
else
    if options.error
        error("All Matrices are singular!")
    else
        warning("All Matrices are singular!")
    end
end

