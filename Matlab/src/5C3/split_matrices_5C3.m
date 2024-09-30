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
p_idx =   [ 4,14,6,16,2,12,18,3,13,19,1,11,17;...
    2,11,9,16,4,12,17,8,15,19,7,14,18;...
    3,11,8,14,6,13,17,9,15,18,10,16,19];
q_idx = [7,8,9,10,5,15;1,3,6,10,5,13;1,2,4,7,5,12];
Q = cell(1,3);
ranks = zeros(3,1);
conds = nan(3,1);
for i=1:3
    Q{i} = -[C(:,q_idx(i,1)),C(:,q_idx(i,2)),C(:,q_idx(i,3)),...
        C(:,q_idx(i,4)),C(:,q_idx(i,5))*x+C(:,q_idx(i,6))];
    ranks(i) = rank(Q{i});
    conds(i) = cond(Q{i});
end
v = [x;y;z];

if any(ranks==5)
    [~,I] = min(conds,[],"omitnan");
    if options.verbose > 0
        fprintf("Using P(%s)\n",string(v(I)))
    end
    p_var = v(I);
    P  = [C(:,p_idx(I,1)).*p_var+C(:,p_idx(I,2)),...
        C(:,p_idx(I,3)).*p_var+C(:,p_idx(I,4)),...
        C(:,p_idx(I,5)).*p_var^2+C(:,p_idx(I,6)).*p_var+C(:,p_idx(I,7)),...
        C(:,p_idx(I,8)).*p_var^2+C(:,p_idx(I,9)).*p_var+C(:,p_idx(I,10)),...
        C(:,p_idx(I,11)).*p_var^3+C(:,p_idx(I,12)).*p_var^2+C(:,p_idx(I,13)).*p_var+C(:,20)];
    lin_vars = [v(setdiff(1:3,I)).^2;v(setdiff(1:3,I));1];
    P2 = Q{I}\P;
else
    if options.error
        error("All Matrices are singular!")
    else
        warning("All Matrices are singular!")
    end
end