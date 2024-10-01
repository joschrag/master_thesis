function [lin_vars, p_var, P2] = split_matrices_E3Q3_Fp(C,prime, x, y, z, opt)
%SPLIT_MATRICES_E3Q3_FP Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_E3Q3_FP(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices Q and P for all variable choices. If Q is invertable,
%   the matrix and the combinations of variables for the approach are returned.
arguments
    C (3,10) {mustBeReal}; % Coefficient matrix
    prime (1,1) {mustBePrimeOrZero,mustBePositive};
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 0;
end
q_idx = [4,5,6;1,3,5;1,2,4];
Q = cell(1,3);
ranks = zeros(1,3);

for i=1:3
    Q{i} = FF(-[C(:,q_idx(i,:))],prime);
    ranks(i) = rank(Q{i});
end

p_idx =   [2,8,3,9,1,7;...
    2,7,5,9,4,8;...
    3,7,5,8,6,9];
v = [x;y;z];

% Find the first invertible matrix in Q and proceed
if any(ranks==3)
    I = find(ranks==3,1);
    if opt.verbose > 0
        fprintf("Using P(%s)\n",string(v(I)))
    end
    % Set the variable for the coefficient space and set the variables in the
    % second group
    p_var = v(I);
    lin_vars = [v(setdiff(1:3,I));1];
    % Build the polynomial matrix P based on the choice of matrix beforehand
    P = FF([...
        C(:,p_idx(I,1)).*p_var+C(:,p_idx(I,2)),...
        C(:,p_idx(I,3)).*p_var+C(:,p_idx(I,4)),...
        C(:,p_idx(I,5)).*p_var^2+C(:,p_idx(I,6)).*p_var+C(:,10)...
        ],prime);
    P2 = inv(Q{I})*P; %#ok<MINV>
else
    error("All Matrices are singular!")
end