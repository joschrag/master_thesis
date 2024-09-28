function [lin_vars, p_var, P2] = split_matrices_E3Q3(C, x, y, z, opt)
%SPLIT_MATRICES_E3Q3 Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_E3Q3(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices Q and P for all variable choices. If Q is invertable,
%   the matrix and the combinations of variables for the approach are returned.
arguments(Input)
    C (3,10) {mustBeReal};
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 0;
end
arguments(Output)
    lin_vars (3,1) sym
    p_var (1,1) sym
    P2 (3,3) sym
end
q_idx = [4,5,6;1,3,5;1,2,4];
Q = cell(1,3);
ranks = zeros(1,3);
conds = zeros(1,3);

for i=1:3
    Q{i} = -[C(:,q_idx(i,:))];
    ranks(i) = rank(Q{i});
    conds(i) = cond(vpa(Q{i}));
end


p_idx =   [2,8,3,9,1,7;...
    2,7,5,9,4,8;...
    3,7,5,8,6,9];
v = [x;y;z];

% Compute condition of matrix candidates and choose the matrix with the lowest
% condition number
if any(ranks==3)
    [~,I] = min(conds,[],"omitnan");
    if opt.verbose > 0
        fprintf("Using P(%s)\n",string(v(I)))
    end
    % Set the variable for the coefficient space and set the variables in the
    % second group
    p_var = v(I);
    lin_vars = [v(setdiff(1:3,I));1];
    % Build the polynomial matrix P based on the choice of matrix beforehand
    P = [...
        C(:,p_idx(I,1)).*p_var+C(:,p_idx(I,2)),...
        C(:,p_idx(I,3)).*p_var+C(:,p_idx(I,4)),...
        C(:,p_idx(I,5)).*p_var^2+C(:,p_idx(I,6)).*p_var+C(:,10)...
        ];
    P2 = Q{I}\P;
    if opt.verbose >= 2
        fprintf("Q:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(Q{I}))
        fprintf("P:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(P))
        fprintf("P2:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(P2))
    end
else
    error("All Matrices are singular!")
end
