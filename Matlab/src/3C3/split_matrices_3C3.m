function [lin_vars, p_var, P2] = split_matrices_3C3(c, x, y, z, opt)
%SPLIT_MATRICES_3C3 Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_3C3(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices for all choices of variable choices. The product of
%   P2 = inv(Q)*P is computed. If it passes the criteria, the matrix and the
%   combinations of variables for the approach are returned.
arguments
    c (3,20) {mustBeReal};
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 1;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
end
vars = [x,y,z];
complete_idx = {9,10,[4,14],[6,16],[2,12,18],[3,13,19],[1,11,17,20];...
    7,8,[4,14],[6,16],[2,12,18],[3,13,19],[1,11,17,20];...
    6,10,[2,11],[9,16],[4,12,17],[8,15,19],[7,14,18,20];...
    1,3,[2,11],[9,16],[4,12,17],[8,15,19],[7,14,18,20];...
    4,7,[3,11],[8,14],[6,13,17],[9,15,18],[10,16,19,20];...
    1,2,[3,11],[8,14],[6,13,17],[9,15,18],[10,16,19,20]};
q_idx = [7,8,5,15;...
    9,10,5,15;...
    1,3,5,13;...
    6,10,5,13;...
    1,2,5,12;...
    4,7,5,12];    %indices of cubic monomials of lin_vars
var_ind_vec = [1:3;1,3,2;2,1,3;2,3,1;3,1,2;3,2,1];
assert(all(c(:,5)==zeros(3,1)))
Q = cell(1,6);
ranks = zeros(1,6);
conds = zeros(1,6);
var_idx = zeros(6,3);
for i=1:6
    Q{i} = -[c(:,q_idx(i,1)),c(:,q_idx(i,2)),c(:,q_idx(i,3))*x+c(:,q_idx(i,4))];
    ranks(i) = rank(Q{i});
    conds(i) = real(cond(Q{i}));
    var_idx(i,:) = var_ind_vec(i,:);
end
valid_P2 = false;

if any(ranks==3)
    for I = find(ranks==3)
        p_var = vars(var_idx(I,1));
        lin_var_idx = var_idx(I,2:3);
        all_lin_vars = [vars(lin_var_idx(1))*vars(lin_var_idx(2)).^2,vars(lin_var_idx(2)).^3,vars(lin_var_idx(1))^2,vars(lin_var_idx(2)).^2,vars(lin_var_idx),1]';
        p_var_pow = [1;p_var;p_var^2;p_var^3];
        P = sym.zeros(3,7);
        for j=1:7
            P(:,j) = c(:,complete_idx{I,j})*p_var_pow(numel(complete_idx{I,j}):-1:1);
        end
        lin_vars = all_lin_vars;
        P2 = Q{I}\P;
        if all(all(P2(:,1:6)==sym(0)))
            if opt.verbose > 0
                fprintf("Using P(%s) with G_1 = {%s,%s,%s}\n",string(p_var),vars(lin_var_idx(1)).^3,vars(lin_var_idx(1)).^2*vars(lin_var_idx(2)),vars(lin_var_idx(1))*vars(lin_var_idx(2)))
            end
            valid_P2 = true;
            break
        end
    end
else
    error("All Matrices are singular!")
end

if ~valid_P2
    error("Substitution not possible!")
end
end