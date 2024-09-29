function [lin_vars, p_var, P2] = split_matrices_4C3_Fp(C, prime, x, y, z, opt)
%SPLIT_MATRICES_4C3_FP Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_4C3_FP(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices for all choices of variable choices. The product of
%   P2 = inv(Q)*P is computed. If it passes the criteria, the matrix and the
%   combinations of variables for the approach are returned.
arguments
    C (4,20) {mustBeReal};
    prime (1,1) {mustBePrime};
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 1;
end
vars = [x,y,z];
complete_idx = {[4,14],[5,15],[6,16],[2,12,18],[3,13,19],[1,11,17,20];...
    [2,11],[5,13],[9,16],[4,12,17],[8,15,19],[7,14,18,20];...
    [3,11],[5,12],[8,14],[6,13,17],[9,15,18],[10,16,19,20]};
q_idx = [7,8,9,10;1,3,6,10;1,2,4,7];    %indices of cubic monomials of lin_vars
Q = cell(1,3);
ranks = zeros(1,3);
var_idx = [1,2,3;2,1,3;3,1,2];
for i=1:3
    Q{i} = FF(-[C(:,q_idx(i,1)),C(:,q_idx(i,2)),C(:,q_idx(i,3)),C(:,q_idx(i,4))],prime);
    ranks(i) = rank(Q{i});
end
if any(ranks==4)
    for I = find(ranks==4)
        if all(C(:,[complete_idx{I,1:3}])==0)
            p_var = vars(var_idx(I,1));
            lin_var_idx = var_idx(I,2:3);
            if opt.verbose > 0
                fprintf("Using P(%s) with G_1 = {%s,%s,%s}\n",string(vars(I)),vars(lin_var_idx(1)).^3,vars(lin_var_idx(1)).^2*vars(lin_var_idx(2)),vars(lin_var_idx(1))*vars(lin_var_idx(2)))
            end
            lin_vars = [vars(lin_var_idx(1))^2,vars(lin_var_idx(1))*vars(lin_var_idx(2)),vars(lin_var_idx(2)).^2,vars(lin_var_idx),1]';
            p_var_pow = [1;p_var;p_var^2;p_var^3];
            P = sym.zeros(4,6);
            for j=1:6
                P(:,j) = C(:,complete_idx{I,j})*p_var_pow(numel(complete_idx{I,j}):-1:1);
            end
            P2 = inv(Q{I})*FF(P,prime); %#ok<MINV>
            return
        end
    end
else
    error("All Matrices are singular!")
end
error("Substitution not possible!")
end