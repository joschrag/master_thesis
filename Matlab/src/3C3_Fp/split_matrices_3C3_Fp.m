function [lin_vars, p_var, P2] = split_matrices_3C3_Fp(C,prime, x, y, z, opt)
%SPLIT_MATRICES_3C3_FP Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_3C3_FP(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices Q and P for all variable choices. If Q is invertable,
%   the matrix and the combinations of variables for the approach are returned.
arguments
    C (3,20) {mustBeReal}; % Coefficient matrix
    prime (1,1) {mustBePrimeOrZero,mustBePositive};
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 0;
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
assert(all(mod(C(:,5),prime)==0))
Q = cell(1,6);
ranks = zeros(1,6);
var_idx = zeros(6,3);

for i=1:6
    Q{i} = FF(-[C(:,q_idx(i,1)),C(:,q_idx(i,2)),C(:,q_idx(i,3))*x+C(:,q_idx(i,4))],prime);
    ranks(i) = rank(Q{i});
    var_idx(i,:) = var_ind_vec(i,:);
end
if any(ranks==3)
    for I = find(ranks==3)
        % Check if coefficient conditions for algorithm are met
        if all(C(:,[complete_idx{I,1:6}])==0)
            % Set the variable for the coefficient space and set the variables in the
            % second group
            p_var = vars(var_idx(I,1));
            lin_var_idx = var_idx(I,2:3);
            if opt.verbose > 0
                fprintf("Using P(%s) with G_1 = {%s,%s,%s}\n",string(vars(I)),vars(lin_var_idx(1)).^3,vars(lin_var_idx(1)).^2*vars(lin_var_idx(2)),vars(lin_var_idx(1))*vars(lin_var_idx(2)))
            end
            lin_vars = [vars(lin_var_idx(1))*vars(lin_var_idx(2))^2,vars(lin_var_idx(2)).^3,vars(lin_var_idx(1))^2,vars(lin_var_idx(2)).^2,vars(lin_var_idx),1]';
            p_var_pow = [1;p_var;p_var^2;p_var^3];
            % Build the polynomial matrix P based on the choice of matrix beforehand
            P = sym.zeros(3,7);
            for j=1:7
                P(:,j) = C(:,complete_idx{I,j})*p_var_pow(numel(complete_idx{I,j}):-1:1);
            end
            if opt.verbose > 1
                disp(lin_vars)
                disp(P)
            end
            P2 = inv(Q{I})*FF(P,prime); %#ok<MINV>
            if opt.verbose > 1
                fprintf("P2:\n")
                disp(P2.value)
            end
            return
        end
    end
else
    error("All matrices singular!")
end
error("No substitution possible!")
end