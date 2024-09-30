function [lin_vars, p_var, P2] = split_matrices_5C3_Fp(C,prime, x, y, z, opt)
%SPLIT_MATRICES_5C3_FP Handles the matrix splitting and multiplication part of the algorithm.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES_5C3_FP(C, X, Y, Z, opt) splits the matrix C
%   into two submatrices Q and P for all variable choices. If Q is invertable,
%   the matrix and the combinations of variables for the approach are returned.
arguments
    C (5,20) {mustBeReal}; % Coefficient matrix
    prime (1,1) {mustBePrime};
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 0;
end

assert(~any(C(:,5)))
v = [x;y;z];
q_idx = [7,8,9,10,5,15;1,3,6,10,5,13;1,2,4,7,5,12];
p_idx =   [ 4,14,6,16,2,12,18,3,13,19,1,11,17;...
    2,11,9,16,4,12,17,8,15,19,7,14,18;...
    3,11,8,14,6,13,17,9,15,18,10,16,19];
% Compute rank of matrix candidates. If candidate is invertible compute other matrices
% and set variable for the coefficient space and set the variables in the second group
for I=1:3
    Q = -FF([C(:,q_idx(I,1)),C(:,q_idx(I,2)),C(:,q_idx(I,3)),C(:,q_idx(I,4)),C(:,q_idx(I,5))*x+C(:,q_idx(I,6))],prime);
    if rank(Q) == 5
        if opt.verbose > 0
            fprintf("Using P(%s)\n",string(v(I)))
        end
        % Set the variable for the coefficient space and set the variables in the
        % second group
        p_var = v(I);
        lin_vars = [v(setdiff(1:3,I)).^2;v(setdiff(1:3,I));1];
        % Build the polynomial matrix P based on the choice of matrix beforehand
        P  = FF([C(:,p_idx(I,1)).*p_var+C(:,p_idx(I,2)),...
            C(:,p_idx(I,3)).*p_var+C(:,p_idx(I,4)),...
            C(:,p_idx(I,5)).*p_var^2+C(:,p_idx(I,6)).*p_var+C(:,p_idx(I,7)),...
            C(:,p_idx(I,8)).*p_var^2+C(:,p_idx(I,9)).*p_var+C(:,p_idx(I,10)),...
            C(:,p_idx(I,11)).*p_var^3+C(:,p_idx(I,12)).*p_var^2+C(:,p_idx(I,13)).*p_var+C(:,20)],prime);
        P2 = inv(Q)*P;%#ok<MINV>
        return
    end
    error("No matrix splitting possible.\n")
end