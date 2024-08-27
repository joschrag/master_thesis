function [lin_vars, p_var, P2] = split_matrices(c, m, x, y, z, options)
%SPLIT_MATRICES Splits a matrix into linear and polynomial parts.
%
%   [LIN_VARS, P_VAR, P2] = SPLIT_MATRICES(C, M, X, Y, Z) splits the matrix C
%   into linear and polynomial parts using one of three possible methods based on
%   the size of C. The output LIN_VARS is a vector of the linear variables used in
%   the split, P_VAR is the polynomial variable used in the split, and P2 is the
%   resulting polynomial matrix. X, Y, and Z are optional symbolic input variables
%   that default to "x", "y", and "z" respectively.
%
arguments
    c (:,:) {mustBeReal}; % Coefficient matrix
    m (1,1) {mustBeInteger,mustBePositive} = size(c,1); % Number of rows in coefficient matrix
    x (1,1) sym = sym("x","real"); % Symbolic variable for x-direction
    y (1,1) sym = sym("y","real"); % Symbolic variable for y-direction
    z (1,1) sym = sym("z","real"); % Symbolic variable for z-direction
    options.verbose (1,1) {mustBeInteger, mustBeInRange(options.verbose,0,2)} = 0;
end
assert(size(c,2)==10 || size(c,2)==20,"Coefficient matrix has invalid size!")
if size(c,2)==10
    A = -[c(:,4),c(:,5),c(:,6)];
    B = -[c(:,1),c(:,3),c(:,5)];
    C = -[c(:,1),c(:,2),c(:,4)];
    min_num_rows = 3;
    i =   [2,8,3,9,1,7;...
        2,7,5,9,4,8;...
        3,7,5,8,6,9];
else
    assert(~any(c(:,5)))
    A = -[c(:,7),c(:,8),c(:,9),c(:,10),c(:,5)*x+c(:,15)];
    B = -[c(:,1),c(:,3),c(:,6),c(:,10),c(:,5)*y+c(:,13)];
    C = -[c(:,1),c(:,2),c(:,4),c(:,7),c(:,5)*z+c(:,12)];
    min_num_rows=5;
    i =   [ 4,14,6,16,2,12,18,3,13,19,1,11,17;...
        2,11,9,16,4,12,17,8,15,19,7,14,18;...
        3,11,8,14,6,13,17,9,15,18,10,16,19];
end
Q_list = {A,B,C};
conds = real([vpa(cond(A)),vpa(cond(B)),vpa(cond(C))]);
[M,I] = min(conds,[],"omitnan");


v = [x;y;z];
if m==min_num_rows
    if ~isinf(M)
        if options.verbose > 0
            fprintf("Using P(%s)\n",string(v(I)))
        end
        Q = Q_list{I};
        p_var = v(I);
        if min_num_rows == 3
            P = [c(:,i(I,1)).*p_var+c(:,i(I,2)),c(:,i(I,3)).*p_var+c(:,i(I,4)),c(:,i(I,5)).*p_var^2+c(:,i(I,6)).*p_var+c(:,10)];
            lin_vars = [v(setdiff(1:3,I));1];
        else
            P  = [c(:,i(I,1)).*p_var+c(:,i(I,2)),...
                c(:,i(I,3)).*p_var+c(:,i(I,4)),...
                c(:,i(I,5)).*p_var^2+c(:,i(I,6)).*p_var+c(:,i(I,7)),...
                c(:,i(I,8)).*p_var^2+c(:,i(I,9)).*p_var+c(:,i(I,10)),...
                c(:,i(I,11)).*p_var^3+c(:,i(I,12)).*p_var^2+c(:,i(I,13)).*p_var+c(:,20)];
            lin_vars = [v(setdiff(1:3,I)).^2;v(setdiff(1:3,I));1];
        end
        P2 = Q\P;
    else
        warning("All Matrices are singular!")
    end
else
    if options.verbose > 0
        fprintf("Using MP-inverse\n")
    end
    if min_num_rows==3
        P_X = [c(:,2).*x+c(:,8),c(:,3).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
        lin_vars = [y;z;1];
    else
        P_X =  [c(:,4).*x+c(:,14),...
            c(:,6).*x+c(:,16),...
            c(:,2).*x^2+c(:,12).*x+c(:,18),...
            c(:,3).*x^2+c(:,13).*x+c(:,19),...
            c(:,1).*x^3+c(:,11).*x^2+c(:,17).*x+c(:,20)];
        lin_vars = [y^2;z^2;y;z;1];
    end
    P2 = pinv(A)*P_X;
    p_var = x;
end
