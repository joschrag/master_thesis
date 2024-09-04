function [lin_vars, p_var, P2] = split_matrices_Fp(c,prime, m, x, y, z, opt)
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
    prime (1,1) {mustBeInteger,mustBePositive}
    m (1,1) {mustBeInteger,mustBePositive} = size(c,1); % Number of rows in coefficient matrix
    x (1,1) sym = sym("x","real"); % Symbolic variable for x-direction
    y (1,1) sym = sym("y","real"); % Symbolic variable for y-direction
    z (1,1) sym = sym("z","real"); % Symbolic variable for z-direction
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 0;
    opt.error (1,1) {mustBeNumericOrLogical} = true;
end
assert(size(c,2)==10 || size(c,2)==20,"Coefficient matrix has invalid size!")
if size(c,2)==10
    A = -FF([c(:,4),c(:,5),c(:,6)],prime);
    B = -FF([c(:,1),c(:,3),c(:,5)],prime);
    C = -FF([c(:,1),c(:,2),c(:,4)],prime);
    min_num_rows = 3;
    i =   [2,8,3,9,1,7;...
        2,7,5,9,4,8;...
        3,7,5,8,6,9];
    ranks = [rank(A),rank(B),rank(C)];
    Q_list = {A,B,C};
else
    min_num_rows=5;
    assert(~any(c(:,5)))
    ranks = zeros(1,3);
    idx = [7,8,9,10,5,15;1,3,6,10,5,13;1,2,4,7,5,12];
    Q_list = cell(1,3);
    i =   [ 4,14,6,16,2,12,18,3,13,19,1,11,17;...
        2,11,9,16,4,12,17,8,15,19,7,14,18;...
        3,11,8,14,6,13,17,9,15,18,10,16,19];
    for j=1:3
        Q_list{j} = -FF([c(:,idx(j,1)),c(:,idx(j,2)),c(:,idx(j,3)),c(:,idx(j,4)),c(:,idx(j,5))*x+c(:,idx(j,6))],prime);
        ranks(j) = rank(Q_list{j});

        if ranks(j) == min_num_rows
            break
        end
    end

end

if any(ranks==min_num_rows)
    M = 1;
    I = find(ranks==min_num_rows,1);
else
    M = Inf;
end
if m > min_num_rows
    inv_func =@(M) l_inv(M);
elseif m == min_num_rows
    inv_func =@(M) inv(M);
end

v = [x;y;z];
if ~isinf(M)
    if opt.verbose > 0
        fprintf("Using P(%s)\n",string(v(I)))
    end
    Q = Q_list{I};
    p_var = v(I);
    if min_num_rows == 3
        P = FF([c(:,i(I,1)).*p_var+c(:,i(I,2)),c(:,i(I,3)).*p_var+c(:,i(I,4)),c(:,i(I,5)).*p_var^2+c(:,i(I,6)).*p_var+c(:,10)],prime);
        lin_vars = FF([v(setdiff(1:3,I));1],prime);
    else
        P  = FF([c(:,i(I,1)).*p_var+c(:,i(I,2)),...
            c(:,i(I,3)).*p_var+c(:,i(I,4)),...
            c(:,i(I,5)).*p_var^2+c(:,i(I,6)).*p_var+c(:,i(I,7)),...
            c(:,i(I,8)).*p_var^2+c(:,i(I,9)).*p_var+c(:,i(I,10)),...
            c(:,i(I,11)).*p_var^3+c(:,i(I,12)).*p_var^2+c(:,i(I,13)).*p_var+c(:,20)],prime);
        lin_vars = FF([v(setdiff(1:3,I)).^2;v(setdiff(1:3,I));1],prime);
    end
    P2 = inv_func(Q)*P;
    if opt.verbose >= 2
        fprintf("Q:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(Q.value))
        fprintf("P:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(P.value))
        fprintf("P2:\t %s, %s, %s\n  \t%s, %s, %s\n  \t%s, %s, %s\n",string(P2.value))
    end
else
    if opt.error
        error("No matrix splitting possible.\n")
    else
        warning("No matrix splitting possible.\n")
    end
end
