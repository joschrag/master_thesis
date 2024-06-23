function [lin_vars, p_var, P2] = split_matrices(c, m, x, y, z)
arguments
    c (:,10) {mustBeReal};
    m (1,1) {mustBeInteger,mustBePositive} = size(c,1);
    x (1,1) sym = sym("x","real");
    y (1,1) sym = sym("y","real");
    z (1,1) sym = sym("z","real");
end
A = -[c(:,2),c(:,3),c(:,6)];
B = -[c(:,1),c(:,3),c(:,5)];
C = -[c(:,1),c(:,2),c(:,4)];
Q_list = {A,B,C};
conds = [vpa(cond(A)),vpa(cond(B)),vpa(cond(C))];
[M,I] = min(conds,[],"omitnan");

i =   [4,8,5,9,1,7;...
    4,7,6,9,2,8;...
    5,7,6,8,3,9];
v = [x;y;z];
if m==3
    if ~isinf(M)
        fprintf("Using P(%s)\n",string(v(I)))
        Q = Q_list{I};
        P = [c(:,i(I,1)).*x+c(:,i(I,2)),c(:,i(I,3)).*x+c(:,i(I,4)),c(:,i(I,5)).*x^2+c(:,i(I,6)).*x+c(:,10)];
        lin_vars = [v(setdiff(1:3,I));1];
        p_var = v(I);
        P2 = Q\P;
    else
        warning("All Matrices are singular!")
    end
else
    fprintf("Using MP-inverse\n")
    P_X = [c(:,4).*x+c(:,8),c(:,5).*x+c(:,9),c(:,1).*x^2+c(:,7).*x+c(:,10)];
    P2 = pinv(A)*P_X; %idk
    lin_vars = [y;z;1];
    p_var = x;
end
