function [ev_count,d0,d1,const] = setup_test(coeff_vec,tolerance)
arguments
    coeff_vec (1,10) {mustBeReal} = [1,1,1,1,1,1,0,0,0,-5];
    tolerance (1,1) {mustBePositive,mustBeReal} = 10^-10;
end
A = [coeff_vec(1),coeff_vec(4)/2,coeff_vec(5)/2;...
    coeff_vec(4)/2,coeff_vec(2),coeff_vec(6)/2;...
    coeff_vec(5)/2,coeff_vec(6)/2,coeff_vec(3)];
B = [coeff_vec(7),coeff_vec(8),coeff_vec(9)];
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
var = [x;y;z];
fun1 = matlabFunction(simplify(expand(var'*A*var)),'Vars',{[x;y;z]});
fun2 = matlabFunction(simplify(expand(dot(B,var))),'Vars',{[x;y;z]});
[V,D] = eig(A);
V = [V(:,1)./norm(V(:,1)),V(:,2)./norm(V(:,2)),V(:,3)./norm(V(:,3))];
eig_vals = diag(D);
eig_vals = (abs(eig_vals) > tolerance).*eig_vals;
if abs(det(V)+1) < tolerance
    V = [V(:,2),V(:,1),V(:,3)];
    eig_vals = [eig_vals(2),eig_vals(1),eig_vals(3)];
end
lambda1 = eig_vals(1);
lambda2 = eig_vals(2);
lambda3 = eig_vals(3);
pos_eig_vals = sum(eig_vals > 0);
neg_eig_vals = sum(eig_vals < 0);
zero_eig_vals = 3-pos_eig_vals-neg_eig_vals;
xi = sym("xi","real");
eta = sym("eta","real");
zeta = sym("zeta","real");
new_vars = [xi;eta;zeta];
const = coeff_vec(10);
if any(B~=0)
    lhs_lin = simplify(fun2(V*new_vars));
    if lambda1 ~= 0
        [tmp1,~] = coeffs(lhs_lin,xi,"All");
        if numel(tmp1) > 1
            factor = double(tmp1(1));
            [~,const] = square_complete(lambda1,factor,const);
        end
    end
    if lambda2 ~= 0
        [tmp1,~] = coeffs(lhs_lin,eta,"All");
        if numel(tmp1) > 1
            factor = double(tmp1(1));
            [~,const] = square_complete(lambda2,factor,const);
        end
    end
    if lambda3 ~= 0
        [tmp1,~] = coeffs(lhs_lin,zeta,"All");
        if numel(tmp1) > 1
            factor = double(tmp1(1));
            [~,const] = square_complete(lambda3,factor,const);
        end
    end
end
const = (abs(const)>tolerance)*const;
[d0, d1] = detect_linear_factors(A,B, V);
ev_count = {pos_eig_vals,neg_eig_vals,zero_eig_vals};
if const > 0
    const = -const;
    d0 = -d0;
    d1 = -d1;
    ev_count = {neg_eig_vals,pos_eig_vals,zero_eig_vals};
end
end

