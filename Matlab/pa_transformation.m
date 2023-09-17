function pa_transformation(coeff_vec,tol)
%PA_TRANSFORMATION Summary of this function goes here
%   Detailed explanation goes here
%   arguments:
%       coeffs: coefficient vector of (x^2,y^2,z^2,xy,xz,yz,x,y,z,1)
arguments
    coeff_vec (1,10) {mustBeReal} = ones(1,10);
    tol (1,1) {mustBeReal,mustBePositive} = 10^-6;
end
% if abs(coeff_vec(10)) >= tol
%     coeff_vec = coeff_vec./(-coeff_vec(10));
%     rhs = 1;
% else
%     rhs = 0;
% end
M = [coeff_vec(1),coeff_vec(4)/2,coeff_vec(5)/2,coeff_vec(7)/2;...
    coeff_vec(4)/2,coeff_vec(2),coeff_vec(6)/2,coeff_vec(8)/2;...
    coeff_vec(5)/2,coeff_vec(6)/2,coeff_vec(3),coeff_vec(9)/2;...
    coeff_vec(7)/2,coeff_vec(8)/2,coeff_vec(9)/2,coeff_vec(10)];
A = [coeff_vec(1),coeff_vec(4)/2,coeff_vec(5)/2;...
    coeff_vec(4)/2,coeff_vec(2),coeff_vec(6)/2;...
    coeff_vec(5)/2,coeff_vec(6)/2,coeff_vec(3)];
B = [coeff_vec(7),coeff_vec(8),coeff_vec(9)];
x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
var = [x,y,z,1];
fun1 = matlabFunction(simplify(expand(var(1:3)*A*var(1:3)')),'Vars',{[x;y;z]});
fun2 = matlabFunction(simplify(expand(dot(B,var(1:3)))),'Vars',{[x;y;z]});

[V,D] = eig(A);
V = [V(:,1)./norm(V(:,1)),V(:,2)./norm(V(:,2)),V(:,3)./norm(V(:,3))];
eig_vals = diag(D);
if round(det(V)) == -1
    V = [V(:,2),V(:,1),V(:,3)];
    eig_vals = [eig_vals(2),eig_vals(1),eig_vals(3)];
end
pos_eig_vals = sum(eig_vals > 0);
neg_eig_vals = sum(eig_vals < 0);
xi = sym("xi","real");
eta = sym("eta","real");
zeta = sym("zeta","real");
lhs_quad = simplify(fun1(V*[xi;eta;zeta]));
const = coeff_vec(10);
offsets = zeros(3,1);
if any(B~=0)
    lhs_lin = simplify(fun2(V*[xi;eta;zeta]));
    [tmp1,~] = coeffs(lhs_lin,[xi,eta,zeta]);
    %Quadratische Ergänzung
    [offsets(1),const] = square_complete(eig_vals(1),tmp1(1),const);
    [offsets(2),const] = square_complete(eig_vals(2),tmp1(2),const);
    [offsets(3),const] = square_complete(eig_vals(3),tmp1(3),const);
end
tmp = 1./(sqrt(abs(eig_vals./(-const))));
a = double(tmp(1));
b = double(tmp(2));
c = double(tmp(3));

if all(eig_vals ~= 0)
    plot_quadric_3d(pos_eig_vals,neg_eig_vals,a,b,c,V,offsets);
end
end