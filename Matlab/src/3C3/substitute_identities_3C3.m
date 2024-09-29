function A = substitute_identities_3C3(P2,lin_vars)
%SUBSTITUTE_IDENTITIES3C3 Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 (3,7) {mustBeReal};
    lin_vars (7,1) sym;
end
cube = P2(1,:)*lin_vars;  %y^3
quad = P2(2,:)*lin_vars;   %y^2*z
mixed = P2(3,:)*lin_vars;   %y*z
%Create identities and directly perform first substitution
identities = [...
    cube*lin_vars(6)-quad*lin_vars(5);...
    cube*lin_vars(6) - mixed*lin_vars(5)^2;...
    quad*lin_vars(6)-mixed^2;...
    cube*lin_vars(6)^3-mixed^3;...
    quad-mixed*lin_vars(5);...
    cube*lin_vars(6)^2-quad*mixed;...
    ];
%Collect the identities terms and perform second and last substitution
sub_new = collect(identities);
old_vars = [lin_vars(5)^3,lin_vars(5)^2*lin_vars(6),lin_vars(5)*lin_vars(6)];
new_vars = collect([cube,quad,mixed],lin_vars(5:6));
sub_new = expand(subs(expand(sub_new),old_vars,new_vars));
%Create new variables for expressions of higher degree to "linearize" system
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
%Substitute the nonlinear expressions with new variables
sub_final = subs(sub_new,lin_vars(2:4),[u;v;w]);
[A_,b] = equationsToMatrix(sub_final,[u;v;w;lin_vars(5:6)]);
A = [A_,-b];
end