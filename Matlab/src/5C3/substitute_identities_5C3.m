function A = substitute_identities_5C3(P2,lin_vars)
%SUBSTITUTE_IDENTITIES_5C3 Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 (5,5) {mustBeReal};
    lin_vars (5,1) sym;
end
cube_1 = P2(1,:)*lin_vars;   %y^3
quad12 = P2(2,:)*lin_vars;   %y^2*z
quad21 = P2(3,:)*lin_vars;   %y*z^2
cube_2 = P2(4,:)*lin_vars;  %z^3
mixed = P2(5,:)*lin_vars;    %y*z
% Create identities and directly perform first substitution
identities = [(cube_1)*lin_vars(4) - (quad12)*lin_vars(3);...      % (y^3*)z = (y^2*z)*y
    (cube_2)*lin_vars(3) - (quad21)*lin_vars(4);...                % (z^3)*y = (z^2*y)*z
    (quad21)*lin_vars(3) - (quad12)*lin_vars(4);...                % (y^2*z)*z = (z^2*y)*y
    (quad12) - (mixed)*lin_vars(3);...                     % (y^2*z) = (y*z)*y
    (mixed)*lin_vars(4) - (quad21)];                      % (y*z)*z = (z^2*y)
old_vars = [lin_vars(3)^3,lin_vars(4)^3,lin_vars(3)^2*lin_vars(4),lin_vars(3)*lin_vars(4)^2,lin_vars(3)*lin_vars(4)];
new_vars = collect([cube_1,cube_2,quad12,quad21,mixed],lin_vars(3:4));
% Collect the identities terms and perform second and last substitution
sub_new = expand(subs(expand(collect(identities)),old_vars,new_vars));
% Subs nonlinear terms to construct matrix-form
u = sym("u","real");
v = sym("v","real");
sub_final = subs(sub_new,lin_vars(1:2),[u;v]);
% Create matrix form from equations
[A_,b] = equationsToMatrix(sub_final,[u;v;lin_vars(3:4)]);
A = [A_,-b];
end