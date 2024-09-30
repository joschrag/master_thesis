function A = substitute_identities_4C3(P2,lin_vars)
%SUBSTITUTE_IDENTITIES_5C3 Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 (4,6) {mustBeReal};
    lin_vars (6,1) sym;
end

cube_1 = P2(1,:)*lin_vars;   %y^3
quad1 = P2(2,:)*lin_vars;   %y^2*z
quad2 = P2(3,:)*lin_vars;   %y*z^2
cube_2 = P2(4,:)*lin_vars;  %z^3
% Create identities and directly perform first substitution
identities = [cube_1*cube_2 - quad1*quad2;...               % (z^3)*y^3 = (z^2*y)*(y^2*z)
    (cube_2)*lin_vars(4)^2 - (quad1)*lin_vars(5)^2;...     % (z^3)*y^2 = (y^2*z)*z^2
    (cube_1)*lin_vars(5)^2 - (quad2)*lin_vars(4)^2;...      % (y^3)*z^2 = (z^2*y)*y^2
    (cube_2)*lin_vars(4) - (quad2)*lin_vars(5);...          % (z^3)*y = (z^2*y)*z
    (quad2)*cube_1 - (quad1)^2;...                           % (z^3)*y^2 = (y^2*z)*z^2
    (quad1)*cube_2 - (quad2)^2];                            % (y^2*z)*z = (z^2*y)*y
old_vars = [lin_vars(4)^3;lin_vars(5)^3;lin_vars(4)^2*lin_vars(5);lin_vars(4)*lin_vars(5)^2];
new_vars = collect([cube_1;cube_2;quad1;quad2],lin_vars(4:5));
% Collect the identities terms and perform second and last substitution
sub_new = expand(subs(expand(collect(identities)),old_vars,new_vars));
% Subs nonlinear terms to construct matrix-form
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
sub_final = subs(sub_new,lin_vars(1:3),[u;v;w]);
% Create matrix form from equations
[A_,b] = equationsToMatrix(sub_final,[u;v;w;lin_vars(4:5)]);
A = [A_,-b];
end