function A = substitute_indentities_5C3_Fp(P2,lin_vars,prime)
%SUBSTITUTE_IDENTITIES_5C3_FP Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 FF {mustBeSizeFF(P2,[5,5])};
    lin_vars (5,1) sym;
    prime (1,1) {mustBePrimeOrZero,mustBePositive};
end
% Compute substitution values
cube_1 = FF(P2.value(1,:)*lin_vars,prime);   %y^3
quad1 = FF(P2.value(2,:)*lin_vars,prime);   %y^2*z
quad2 = FF(P2.value(3,:)*lin_vars,prime);   %y*z^2
cube_2 = FF(P2.value(4,:)*lin_vars,prime);   %z^3
mixed = FF(P2.value(5,:)*lin_vars,prime);    %y*z
% Create identities and directly perform first substitution
identities = [(cube_1)*FF(lin_vars(4),prime) - (quad1)*FF(lin_vars(3),prime);...      % (y^3*)z = (y^2*z)*y
    (cube_2)*FF(lin_vars(3),prime) - (quad2)*FF(lin_vars(4),prime);...                % (z^3)*y = (z^2*y)*z
    (quad2)*FF(lin_vars(3),prime) - (quad1)*FF(lin_vars(4),prime);...                % (y^2*z)*z = (z^2*y)*y
    (quad1) - (mixed)*FF(lin_vars(3),prime);...                     % (y^2*z) = (y*z)*y
    (mixed)*FF(lin_vars(4),prime) - (quad2)];                      % (y*z)*z = (z^2*y)

old_vars = [lin_vars(3)^3,...
    lin_vars(4)^3,...
    lin_vars(3)^2*lin_vars(4),...
    lin_vars(3)*lin_vars(4)^2,...
    lin_vars(3)*lin_vars(4)];
new_vars = [cube_1.value,cube_2.value,quad1.value,quad2.value,mixed.value];
% Collect the identities terms and perform second and last substitution
for i = 1:5
    identities(i) = subs(collect(identities(i)),old_vars,new_vars);
end
% Subs nonlinear terms to construct matrix-form
u = sym("u","real");
v = sym("v","real");
sub_final = subs(get_value(identities),lin_vars(1:2),[u;v]);
% Create matrix form from equations
[A_,b] = equationsToMatrix(sub_final,[u;v;lin_vars(3:4)]);
A = FF([A_,-b],prime);
end