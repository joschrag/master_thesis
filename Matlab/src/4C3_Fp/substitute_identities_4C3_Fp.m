function A = substitute_identities_4C3_Fp(P2,lin_vars,prime)
%substitute_identities_4C3_Fp Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 FF {mustBeSizeFF(P2,[4,6])};
    lin_vars (6,1) sym;
    prime (1,1) {mustBePrime};
end
cube_1 = FF(P2.value(1,:)*lin_vars,prime);   %y^3
quad1 = FF(P2.value(2,:)*lin_vars,prime);   %y^2*z
quad2 = FF(P2.value(3,:)*lin_vars,prime);   %y*z^2
cube_2 = FF(P2.value(4,:)*lin_vars,prime);  %z^3
% Create identities and directly perform first substitution
identities = [cube_1*cube_2 - quad1*quad2;...
    (cube_1)*FF(lin_vars(5)^2,prime) - (quad2)*FF(lin_vars(4)^2,prime);...                % (z^3)*y = (z^2*y)*z
    (cube_2)*FF(lin_vars(4)^2,prime) - (quad1)*FF(lin_vars(5)^2,prime);...                % (y^2*z)*z = (z^2*y)*y
    (cube_2)*FF(lin_vars(4),prime) - (quad2)*FF(lin_vars(5),prime);...
    (quad2)*cube_1 - (quad1)^2;...
    (quad1)*cube_2 - (quad2)^2];% (y*z)*z = (z^2*y)
old_vars = [lin_vars(4)^3;lin_vars(5)^3;lin_vars(4)^2*lin_vars(5);lin_vars(4)*lin_vars(5)^2];
new_vars = collect([cube_1.value;cube_2.value;quad1.value;quad2.value],lin_vars(4:5));
for i = 1:5
    identities(i) = subs(collect(identities(i)),old_vars,new_vars);
end
% Subs nonlinear terms to construct matrix-form
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
sub_final = subs(get_value(identities),lin_vars(1:3),[u;v;w]);
% Create matrix form from equations
[A_,b] = equationsToMatrix(sub_final,[u;v;w;lin_vars(4:5)]);
A = FF([A_,-b],prime);
end