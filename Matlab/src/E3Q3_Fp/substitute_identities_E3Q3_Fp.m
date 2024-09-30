function A = substitute_identities_E3Q3_Fp(P2,lin_vars,prime)
%SUBSTITUTE_IDENTITIES_5C3_FP Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 FF {mustBeSizeFF(P2,[3,3])};
    lin_vars (3,1) sym;
    prime (1,1) {mustBePrime};
end
% Compute substitution values
quad_1 = FF(P2.value(1,:),prime)*FF(lin_vars,prime);
quad_2 = FF(P2.value(3,:),prime)*FF(lin_vars,prime);
mixed = FF(P2.value(2,:),prime)*FF(lin_vars,prime);
% Create identities and directly perform first substitution
identities = [(quad_1)*FF(lin_vars(2),prime) - (mixed)*FF(lin_vars(1),prime);...
    (mixed)*FF(lin_vars(2),prime) - (quad_2)*FF(lin_vars(1),prime);...
    (mixed)*(mixed) - quad_1*quad_2];
% Collect the identities terms and perform second and last substitution
old_vars = [lin_vars(1)^2,...
    lin_vars(2)^2,...
    lin_vars(1)*lin_vars(2)];
new_vars = [quad_1.value,quad_2.value,mixed.value];
for i = 1:3
    identities(i) = subs(collect(identities(i)),old_vars,new_vars);
end
% Create matrix form from equations
[A_,b] = equationsToMatrix(get_value(identities),lin_vars(1:2));
A = FF([A_,-b],prime);
end