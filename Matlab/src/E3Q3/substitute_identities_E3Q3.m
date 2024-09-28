function A = substitute_identities_E3Q3(P2,lin_vars)
%SUBSTITUTE_IDENTITIES_E3Q3 Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in expressions from lin_vars. After two substitutions
%   the new equations are written into matrix form and the coefficient matrix A is returned.
arguments
    P2 (3,3) {mustBeReal};
    lin_vars (3,1) sym;
end
% Compute substitution values
quad_1 = P2(1,:)*lin_vars;  %y^2
mixed = P2(2,:)*lin_vars;   %y*z
quad_2 = P2(3,:)*lin_vars;  %z^2
% Create identities and directly perform first substitution
identities = [(quad_1)*lin_vars(2) == (mixed)*lin_vars(1);...
    (mixed)*lin_vars(2) == (quad_2)*lin_vars(1);...
    (mixed)*(mixed) == quad_1*quad_2];
% Collect the identities terms and perform second and last substitution
old_vars = [lin_vars(1)^2,lin_vars(2)^2,lin_vars(1)*lin_vars(2)];
new_vars = [quad_1,quad_2,mixed];
identities = subs(expand(collect(identities)),old_vars,new_vars);
% Create matrix form from equations
[A_,b] = equationsToMatrix(identities,lin_vars(1:2));
A = [A_,-b];
assert(rank(A)==3)
end