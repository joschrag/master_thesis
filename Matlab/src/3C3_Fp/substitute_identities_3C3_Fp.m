function A = substitute_identities_3C3_Fp(P2,lin_vars,prime)
%SUBSTITUTE_IDENTITIES_5C3_FP Substitute expressions into identities to eliminate variables.
%   This function uses a set of identities to substitute monomials with linear
%   combinations in lin_vars. After two substitutions the new equations are
%   written into matrix form and the coefficient matrix A is returned.
arguments
    P2 FF {mustBeSizeFF(P2,[3,7])};
    lin_vars (7,1);
    prime (1,1) {mustBePrime};
end

cube = FF(P2.value(1,:)*lin_vars,prime);  %y^3
quad = FF(P2.value(2,:)*lin_vars,prime);   %y^3
mixed = FF(P2.value(3,:)*lin_vars,prime);   %y^3

identities = [...
    cube*FF(lin_vars(6),prime)-quad*FF(lin_vars(5),prime);...
    cube*FF(lin_vars(6),prime) - mixed*FF(lin_vars(5)^2,prime);...
    quad*FF(lin_vars(6),prime)-mixed^2;...
    cube*FF(lin_vars(6)^3,prime)-mixed^3;...
    quad-mixed*FF(lin_vars(5),prime);...
    cube*FF(lin_vars(6)^2,prime)-quad*mixed;...
    ];
for i=1:6
    identities(i) = collect(identities(i));
end
old_vars = [lin_vars(5)^3,lin_vars(5)^2*lin_vars(6),lin_vars(5)*lin_vars(6)];
new_vars = collect([cube.value,quad.value,mixed.value],lin_vars(5:6));
for i=1:6
    identities(i) = expand(subs(expand(identities(i)),old_vars,new_vars));
end
u = sym("u","real");
v = sym("v","real");
w = sym("w","real");
for i=1:6
    identities(i) = subs(identities(i),lin_vars(2:4),[u;v;w]);
end
eq = [identities(1).value,identities(2).value,identities(3).value,identities(4).value,identities(5).value,identities(6).value];
[A_,b] = equationsToMatrix(eq,[u;v;w;lin_vars(5:6)]);
A = [A_,-b];
end