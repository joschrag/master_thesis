x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';

%%
prime = 5;
c = [...
4.0,   0, 2.0, 2.0, 4.0, 1.0, 1.0, 5.0, 5.0,   0;...
2.0, 3.0, 2.0, 4.0, 3.0, 1.0, 3.0, 4.0, 4.0, 1.0;...
  0, 1.0, 1.0, 1.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0;...
        ];
E3Q3_Fp(c,prime);
latex(simplify(FF(c*v,prime)).value)
sol_check(c,prime)