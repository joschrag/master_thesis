x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
v = [x^2,x*y,x*z,y^2,y*z,z^2,x,y,z,1]';
%%
prime = 17;
c = [...
        1,1,1,2,0,1,1,1,1,-10;...
        1,0,1,1,-1,3,7,0,0,-10;...
        1,0,0,1,1,1,-15,0,0,-10;...
        ];
E3Q3_Fp(c,prime);
latex(simplify(FF(c*v,prime)).value)
%%
prime=69427;
c = [...
29236, 0,    0,     0, 52080,     0,  5262, 13670, 29511, 64383;...
    0, 0, 9520,     0, 16456, 36592, 63338,  6437, 65808, 20861;...
    0, 0,    0, 40362,     0,     0, 17093, 18928,  1294, 29534;...
        ];
E3Q3_Fp(c,prime);
latex(simplify(FF(c*v,prime)).value)
%%
prime=nextprime(10^6);
c = [...
29236, 0,    0,     0, 52080,     0,  5262, 13670, 29511, 64383;...
    0, 0, 9520,     0, 16456, 36592, 63338,  6437, 65808, 20861;...
    0, 0,    0, 40362,     0,     0, 17093, 18928,  1294, 29534;...
        ];
E3Q3_Fp(c,prime);
latex(simplify(FF(c*v,prime)).value)