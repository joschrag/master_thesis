x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
%%
c  = [1,-3,-2,0,0,0,-1,4,-1,-2,2,6,2,-1,-4,5,-1,-5,2,-2;...
        1,2,-1,-3,0,-1,0,1,-2,2,-4,2,4,3,4,-5,2,-8,2,1;...
        -2,2,1,-4,0,2,-5,0,2,2,1,2,-8,22,-6,-8,9,-20,13,-3;...
        -2,2,-2,-1,0,2,2,-2,-2,2,6,-2,1,-4,9,-5,-4,-4,-5,9;...
        0,-2,0,-1,0,0,-3,0,4,0,-2,6,1,6,-8,-2,3,2,3,-7;...
        ];
E5C3(c);
latex(c*var_vec==0)