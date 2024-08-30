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
%%
c  = [   0,   0,   0, -2.0, 0,    0,    0,   0,   0,    0, -1.0,  12.571428571428571428571428571429,   1.0,                            -183.25,                               -2.0,                                  0, -206.25510204081632653061224489796,  1162.3571428571428571428571428571, 96.910714285714285714285714285714, -10531.733099489796586567535996437;...
   0,   0, 2.0,  2.0, 0,    0, -1.0,   0,   0,    0, -10.5, -12.571428571428571428571428571429, 362.5,  190.67857142857142857142857142857,                                  0,                               -1.0, -1882.3698979591836177860386669636, -1170.9183673469387755102040816327,                       16438.28125, -84355.153493986887042410671710968;...
   0, 1.0,   0,  1.0, 0, -2.0,    0, 2.0, 1.0, -1.0, -3.1428571428571428571428571428571,  174.96428571428571428571428571429,  21.0,                             80.125, -24.071428571428571428571428571429, -168.64285714285714285714285714286, -614.89030612244903295504627749324,  7744.0602678571431169984862208366, 1876.3354591836734925891505554318, -29980.913903061224118573591113091;...
-2.0,   0,   0,    0, 0,    0,  2.0,   0,   0,    0,                            -542.75,                                  0,     0, -18.857142857142857142857142857143,                                  0,                                  0,                       -49096.09375,  58.265306122448979591836734693878,                                 0, -1480432.4797626640647649765014648;...
   0, 1.0,   0,    0, 0,    0,    0,   0, 1.0,    0, -3.1428571428571428571428571428571,                             180.25,     0,                                  0,                              -10.5, -3.1428571428571428571428571428571,                             -565.5,                        8149.828125,                              32.0, -25517.870535714286233996972441673];
E5C3(c);