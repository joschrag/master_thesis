x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
%%
prime = 7;
c = [...
    3, 4, 4, 6, 0, 5, 2, 6, 5,   0, 5,   0, 2, 3, 6,   0, 3, 3,   0, 4;...
5, 6, 1, 4, 0, 4, 5, 6, 2, 1, 2, 6, 4, 6,   0,   0, 6, 5,   0, 1;...
2, 5, 5, 1, 0, 1, 3, 2, 5, 1, 1, 5, 3, 3, 2, 4,   0, 3, 3, 2;...
2, 5, 5, 4, 0,   0, 4,   0, 1, 6, 5,   0, 1, 2, 1, 2, 6,   0, 5, 2;...
3, 2, 4, 5, 0, 4, 4,   0, 4, 1,   0, 2, 4, 2, 1, 2, 2,   0,   0, 1];
E5C3_Fp(c,prime)
%%
prime = 5237;
c = [...
    2712, 3031, 3379, 5185, 0, 4293, 2163, 4589, 4314,  285, 1587, 2444, 3977, 3153, 1340, 3788, 5069, 3830, 1696, 1129;...
3756, 5206, 1338, 3515, 0, 3136, 3756, 4908, 1842, 1328, 3130, 4428, 2496, 4165, 4774, 1374, 3070, 1141, 1718, 4182;...
2086, 4037, 4005, 1486, 0,  992, 2467, 1753, 3847,  990, 4125,  636,   32, 1683, 4085, 1182, 1958, 1242,  823, 1007;...
2116, 3885, 4325, 3578, 0,  729, 3680,  308, 1008, 4843, 3265, 4943, 5035, 2183, 4650, 3669, 1831, 3705, 2656, 4034;...
2686, 1787, 3456, 4422, 0, 3162, 3075,  612, 3720, 1071,  988,  697, 2548, 3060, 4763, 3688, 3736,  717, 2547, 3014];
E5C3_Fp(c,prime)