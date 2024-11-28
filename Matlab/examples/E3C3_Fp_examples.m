x = sym("x","real");
y = sym("y","real");
z = sym("z","real");
var_vec = [x.^3,x.^2.*y,x.^2.*z,x.*y.^2,...
    x.*y.*z,x.*z.^2,y.^3,y.^2.*z,y.*z.^2,z.^3,x.^2,x.*y,x.*z,y.^2,y.*z,z.^2,x,y,z,1]';
%%
prime = 5;
c = [...
    4,	0,	0,	0,	0,	0,	4,	0,	0,	0,	4,	0,	0,	0,	0,	0,	1,	0,	0,	2;...
    4,	0,	0,	0,	0,	0,	0,	4,	0,	0,	1,	0,	0,	0,	0,	0,	2,	0,	0,	3;...
    3,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	4,	0,	2,	0,	0,	2;...
    ];
E3C3_Fp(c,prime);
latex(c*var_vec)
%%
prime = 13;
c = [...
    0,	0,	4,	0,	0,	0,	7,	0,	0,	0,	0,	0,	6,	2,	0,	0,	0,	4,	0,	7;...
    -1,	0,	2,	0,	0,	0,	4,	0,	0,	0,	0,	0,	-4,	10,	0,	0,0,	5,	0,	-2;...
    0,	0,	7,	0,	0,	0,	8,	0,	0,	0,	0,	0,	2,	-3,	0,	0,	0,	0,	0,	-3;...
    ];
E3C3_Fp(c,prime);
latex(c*var_vec==0)
%%
prime = 67;
c = [39	0	17	0	0	0	12	0	0	0	0	0	66	48	0	0	0	23	0	60;...
    27	0	50	0	0	0	52	0	0	0	0	0	53	33	0	0	0	4	0	12;...
    20	0	66	0	0	0	13	0	0	0	0	0	28	54	0	0	0	39	0	28];
E3C3_Fp(c,prime);
latex(c*var_vec==0)
%%
prime = nextprime(280);
c = [121	0	136	0	0	0	41	0	0	0	0	0	143	142	0	0	0	21	0	118;...
    134	0	94	0	0	0	81	0	0	0	0	0	23	72	0	0	0	62	0	142;...
    18	0	14	0	0	0	142	0	0	0	0	0	144	119	0	0	0	136	0	97;...
    ];
E3C3_Fp(c,prime);
latex(c*var_vec==0)
%%
prime = 5;
c = [2	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	1	0	0	2;...
    2	0	0	0	0	0	0	-1	0	0	2	0	0	0	0	0	2	0	0	2;...
    1	0	0	0	0	0	0	0	0	0	2	0	0	0	-1	0	3	0	0	3];
E3C3_Fp(c,prime);
latex(c*var_vec==0)
%%
prime = 100003;
c = [1,	0,	0,	0,	0,	0,	5,	99992,	0,	0,	55588,	0,	0,	0,	92749,	0,	77685,	0,	0,	90961;...
7,	0,	0,	0,	0,	0,	1,	1,	0,	0,	89101,	0,	0,	0,	18839,	0,	21432,	0,	0,	51045;...
2,	0,	0,	0,	0,	0,	100001,	5,	0,	0,	11172,	0,	0,	0,	94208,	0,	18308,	0,	0,	73380];
E3C3_Fp(c,prime);
latex(c*var_vec==0)