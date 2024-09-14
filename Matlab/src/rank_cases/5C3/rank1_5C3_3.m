function        [u_sol,v_sol] = rank1_5C3_3(r)
arguments
    r (1,4) {mustBeReal}
 end
v_sol = (-5:0.1:5)';
u_sol = -r(1,3).*v_sol-r(1,4);

end