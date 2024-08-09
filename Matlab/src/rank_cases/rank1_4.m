function        [u_sol,v_sol] = rank1_4(r)
arguments
    r (1,4) {mustBeReal}
 end
u_sol = (-5:0.1:5)';
v_sol = repmat(-r(1,4),numel(u_sol),1);

end