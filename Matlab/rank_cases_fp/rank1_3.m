function        [u_sol,v_sol] = rank1_3(r,prime)
arguments
    r (1,4) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
v_sol = (0:prime-1)';
u_sol = FF(-r(1,3).*v_sol-r(1,4),prime).value;

end