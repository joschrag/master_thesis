function        [u_sol,v_sol] = rank1_4(r,prime)
arguments
    r (1,4) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
u_sol = (0:prime-1)';
v_sol = repmat(FF(-r(1,4),prime).value,prime,1);

end