function        [u_sol,v_sol] = rank1_1(r,prime)
arguments
    r (1,4) {mustBeReal}
    prime (1,1) {mustBeInteger,mustBePositive}
end
u = 0:prime-1;
v = 0:prime-1;
u_sol = [];
v_sol = [];

for u0=u
    res = FF(-r(1,1).*v.^2-r(1,2)*u0-r(1,3).*v-r(1,4)-u0^2,prime).value;
    if any(res == 0)
        v0 = v(res == 0)';
        if ~isempty(v0)
            v_sol = [v_sol;v0];
            u_sol = [u_sol;repmat(u0,size(v0))];
        end
    end
end


end