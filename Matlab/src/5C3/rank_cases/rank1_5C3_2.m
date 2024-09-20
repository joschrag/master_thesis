function        [u_sol,v_sol] = rank1_5C3_2(r)
arguments
    r (1,4) {mustBeReal}
end
u = sym("u","real");
u_0 = roots([-4*r(1,2),r(1,3)^2-4*r(1,4)]);
if ~isempty(u_0)
    u_0 = linspace(u_0,u_0-sign(r(1,2))*5,50);
v_0 = roots([-1,-r(1,3),-r(1,2)*u-r(1,4)]);

v_sol = subs(v_0,u,u_0);
v_sol = vpa(real(reshape(v_sol',[numel(v_sol),1])));
u_sol = vpa(real(repmat(u_0,1,2)))';
else
    u_sol = [];
    v_sol = [];
end


end