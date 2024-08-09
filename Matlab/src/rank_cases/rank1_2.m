function        [u_sol,v_sol] = rank1_2(r)
arguments
    r (1,4) {mustBeReal}
end
u = sym("u","real");
u_0 = roots([-4*r(1,2),r(1,3)^2-4*r(1,4)]);
u_0 = linspace(u_0,u_0-sign(r(1,2))*5,50)';
v_0 = roots([-1,-r(1,3),-r(1,2)*u-r(1,4)]);

v_sol = subs(v_0,u,u_0);
v_sol = reshape(v_sol',[numel(v_sol),1]);
u_sol = [u_0;u_0];

end