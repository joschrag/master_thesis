function        [u_sol,v_sol] = rank1_2(r)
arguments
    r (1,4) {mustBeReal}
end
u = sym("u","real");
u_0 = roots([-4*r(1,2),r(1,3)^2-4*r(1,4)]);
u_0 = sort(u_0(imag(u_0)==0));
if numel(u_0) == 2
    if r(1,1) > 0
        u_0 = linspace(u_0(1),u_0(2),50);
    else
        u_0 = [linspace(u_0(1)-2,u_0(1),25);linspace(u_0(2),u_0(2)+2,25)];
    end
end
v_0 = roots([-1,-r(1,3),-r(1,2)*u-r(1,4)]);

v_sol = subs(v_0,u,u_0);
v_sol = vpa(real(reshape(v_sol',[numel(v_sol),1])));
u_sol = vpa(real([u_0;u_0]));

end