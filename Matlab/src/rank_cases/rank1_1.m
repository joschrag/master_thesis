function        [u_sol,v_sol] = rank1_1(r)
%RANK1_1 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (1,4) {mustBeReal}
end
u_sol = [];
v_sol = [];
u = sym("u","real");
v = sym("u","real");
v_0 = roots([-4*r(1,1),-4*r(1,3),-4*r(1,4)+r(1,2)^2]);
v_0 = v_0(imag(v_0)==0);
if numel(v_0) == 2
    if r(1,1) > 0
        v_0 = linspace(v_0(1),v_0(2),50);
    else
        v_0 = [linspace(v_0(1)-2,v_0(1),25);linspace(v_0(2),v_0(2)+2,25)];
    end
end
if ~isempty(v_0)
    u_root = roots([-1,-r(1,2),-r(1,1).*v.^2-r(1,3).*v-r(1,4)]);
    u_0 = subs(u_root,v,v_0);
    if size(v_0,1) == 2
    u_0 = [reshape(u_0(1:2:end,:)',1,[]);reshape(u_0(2:2:end,:)',1,[])];
    else
        u_0 = reshape(u_0',1,[]);
    end
    u_sol = [u_sol;u_0'];
    v_sol = [v_sol;repmat(v_0',size(u_0,2)/size(v_0,2),size(u_0,1)/size(v_0,1))];
end
u_0 = roots([-4*r(1,1),-4*r(1,1)*r(1,2),-4*r(1,1)*r(1,4)+r(1,3)^2]);
u_0 = u_0(imag(u_0)==0);
if numel(u_0) == 2
    if r(1,1) > 0
        u_0 = linspace(u_0(1),u_0(2),50);
    else
        u_0 = [linspace(u_0(1)-2,u_0(1),25);linspace(u_0(2),u_0(2)+2,25)];
    end
end
if ~isempty(u_0)
    v_root = roots([-r(1,1),-r(1,3),-u.^2-r(1,2).*u-r(1,4)]);
    v_0 = subs(v_root,u,u_0);
    if size(u_0,1)==2
        v_0 = [reshape(v_0(1:2:end,:)',1,[]);reshape(v_0(2:2:end,:)',1,[])]; 
    else
        v_0 = reshape(v_0',1,[]);
    end
    v_sol = [v_sol;v_0'];
    u_sol = [u_sol;repmat(u_0',size(v_0,2)/size(u_0,2),size(v_0,1)/size(u_0,1))];
end
end