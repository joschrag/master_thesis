function        [t_sol,u_sol] = rank1_1(r)
%RANK1_1 Summary of this function goes here
%   Detailed explanation goes here
arguments
    r (1,4) {mustBeReal}
end
t_sol = [];
u_sol = [];
t = sym("t","real");
u = sym("u","real");
u_0 = roots([-4*r(1,1),-4*r(1,3),-4*r(1,4)+r(1,2)^2]);
u_0 = sort(u_0(imag(u_0)==0));
if numel(u_0) == 2
    if r(1,1) > 0
        u_0 = linspace(u_0(1),u_0(2),50);
    else
        u_0 = [linspace(u_0(1)-2,u_0(1),25);linspace(u_0(2),u_0(2)+2,25)];
    end
end
if ~isempty(u_0)
    t_root = roots([-1,-r(1,2),-r(1,1).*u.^2-r(1,3).*u-r(1,4)]);
    t_0 = subs(t_root,u,u_0);
    if size(u_0,1) == 2
        t_0 = [reshape(t_0(1:2:end,:)',1,[]);reshape(t_0(2:2:end,:)',1,[])];
    else
        t_0 = reshape(t_0',1,[]);
    end
    t_sol = [t_sol;t_0'];
    u_sol = [u_sol;repmat(u_0',size(t_0,2)/size(u_0,2),size(t_0,1)/size(u_0,1))];
end
t_0 = roots([-4*r(1,1),-4*r(1,1)*r(1,2),-4*r(1,1)*r(1,4)+r(1,3)^2]);
t_0 = sort(t_0(imag(t_0)==0));
if numel(t_0) == 2
    if r(1,1) > 0
        t_0 = linspace(t_0(1),t_0(2),50);
    else
        t_0 = [linspace(t_0(1)-2,t_0(1),25);linspace(t_0(2),t_0(2)+2,25)];
    end
end
if ~isempty(t_0)
    u_root = roots([-r(1,1),-r(1,3),-t.^2-r(1,2).*t-r(1,4)]);
    u_0 = subs(u_root,t,t_0);
    if size(t_0,1)==2
        u_0 = [reshape(u_0(1:2:end,:)',1,[]);reshape(u_0(2:2:end,:)',1,[])]; 
    else
        u_0 = reshape(u_0',1,[]);
    end
    u_sol = [u_sol;u_0'];
    t_sol = [t_sol;repmat(t_0',size(u_0,2)/size(t_0,2),size(u_0,1)/size(t_0,1))];
end
if ~isempty(u_sol) && ~isempty(t_sol)
    u_sol = vpa(real(u_sol));
    t_sol = vpa(real(t_sol));
end