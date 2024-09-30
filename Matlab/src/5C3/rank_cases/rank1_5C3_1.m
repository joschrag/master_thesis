function [v_sol,w_sol] = rank1_5C3_1(r)
%RANK1_5C3_1 Solve the resulting subsystem of equations for the case R1.
arguments
    r (1,4) {mustBeReal}
end
% Obtain solution from equations through abc formula and substitution
v_sol = [];
w_sol = [];
v = sym("v","real");
w = sym("w","real");
w_0 = roots([-4*r(1,1),-4*r(1,3),-4*r(1,4)+r(1,2)^2]);
w_0 = sort(w_0(imag(w_0)==0));
if numel(w_0) == 2
    if r(1,1) > 0
        w_0 = linspace(w_0(1),w_0(2),50);
    else
        w_0 = [linspace(w_0(1)-2,w_0(1),25);linspace(w_0(2),w_0(2)+2,25)];
    end
end
if ~isempty(w_0)
    v_root = roots([-1,-r(1,2),-r(1,1).*w.^2-r(1,3).*w-r(1,4)]);
    v_0 = subs(v_root,w,w_0);
    if size(w_0,1) == 2
        v_0 = [reshape(v_0(1:2:end,:)',1,[]);reshape(v_0(2:2:end,:)',1,[])];
    else
        v_0 = reshape(v_0',1,[]);
    end
    v_sol = [v_sol;v_0'];
    w_sol = [w_sol;repmat(w_0',size(v_0,2)/size(w_0,2),size(v_0,1)/size(w_0,1))];
end
v_0 = roots([-4*r(1,1),-4*r(1,1)*r(1,2),-4*r(1,1)*r(1,4)+r(1,3)^2]);
v_0 = sort(v_0(imag(v_0)==0));
if numel(v_0) == 2
    if r(1,1) > 0
        v_0 = linspace(v_0(1),v_0(2),50);
    else
        v_0 = [linspace(v_0(1)-2,v_0(1),25);linspace(v_0(2),v_0(2)+2,25)];
    end
end
if ~isempty(v_0)
    u_root = roots([-r(1,1),-r(1,3),-v.^2-r(1,2).*v-r(1,4)]);
    w_0 = subs(u_root,v,v_0);
    if size(v_0,1)==2
        w_0 = [reshape(w_0(1:2:end,:)',1,[]);reshape(w_0(2:2:end,:)',1,[])];
    else
        w_0 = reshape(w_0',1,[]);
    end
    w_sol = [w_sol;w_0'];
    v_sol = [v_sol;repmat(v_0',size(w_0,2)/size(v_0,2),size(w_0,1)/size(v_0,1))];
end
if ~isempty(w_sol) && ~isempty(v_sol)
    w_sol = vpa(real(w_sol));
    v_sol = vpa(real(v_sol));
end