function [v_sol,w_sol] = rank3_3C3_123(r,complex)
%RANK_3C3_123 Solve the resulting subsystem of equations for the case R123.
arguments
    r (3,3) 
    complex (1,1) {mustBeNumericOrLogical} = false;
end
v = sym("v","real");
w = sym("w","real");
v_sol = [];
w_sol = [];
w_root = [];
v_root = [];
eq = -r(3,1)*v-r(3,2)*w-r(3,3)-w^2;
% Obtain root candidates from equations
if r(2,2) ~= 0
    w_0 = -(v^2+r(2,1)*v+r(2,3))/r(2,2);
    v_root = roots(coeffs(subs(eq,w,w_0),v,"All"));
    if ~complex
    v_root =  v_root(abs(imag(v_root))<10^-10);
    end
    w_root = subs(w_0,v,v_root);
else
    v_0 = roots([-1,-r(2,1),-r(2,3)]);
    if ~complex
    v_0 =  v_0(abs(imag(v_0))<10^-10);
    end
    if ~isempty(v_0)
        sub_w = subs(eq,v,v_0);
        for i=1:size(sub_w,1)
            coefs = coeffs(sub_w(i),w,"All");
            w_0 = roots(coefs);
            if ~complex
            w_0 =  w_0(abs(imag(w_0))<10^-10);
            end
            w_root = [w_root;w_0];
            v_root = [v_root;repmat(v_0(i),size(w_0))];
        end
    end
end
if ~isempty(v_root) && ~isempty(w_root)
    % Remove roots not satysfying the conditions
    conds = abs(-r(1,1).*v_root-r(1,2).*w_root-r(1,3)-w_root.^3) <10^-10;
    v_sol = v_root(conds);
    w_sol = w_root(conds);
end
end