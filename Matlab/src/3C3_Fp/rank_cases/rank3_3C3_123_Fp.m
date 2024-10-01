function [v_sol,w_sol] = rank3_3C3_123_Fp(r,prime)
%RANK_3C3_123_FP Solve the resulting subsystem of equations for the case R123.
arguments
    r (3,3) {mustBeReal}
    prime (1,1) {mustBePrime}
end
v = sym("v","integer");
w = sym("w","integer");
v_sol = [];
w_sol = [];
w_root = [];
v_root = [];
eq = FF(-r(3,1)*v-r(3,2)*w-r(3,3)-w^2,prime).value;
% Obtain root candidates from equations
if r(2,2) ~= 0
    w_0 = FF(-r(2,2),prime)^-1*FF(v^2+r(2,1)*v+r(2,3),prime);
    v_root = get_gf_root(coeffs(FF(subs(eq,w,w_0.value),prime).value,v,"All"),prime);
    w_root = FF(subs(w_0,v,v_root),prime).value;
else
    v_0 = get_gf_root([-1,-r(2,1),-r(2,3)],prime);
    if ~isempty(v_0)
        sub_w = FF(subs(eq,v,v_0),prime).value;
        for i=1:size(sub_w,1)
            coefs = coeffs(sub_w(i),w,"All");
            w_0 = get_gf_root(coefs,prime);
            w_root = [w_root;w_0];
            v_root = [v_root;repmat(v_0(i),size(w_0))];
        end
    end
end
if ~isempty(v_root) && ~isempty(w_root)
    % Remove roots not satysfying the conditions
    conds = FF(-r(1,1).*v_root-r(1,2).*w_root-r(1,3)-w_root.^3,prime).value ==0;
    v_sol = v_root(conds);
    w_sol = w_root(conds);
end
end