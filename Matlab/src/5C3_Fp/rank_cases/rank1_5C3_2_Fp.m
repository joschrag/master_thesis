function [v_sol,w_sol] = rank1_5C3_2_Fp(r,prime)
%RANK1_5C3_2_FP Solve the resulting subsystem of equations for the case R2.
arguments
    r (1,4) {mustBeReal}
    prime (1,1) {mustBePrime}
end
v = 0:prime-1;
w = 0:prime-1;
v_sol = [];
w_sol = [];
% Solve equation by brute force.
for v0=v
    res = FF(-w.^2-r(1,2)*v0-r(1,3).*w-r(1,4),prime).value;
    if any(res == 0)
        w0 = w(res == 0);
        if ~isempty(w0)
            w_sol = [w_sol;w0'];
            v_sol = [v_sol;repmat(v0,numel(w0),1)];
        end
    end
end
end