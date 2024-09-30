function [v_sol,w_sol] = rank1_5C3_2(r)
%RANK1_5C3_2 Solve the resulting subsystem of equations for the case R2.
arguments
    r (1,4) {mustBeReal}
end
% Obtain solution from equations through abc formula and substitution
v = sym("v","real");
v_0 = roots([-4*r(1,2),r(1,3)^2-4*r(1,4)]);
if ~isempty(v_0)
    v_0 = linspace(v_0,v_0-sign(r(1,2))*5,50);
    w_0 = roots([-1,-r(1,3),-r(1,2)*v-r(1,4)]);
    
    w_sol = subs(w_0,v,v_0);
    w_sol = vpa(real(reshape(w_sol',[numel(w_sol),1])));
    v_sol = vpa(real(repmat(v_0,1,2)))';
else
    v_sol = [];
    w_sol = [];
end


end