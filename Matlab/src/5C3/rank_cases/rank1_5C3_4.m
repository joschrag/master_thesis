function [v_sol,w_sol] = rank1_5C3_4(r)
%RANK1_5C3_4 Solve the resulting subsystem of equations for the case R4.
arguments
    r (1,4) {mustBeReal}
end
%Obtain solutions from equations
v_sol = (-5:0.1:5)';
w_sol = repmat(-r(1,4),numel(v_sol),1);
end