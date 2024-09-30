function [v_sol,w_sol] = rank1_5C3_3(r)
%RANK1_5C3_3 Solve the resulting subsystem of equations for the case R3.
arguments
    r (1,4) {mustBeReal}
end
%Obtain solutions from equations
w_sol = (-5:0.1:5)';
v_sol = -r(1,3).*w_sol-r(1,4);
end