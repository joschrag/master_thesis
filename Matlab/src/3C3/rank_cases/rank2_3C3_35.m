function [v_sol,w_sol] = rank2_3C3_35(r)
    %RANK_3C3_23 Solve the resulting subsystem of equations for the case R23.
    arguments
        r (2,4) {mustBeReal}
    end
    v_sol = [];
    w_sol = [];
    w_root = -r(2,4);
    if r(1,3) ~= 0
        v_sol = -(r(1,4)+r(2,4)^2)/r(1,3);
        w_sol = w_root;
    elseif abs(r(2,4)^2+r(1,4)) < 10^-10
        v_sol = (-5:0.1:5)';
        w_sol = repmat(w_root,size(v_sol));
    end
    end