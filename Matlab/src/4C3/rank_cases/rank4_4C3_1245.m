function [v_sol,w_sol] = rank4_4C3_1245(r)
    %RANK5_12345 Summary of this function goes here
    %   Detailed explanation goes here
    arguments
        r (4,2) {mustBeReal}
    end
    
    v_root = -r(3,2);
    w_root = -r(4,2);
    
    conds = abs(-r(1,1).*w_root.^2-r(1,2)-v_root.^2) < 10^-10 & abs(-r(2,1).*w_root.^2-r(2,2)-w_root.*v_root) < 10^-10;
    v_sol = v_root(conds);
    w_sol = w_root(conds);
    end
    
    