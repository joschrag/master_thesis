function [v_sol,w_sol] = rank4_3C3_2345(r)
    %RANK5_12345 Summary of this function goes here
    %   Detailed explanation goes here 
    arguments
        r (4,2) {mustBeReal}
    end
    
    v_root = -r(3,2);
    w_root = -r(4,2);
    
    conds = abs(-r(1,2)-v_root.^2)< 10^-10 & abs(-r(2,2)-w_root.^2)< 10^-10;
    v_sol = v_root(conds);
    w_sol = w_root(conds);
    end
    
    