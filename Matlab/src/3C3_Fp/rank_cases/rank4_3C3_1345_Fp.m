function [v_sol,w_sol] = rank4_3C3_1345_Fp(r,prime)
    %RANK5_12345 Summary of this function goes here
    %   Detailed explanation goes here
    arguments
        r (4,2) {mustBeReal}
        prime (1,1) {mustBeInteger,mustBePositive}
    end
    
    v_root = FF(-r(3,2),prime).value;
    w_root = FF(-r(4,2),prime).value;
    
    conds = FF(-r(1,1).*v_root.^2-r(1,2)-w_root.^3,prime).value == 0 & FF(-r(2,2)-w_root.^2,prime).value == 0;
    v_sol = v_root(conds);
    w_sol = w_root(conds);
    end
    
    