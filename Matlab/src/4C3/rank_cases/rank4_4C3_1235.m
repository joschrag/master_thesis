function [v_sol,w_sol] = rank4_4C3_1235(r)
    %RANK5_12345 Summary of this function goes here
    %   Detailed explanation goes here
    arguments
        r (4,2) {mustBeReal}
    end
    v_root = roots([-1,-r(1,1),-r(1,2)]);
    v_root = v_root(abs(imag(v_root))<10^-10);
    w_root = repmat(-r(4,2),size(v_root));

    
    conds = abs(-r(3,1).*v_root-r(3,2)-w_root.^2) < 10^-10 & abs(-r(2,1).*v_root-r(2,2)-v_root.*w_root) < 10^-10;
    v_sol = v_root(conds);
    w_sol = w_root(conds);
    end
    
    