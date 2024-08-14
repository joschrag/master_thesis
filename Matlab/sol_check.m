function [] = sol_check(c,p)
for x=0:p-1
    for y=0:p-1
        for z=0:p-1
            m=mod(c*[x^2;x*y;x*z;y^2;y*z;z^2;x;y;z;1],p);
            if  m == zeros(3,1)
                disp([x,y,z])
            end
        end
    end
end
end

