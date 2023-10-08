function [result] = classify_quadrics(eigs,r,g,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pos_eig = eigs{1};
neg_eig = eigs{2};
zero_eig = eigs{3};
if zero_eig == 0
    if r < 0
        if pos_eig == 3 || neg_eig == 3
            %point solution
            result = "point solution";
        else
            %ell cone
            result = "ell cone";
        end
    elseif r == 0
        if pos_eig == 3
            %Ellipsoid
            result = "ellipsoid";
        elseif pos_eig == 2
            %Hyperboloid 1
            result = "hyperboloid 1";
        elseif pos_eig == 1
            %Hyperboloid 2
            result = "hyperboloid 2";
        elseif pos_eig == 0
            error("Equation contains no real solutions.")
        else
            error("r>0 should not happen!")
        end
    end
elseif zero_eig == 1
    if g ~= 0
        if pos_eig == 2 || neg_eig == 2
            % ell paraboloid
            result = "ell paraboloid";
        else
            % hyperb paraboloid
            result = "hyp paraboloid";
        end
    elseif g == 0
        if pos_eig == 2
            if r < 0
                % ell cylinder
                result = "ell cylinder";
            elseif r == 0
                % line
                result = "line";
            elseif r >0
                error("Equation contains no real solutions.")
            end
        elseif neg_eig == 2
            if r < 0
                error("Equation contains no real solutions.")
            elseif r == 0
                % line
                result = "line";
            elseif r >0
                %ell cylinder
                result = "ell cylinder";
            end
        elseif neg_eig == 1 && pos_eig == 1
            if r < 0
                % hyberb cylinder
                result = "hyp cylinder";
            else    %if r == 0 || r <0
                % two planes
                result = "two planes";
            end
        end
    end
elseif zero_eig == 2
    if (g ~=0) ||  (h~= 0)
        % parabol cylinder
        result = "par cylinder";
    elseif pos_eig == 1
        if r < 0
            %two planes
            result = "two planes";
        elseif r == 0
            % one plane
            result = "one plane";
        else    %r>0
            error("Equation contains no real solutions.")
        end
    else %if neg_eig == 1
        if r > 0
            %two planes
            result = "two planes";
        elseif r == 0
            % one plane
            result = "one plane";
        else    %r<0
            error("Equation contains no real solutions.")
        end
    end
else
    error("No quadratic factors included.")
end
end

