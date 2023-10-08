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
            fprintf("point solution\n")
            result = 1;
        else
            %ell cone
            fprintf("ell cone\n")
            result = 2;
        end
    elseif r == 0
        if pos_eig == 3
            %Ellipsoid
            fprintf("ellipsoid\n")
            result = 3;
        elseif pos_eig == 2
            %Hyperboloid 1
            fprintf("hyperboloid1\n")
            result = 4;
        elseif pos_eig == 1
            %Hyperboloid 2
            fprintf("hyperboloid2\n")
            result = 5;
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
            fprintf("ell paraboloid\n")
            result = 6;
        else
            % hyperb paraboloid
            fprintf("hyp paraboloid\n")
            result = 7;
        end
    elseif g == 0
        if pos_eig == 2
            if r < 0
                % ell cylinder
                fprintf("ell cylinder\n")
                result = 8;
            elseif r == 0
                % line
                fprintf("line\n")
                result = 9;
            elseif r >0
                error("Equation contains no real solutions.")
            end
        elseif neg_eig == 2
            if r < 0
                error("Equation contains no real solutions.")
            elseif r == 0
                % line
                fprintf("line\n")
                result = 9;
            elseif r >0
                %ell cylinder
                fprintf("ell cylinder\n")
                result = 8;
            end
        elseif neg_eig == 1 && pos_eig == 1
            if r < 0
                % hyberb cylinder
                fprintf("hyp cylinder\n")
                result = 10;
            else    %if r == 0 || r <0
                % two planes crossing
                fprintf("two planes crossing\n")
                result = 14;
            end
        end
    end
elseif zero_eig == 2
    if (g ~=0) ||  (h~= 0)
        % parabol cylinder
        fprintf("parabol cylinder\n")
        result = 11;
    elseif pos_eig == 1
        if r < 0
            %two planes parallel
            fprintf("two planes parallel\n")
            result = 12;
        elseif r == 0
            % one plane
            fprintf("one plane\n")
            result = 13;
        else    %r>0
            error("Equation contains no real solutions.")
        end
    else %if neg_eig == 1
        if r > 0
            %two planes
            fprintf("two planes parallel\n")
            result = 12;
        elseif r == 0
            % one plane
            fprintf("one plane\n")
            result = 13;
        else    %r<0
            error("Equation contains no real solutions.")
        end
    end
else
    error("No quadratic factors included.")
end
end

