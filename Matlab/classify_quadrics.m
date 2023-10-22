function [result] = classify_quadrics(eigs,d0,d1,r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    eigs (3,1)
    d0 (1,1) {mustBeReal}
    d1 (1,1) {mustBeReal}
    r (1,1) {mustBeReal}
end
pos_eig = eigs{1};
neg_eig = eigs{2};
zero_eig = eigs{3};

switch zero_eig
    case 0
        if r<0
            switch pos_eig
                case 3
                    fprintf("ellipsoid\n")
                    result = 3;
                case 2
                    fprintf("hyperboloid1\n")
                    result = 4;
                case 1
                    fprintf("hyperboloid2\n")
                    result = 5;
                otherwise
                    error("Equation contains no real solutions.")
            end
        elseif r == 0
            if pos_eig == 3 || neg_eig == 3
                %point solution
                fprintf("point solution\n")
                result = 1;
            else
                %ell cone
                fprintf("ell cone\n")
                result = 2;
            end
        else
            error("r should be <= 0")
        end
    case 1
        if d0 == 0
            if r < 0
                switch pos_eig
                    case {2,0}
                        fprintf("ell cylinder\n")
                        result = 8;
                    case 1
                        fprintf("hyp cylinder\n")
                        result = 10;
                    otherwise
                        error("Equation contains no real solutions.")
                end
            elseif r == 0
                switch pos_eig
                    case {2,0}
                        fprintf("line\n")
                        result = 9;
                    case 1
                        fprintf("two planes crossing\n")
                        result = 14;
                    otherwise
                        error("Equation contains no real solutions.")
                end
            else
                error("r should be <= 0")
            end
        else
            switch pos_eig
                case 2
                    fprintf("ell paraboloid\n")
                    result = 6;

                case 1
                    fprintf("hyp paraboloid\n")
                    result = 7;
                otherwise
                    error("Equation contains no real solutions.")
            end
        end
    case 2
        if r == 0
            if d0 ~= 0 || d1 ~= 0
                fprintf("parabol cylinder\n")
                result = 11;
            else
                fprintf("one plane\n")
                result = 13;
            end
        elseif r < 0
            if pos_eig == 1
                fprintf("two planes parallel\n")
                result = 12;
            else
                error("Equation contains no real solutions.")
            end
        else
            error("r should be <= 0")
        end
    otherwise
        error("No quadratic factors included.")
end