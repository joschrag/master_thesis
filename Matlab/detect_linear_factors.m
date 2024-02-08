function [g, h] = detect_linear_factors(A, B, P)
g = 0;
h = 0;

if numel(single_vars) >= 1
    tmp = coeffs(dot(B,P*sub_vars),single_vars(1),"All");
    if numel(tmp) == 2
        g = double(tmp(1));
    end
end
if numel(single_vars) >= 2
    tmp = coeffs(dot(B,P*sub_vars),single_vars(2),"All");
    if numel(tmp) == 2
        h = double(tmp(1));
    end
end
end