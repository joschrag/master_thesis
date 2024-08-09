function [limits] = get_func_limits(p)
%GET_FUNC_LIMITS Computes the limits of a polynomial function.
%
%   LIMITS = GET_FUNC_LIMITS(P) computes the limits of the polynomial function defined by
%   the matrix P, where each row represents a polynomial term and its corresponding
%   coefficients. The output LIMITS is a 2-by-n matrix containing the lower and upper
%   limits for each polynomial term.
%
arguments
    p (:,:) {mustBeReal}; % Coefficients of the polynomial terms
end
[n,m] = size(p);
c = zeros(n,1);
s = zeros(n,1);
for i=1:n
c(i) = find(p(i,:)~=0,1,"first"); % Find the index of the first non-zero coefficient for each polynomial term
s(i) = sign(p(i,c(i))); % Determine the sign of the corresponding coefficient
end
is_even_degree = logical(mod(m-c,2)); % Check if the degree of each polynomial term is even or odd
limits = [(-1).^(is_even_degree).*s, s]'; % Compute the limits based on the degree and sign of each coefficient
end