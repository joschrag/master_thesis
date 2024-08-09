function [roots] = get_gf_root(poly, order)
%GET_GF_ROOT Computes the roots of a polynomial over a Galois field.
%
%   ROOTS = GET_GF_ROOT(POLY, ORDER) computes the roots of the polynomial POLY over
%   the Galois field with prime order ORDER. If POLY is empty, then all possible
%   roots are returned as a vector from 0 to ORDER-1.
%
arguments
    poly (1,:); % Coefficients of the polynomial
    order (1,1) {mustBeInteger, mustBeNonnegative}; % Order of the Galois field
end
if isa(poly, "sym")
    poly = double(poly);
end
if ~isprime(order)
    error("FF:order", "Order must be prime.");
end
if isempty(poly)
    roots = 0:order-1;
else
    prim_pol = gfprimfd(1, "min", order); % Find a primitive polynomial for the Galois field
    q = fliplr(mod(poly, order)); % Reverse and reduce the coefficients of the polynomial modulo ORDER
    [~, roots] = gfroots(q, prim_pol, order); % Compute the roots of the polynomial over the Galois field
end
end