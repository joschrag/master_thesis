function [offset, const] = square_complete(a, b, c)
%SQUARE_COMPLETE Completes a quadratic equation into square form.
%
%   [OFFSET, CONST] = SQUARE_COMPLETE(A, B, C) completes the quadratic
%   equation Ax^2 + Bx + C into square form and returns its offset and constant.
%
%   Example:
%       a = 1; b = -3; c = 2; % Quadratic equation is x^2 - 3x + 2
%       [offset, const] = square_complete(a, b, c);
%       fprintf('The square form of the quadratic equation is (x - %.1f)^2 + %.1f\n', offset, const);
%
arguments
    a (1,1) {mustBeReal, mustBeNonzero} % Coefficient of x^2.
    b (1,1) {mustBeReal}                % Coefficient of x.
    c (1,1) {mustBeReal}                % Constant term.
end
offset = (b/(2*a));
const = c - (b^2/(4*a));
end