function [p_root] = sturm(p,max_iter,epsilon)
%STURM Finds roots of a polynomial using Sturm's theorem.
%
%   p_root = sturm(p) returns the real roots of the polynomial defined by
%   its coefficients in p. The default maximum number of iterations is 5000,
%   and the default tolerance epsilon is 10^-5.
%
%   p_root = sturm(p,max_iter) allows the user to specify the maximum number
%   of iterations used in finding the roots. max_iter must be a positive integer.
%
%   p_root = sturm(p,max_iter,epsilon) also allows the user to specify the tolerance
%   epsilon for the minimum size of the intervals. epsilon must be a positive real number.
%
%   Example:
%       p = [1,-3,2]; % Polynomial is x^2 - 3x + 2
%       r = sturm(p);
%       disp('The roots are: ');
%       disp(r);
%
arguments
    p (1,:);                              % Coefficients of the polynomial.
    max_iter (1,1) {mustBePositive,mustBeInteger} = 5000; % Maximum number of iterations.
    epsilon (1,1) {mustBePositive,mustBeReal} = 10^-5;    % Tolerance for function values at roots.
end
n = numel(p);
p_vec = zeros(n,n);
p_vec(1,:) = p;
p_vec(2,2:end) = p(1:end-1).*(numel(p)-1:-1:1);
for j=3:n
    p1 = p_vec(j-2,j-2:end);
    p2 = p_vec(j-1,j-1:end);
    [~,r] = deconv(p1(find(p1~=0,1,'first'):end),p2(find(p2~=0,1,'first'):end));
    r=r.*(abs(r) > 10^(-5));
    p_vec(j,n+1-numel(r):end) = -r;
    if numel(r) == 1 || ~any(r)
        p_vec = p_vec(1:j,:);
        break
    end
end
p_vec=p_vec(any(p_vec,2),:);
m = size(p_vec,1);
lims = get_func_limits(p_vec);
intervals = [-inf,inf];
roots_per_interval = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1));
iter = 0;
%     less than 1 root per interval || ±inf as interval boundary
while (~all(roots_per_interval <= 1) || any(abs(intervals)==inf)) && iter <= max_iter
    [intervals,insert_indices] = insert_intervals(intervals,roots_per_interval);
    p_vals = polyval(double(p),intervals);
    if any(abs(p_vals) < epsilon)
        intervals(abs(p_vals) < epsilon) = intervals(abs(p_vals) < epsilon) + 0.1;
    end
    for i = insert_indices
        lims = [lims(1:i-1, :); zeros(1,m); lims(i:end, :)];
    end
    for i = 1:m
        lims(insert_indices,i) = sign(polyval(p_vec(i,:),intervals(insert_indices)));
    end
    roots_per_interval = -diff(sum(diff(lims, 1, 2) ~= 0, 2),1);
    padded = [0;roots_per_interval;0];
    running_sum = padded(1:end-1) + padded(2:end);
    to_be_trimmed = logical(running_sum);
    intervals = intervals(to_be_trimmed);
    lims = lims(to_be_trimmed,:);
    roots_per_interval = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1));
    iter = iter+1;
end
p_root = nan(sum(roots_per_interval~=0),1);
j = 1;
options = optimset('TolFun',1e-20,'TolX',1e-20);
for i = 1:numel(roots_per_interval)
    if roots_per_interval(i) ~= 0
        if sign(polyval(double(p),intervals(i))) ~= sign(polyval(double(p),intervals(i+1)))
            polynomial =@(x) polyval(double(p),x);
            p_root(j) = fzero(polynomial,intervals(i:i+1),options);
            j = j+1;
        end
    end
end
p_root = p_root(~isnan(p_root));
end