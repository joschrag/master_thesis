function [p_root] = sturm(p,tolerance,max_iter,epsilon)
arguments
    p;
    tolerance (1,1) {mustBePositive,mustBeReal} = 10^-1;
    max_iter (1,1) {mustBePositive,mustBeInteger} = 5000;
    epsilon (1,1) {mustBePositive,mustBeReal} = 10^-5;
end
n = numel(p);
p_vec = zeros(n,n);
p_vec(1,:) = p;
p_vec(2,2:end) = polyder(p);
for j=3:n
    p1 = p_vec(j-2,j-2:end);
    p2 = p_vec(j-1,j-1:end);
    [~,r] = deconv(p1(find(p1~=0,1,'first'):end),p2(find(p2~=0,1,'first'):end));
    p_vec(j,n+1-numel(r):end) = -r;
    if numel(r) == 1
        p_vec = p_vec(1:j,:);
        break
    end
end
p_vec=p_vec(any(p_vec,2),:);
m = size(p_vec,1);
lims = zeros(2,m);
for j = 1:size(p_vec,1)
    lims(:,j) = get_func_limits(p_vec(j,:));
end
intervals = [-inf,inf];
roots_per_interval = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1));
iter = 0;
%     less than 1 root per interval || ± as interval boundary
while (~all(roots_per_interval <= 1) || any(ismember(intervals,[-inf,inf]))) && iter <= max_iter
    %         || ~all(abs(diff(intervals)) < tolerance)...
    [intervals,insert_indices] = insert_intervals2(intervals,roots_per_interval);
    p_vals = polyval(p,intervals);
    if any(abs(p_vals) < epsilon)
        intervals(abs(p_vals) < epsilon) = intervals(abs(p_vals) < epsilon) + 0.1;
    end
    for i = insert_indices
        lims = [lims(1:i-1, :); zeros(1,m); lims(i:end, :)];
    end
    for i = 1:size(p_vec,1)
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
p_root = zeros(sum(roots_per_interval~=0),1);
j = 1;
for i = 1:numel(roots_per_interval)
    if roots_per_interval(i) ~= 0
        polynomial =@(x) sum(p .* x.^(numel(p)-1:-1:0));
        p_root(j) = fzero(polynomial,intervals(i:i+1));
        j = j+1;
    end
end
end