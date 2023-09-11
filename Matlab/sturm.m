function [p_root] = sturm(p)
%STURM Summary of this function goes here
%   Detailed explanation goes here
n = numel(p);
p_vec = zeros(n,n);
p_vec(1,:) = p;
p_vec(2,2:end) = polyder(p);
for j=3:n
    [~,r] = deconv(p_vec(j-2,j-2:end),p_vec(j-1,j-1:end));
    p_vec(j,n+1-numel(r):end) = -r;
    if numel(r) == 1
        p_vec = p_vec(1:j,:);
        break
    end
end
m = size(p_vec,1);
lims = zeros(2,m);
for j = 1:size(p_vec,1)
    lims(:,j) = get_func_limits(p_vec(j,:));
end
intervals = [-inf,inf];

sc_intervals = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1));
root_indices = find(sc_intervals~=0);

while ~all(sc_intervals <= 1) || any(intervals == inf) || any(intervals == -inf)
    %TODO check if on root -> handle root accordingly
    [intervals,insert_indices] = insert_intervals(intervals,root_indices,sc_intervals);
    for i = insert_indices
        lims = [lims(1:i-1, :); zeros(1,m); lims(i:end, :)];
    end
    for i = 1:size(p_vec,1)
        lims(insert_indices,i) = sign(polyval(p_vec(i,:),intervals(insert_indices)));
    end
    sc_intervals = -diff(sum(diff(lims, 1, 2) ~= 0, 2),1);

    %interval trimming
    padded = [0;sc_intervals;0];
    running_sum = padded(1:end-1) + padded(2:end);
    to_be_trimmed = logical(running_sum);
    intervals = intervals(to_be_trimmed);
    lims = lims(to_be_trimmed,:);
    sc_intervals = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1));
    root_indices = find(sc_intervals~=0);
end
p_root = zeros(sum(sc_intervals~=0),1);
j = 1;
for i = 1:numel(sc_intervals)
    if sc_intervals(i) ~= 0
        polynomial =@(x) sum(p .* x.^(numel(p)-1:-1:0));
        p_root(j) = fzero(polynomial,intervals(i:i+1));
        j = j+1;
    end
end
end

