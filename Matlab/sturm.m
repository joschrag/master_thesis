function [outputArg1] = sturm(p)
%STURM Summary of this function goes here
%   Detailed explanation goes here
n = numel(p) -1 ;
p_vec = zeros([n+1,n+1]);
p_vec(1,:) = p;
p_vec(2,2:end) = polyder(p);
for j=3:n+1
    [~,r] = deconv(p_vec(j-2,j-2:end),p_vec(j-1,j-1:end));
    p_vec(j,n+2-numel(r):end) = -r;
    if numel(r) == 1
        p_vec = p_vec(1:j,:);
        break
    end
end
p_vec
x = sym("x","real");
m = size(p_vec,1);
lims = zeros(2,m);
for j = 1:size(p_vec,1)
    q = poly2sym(p_vec(j,:), x);
    lims(:,j) = sign([limit(q, x, -inf),limit(q, x, inf)]);
end
intervals = [-inf,inf];
lims
sc_intervals = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1))
root_indices = find(sc_intervals~=0)

while ~all(sc_intervals ==1) || any(intervals == inf) || any(intervals == -inf)
    %TODO check if on root -> handle root accordingly
    [intervals,insert_indices] = insert_intervals(intervals,root_indices,sc_intervals)
    for i = insert_indices
        lims = [lims(1:i-1, :); zeros(1,m); lims(i:end, :)]
    end
    for i = 1:size(p_vec,1)
        lims(insert_indices,i) = sign(polyval(p_vec(i,:),intervals(insert_indices)))
    end
    lims
    sc_intervals = -diff(sum(diff(lims, 1, 2) ~= 0, 2),1)
    root_indices = find(sc_intervals~=0)

    %interval trimming
    padded = [0;sc_intervals;0]
    sum_of_following_elements = padded(1:end-1) + padded(2:end)
    to_be_trimmed = logical([0;sc_intervals] + [sc_intervals;0])
    intervals = intervals(to_be_trimmed)
    lims = lims(to_be_trimmed,:);
    sc_intervals = abs(diff(sum(diff(lims, 1, 2) ~= 0, 2),1))
    root_indices = find(sc_intervals~=0)

end
disp("DONE")
intervals
end

