function     [new_intervals,insert_indices] = insert_intervals(intervals,roots_per_interval,step,tolerance)
%INSERT_INTERVALS Summary of this function goes here
%   Detailed explanation goes here
arguments
    intervals (1,:) {mustBeReal};
    roots_per_interval (1,:) {mustBeNonnegative,mustBeInteger};
    step (1,1) {mustBeNonnegative,mustBeInteger} = 50;
    tolerance (1,1) {mustBeReal,mustBePositive} = 10^-1;
end
if numel(intervals) == 2 && all(intervals == [-inf,inf])
    new_intervals = [-inf,-100,100,inf];
    insert_indices = [2,3];
else
    valid_intervals = roots_per_interval & diff(intervals) >= tolerance;
    valid_indices = find(valid_intervals);

    roots_finite_interval = (isfinite(diff(intervals)).*roots_per_interval);
    valid_roots_finite_interval = roots_finite_interval(valid_indices);
    num_finite_intervals_add = sum(valid_roots_finite_interval);

    roots_infinite_interval = ~isfinite(diff(intervals));
    valid_roots_infinite_interval = roots_infinite_interval(valid_indices);
    num_infinite_intervals_add = sum(valid_roots_infinite_interval);

    roots_interval = roots_infinite_interval + roots_finite_interval;
    roots_interval_valid = roots_interval(valid_indices);
    num_intervals_add = num_finite_intervals_add + num_infinite_intervals_add;
    new_intervals = NaN(1,numel(intervals)+num_intervals_add);
    insert_indices = zeros(1, sum(roots_interval_valid));
    csum = cumsum(roots_interval_valid(1:end-1));
    new_indices = valid_indices + [1,csum+1];
    index = [1,csum+1];
    
    for i = 1:numel(roots_interval_valid)
        offset = roots_interval_valid(i)-1;
        insert_indices(index(i):index(i)+offset) = new_indices(i):new_indices(i)+offset;
    end

    for i = 1:numel(valid_indices)
        index = valid_indices(i);
        lower_bound = intervals(index);
        upper_bound = intervals(index+1);
        if ~isfinite(lower_bound)
            new_intervals(2) = upper_bound - step;
        elseif ~isfinite(upper_bound)
            new_intervals(end-1) = lower_bound + step;
        else
            parts = roots_per_interval(index)+1;
            step_length = (upper_bound - lower_bound)/parts;
            if parts == 2
                new_intervals(new_indices(i):new_indices(i)+roots_per_interval(index)-1) = lower_bound+step_length;
            else
                new_intervals(new_indices(i):new_indices(i)+roots_per_interval(index)-1) = lower_bound+step_length:step_length:upper_bound-step_length;
            end
        end
    end
    new_intervals(isnan(new_intervals)) = intervals;
end
end