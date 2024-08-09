function [new_intervals, insert_indices] = insert_intervals(intervals, roots_per_interval, step, tolerance)
%INSERT_INTERVALS Insert new intervals into an existing array of intervals.
%   This function takes in an array of intervals and the number of additional
%   intervals to insert within each interval. The function returns a new array
%   of intervals with the inserted intervals, as well as an array of indices
%   where the new intervals were inserted.
%
%   Parameters:
%       intervals (1,:) - Array of real-valued numbers representing the existing
%                         intervals. Must have at least two elements and be sorted
%                         in ascending order.
%       roots_per_interval (1,:) - Array of non-negative integers representing the
%                                  number of additional intervals to insert within
%                                  each existing interval. Must have the same length
%                                  as the `intervals` array.
%       step (1,1) - Positive integer representing the size of the new intervals
%                    inserted at either end of an infinite interval. Default value is 50.
%       tolerance (1,1) - Real-valued number greater than zero representing the minimum
%                         width required for a finite interval to have additional
%                         intervals inserted within it. Default value is 1e-6.
%
%   Returns:
%       new_intervals (1,:) - Array of real-valued numbers representing the new array
%                             of intervals with the inserted intervals. Will have the
%                             same length as the `intervals` array plus any additional
%                             intervals that were inserted.
%       insert_indices (1,:) - Array of integers representing the indices where the new
%                              intervals were inserted within the `new_intervals` array.
%                              Will have a length equal to the total number of new
%                              intervals that were inserted.
%
%   Example usage:
%       >> intervals = [-inf, -10, 10, inf];
%       >> roots_per_interval = [2, 3, 0];
%       >> step = 5;
%       >> tolerance = 1e-6;
%       >> [new_intervals, insert_indices] = insert_intervals(intervals, roots_per_interval, step, tolerance);
%       new_intervals = [-Inf, -15, -10, -5, 0, 5, 10, Inf];
%       insert_indices = [2, 4, 5, 6]
arguments
    intervals (1,:) {mustBeReal};
    roots_per_interval (1,:) {mustBeNonnegative,mustBeInteger};
    step (1,1) {mustBeNonnegative,mustBeInteger} = 50;
    tolerance (1,1) {mustBeReal,mustBePositive} = 10^-6;
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