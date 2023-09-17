function [new_intervals,new_indices] = insert_intervals(intervals,root_indices,num_roots)
%INSERT_INTERVALS Summary of this function goes here
%   Detailed explanation goes here
arguments
    intervals;
    root_indices;
    num_roots;
end
step = 10;
if numel(intervals) == 2 && all(intervals == [-inf,inf])
    new_intervals = [-inf,-100,100,inf];
    new_indices = [2,3];
else
    max_size = sum(num_roots);
    new_intervals = intervals;
    new_indices = zeros(1,max_size);

    for i = numel(root_indices):-1:1
        if i ~=1
            ind_before = sum(num_roots(1:i-1));
        else
            ind_before = 0;
        end
   
        index = root_indices(i);
        parts = num_roots(i)+1;
        if intervals(index) == -inf
            new_intervals = [-inf,new_intervals(2)-step,new_intervals(2:end)];
            insertions = 1;
        elseif intervals(index+1) == inf
            new_intervals = [new_intervals(1:end-1),new_intervals(end-1)+step,inf];
            insertions = 1;
        else
            step_length = (intervals(index+1) - intervals(index))/parts;
            in_intervals = intervals(index)+step_length:step_length:intervals(index+1)-step_length;
            new_intervals = [new_intervals(1:index),in_intervals,new_intervals(index+1:end)];
            insertions = num_roots(i);
        end
        new_indices(ind_before+1:ind_before+insertions) = root_indices(i)+1+ind_before:root_indices(i)+insertions+ind_before;
    end
    new_indices = new_indices(new_indices~=0);
end
end

