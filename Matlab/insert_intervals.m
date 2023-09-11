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
    new_intervals = intervals;
    new_indices = [];

    [sorted_indices,order] = sort(root_indices,'descend');
    sorted_num_roots = num_roots(order);

    for i = 1:numel(root_indices)
        index = sorted_indices(i);
        parts = sorted_num_roots(i)+1;
        if intervals(index) == -inf
            new_intervals = [-inf,new_intervals(2)-step,new_intervals(2:end)];
            new_indices = [new_indices+1,2]; %consider added
        elseif intervals(index+1) == inf
            new_intervals = [new_intervals(1:end-1),new_intervals(end-1)+step,inf];
            new_indices = [new_indices,index+1]; %consider added
        else
            step_length = (intervals(index+1) - intervals(index))/parts;
            in_intervals = intervals(index)+step_length:step_length:intervals(index+1)-step_length;
            new_intervals = [new_intervals(1:index),in_intervals,new_intervals(index+1:end)];
            new_indices = [new_indices+numel(in_intervals),index+1:index+numel(in_intervals)]; %consider added
        end
    end
end
new_indices = sort(new_indices);

end

