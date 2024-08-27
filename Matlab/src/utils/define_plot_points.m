function [points] = define_plot_points(range)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
arguments
    range (1,2) {mustBeReal};
end
lower_bound = -5;
upper_bound = 5;
points = horzcat(range(1):0.5:lower_bound,...
    max(lower_bound,range(1)):0.1:min(upper_bound,range(2)),...
    upper_bound:0.5:range(2));
end