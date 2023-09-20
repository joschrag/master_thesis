function plot_and_color(c)
%PLOT_AND_COLOR Summary of this function goes here
%   Detailed explanation goes here
f = gcf;
pa_transformation(c(1,:))
pa_transformation(c(2,:))
pa_transformation(c(3,:))
plt_list = f.Children.Children;
color_list = [[1,0.411764705882353,0.160784313725490];...
    [0.058823529411765,1,1];...
    [1,0,1];...
    [0.717647058823529,0.274509803921569,1];...
    [0.650980392156863,0.650980392156863,0.650980392156863];...
    [0,0,0]];
j =1;
skip = false;
for i=1:numel(plt_list)
    if ~isa(plt_list(i),"matlab.graphics.chart.primitive.Scatter")
        %fprintf("%i : %s\n",i,class(plt_list(i)))
        switch plt_list(i).DisplayName
            case {"2 Hyperboloid","Ell Cone"}
                plt_list(i).FaceColor = color_list(j,:);
                plt_list(i).FaceAlpha = 0.7;
                if skip
                    j = j+1;
                    continue
                end
                skip = true;
            otherwise
                plt_list(i).FaceColor = color_list(j,:);
                plt_list(i).FaceAlpha = 0.7;
                skip = false;
                j = j+1;
        end
    end
end
end

