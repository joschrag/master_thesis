function plot_and_color(c)
%PLOT_AND_COLOR Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (:,10) {mustBeReal}
end
f = gcf;
for cc =c'
    pa_transformation(cc)
end
plt_list = f.Children.Children;
color_list = [[1,0.411764705882353,0.160784313725490];...
    [0.058823529411765,1,1];...
    [1,0,1];...
    [0.960784313725490   1.000000000000000   0.721568627450980];...
    [0,     1,     0];...
    [0.650980392156863,0.650980392156863,0.650980392156863];...
    [0.717647058823529,0.274509803921569,1];...
    [0.494117647058824  , 0.184313725490196  , 0.556862745098039];...
    [1.000000000000000,   0.074509803921569,   0.650980392156863];...
    [0,0,0]];
j =1;
skip = false;
for i=1:numel(plt_list)
    fprintf("%i : %s | %s\n",i,class(plt_list(i)),plt_list(i).DisplayName)
    switch plt_list(i).DisplayName
        case {"hyperboloid2","ell cone","cross planes","par planes","hyp cylinder"}
            plt_list(i).FaceColor = color_list(j,:);
            plt_list(i).FaceAlpha = 0.7;
            if skip
                j = j+1;
                skip = false;
                continue
            end
            skip = true;
        case "solution"
            continue
        case "line"
            plt_list(i).Color = color_list(j,:);
            j = j+1;
        case "point"
            plt_list(i).MarkerEdgeColor = color_list(j,:);
            plt_list(i).MarkerFaceColor = color_list(j,:);
            j = j+1;
        otherwise
            plt_list(i).FaceColor = color_list(j,:);
            plt_list(i).FaceAlpha = 0.7;
            skip = false;
            j = j+1;
    end
end
xlabel("x")
ylabel("y")
zlabel("z")
end

