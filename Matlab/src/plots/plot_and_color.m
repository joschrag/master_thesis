function plot_and_color(c,opt)
%PLOT_AND_COLOR Summary of this function goes here
%   Detailed explanation goes here
arguments
    c (:,10) {mustBeReal};
    opt.verbose (1,1) {mustBeInteger, mustBeInRange(opt.verbose,0,2)} = 0;
    opt.plotRange (1,2) {mustBeReal} = [-5,5];
    opt.tolerance (1,1) {mustBeReal,mustBePositive} = 10^-10;
end
f = gcf;
for cc =c'
    pa_transformation(cc,plotRange=opt.plotRange,tolerance=opt.tolerance)
end
plt_list = f.Children.Children;
color_list = ["#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"];
j =1;
skip = false;
for i=numel(plt_list):-1:1
    if opt.verbose > 0
        fprintf("%i : %s | %s\n",i,class(plt_list(i)),plt_list(i).DisplayName)
    end
    switch plt_list(i).DisplayName
        case {"hyperboloid2","ell cone","cross planes","par planes","hyp cylinder"}
            plt_list(i).FaceColor = color_list(j);
            plt_list(i).FaceAlpha = 0.4;
            if skip
                j = j+1;
                skip = false;
                continue
            end
            skip = true;
        case "solution"
            continue
        case "line"
            plt_list(i).Color = color_list(j);
            j = j+1;
        case "point"
            plt_list(i).MarkerEdgeColor = color_list(j);
            plt_list(i).MarkerFaceColor = color_list(j);
            j = j+1;
        otherwise
            plt_list(i).FaceColor = color_list(j);
            plt_list(i).FaceAlpha = 0.4;
            skip = false;
            j = j+1;
    end
end
xlabel("x")
ylabel("y")
zlabel("z")
end

