function [M,D] = rref_matrix(ind_cols,len,prime,elements,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    ind_cols (1,:) {mustBeReal,mustBeVector}
    len      (1,1) {mustBeInteger,mustBePositive}
    prime    (1,1) {mustBeInteger,mustBeNonnegative} = 0
    elements (1,:) {mustBeInteger,mustBePositive,mustBeVector} = []
    opt.mode (1,1) string {mustBeMember(opt.mode,["real","integer"])} = "real";
end
dep_cols = setdiff(1:len,ind_cols);
M = zeros(numel(ind_cols),len);
if prime ~= 0
    rand_func =@(sz) randi(prime,1,sz)-1;
elseif opt.mode == "integer"
    rand_func =@(sz) randi(20,1,sz)-10;
else
    rand_func =@(sz) 10.*rand(1,sz)-5;
end

for d=dep_cols
    sz = sum(ind_cols<d);
    if sz > 0
        if ~isempty(elements)
            M(1:sz,d) = elements(1:sz);
            elements = elements(sz:end);
        else
            M(1:sz,d) = rand_func(sz);
        end
    end
end
M(:,ind_cols) = eye(numel(ind_cols));
D = M(:,dep_cols);
end