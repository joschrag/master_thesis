function mustBeSizeFF(obj,obj_size)
% Test for equal size
if size(obj.value)~=obj_size
    eid = "FF:validate:sizeEquals";
    msg = "Size of FF value does not match specified size!\nSize of obj: %s\nSize specified: %s";
    error(eid,msg,mat2str(size(obj.value)),mat2str(obj_size))
end
end