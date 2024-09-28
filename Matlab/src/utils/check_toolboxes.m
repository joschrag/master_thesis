function opt = check_toolboxes(opt)
arguments
    opt struct;
end
S = dbstack(1);
caller = S.name;
required_toolboxes = "Symbolic Math Toolbox";
if endsWith(caller,"_Fp")
    required_toolboxes = [required_toolboxes,"Communications Toolbox"];
end
version = ver;
installed_toolboxes = {version.Name};
missing_toolboxes = string.empty(0,numel(required_toolboxes));
for i =1:numel( required_toolboxes)
    if ~ismember(required_toolboxes(i),installed_toolboxes)
        missing_toolboxes(i) = required_toolboxes(i);
    end
end
if ~isempty(missing_toolboxes)
    error("Toolbox %s is required to run the program!",join(missing_toolboxes,","))
end

if ~ismember("Database Toolbox",installed_toolboxes)
    opt.log_db = 0;
end

end