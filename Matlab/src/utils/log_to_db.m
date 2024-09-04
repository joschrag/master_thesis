function [] = log_to_db(c,result,completion_time,tolerance,modulus)
%LOG_TO_DB Summary of this function goes here
%   Detailed explanation goes here
S = dbstack(1);
caller = S.name;
dbfile = "db.sqlite";
conn = sqlite(dbfile);
time = uint64(posixtime(datetime("now",TimeZone="UTC"))*10^6);
data = table(str_to_hash(time),string(caller),modulus,tolerance,completion_time,'VariableNames',["ID", ...
    "mode","modulus","tolerance","completion_time"]);
sqlwrite(conn,"Runs",data);
for row = double(result)'
    data = table(str_to_hash(time),row(1),row(2),row(3),'VariableNames',["ID", ...
        "solution_x","solution_y","solution_z"]);
    sqlwrite(conn,"Solutions",data);
end
for row = c'
    data = table(str_to_hash(time),string(mat2str(row')),'VariableNames',["ID", ...
        "coefficient_list",]);
    sqlwrite(conn,"Coefficients",data);
end
end

