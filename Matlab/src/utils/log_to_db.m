function log_to_db(C,result,completion_time,tolerance,modulus)
%LOG_TO_DB Log run results to sqlite database.
arguments
    C (:,:) double {mustBeReal};
    result (:,3);
    completion_time (1,1) double;
    tolerance (1,1) {mustBeReal};
    modulus (1,1) {mustBeInteger,mustBeNonnegative};
end
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
for row = C'
    data = table(str_to_hash(time),string(mat2str(row')),'VariableNames',["ID", ...
        "coefficient_list",]);
    sqlwrite(conn,"Coefficients",data);
end
end

