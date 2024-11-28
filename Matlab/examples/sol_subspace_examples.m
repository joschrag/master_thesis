%% 1sol
M = [1,-2,-3,4,5;zeros(4,5)];
solve_subsystem(M,1,1:3);
%% 1no sol
M = [1,2,3,4,5;zeros(4,5)];
solve_subsystem(M,1,1:3);
%% 2sol
M = [0,1,-3,4,5;zeros(4,5)];
solve_subsystem(M,1,1:3);
%% 3sol     WRONG!
M = [0,0,1,4,5;zeros(4,5)];
solve_subsystem(M,1,1:3);
%% 3no sol  WRONG!
M = [0,0,1,4,5;zeros(4,5)];
solve_subsystem(M,1,1:3)
%% 4sol     WRONG!
M = [0,0,0,1,5;zeros(4,5)];
solve_subsystem(M,1,1:3)
%% 4no sol  WRONG!
M = [0,0,0,1,-2;zeros(4,5)];
solve_subsystem(M,1,1:3)
%% 12two sol
M = [1,0,1,2,-1;...
    0,1,2,1,-1;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 12four sol
M = [1,0,1,2,-1;...
    0,1,-2,5,-4;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 12no sol
M = [1,0,2,0,2;...
    0,1,-2,3,3;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 13two sol    2nd Look!
M = [1,-1,0,2,-4;...
    0,0,1,-2,-1;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 13no sol
M = [1,3,0,5,-1;...
    0,0,1,1,3;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 14two sol    2nd Look!
M = [1,-1,1,0,-4;...
    0,0,0,1,-2;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 14no sol
M = [1,3,0,0,-1;...
    0,0,0,1,3;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 23two sol    2nd Look!
M = [0,1,0,1,-4;...
    0,0,1,1,-2;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 23no sol
M = [0,1,0,1,-1;...
    0,0,1,1,3;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 24
M = [0,1,1,0,-4;...
    0,0,0,1,-2;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 24
M = [0,1,1,0,-1;...
    0,0,0,1,3;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 34
M = [0,0,1,0,-4;...
    0,0,0,1,-2;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 34
M = [0,0,1,0,-1;...
    0,0,0,1,3;zeros(3,5)];
solve_subsystem(M,1,1:3)
%% 123 one sol
M = [1,0,0,1,-3;...
    0,1,0,0,-1;...
    0,0,1,1,-1;zeros(2,5)];
solve_subsystem(M,1,1:3)
%% 124 one sol
M = [1, 0,-6, 0,8;...
     0, 1,-2, 0, 5;...
     0, 0, 0, 1,-1;zeros(2,5)];
solve_subsystem(M,1,1:3)
%% 234
while true
M = [rref_matrix([2,3,4]);zeros(2,5)];
res = solve_subsystem(M,1,1:3,plot_subspace=0);
if ~isempty(res)
    disp(M)
    break
end
end
solve_subsystem(M,1,1:3,plot_subspace=1);
%% 1234
M = [1,0,0,0,-4;...
    0,1,0,0,-9;...
    0,0,1,0,2;
    0,0,0,1,-3;zeros(1,5)];
solve_subsystem(M,1,1:3,plot_subspace=1);