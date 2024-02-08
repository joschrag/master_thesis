function tests = test_square_complete
tests = functiontests(localfunctions);
end

function testSquareComplete(testCase)
a = 1;
b = 3;
c = 0;
[actSolution(1),actSolution(2)] = square_complete(a,b,c);
expSolution = [3/2,-9/4];
verifyEqual(testCase,actSolution,expSolution)
end

function testSquareComplete2(testCase)
a = 2;
b = 6;
c = 2;
[actSolution(1),actSolution(2)] = square_complete(a,b,c);
expSolution = [3/2,-5/2];
verifyEqual(testCase,actSolution,expSolution)
end

function testSquareComplete3(testCase)
a = 1;
b = -1;
c = -6;
[actSolution(1),actSolution(2)] = square_complete(a,b,c);
expSolution = [-1/2,-25/4];
verifyEqual(testCase,actSolution,expSolution)
end