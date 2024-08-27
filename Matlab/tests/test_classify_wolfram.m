function tests = test_classify_wolfram
tests = functiontests(localfunctions);
end

function testEllipsoids(testCase)
c = [1,0,0,1,0,1,0,0,0,-10];
actSolution = classify_wolfram(c);
expSolution = 3;
verifyEqual(testCase,actSolution,expSolution)
c = [2.15707, 1.31901, 0.727695, 4.55324, 0.680343, 0.909235, 4.34646, 2.89852, 2.7493, 0.724774];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [1.0, -3.0, -1.0, 4.0, 0, 2.0, 5.0, -3.0, 4.0, 2.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
end

function testHyperboloids1(testCase)
c = [1.0, 0, 0, 1.0, 0, -1.0, 0, 0, 0, -10.0];
actSolution = classify_wolfram(c);
expSolution = 4;
verifyEqual(testCase,actSolution,expSolution)
c = [0.379271, 3.89584, 4.67005, 0.269751, 0.649531, 2.65399, 2.84412, 2.34695, 0.0595103, 1.68561];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [-1.0, 0, 6.0, 3.0, 0, -1.0, 0, 0, 0, -1.0];

actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [-1.0, -4.0, -3.0, 4.0, 2.0, -4.0, 3.0, 2.0, 0, 1.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [0, -3.0, -3.0, -4.0, -1.0, 5.0, -3.0, 0, -1.0, 5.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [5,-4,3,-2,0,1,5,0,5,-1];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [5.0, -2.0, 0, -4.0, 1.0, 3.0, 5.0, 0, 5.0, -1.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
end

function testHyperboloids2(testCase)
c = [-1.0, 0, 0, 1.0, 0, -1.0, 0, 0, 0, -10.0];
actSolution = classify_wolfram(c);
expSolution = 5;
verifyEqual(testCase,actSolution,expSolution)
c = [1.7583, 2.74862, 4.58597, 4.15414, 1.4292, 2.92632, 3.786, 3.76865, 1.90223, 2.83911];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [2.25271, 4.56669, 0.76189, 0.419107, 4.12908, 1.14488, 2.69171, 4.98067, 0.390878, 2.21339];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [-1.0, 0, -3.0, -3.0, 1.0, 0, -2.0, -1.0, 1.0, -2.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [4.0, 1.0, -4.0, 4.0, 0, -2.0, -1.0, -3.0, -3.0, 0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [-4.0, 2.0, 2.0, 1.0, 2.0, 0, -4.0, -4.0, -1.0, 1.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [2.0, 3.0, 5.0, 0, 1.0, 4.0, -1.0, -3.0, 2.0, 3.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
end

function testEllParaboloid(testCase)
c = [1.0, 0, 0, 1.0, 0, 0, 0, 0, -4.0, 0];
actSolution = classify_wolfram(c);
expSolution = 6;
verifyEqual(testCase,actSolution,expSolution)
c = [-1.0, 0, 0, -1.0, 0, 0, 0, 0, -1.0, 6.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [9.0, 0, 0, 4.0, 0, 0, -18.0, 16.0, -36.0, 25.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
end

function testHypParaboloid(testCase)
c = [1.0, 0, 0, -1.0, 0, 0, 0, 0, -4.0, 0];
actSolution = classify_wolfram(c);
expSolution = 7;
verifyEqual(testCase,actSolution,expSolution)
c = [0, 0, 0, -5.0, 0, 2.0, 1.0, 0, 0, 0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [2.0, 0, -4.0, -1.0, -2.0, 1.0, -3.0, -3.0, -2.0, 0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
end

function testEllCone(testCase)
c = [1.0, 0, 0, 1.0, 0, -1.0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
expSolution = 2;
verifyEqual(testCase,actSolution,expSolution)
c = [1.0, -2.0, -2.0, 1.0, -2.0, 1.0, 2.0, 2.0, 2.0, -3.0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [2.0, -2.0, -2.0, 2.0, 2.0, -4.0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [4.0, 2.0, -4.0, -1.0, 2.0, 3.0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
c = [3.0, 0, -4.0, -1.0, -3.0, 3.0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
verifyEqual(testCase,actSolution,expSolution)
end

function testEllCylinder(testCase)
c = [1.0, 0, 0, 1.0, 0, 0, 0, 0, 0, -10.0];
actSolution = classify_wolfram(c);
expSolution = 8;
verifyEqual(testCase,actSolution,expSolution)
end

function testHypCylinder(testCase)
c = [1.0, 0, 0, -1.0, 0, 0, 0, 0, 0, -10.0];
actSolution = classify_wolfram(c);
expSolution = 10;
verifyEqual(testCase,actSolution,expSolution)
end

function testParCylinder(testCase)
c = [1.0, 0, 0, 0, 0, 0, 0, -6.0, 0, 0];
actSolution = classify_wolfram(c);
expSolution = 11;
verifyEqual(testCase,actSolution,expSolution)
end

function testCrossPlanes(testCase)
c = [1.0, 0, 0, -2.0, 0, 0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
expSolution = 14;
verifyEqual(testCase,actSolution,expSolution)
end

function testParPlanes(testCase)
c = [4,0,0,0,0,0,0,0,0,-10];
actSolution = classify_wolfram(c);
expSolution = 12;
verifyEqual(testCase,actSolution,expSolution)
end

function testOnePlane(testCase)
c = [14,0,0,0,0,0,0,0,0,0];
actSolution = classify_wolfram(c);
expSolution = 13;
verifyEqual(testCase,actSolution,expSolution)
end

function testLine(testCase)
c = [1.0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
expSolution = 9;
verifyEqual(testCase,actSolution,expSolution)
end

function testSinglePoint(testCase)
c = [1.0, 0, 0, 1.0, 0, 1.0, 0, 0, 0, 0];
actSolution = classify_wolfram(c);
expSolution = 1;
verifyEqual(testCase,actSolution,expSolution)
end
