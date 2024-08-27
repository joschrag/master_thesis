classdef test_pa_transformation < matlab.unittest.TestCase

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        function testEllipsoids(testCase)
            expSolution = 3;
            c = [1.0, -3.0, -1.0, 4.0, 0, 2.0, 5.0, -3.0, 4.0, 2.0];
            actSolution = classify_wolfram(c);
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testHyperboloids1(testCase)
            expSolution = 4;
            c = [5.0, -2.0, 0, -4.0, 1.0, 3.0, 5.0, 0, 5.0, -1.0];
            actSolution = classify_wolfram(c);
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testHyperboloids2(testCase)
            expSolution = 5;
            c = [2.0, 3.0, 5.0, 0, 1.0, 4.0, -1.0, -3.0, 2.0, 3.0];
            actSolution = classify_wolfram(c);
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testEllParaboloid(testCase)
            expSolution = 6;
            c = [9.0, 0, 0, 4.0, 0, 0, -18.0, 16.0, -36.0, 25.0];
            actSolution = classify_wolfram(c);
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testHypParaboloid(testCase)
            expSolution = 7;
            c = [2.0, 0, -4.0, -1.0, -2.0, 1.0, -3.0, -3.0, -2.0, 0];
            actSolution = classify_wolfram(c);
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testEllCone(testCase)
            expSolution = 2;
            c = [3.0, 0, -4.0, -1.0, -3.0, 3.0, 0, 0, 0, 0];
            actSolution = classify_wolfram(c);
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testEllCylinder(testCase)
            c = [1.0, 0, 0, 1.0, 0, 0, 0, 0, 0, -10.0];
            actSolution = classify_wolfram(c);
            expSolution = 8;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testHypCylinder(testCase)
            c = [1.0, 0, 0, -1.0, 0, 0, 0, 0, 0, -10.0];
            actSolution = classify_wolfram(c);
            expSolution = 10;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testParCylinder(testCase)
            c = [1.0, 0, 0, 0, 0, 0, 0, -6.0, 0, 0];
            actSolution = classify_wolfram(c);
            expSolution = 11;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testCrossPlanes(testCase)
            c = [1.0, 0, 0, -2.0, 0, 0, 0, 0, 0, 0];
            actSolution = classify_wolfram(c);
            expSolution = 14;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testParPlanes(testCase)
            c = [4,0,0,0,0,0,0,0,0,-10];
            actSolution = classify_wolfram(c);
            expSolution = 12;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testOnePlane(testCase)
            c = [14,0,0,0,0,0,0,0,0,0];
            actSolution = classify_wolfram(c);
            expSolution = 13;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testLine(testCase)
            c = [1.0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0];
            actSolution = classify_wolfram(c);
            expSolution = 9;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end

        function testSinglePoint(testCase)
            c = [1.0, 0, 0, 1.0, 0, 1.0, 0, 0, 0, 0];
            actSolution = classify_wolfram(c);
            expSolution = 1;
            verifyEqual(testCase,actSolution,expSolution)
            pa_transformation(c);
        end
    end
end