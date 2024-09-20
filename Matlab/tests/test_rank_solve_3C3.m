classdef test_rank_solve_3C3 < matlab.unittest.TestCase

    properties (TestParameter)
        testParameter = struct("num_runs",10^5);
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods
        function test_r1234(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,3,4];
                [m,r] = rref_matrix(c,6,mode="integer");
                M = [m;zeros(6-numel(c),6)];
                [u,v] = rank4_3C3_1234(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0;0],AbsTol=ones(6,1).*10^-10)
                    end
                end
            end
        end
        function test_r1235(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,3,5];
                [m,r] = rref_matrix(c,6,mode="integer");
                M = [m;zeros(6-numel(c),6)];
                [u,v] = rank4_3C3_1235(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0;0],AbsTol=ones(6,1).*10^-10)
                    end
                end
            end
        end
        function test_r1245(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,4,5];
                [m,r] = rref_matrix(c,6,mode="integer");
                M = [m;zeros(6-numel(c),6)];
                [u,v] = rank4_3C3_1245(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0;0],AbsTol=ones(6,1).*10^-10)
                    end
                end
            end
        end
        function test_r1345(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,3,4,5];
                [m,r] = rref_matrix(c,6,mode="integer");
                M = [m;zeros(6-numel(c),6)];
                [u,v] = rank4_3C3_1345(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0;0],AbsTol=ones(6,1).*10^-10)
                    end
                end
            end
        end
        function test_r2345(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[2,3,4,5];
                [m,r] = rref_matrix(c,6,mode="integer");
                M = [m;zeros(6-numel(c),6)];
                [u,v] = rank4_3C3_2345(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0;0],AbsTol=ones(6,1).*10^-10)
                    end
                end
            end
        end
        function test_r12345(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,3,4,5];
                [m,r] = rref_matrix(c,6,mode="integer");
                M = [m;zeros(6-numel(c),6)];
                [u,v] = rank5_3C3_12345(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0;0],AbsTol=ones(6,1).*10^-10)
                    end
                end
            end
        end
    end
end