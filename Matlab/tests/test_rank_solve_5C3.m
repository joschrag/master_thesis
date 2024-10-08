classdef test_rank_solve_5C3 < matlab.unittest.TestCase

    properties (TestParameter)
        testParameter = struct("num_runs",10^3);
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods
        function test_r1(testCase)
            rng(0)
            for i=1:ceil(testCase.testParameter.num_runs/10)
                c=1;
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [t,u] = rank1_5C3_1(r);
                t = reshape(t,[],1);
                u = reshape(u,[],1);
                if ~isempty(t)
                    for vec = [t,u]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(abs(res)>10^-10)
                            fprintf("M:\n")
                            disp(M)
                            fprintf("vec:\n")
                            disp(vec)
                            fprintf("t:\n")
                            disp(t)
                            fprintf("u:\n")
                            disp(u)
                            fprintf("res:\n")
                            disp(res)                                
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r2(testCase)
            rng(0)
            for i=1:ceil(testCase.testParameter.num_runs/10)
                c=2;
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank1_5C3_2(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(abs(res)>10^-10)
                            disp(M)
                            disp(vec)
                            disp(res)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r3(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=3;
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank1_5C3_3(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(abs(res)>10^-10)
                            disp(M)
                            disp(vec)
                            disp(res)
                            
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r4(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=4;
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank1_5C3_4(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(abs(res)>10^-10)
                            disp(M)
                            disp(vec)
                            disp(res)
                            
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r12(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank2_5C3_12(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r13(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,3];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank2_5C3_13(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r14(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank2_5C3_14(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r23(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[2,3];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank2_5C3_23(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r24(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[2,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank2_5C3_24(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r34(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[3,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank2_5C3_34(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r123(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,3];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank3_5C3_123(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)
                    end
                end
            end
        end
        function test_r124(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank3_5C3_124(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(M)
                            disp(vec)
                            res
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r134(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,3,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank3_5C3_134(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(M)
                            disp(vec)
                            res
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r234(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[2,3,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank3_5C3_234(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                            res
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
        function test_r1234(testCase)
            rng(0)
            for i=1:testCase.testParameter.num_runs
                c=[1,2,3,4];
                [m,r] = rref_matrix(c,5);
                M = [m;zeros(5-numel(c),5)];
                [u,v] = rank4_5C3_1234(r);
                if ~isempty(u)
                    for vec = [u,v]'
                        res = M*[vec(1)^2;vec(2)^2;vec(1);vec(2);1];
                        if any(res)
                            disp(r)
                            disp(M)
                            disp(vec)
                        end
                        testCase.verifyEqual(double(res),[0;0;0;0;0],AbsTol=ones(5,1).*10^-10)

                    end
                end
            end
        end
    end
end