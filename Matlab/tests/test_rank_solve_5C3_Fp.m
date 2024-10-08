classdef test_rank_solve_5C3_Fp < matlab.unittest.TestCase
    
    properties (TestParameter)
        testParameter = struct("num_runs",5);
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
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=1;
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank1_5C3_1_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                fprintf("M:\n")
                                disp(M.value)
                                fprintf("vec:\n")
                                disp(vec)
                                fprintf("t:\n")
                                disp(u)
                                fprintf("u:\n")
                                disp(v)
                                fprintf("res:\n")
                                disp(res.value)
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r2(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=2;
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank1_5C3_2_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r3(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=3;
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank1_5C3_3_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r4(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=4;
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank1_5C3_4_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r12(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank2_5C3_12_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r13(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,3];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank2_5C3_13_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r14(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank2_5C3_14_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r23(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[2,3];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank2_5C3_23_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r24(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[2,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank2_5C3_24_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r34(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[3,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank2_5C3_34_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r123(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,3];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank3_5C3_123_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r124(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank3_5C3_124_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r134(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,3,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank3_5C3_134_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r234(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[2,3,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank3_5C3_234_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r1234(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,3,4];
                    [m,r] = rref_matrix(c,5,p);
                    M = FF([m;zeros(5-numel(c),5)],p);
                    [u,v] = rank4_5C3_1234_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
    end
end