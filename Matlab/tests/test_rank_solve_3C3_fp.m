classdef test_rank_solve_3C3_fp < matlab.unittest.TestCase

    properties (TestParameter)
        testParameter = struct("num_runs",10^2);
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods
        function test_r123(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,3];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_123_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
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
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_124_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r125(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_125_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
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
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_134_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end

        function test_r135(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,3,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_135_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r145(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_145_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
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
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_234_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r235(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[2,3,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_235_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r245(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[2,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_245_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r345(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[3,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank3_3C3_345_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
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
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank4_3C3_1234_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r1235(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,3,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank4_3C3_1235_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r1245(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank4_3C3_1245_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r1345(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,3,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank4_3C3_1345_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r2345(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[2,3,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank4_3C3_2345_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(r)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
        function test_r12345(testCase)
            for p =[3,5,7,11]
                rng(p)
                for i=1:testCase.testParameter.num_runs
                    c=[1,2,3,4,5];
                    [m,r] = rref_matrix(c,6,p);
                    M = FF([m;zeros(6-numel(c),6)],p);
                    [u,v] = rank5_3C3_12345_Fp(r,p);
                    if ~isempty(u)
                        for vec = [u,v]'
                            res = M*FF([vec(2)^3;vec(1)^2;vec(2)^2;vec(1);vec(2);1],p);
                            if any(res.value)
                                disp(p)
                                disp(FF(-r,p).value)
                                disp(M.value)
                                disp(vec)
                                res.value
                                for t0=1:p
                                    for u0=1:p
                                        cntrl = M*FF([u0^3;t0^2;u0^2;t0;u0;1],p);
                                        if ~any(cntrl.value)
                                            fprintf("%i,%i;%i,%i | %i\n",t0,u0,vec,p)
                                        end
                                    end
                                end
                            end
                            testCase.verifyEqual(res.value,cast([0;0;0;0;0;0],class(res.value)))
                            testCase.verifyEqual(res.order,p)
                        end
                    end
                end
            end
        end
    end
end