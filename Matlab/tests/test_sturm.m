classdef test_sturm < matlab.unittest.TestCase

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
        
        function test_sturm_fun(testCase)
            rng(0)
            for d=1:20
                for i=1:testCase.testParameter.num_runs
                    p = 5.*rand(1,d)-2.5;
                    root = roots(p);
                    root = sort(root(imag(root)==0));
                    roots_sturm = sturm(p);
                    if isempty(roots_sturm)
                        testCase.assertTrue(isempty(root))
                    else
                        testCase.verifyEqual(size(roots_sturm),size(root))
                        testCase.verifyEqual(roots_sturm,root,AbsTol=ones(numel(roots_sturm),1).*10^-10)
                    end
                end
            end
        end
    end
    
end