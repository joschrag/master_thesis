classdef test_FF < matlab.unittest.TestCase

    properties (TestParameter)
        testParameter = struct("num_primes",10^3);
        testParameter2 = struct("scalar",1,"vector",[1 1]);
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods
        function testFF(testCase)
            v = -10^3:10^3;

            for p = primes(testCase.testParameter.num_primes)
                f = FF(v,p);
                testCase.verifyEqual(f.value,mod(v,p))
                testCase.verifyEqual(f.order,p)
            end
            testCase.verifyError(@() FF(v,4),"FF:order")
        end
        function testFFplus(testCase)
            v = -10^3:10^3;
            w = [-10^3/2:10^3/2,1:10^3];

            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = FF(w,p);
                f3 = f1+f2;
                f4 = f2+f1;
                testCase.verifyEqual(f3.value,mod(v+w,p))
                testCase.verifyEqual(f3.order,p)
                testCase.verifyEqual(f4.value,mod(v+w,p))
                testCase.verifyEqual(f4.order,p)
            end
        end
        function testFFminus(testCase)
            v = -10^3:10^3;
            w = [-10^3/2:10^3/2,1:10^3];

            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = FF(w,p);
                f3 = f1-f2;
                f4 = f2-f1;
                testCase.verifyEqual(f3.value,mod(v-w,p))
                testCase.verifyEqual(f3.order,p)
                testCase.verifyEqual(f4.value,mod(w-v,p))
                testCase.verifyEqual(f4.order,p)
            end
        end
        function testFFuplus(testCase)
            v = -10^3:10^3;

            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = +f1;
                testCase.verifyEqual(f2.value,f1.value)
                testCase.verifyEqual(f2.order,f1.order)
            end
        end
        function testFFuminus(testCase)
            v = -10^3:10^3;
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = -f1;
                testCase.verifyEqual(f2.value,mod(-v,p))
                testCase.verifyEqual(f2.order,p)
            end
        end
        function testFFmtimes(testCase)
            u = [31473,-37301,13236;40580,41338,-40245];
            v = [2,23425;1312,124214;123123,23435];
            w = [-1232,243,12323;-1232,-243,766;-4,4536,-9254];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = FF(w,p);
                f3 = FF(u,p);
                f4 = f2*f1;
                f5 = f3*f2;
                f6 = f1*f3;
                testCase.verifyEqual(f4.value,mod(w*v,p))
                testCase.verifyEqual(f4.order,p)
                testCase.verifyEqual(f5.value,mod(u*w,p))
                testCase.verifyEqual(f5.order,p)
                testCase.verifyEqual(f6.value,mod(v*u,p))
                testCase.verifyEqual(f6.order,p)
            end
        end
    end

end