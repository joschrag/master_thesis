classdef test_FF < matlab.unittest.TestCase

    properties (TestParameter)
        testParameter = struct("num_primes",10^2);
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
                testCase.assertEqual(f.value,sym(mod(v,p)))
                testCase.assertEqual(f.order,p)
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
                testCase.assertEqual(f3.value,sym(mod(v+w,p)))
                testCase.assertEqual(f3.order,p)
                testCase.assertEqual(f4.value,sym(mod(v+w,p)))
                testCase.assertEqual(f4.order,p)
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
                testCase.assertEqual(f3.value,sym(mod(v-w,p)))
                testCase.assertEqual(f3.order,p)
                testCase.assertEqual(f4.value,sym(mod(w-v,p)))
                testCase.assertEqual(f4.order,p)
            end
        end
        function testFFuplus(testCase)
            v = -10^3:10^3;

            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = +f1;
                testCase.assertEqual(f2.value,f1.value)
                testCase.assertEqual(f2.order,f1.order)
            end
        end
        function testFFuminus(testCase)
            v = -10^3:10^3;
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = -f1;
                testCase.assertEqual(f2.value,sym(mod(-v,p)))
                testCase.assertEqual(f2.order,p)
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
                testCase.assertEqual(f4.value,sym(mod(w*v,p)))
                testCase.assertEqual(f4.order,p)
                testCase.assertEqual(f5.value,sym(mod(u*w,p)))
                testCase.assertEqual(f5.order,p)
                testCase.assertEqual(f6.value,sym(mod(v*u,p)))
                testCase.assertEqual(f6.order,p)
            end
        end
        function testFFrdivide(testCase)
            u = [31473,-37301,13236;40580,41338,-40245];
            v = [2,2325,1317;12413,12123,23435];
            w = [2;25];

            for p = primes(testCase.testParameter.num_primes)
                r = mod(312312,p);
                u = mod(u,p);
                v = mod(v,p);
                w = mod(w,p);
                u(u==0) = 1;
                v(v==0) = 1;
                w(w==0) = 1;
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(w,p);
                f4 = f1./f2;
                f5 = f1./f3;
                f6 = f1./r;
                testCase.assertEqual(f4.value,sym(mod(u./v,p)))
                testCase.assertEqual(f4.order,p)
                testCase.assertEqual(f5.value,sym(mod(u./w,p)))
                testCase.assertEqual(f5.order,p)
                testCase.assertEqual(f6.value,sym(mod(u./r,p)))
                testCase.assertEqual(f6.order,p)
            end
        end

        function testFFldivide(testCase)
            u = [31473,-37301,13236;40580,41338,-40245];
            v = [2,2325,1317;12413,12123,23435];
            w = [2;25];

            for p = primes(testCase.testParameter.num_primes)
                r = mod(312312,p);
                u = mod(u,p);
                v = mod(v,p);
                w = mod(w,p);
                u(u==0) = 1;
                v(v==0) = 1;
                w(w==0) = 1;
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(w,p);
                f4 = f1.\f2;
                f5 = f1.\f3;
                f6 = f1.\r;
                testCase.assertEqual(f4.value,sym(mod(u.\v,p)))
                testCase.assertEqual(f4.order,p)
                testCase.assertEqual(f5.value,sym(mod(u.\w,p)))
                testCase.assertEqual(f5.order,p)
                testCase.assertEqual(f6.value,sym(mod(u.\r,p)))
                testCase.assertEqual(f6.order,p)
            end
        end

        function testFFpower(testCase)
            u = [31473,-37301,13236;40580,41338,-40245];
            v = [2,2325,1317;12413,12123,23435];
            w = [2;25];
            powers = -10:10;
            for p = primes(testCase.testParameter.num_primes)
                for pow = powers
                    u = mod(u,p);
                    v = mod(v,p);
                    w = mod(w,p);
                    if pow < 0
                        u(u==0) = 1;
                        v(v==0) = 1;
                        w(w==0) = 1;
                    end
                    f1 = FF(u,p);
                    f2 = FF(v,p);
                    f3 = FF(w,p);
                    f4 = f1.^pow;
                    f5 = f2.^pow;
                    f6 = f3.^pow;
                    if pow >= 0
                        testCase.assertEqual(f4.value,sym(mod(mod(u,p).^pow,p)))
                        testCase.assertEqual(f5.value,sym(mod(mod(v,p).^pow,p)))
                        testCase.assertEqual(f6.value,sym(mod(mod(w,p).^pow,p)))
                    end
                    testCase.assertEqual(f4.*f1.^(-pow),FF(ones(size(u)),p))
                    testCase.assertEqual(f5.*f2.^(-pow),FF(ones(size(v)),p))
                    testCase.assertEqual(f6.*f3.^(-pow),FF(ones(size(w)),p))
                    testCase.assertEqual(f4.order,p)
                    testCase.assertEqual(f5.order,p)
                    testCase.assertEqual(f6.order,p)
                end
            end
        end

        function testFFmldivide(testCase)
            u = [31473,-37301,13236;40580,41338,-40245;123,5334,34];
            v = [2,2325;1317,12413;12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = f1\f2;
                if any(size(f3.value))
                    testCase.assertEqual(f2,f1*f3)
                end
            end
        end

        function testFFmrdivide(testCase)
            u = [31473,-37301,13236;40580,41338,-40245;123,5334,34];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = f2/f1;
                if any(size(f3.value))
                    testCase.assertEqual(f2,f3*f1)
                end
            end
        end
        function testFFctranspose(testCase)
            u = [-37301,13236;40580,-40245;123,5334];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = f1.';
                f4 = f2.';
                testCase.assertEqual(f3.value,f1.value.')
                testCase.assertEqual(f3.order,f1.order)
                testCase.assertEqual(f4.value,f2.value.')
                testCase.assertEqual(f4.order,f2.order)

            end
        end
        function testFFtranspose(testCase)
            u = [-37301,13236;40580,-40245;123,5334];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = f1';
                f4 = f2';
                testCase.assertEqual(f3.value,f1.value')
                testCase.assertEqual(f3.order,f1.order)
                testCase.assertEqual(f4.value,f2.value')
                testCase.assertEqual(f4.order,f2.order)
            end
        end
        function testFFeq(testCase)
            u = [-37301,13236;40580,-40245;123,5334];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(u,p);
                f4 = FF(u,nextprime(p+1));
                testCase.verifyFalse(f1==f2)
                testCase.verifyFalse(f1==f4)
                testCase.verifyTrue(f1==f3)
            end
        end
        function testFFne(testCase)
            u = [-37301,13236;40580,-40245;123,5334];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(u,p);
                f4 = FF(u,nextprime(p+1));
                testCase.assertTrue(f1~=f2)
                testCase.assertTrue(f1~=f4)
                testCase.assertFalse(f1~=f3)
            end
        end
        function testFFhorzcat(testCase)
            u = [-37301,13236;40580,-40245;123,5334];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                v1 = [f1,f2,f1];
                v2 = [f1,f2',f1,f2];
                testCase.assertEqual(size(v1),[1,3])
                testCase.assertEqual(size(v2),[1,4])
            end
        end
        function testFFvertcat(testCase)
            u = [-37301,13236;40580,-40245;123,5334];
            v = [2,2325,1317;12413,12123,23435];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                v1 = [f1;f2;f1];
                v2 = [f1;f2';f1;f2];
                testCase.assertEqual(size(v1),[3,1])
                testCase.assertEqual(size(v2),[4,1])
            end
        end

        function testFFmpower(testCase)
            u = [31473,-37301,13236;40580,41338,-40245;424,355,123];
            v = [2,2325;12413,23435];
            w = [2;25];
            x = [0,0;1,1];
            powers = -10:10;
            for p = primes(testCase.testParameter.num_primes)
                for pow = powers

                    f1 = FF(u,p);
                    f2 = FF(v,p);
                    f3 = FF(w,p);
                    f4 = FF(x,p);
                    if (pow < 0) && (det(f1)~=0)
                        f5 = f1^pow;
                        if pow >= 0
                            testCase.assertEqual(f5.value,mod(mod(u,p).^pow,p))
                        end
                        testCase.assertEqual(f5*f1^(-pow),FF(eye(min(size(u))),p))
                        testCase.assertEqual(f5.order,p)
                    end
                    if (pow < 0) && (rank(f2)==size(v,1))
                        f6 = f2^pow;
                        if pow >= 0
                            testCase.assertEqual(f6.value,mod(mod(v,p).^pow,p))
                        end
                        testCase.assertEqual(f6*f2^(-pow),FF(eye(min(size(v))),p))
                        testCase.assertEqual(f6.order,p)
                    end
                    testCase.verifyError(@() f3^pow,"FF:mpower:NonSquareMatrix")
                    if (pow < 0)
                        testCase.verifyError(@() f4^pow,"FF:inv:SingularMatrix")
                    end
                end
            end
        end
        function testFFinv(testCase)
            u = [31473,-37301,13236;40580,41338,-40245;424,355,123];
            v = [2,2325,1317;12413,12123,23435];
            w = [0,0;1,2];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(w,p);
                if det(f1) ~= 0
                    f4 = inv(f1);
                    testCase.assertEqual(f1*f4,FF(eye(min(size(u))),p))
                    testCase.assertEqual(f4*f1,FF(eye(min(size(u))),p))
                end
                testCase.verifyError(@() inv(f2),"FF:inv:NonSquareMatrix")
                testCase.verifyError(@() inv(f3),"FF:inv:SingularMatrix")
            end
        end
        function testFFadj(testCase)
            u = [31473,-37301,13236;40580,41338,-40245;424,355,123];
            v = [2,2325,1317;12413,12123,23435];
            w = [0,3;1,2];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(w,p);
                f4 = adj(f1);
                testCase.assertEqual(size(f1.value),size(f4.value))
                testCase.assertEqual(f1.order,f4.order)
                testCase.verifyError(@() adj(f2),"FF:adj:NonSquareMatrix")
                f6 = adj(f3);
                testCase.assertEqual(size(f3.value),size(f6.value))
                testCase.assertEqual(f3.order,f6.order)
            end
        end
        function testFFrank(testCase)
            u = [31473,-37301,13236;40580,41338,-40245;424,355,123];
            v = [2,2325,1317;12413,12123,23435];
            w = [0,3;1,2];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(u,p);
                f2 = FF(v,p);
                f3 = FF(w,p);
                for f = [f1,f2,f3]
                    if size(f.value,1) == size(f.value,2)
                        if det(f)==0
                            testCase.assertLessThan(rank(f),min(size(f.value)))
                        else
                            testCase.assertEqual(rank(f),min(size(f.value)))
                        end
                    else
                        testCase.assertLessThanOrEqual(rank(f),min(size(f.value)))
                    end
                end
            end
        end
        function testFFsum(testCase)
            v = -10^3:10^3;
            w = -5^3:15^3;
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = FF(w,p);
                s1 = sum(f1);
                s2 = sum(f2);
                testCase.assertEqual(s1.value,mod(sum(f1.value),p))
                testCase.assertEqual(s2.value,mod(sum(f2.value),p))
                testCase.assertEqual(s1.order,p)
                testCase.assertEqual(s2.order,p)
            end
        end
        function testFFl_inv(testCase)
            v = [3,453465,3453;6768,234235,234243;1234,123445,12323];
            w = [7,18;234324,1233;134,987;12313,76457;1357,32424];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = FF(w,p);
                if rank(f1) == size(v,2)
                    l1 = l_inv(f1);
                    testCase.verifyEqual(l1*f1,FF(eye(size(v,2)),p))
                    testCase.assertEqual(l1.order,p)
                end
                if rank(f2) == size(w,2)
                    l2 = l_inv(f2);
                    testCase.verifyEqual(l2*f2,FF(eye(size(w,2)),p))
                    testCase.assertEqual(l2.order,p)
                end

            end
        end
                function testFFr_inv(testCase)
            v = [3,453465,3453;6768,234235,234243;1234,123445,12323];
            w = [7,18,2344,1233;134,987,12313,76457;1357,32424,1314,89347];
            for p = primes(testCase.testParameter.num_primes)
                f1 = FF(v,p);
                f2 = FF(w,p);
                if rank(f1) == size(v,1)
                    r1 = r_inv(f1);
                    testCase.verifyEqual(f1*r1,FF(eye(size(v,1)),p))
                    testCase.assertEqual(r1.order,p)
                end
                if rank(f2) == size(w,1)
                    r2 = r_inv(f2);
                    testCase.verifyEqual(f2*r2,FF(eye(size(w,1)),p))
                    testCase.assertEqual(r2.order,p)
                end

            end
        end
    end

end