function mustBePrimeOrZero(num)
%mustBePrimeOrZero Validate that input is prime or zero.
if ~isprime(num) && num~=0
    eid = "Validate:mustBePrime";
    msg = "Argument must be prime!\n%i is not prime.\nNext prime would be %i.";
    error(eid,msg,num,nextprime(num))
end
end