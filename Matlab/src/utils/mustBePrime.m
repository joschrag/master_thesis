function mustBePrime(num)
%MUSTBEPRIME Validate that input is prime.
if ~isprime(num)
    eid = "Validate:mustBePrime";
    msg = "Argument must be prime!\n%i is not prime.\nNext prime would be %i.";
    error(eid,msg,num,nextprime(num))
end
end