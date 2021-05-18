function boolVal=isProbablePrime(n,numRuns)
%%ISPROBABLEPRIME Determine whether a positive integer is prime, with an
%                 algorithm that has a small probablity of falsely
%                 declaring a non-prime to be a prime, but which is faster
%                 than dbrute force division.
%
%INPUTS: n A real positive integer. This value must be <=sqrt(flintmax).
%          The algorithm gets slower for smaller numbers.
%  numRuns An optional parameter specifying how many times the algorithm is
%          rerun when it initially detects a probable prime. The default if
%          omitted or an empty matrix is passed is 25. The more runs, the
%          more conficence one can have in the solution. In Chapter 4.5.4
%          of [1], it is noted that the algorithm might have a pase
%          positive up to 1/4 of the time, so the likelihood of having 25
%          of them is about (1/4)^25=8.8818e-16 is very slow.
%
%OUTPUTS: boolVal This is true if n is a probable prime and false if n is
%                 definitely not a prime.
%
%This function implements Algorithm P in Section 4.5.4 of [1].
%
%EXAMPLE:
%Here, we test a known Mersenne prime and a nearby number that is not
%prime.
% boolVal0=isProbablePrime(2^19-1)
% boolVal1=isProbablePrime(2^19-3)
%One should get boolVal0=true and boolVal1=false.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 2: Seminumerical
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%
%October 2018 David Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n<0||~isreal(n)||n~=fix(n))
    error('n must be a real non-negative integer.')
end

if(n>sqrt(flintmax))
   error('The algorithm will not work with n>sqrt(flintmax).') 
end

if(nargin<2||isempty(numRuns))
    numRuns=25;
end

%If divisible by 2.
if(mod(n,2)==0)
   boolVal=false;
   return
end

%We now want to find an odd q such that n-1=2^k*q
q=n-1;
k=0;
pow2Val=1;
while(mod(q,2)==0)
    pow2Val=pow2Val*2;
    q=q/2;
    k=k+1;
end

for curRun=1:numRuns
    %Step P1
    x=randi(n-2)+1;%1<x<n
    %Step P2
    %This calculation method helps avoid overflow.
    j=0;
    %This does y=mod(x^q,n), but in a manner that avoids numerical
    %overflows unless x is very big.
    y=powModuloC(x,q,n);

    while(1)
        %Step P3
        if(y==n-1||y==1&&j==0)
            boolVal=true;
            break;
        end

        %Step P4
        if(y==1&&j>0)
            boolVal=false;
            return
        end

        j=j+1;
        if(j<k)
            y=mod(y^2,n);
        else
            boolVal=false;
            return;
        end
    end
end
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
