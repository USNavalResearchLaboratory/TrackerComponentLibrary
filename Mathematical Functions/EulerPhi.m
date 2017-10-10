function val=EulerPhi(n,method)
%%EULERPHI Compute the Euler toitent function phi(n). This is the number of
%          positive integers n that are relatively prime to n, where 1 is
%          considered to be relatively prime to all numbers. Relatively
%          prime numbers play a role in range and Doppler disambiguation in
%          pulse Doppler radar.
%
%INPUTS: n A nonnegative integer or a matrix of nonnegative integers.
%   method An optional parameter that specifies the algorithm to use.
%          Both approaches will be slow for large value sof n. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use a
%             brute-force apprinvolving the unique factors of a prime
%             factorization of n as discussed in [2].
%          1  Use the algorithm given in Chapter 8.1 of [1].
%
%OUTPUTS: val A matrix having the same dimensions as n where each entry if
%             the value of the Euler toitent function for the corresponding
%             n. For the special case of n=0, this is 0; for n=1, it is 1.
%
%REFERENCES:
%[1] R. P. Grimaldi, Discrete and Combinatorial Mathematics: An Applied
%    Introduction, 2nd ed. Reading, MA: Addison-Weslet, 1989.
%[2] Weisstein, Eric W. "Totient Function." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/TotientFunction.html
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(method))
   method=0; 
end

numN=numel(n);
val=zeros(size(n));

if(any(n(:)~=fix(n(:)))||any(n(:)<0)||any(~isfinite(n)))
   error('n Must be composed of finite positive integers.') 
end

switch(method)
    case 0
        for curN=1:numN
            val(curN)=EulerToitentBF(n(curN));
        end
    case 1
        for curN=1:numN
            val(curN)=EulerTotientGrimaldi(n(curN));
        end
    otherwise
        error('Unknown method specified.')
end
end

function phi=EulerToitentBF(n)
%%EULERTOTIENTBF This is a brute-force implementation fo the Euler totient
%%function.
%
%The function is implemented as a product involving the unique factors of a
%prime factorization of n as discussed in [1]. This is essentially a
%brute-force implementation.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Totient Function." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/TotientFunction.html
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(n==0)
        phi=0;
    elseif(n==1)
        phi=1;
    else
        nCur=n;
        p=unique(factor(nCur));

        %The fix function deals with finite precision errors. The result should
        %always be an integer.
        phi=fix(nCur*prod(1-1./p));
    end
end

function phi=EulerTotientGrimaldi(n)
%%EULERTOTIENTGRIMALDI This is the implementation fo the Euler totient
%                      function that is given in Chapter 8.1 of [1].
%
%REFERENCES:
%[1] R. P. Grimaldi, Discrete and Combinatorial Mathematics: An Applied
%    Introduction, 2nd ed. Reading, MA: Addison-Weslet, 1989.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==0)
    phi=0;
    return
end

phi=n;
if(mod(n,2)==0)
    phi=fix(phi/2);
    while(mod(n,2)==0)
       n=n/2; 
    end
end

if(mod(n,3)==0)
    phi=fix(2*phi/3);
    while(mod(n,3)==0)
        n=n/3; 
    end
end

i=5;
while(n>=5)
    j=1;
    while(1)
        j=j+1;
        k=mod(i,j);
        if(k==0||j==fix(sqrt(i)))
            break;
        end
    end
    
    if(k~=0&&mod(n,i)==0)
        phi=fix(phi*(i-1)/i);
        while(mod(n,i)==0)
           n=n/i; 
        end
    end
    i=i+2;
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
