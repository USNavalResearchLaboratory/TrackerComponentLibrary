function val=getFactor(n,algorithm)
%%GETFACTOR Given a real positive integer n, get a factor of n >1, if n is
%           not prime, or 1 if n is prime.
%
%INPUTS: n A real positive integer. This value must be <=flintmax.
% algorithm An optional parameter selecting the algorithm to use. Possibel
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Algorithm C in Section 4.5.4 of [1]. The algorithm returns the
%            largets factor of n <=sqrt(n).
%          1 Use the algorithm in Section 2 of [2].
%
%OUTPUTS: val A factor of n or 1 if n is prime.
%
%Faster algorithms for large numbers exist. However, the algorithms used
%here are fairely simple. Algorithm 1 is theoretically faster. However, on
%certain values such as 59138501, algorithm 1 is be much slower than
%algorithm 0.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 2: Seminumerical
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%[2] W. B. Hart, "A one line factoring algorithm," Journal of the
%    Australian Mathematical Society, vol. 92, no. 1, 2012, pp. 61-69.
%
%October 2018 David Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
   algorithm=1; 
end

if(n<0||~isreal(n)||n~=fix(n))
    error('n must be a real non-negative integer.')
end

if(n>flintmax)
   error('The algorithm will not work with n>flintmax.') 
end

switch(algorithm)
    case 0%Algorithm C in [2].
        sqrtN=floor(realsqrt(n));
        
        %Step C1
        a=2*sqrtN+1;
        b=1;
        r=sqrtN^2-n;
        while(r~=0)
            %Step C3
            r=r+a;
            a=a+2;
            %Step C4
            r=r-b;
            b=b+2;
            %Step C5
            while(r>0)
                r=r-b;
                b=b+2;
            end
        end
        
        %Here, r==0.
        val=(a-b)/2;
    case 1%The "one-line" algorithm of [2].
        i=1;
        while(1)
            s=ceil(realsqrt(n*i));
            m=mod(s^2,n);
            sqrtM=fix(realsqrt(m));
            if(sqrtM^2==m)
                val=gcd(n,s-sqrtM);
                return
            end
            i=i+1;
        end
    otherwise
        error('Unknown algorithm specified.')
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
