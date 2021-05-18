function val=greatestCommonDivisor(m,n)
%%GREATESTCOMMONDIVISOR Given two integers m and n, determine the integer
%           that is the greatest common divisor. This is like the gcd
%           function that is built into Matlab, but it can be useful to
%           have an explcit function if one wishes to use a different data
%           type with appropriately overloaded operators, such as for
%           extended precision arithmetic.
%
%INPUTS: m, n Two positive real integers.
%
%OUTPUTS: val The greatest common divisor of the integers.
%
%This function implements Algorithm F in Problem 3 of Section 1.1 of [1].
%
%EXAMPLE:
% val=greatestCommonDivisor(24,32)
%One will get val=8.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 1: Fundamental
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

while(1)
    %Step F1
    m=rem(m,n);

    %Step F2
    if(m==0)
        val=n;
        return
    end

    %Step F3
    n=rem(n,m);

    %Step F4
    if(n==0)
        val=m;
        return
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
