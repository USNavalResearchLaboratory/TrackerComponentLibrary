function val=EulerPhi(n)
%%EULERPHI Compute the Euler toitent function phi(n). This is the number of
%          positive integers n that are relatively prime to n, where 1 is
%          considered to be relatively prime to all numbers. Relatively
%          prime numbers play a role in range and Doppler disambiguation in
%          pulse Doppler radar.
%
%INPUTS: n A nonnegative integer or a matrix of nonnegative integers.
%
%OUTPUTS: val A matrix having the same dimensions as n where each entry if
%             the value of the Euler toitent function for the corresponding
%             n. For the special case of n=0, this is 0; for n=1, it is 1.
%
%The function is implemented as a product involving the unique factors of a
%prime factorization of n as discussed in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Totient Function." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/TotientFunction.html
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numN=numel(n);
val=zeros(size(n));

for curN=1:numN
    nCur=n(curN);
    p=unique(factor(nCur));

    %The fix function deals with finite precision errors. The result should
    %always be an integer.
    val(curN)=fix(nCur*prod(1-1./p));
end

%Deal with the special cases.
val(n==1)=1;
val(n==0)=0;

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
