function binomTable=makeBinomTable(maxDeg)
%%MAKEBINOMTABLE  Create a lookup table of binomial coefficients.
%                 binomTable(n+1,k+1) is the value of binomial(n,k). This
%                 just build's Pascal's triangle. Consider, for example [1].
%                 The table uses about twice the memory it actually needs
%                 just so that it can be easily accessed as a matrix. This
%                 function is easier to use than the output of the pascal
%                 command in Matlab.
%
%INPUTS: maxDeg The maximum degree of the table. This is a value >=0.
%
%OUTPUTS: binomTable A (maxDeg+1)X(maxDeg+1) matrix where
%               binomTable(n+1,k+1) is the value of binomial(n,k). Unused
%               components are set to zero.
%
%This function can be useful when a large number of low-order binomial
%values are needed. The recursion used to fill in the table is from [1].
%
%REFERENCES:
%[1] Stover, Christopher and Weisstein, Eric W. "Pascal's Triangle." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/PascalsTriangle.html
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    binomTable=zeros(maxDeg+1,maxDeg+1);
    binomTable(:,0+1)=1;%Boundary values.
    for n=1:maxDeg
        for k=1:(n-1)
            binomTable(n+1,k+1)=binomTable(n-1+1,k-1+1)+binomTable(n-1+1,k+1);
        end
        binomTable(n+1,n+1)=1;%Boundary value
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
