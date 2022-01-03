function s=accurateSum(a,algorithm,sortVals)
%%ACCURATESUM Sum the elements in a vector in a manner that reduces finite
%             precision errors when summing large and small values as
%             compared to sequentially summing the terms.
%
%INPUTS: a An nX1 or 1Xn vector of real values.
% algorithm An optional parameter that specifies the algorithm to use for
%          the summation. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the 
%            improved Kahan-Babuska procedure as described in [1].
%          1 Use the Kahan-Babuska procedure, which is in [1] and shown to
%            be inferior to the improved procedure.
%   doSort This parameter indicates whether the values in a should be
%          sorted in ascending order by magnitude prior to summing. This
%          also improves the accuracy when given values of highly differing
%          magnitudes. The default if omitted or an empty matrix is passed
%          is false.
%
%OUTPUTS: s The sum of the values in a.
%
%EXAMPLE 1:
%Consider the following summation, whose correct value is 3.
% a=[1,1e100,2,-1e100];
% sMatlab=sum(a)
% s=accurateSum(a)
%As of Matlab 2018b, one will get sMatlab=0 (incorrect) and s=3 (correct).
%Note that Matlab's sum function still has the same result even if one
%sorts the values.
%
%EXAMPLE 2:
%Here, we demonstrate how sorting can improve the accuracy of the results.
% a=[1e100;2;1;eps();eps();eps();-1e100];
% sUnsorted=accurateSum(a,[],false)
% sSorted=accurateSum(a,[],true)
% sMatlab=sum(a)
%One will get sUnsorted=3, sSorted=3.000000000000001, which is more
%correct, and sMatlab=0 (as of Matlab 2018b), which is just wrong.
%
%REFERENCES
%[1] A. Neumaier, "Rundungsfehleranalyse einiger Verfahren zur Summation
%    endlicher Summen," Zeitschrift für angewandte Mathematik und Mechanik,
%    vol. 54, no. 1, 1974, pp. 39-51.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

if(nargin<3||isempty(sortVals))
    sortVals=false;
end

if(sortVals)
   [~,idx]=sort(abs(a),'ascend'); 
   a=a(idx);
end

n=length(a);
switch(algorithm)
    case 0%The improved Kahan-Babuska procedure:
        s=0;
        w=0;
        for m=1:n
            sPrev=s;
            s=a(m)+s;
            if(abs(a(m))<=abs(sPrev))
                w=w+(a(m)+(sPrev-s));%Parentheses are important. 
            else
                w=w+(sPrev+(a(m)-s));%Parentheses are important.
            end
        end
        s=s+w;
    case 1%The Kahan-Babuska procedure
        s=0;
        w=0;
        for m=1:n
           sPrev=s;
           s=a(m)+s;
           w=w+(a(m)+(sPrev-s));%Parentheses are important.
        end
        s=s+w;
    otherwise
        error('Unknown algorithm specified.')
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
