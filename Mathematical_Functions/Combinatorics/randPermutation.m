function thePerm=randPermutation(n)
%%RANDPERMUTATION Generate a random permutation of n values from 0 to n. 
%
%INPUTS: n The length of the permutation.
%
%OUTPUTS: thePerm The nX1 permutation of values from 1 to n.
%
%Floyd's algorithm in [1] is used.
%
%REFERENCES:
%[1] J. Bentley and B. Floyd, "Programming pearls: A sample of brilliance,"
%   Communications of the ACM, vol. 30, no. 9, pp. 754-757, Sep. 1987.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==0)
    thePerm=[];
    return;
end

s=zeros(2*n-1,1);
s(n)=1;
startIdx=n-1;
endIdx=n+1;
for j=2:n
    t=randi(j);
    tIdx=find(t==s((startIdx+1):(endIdx-1)));
    if(isempty(tIdx))
        %If t is not already in s, prefix t to s.
        s(startIdx)=t;
        startIdx=startIdx-1;
    else
        %t is in s, so we add j to s after t.
        s((startIdx+tIdx+2):(endIdx))=s((startIdx+tIdx+1):(endIdx-1));
        
        s(startIdx+tIdx+1)=j;
        endIdx=endIdx+1;
    end
end

thePerm=s((startIdx+1):(endIdx-1));

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
