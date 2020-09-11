function theArr=randArrangement(n,m)
%%RANDARRANGEMENT Generate a random arrangement of n items put into m
%                 spaces. If n=m, then this is a permutation.
%
%INPUTS: n The number of items to place into spaced.
%        m The number of spaces that can take items. m<=n.
%
%
%OUTPUTS: theArr The mX1 arrangement, which can hold values from 1 to n.

%This function implements the algorithm of [1].
%
%EXAMPLE:
%Here, we generate random arrangement of length 3 from 10 items and
%demonstrate that the distribution of values in each index is apprximately
%uniform.
% numRuns=10000;
% n=10;
% m=3;
% firstIdx=zeros(m,numRuns);
% for curRun=1:numRuns
%     theArr=randArrangement(n,m);
%     firstIdx(:,curRun)=theArr;
% end
% figure(1)
% clf
% histogram(firstIdx(1,:))
% title('First Index')
% figure(2)
% clf
% histogram(firstIdx(2,:))
% title('Second Index')
% figure(3)
% clf
% histogram(firstIdx(3,:))
% title('Third Index')
%
%REFERENCES:
%[1] J. Bentley and B. Floyd, "Programming pearls: A sample of brilliance,"
%   Communications of the ACM, vol. 30, no. 9, pp. 754-757, Sep. 1987.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==0)
    theArr=[];
    return;
end

s=zeros(2*m-1,1);

s(m)=randi(n-m+1);
startIdx=m-1;
endIdx=m+1;
for j=(n-m+2):n
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

theArr=s((startIdx+1):(endIdx-1));

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
