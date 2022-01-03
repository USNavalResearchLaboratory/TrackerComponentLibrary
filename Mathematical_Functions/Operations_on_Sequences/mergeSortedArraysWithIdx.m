function [z,idxX,idxY]=mergeSortedArraysWithIdx(x,y)
%MERGESORTEDARRAYSWITHIDX Given two arrays that are each sorted in
%           ascending order, merge them into one large array that is in
%           ascending order and keep track of the indices of z into which
%           each element of x and y is placed.
%
%INPUTS: x An xDimX1 or 1XxDim array that has been sorted in ascending
%          order.
%        y A yDimX1 or 1XyDim array that has been sorted in ascending
%          order.
%
%OUTPUTS: z An (xDim+yDim)X1 array containing all of the elements of x and
%           y, in ascending order.
%      idxX An xDimX1 vector such that z(i)=x(idxX(i)).
%      idxY A yDimX1 vector such that z(i)=y(idxY(i)).
%
%The algorthm is taken from Chapter 5.2.4 of [1] and modified to keep track
%of the indices where the merged items land.
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
%    ed. Reading, MA: Addison-Wesley, 1998, vol. 3.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=length(x);
n=length(y);
zLen=n+m;
z=zeros(zLen,1);
idxX=zeros(m,1);
idxY=zeros(n,1);

%Initialize
i=1;
j=1;
k=1;

%Find smaller
while(1)
    if(x(i)<=y(j))
        z(k)=x(i);
        idxX(i)=k;
        k=k+1;
        i=i+1;
        if(i<=m)
            continue;
        else
            z(k:zLen)=y(j:n);
            idxY(j:n)=k:zLen;
            break;
        end
    else
        z(k)=y(j);
        idxY(j)=k;
        k=k+1;
        j=j+1;

        if(j<=n)
            continue;
        else
           z(k:zLen)=x(i:m);
           idxX(i:m)=k:zLen;
           break;
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
