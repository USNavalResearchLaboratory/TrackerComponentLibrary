function [z,idxX,idxY]=mergeSortedArraysUnique(x,y)
%MERGESORTEDARRAYSUNIQUE Given two arrays that are each sorted in ascending
%           order, merge them into one large array that is in ascending
%           order, removing duplicates.
%
%INPUTS: x An xDimX1 or 1XxDim array that has been sorted in ascending
%          order.
%        y A yDimX1 or 1XyDim array that has been sorted in ascending
%          order.
%
%OUTPUTS: z An NX1 array containing all of the elements of x and y, in
%           ascending order. N<=xDim+yDim.
%
%The algorthm is taken from Chapter 5.2.4 of [1] and modified to remove
%duplicates.
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

%Initialize
i=1;
j=1;
k=1;

%Find smaller
zPrev=NaN;
while(1)
    if(x(i)<=y(j))
        if(x(i)==zPrev)
            i=i+1;
        else
            z(k)=x(i);
            zPrev=z(k);
            k=k+1;
            i=i+1;
        end

        if(i<=m)
            continue;
        else
            %Copy in the remaining values, but make sure that no duplicates
            %are added.
            while(j<=n)
                if(y(j)==zPrev)
                    j=j+1;
                else
                    z(k)=y(j);
                    zPrev=z(k);
                    k=k+1;
                    j=j+1;
                end
            end
            break;
        end
    else
        if(y(j)==zPrev)
            j=j+1;
        else
            z(k)=y(j);
            zPrev=z(k);
            k=k+1;
            j=j+1;
        end

        if(j<=n)
            continue;
        else
            %Copy in the remaining values, but make sure that no duplicates
            %are added.
            while(i<=m)
                if(x(i)==zPrev)
                    i=i+1;
                else
                    z(k)=x(i);
                    zPrev=z(k);
                    k=k+1;
                    i=i+1;
                end
            end
            break;
        end
    end
end
%Size to fit.
z=z(1:(k-1));
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
